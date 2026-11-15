# Remote backend for Snakemake's between-workflow output cache (--cache).
#
# Snakemake can keep its output cache in remote storage only when a default
# storage provider is set (--default-storage-provider), which would also move
# every other input and output of the workflow to that storage. This file
# instead selects the remote cache from $SNAKEMAKE_OUTPUT_CACHE alone:
#
#   SNAKEMAKE_OUTPUT_CACHE=pelican://<federation>/<namespace>/<dir>  -> pelican plugin
#   SNAKEMAKE_OUTPUT_CACHE=osdf:///<namespace>/<dir>                -> pelican plugin
#   SNAKEMAKE_OUTPUT_CACHE=root[s]://<host>[:<port>]//<path>/<dir>  -> xrootd plugin
#   SNAKEMAKE_OUTPUT_CACHE=<local path>                             -> unchanged (local cache)
#
# BENCHMARKS_CACHE_STORAGE_PROVIDER overrides the plugin (e.g. "fs" for tests).
# The plugin (snakemake-storage-plugin-pelican or -xrootd) must be installed.
#
# Cache entries are stored as <SNAKEMAKE_OUTPUT_CACHE>/<provenance hash><suffix>,
# the same names as in the local cache (doc/cache-policy.md). By default a cache
# hit downloads the entry to the job's output path.
#
# Remote-only mode (BENCHMARKS_CACHE_REMOTE_ONLY=1): for rules listed in
# REMOTE_CACHE_READY_RULES, the job's output path holds a small stub file with
# the entry's read URL instead of the data, both on a cache hit and after a job
# ran and uploaded its output. Consumers turn the path into something ROOT,
# podio or uproot can open with cache_resolve (shell, see Snakefile) or
# remote_cache_resolve (Python). The read URL is the entry's location with the
# prefix $SNAKEMAKE_OUTPUT_CACHE replaced by $BENCHMARKS_CACHE_READ_PREFIX
# (default: $SNAKEMAKE_OUTPUT_CACHE), so that e.g. entries written through
# pelican:// are read over root://.
#
# A rule may be added to REMOTE_CACHE_READY_RULES only when every consumer of
# its output resolves the path first.

import os

from snakemake.exceptions import CacheMissException
from snakemake.logging import logger

REMOTE_CACHE_STUB_MAGIC = "#snakemake-remote-cache-entry\n"

# Cached rules whose consumers all read their output through cache_resolve
REMOTE_CACHE_READY_RULES = set()


def remote_cache_resolve(path):
    """Return the read URL stored in a remote-only cache stub, or path itself."""
    try:
        with open(path, "r", encoding="ascii") as fp:
            if fp.read(len(REMOTE_CACHE_STUB_MAGIC)) == REMOTE_CACHE_STUB_MAGIC:
                return fp.readline().strip()
    except (UnicodeDecodeError, IsADirectoryError):
        pass
    return str(path)


def _remote_cache_install():
    import inspect
    import snakemake.workflow

    location = os.environ.get("SNAKEMAKE_OUTPUT_CACHE", "")
    provider_name = os.environ.get("BENCHMARKS_CACHE_STORAGE_PROVIDER")
    if not provider_name:
        scheme = location.split("://", 1)[0] if "://" in location else None
        provider_name = {
            "pelican": "pelican",
            "osdf": "pelican",
            "root": "xrootd",
            "roots": "xrootd",
        }.get(scheme)
    if not provider_name:
        return  # local cache, nothing to change

    read_prefix = os.environ.get("BENCHMARKS_CACHE_READ_PREFIX", location)
    remote_only = os.environ.get("BENCHMARKS_CACHE_REMOTE_ONLY", "0") not in ("", "0")
    provider = workflow.storage_registry.register_storage(
        provider=provider_name, tag="output_cache"
    )
    StorageOutputFileCache = snakemake.workflow.StorageOutputFileCache

    async def _maybe_await(x):
        return await x if inspect.isawaitable(x) else x

    class RemoteOutputFileCache:
        """Wraps Snakemake's storage cache; method signatures differ between
        Snakemake versions (9.6 passes an extra cache_mode), hence *args."""

        def __init__(self):
            self.inner = StorageOutputFileCache(provider)
            # set on the first failed remote operation; later jobs skip the
            # remote instead of each waiting through the client's retries
            self.unreachable = False

        def __getattr__(self, name):
            return getattr(self.inner, name)

        async def _entries(self, job, *args, **kwargs):
            gen = self.inner._get_storage_objects(job, *args, **kwargs)
            if hasattr(gen, "__aiter__"):
                return [entry async for entry in gen]
            return list(gen)

        def _remote_only(self, job):
            return remote_only and job.rule.name in REMOTE_CACHE_READY_RULES

        def _write_stub(self, entry):
            url = read_prefix + entry.query[len(location):]
            path = str(entry.local_path())
            if os.path.lexists(path):
                os.remove(path)
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
            with open(path, "w", encoding="ascii") as fp:
                fp.write(REMOTE_CACHE_STUB_MAGIC + url + "\n")

        # An unreachable remote cache is treated as a cache miss: the job runs
        # and keeps its output locally, instead of failing the pipeline.

        def _disable(self, what, e):
            self.unreachable = True
            logger.warning(f"{what} failed, not using the remote cache for the rest of this run: {e}")

        async def exists(self, job, *args):
            if self.unreachable:
                return False
            try:
                return await self.inner.exists(job, *args)
            except Exception as e:
                self._disable(f"Remote cache lookup for {job}", e)
                return False

        async def fetch(self, job, *args):
            if self.unreachable:
                self.inner.raise_cache_miss_exception(job)
            try:
                if not self._remote_only(job):
                    return await self.inner.fetch(job, *args)
                entries = await self._entries(job, *args)
                for entry in entries:
                    if not await entry.managed_exists():
                        self.inner.raise_cache_miss_exception(job)
            except CacheMissException:
                raise
            except Exception as e:
                self._disable(f"Remote cache fetch for {job}", e)
                raise CacheMissException(str(e))
            for entry in entries:
                self._write_stub(entry)
            logger.info(f"Remote cache hit for {job}: outputs are stubs pointing to {read_prefix}")

        async def store(self, job, *args):
            if self.unreachable:
                return
            entries = await self._entries(job, *args, check_output_exists=True)
            for entry in entries:
                try:
                    await entry.managed_store()
                except Exception as e:
                    # another pipeline may have stored the same entry concurrently
                    try:
                        stored = await entry.managed_exists()
                    except Exception:
                        stored = False
                    if not stored:
                        self._disable(f"Upload of {entry.query} (output kept local)", e)
                        return
                    logger.warning(f"Cache entry {entry.query} already exists: {e}")
            if self._remote_only(job):
                for entry in entries:
                    self._write_stub(entry)

    # Workflow.execute() constructs the cache via this module-level name when no
    # default storage provider is set.
    snakemake.workflow.LocalOutputFileCache = RemoteOutputFileCache
    logger.info(
        f"Output cache: {location} via storage plugin '{provider_name}'"
        + (f", remote-only reads from {read_prefix}" if remote_only else "")
    )


_remote_cache_install()
