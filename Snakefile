configfile: "snakemake.yml"

import functools
import os
import re
from snakemake.logging import logger


_rucio_upload_enabled = bool(config.get("rucio", {}).get("host"))
if not _rucio_upload_enabled:
    logger.warning(
        "Rucio host not configured (config.rucio.host); "
        "benchmark outputs will not be uploaded to Rucio."
    )

_rucio_scope     = config.get("rucio", {}).get("scope", "epic")
_rucio_job_id    = config.get("job_id", os.environ.get("CI_PIPELINE_ID", "local"))
_RUCIO_DID_PREFIX = f"rucio://{_rucio_scope}/VALIDATION/benchmarks/detector/{_rucio_job_id}"


def _rucio_output(path):
    """Return a Rucio storage() output, or a local temp() placeholder when
    Rucio is not configured (e.g. in local test runs).

    When storage() is used, {output.rucio} in shell commands resolves to a
    *local* staging path managed by Snakemake — a plain cp is all that is
    needed.  Snakemake uploads the staged file to the remote Rucio server
    after the job completes successfully.
    """
    if _rucio_upload_enabled:
        return storage(path)
    local = re.sub(r"^rucio://[^/]*/", ".rucio_staging/", path)
    return temp(local)


rule compile_analysis:
    input:
        "{path}/{filename}.cxx",
    output:
        "{path}/{filename}_cxx.d",
        "{path}/{filename}_cxx.so",
        "{path}/{filename}_cxx_ACLiC_dict_rdict.pcm",
    shell:
        """
root -l -b -q -e '.L {input}+'
"""


@functools.cache
def get_spack_package_hash(package_name):
    import json
    try:
        ver_info = json.loads(subprocess.check_output(["spack", "find", "--json", package_name]))
        return ver_info[0]["hash"]
    except FileNotFoundError as e:
        logger.warning("Spack is not installed")
        return ""
    except subprocess.CalledProcessError as e:
        print(e)
        return ""


@functools.cache
def find_epic_libraries():
    import ctypes.util
    # if library is not found (not avaliable) we return an empty list to let DAG still evaluate
    libs = []
    lib = ctypes.util.find_library("epic")
    if lib is not None:
        libs.append(os.environ["DETECTOR_PATH"] + "/../../lib/" + lib)
    return libs


include: "benchmarks/backgrounds/Snakefile"
include: "benchmarks/backwards_ecal/Snakefile"
include: "benchmarks/barrel_ecal/Snakefile"
include: "benchmarks/beamline/Snakefile"
include: "benchmarks/calo_pid/Snakefile"
include: "benchmarks/campaign/Snakefile"
include: "benchmarks/ecal_gaps/Snakefile"
include: "benchmarks/far_forward_dvcs/Snakefile"
include: "benchmarks/lowq2_reconstruction/Snakefile"
include: "benchmarks/material_scan/Snakefile"
include: "benchmarks/secondary_vertexing_dis/Snakefile"
include: "benchmarks/tracking_performances/Snakefile"
include: "benchmarks/tracking_performances_dis/Snakefile"
include: "benchmarks/lfhcal/Snakefile"
include: "benchmarks/zdc_lyso/Snakefile"
include: "benchmarks/zdc_neutron/Snakefile"
include: "benchmarks/insert_muon/Snakefile"
include: "benchmarks/zdc_lambda/Snakefile"
include: "benchmarks/zdc_photon/Snakefile"
include: "benchmarks/zdc_pi0/Snakefile"
include: "benchmarks/zdc_sigma/Snakefile"
include: "benchmarks/insert_neutron/Snakefile"
include: "benchmarks/insert_tau/Snakefile"
include: "benchmarks/femc_electron/Snakefile"
include: "benchmarks/femc_photon/Snakefile"
include: "benchmarks/femc_pi0/Snakefile"
include: "benchmarks/nhcal_acceptance/Snakefile"
include: "benchmarks/nhcal_basic_distribution/Snakefile"

use_s3 = config["remote_provider"].lower() == "s3"
use_xrootd = config["remote_provider"].lower() == "xrootd"
use_rucio = config["remote_provider"].lower() == "rucio"

_xrootd_cfg = config.get("xrootd", {})
_xrootd_server = _xrootd_cfg.get("server", "root://dtn-eic.jlab.org")
_xrootd_base = _xrootd_cfg.get("base_path", "/volatile/eic").strip("/")
# Full base URL used by get_remote_path() and the fetch_epic rule.
_xrootd_base_url = f"{_xrootd_server}//{_xrootd_base}"


def get_remote_path(path):
    if use_s3:
        return f"s3https://eics3.sdcc.bnl.gov:9000/eictest/{path}"
    elif use_xrootd or use_rucio:
        # EIC-XRD is an XRootD-backed RSE; the configured XRootD server is
        # used for both xrootd and rucio remote_provider modes.
        return f"{_xrootd_base_url}/{path}"
    else:
        raise RuntimeError(f'Unexpected value for config["remote_provider"]: {config["remote_provider"]}')


rule fetch_epic:
    output:
        filepath="EPIC/{PATH}"
    params:
        # wildcards are not included in hash for caching, we need to add them as params
        PATH=lambda wildcards: wildcards.PATH,
        xrd_base=_xrootd_base_url,
    cache: True
    retries: 3
    shell: """
xrdcp --debug 2 {params.xrd_base}/{output.filepath} {output.filepath}
""" if use_xrootd or use_rucio else """
mc cp S3/eictest/{output.filepath} {output.filepath}
""" if use_s3 else f"""
echo 'Unexpected value for config["remote_provider"]: {config["remote_provider"]}'
exit 1
"""


rule warmup_run:
    output:
        "warmup.edm4hep.root",
    message: "Ensuring that calibrations/fieldmaps are available"
    shell: """
set -m # monitor mode to prevent lingering processes
exec ddsim \
  --runType batch \
  --numberOfEvents 1 \
  --compactFile "$DETECTOR_PATH/epic_ip6.xml" \
  --outputFile "{output}" \
  --enableGun
"""


rule matplotlibrc:
    output:
        ".matplotlibrc",
    run:
        with open(output[0], "wt") as fp:
            fp.write("backend: Agg\n")
            # interactive mode prevents plt.show() from blocking
            fp.write("interactive : True\n")


rule org2py:
    input:
        notebook=workflow.basedir + "/{NOTEBOOK}.org",
        converter=workflow.source_path("benchmarks/common/org2py.awk"),
    output:
        "{NOTEBOOK}.py"
    shell:
        """
awk -f {input.converter} {input.notebook} > {output}
"""


rule metadata:
    output:
        "results/metadata.json"
    shell:
        """
cat > {output} <<EOF
{{
  "CI_COMMIT_REF_NAME": "${{CI_COMMIT_REF_NAME:-}}",
  "CI_COMMIT_SHA": "${{CI_COMMIT_SHA:-}}",
  "CI_PIPELINE_ID": "${{CI_PIPELINE_ID:-}}",
  "CI_PIPELINE_SOURCE": "${{CI_PIPELINE_SOURCE:-}}",
  "CI_PROJECT_ID": "${{CI_PROJECT_ID:-}}",
  "GITHUB_REPOSITORY": "${{GITHUB_REPOSITORY:-}}",
  "GITHUB_SHA": "${{GITHUB_SHA:-}}",
  "GITHUB_PR": "${{GITHUB_PR:-}}",
  "PIPELINE_NAME": $(echo "${{PIPELINE_NAME:-}}" | jq -Rs .)
}}
EOF
# validate JSON
jq '.' {output}
"""
