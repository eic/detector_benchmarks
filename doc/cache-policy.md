Snakemake output cache policy
=============================

This document specifies what the detector benchmarks store in Snakemake's
between-workflow output cache, how a cache entry is named, and under which
conditions the cache is not used. It describes the behaviour of the Snakemake
version in the CI container (9.6.3) together with this repository's
`Snakefile`, `.gitlab-ci.yml` and `garbage-collect.py`.

Update this document whenever a rule gains or loses `cache: True`, a cached
rule changes its `params`/`input`, or the CI cache handling changes.


1. Notation
-----------

The grammar uses EBNF (ISO/IEC 14977):

| Notation | Meaning |
|---|---|
| `a = b ;` | `a` is defined as `b` |
| `a , b` | `a` followed by `b` (concatenation) |
| `a \| b` | `a` or `b` |
| `{ a }` | zero or more `a`, in the order stated in the comment |
| `[ a ]` | optional `a` |
| `"x"` | the literal string `x` |
| `? text ?` | a value defined in prose |
| `(* text *)` | comment |

`sha256hex(x)` is the lowercase hex SHA-256 digest of the UTF-8 bytes of `x`.


2. Cache entries
----------------

```ebnf
cache entry       = cache directory , "/" , provenance hash , output suffix ;

cache directory   = shared cache | pipeline cache ;
shared cache      = ? $SNAKEMAKE_OUTPUT_CACHE as set by the eicweb CI/CD variable,
                      i.e. /scratch/snakemake-cache; shared with physics_benchmarks ? ;
pipeline cache    = ? $LOCAL_DATA_PATH ? , "/snakemake-cache" ;
                    (* /scratch/detector_benchmarks_<CI_PIPELINE_ID>/snakemake-cache *)

output suffix     = ? the last extension of the job's single output file
                      (Python Path(output).suffix), e.g. ".root" for
                      "x.edm4hep.rnt.root", ".hepmc" for "x.hepmc" ? ;
```

Every cached rule has exactly one output file, so a job has exactly one cache
entry. The job's output path in the pipeline (`sim_output/...`, `EPIC/...`)
becomes a symbolic link to the cache entry.


3. Provenance hash (cache key)
------------------------------

```ebnf
provenance hash   = sha256hex( job hash , "0.1" ) ;
                    (* "0.1" is the version of Snakemake's hashing scheme *)

job hash          = sha256hex( code ,
                               { param } ,              (* sorted by param name *)
                               { input content hash } , (* sorted by hash value *)
                               software ,
                               { upstream hash } ) ;    (* sorted by hash value *)

code              = shell template | script source | notebook source | wrapper source ;
shell template    = ? the rule's shell: string before formatting: placeholders such as
                      {wildcards.X}, {params.X}, {input}, {output}, {threads} appear
                      literally, never their values ? ;

param             = [ param name ] , json value ;
param name        = ? name of a params: entry ? ;      (* absent for positional params *)
json value        = ? json.dumps(value, sort_keys=True) of the param's value
                      for this job, after lambdas are evaluated ? ;

input content hash = sha256hex( ? the bytes of an input file, read through
                                  symbolic links ? ) ;
                     (* only for inputs that no job in the DAG produces *)

software          = "" ;
                    (* Snakemake adds container or conda identifiers only with
                       --software-deployment-method apptainer/conda; CI runs without
                       it, so the singularity: directives of the rules are ignored *)

upstream hash     = job hash ;
                    (* of each job that produces one of this job's inputs, cached or
                       not, computed with this same definition, recursively *)
```

Consequences:

 - Only `params` values and the contents of non-generated input files carry
   per-job information into the key. A wildcard that is used only in the
   `shell:` command or in the `input:`/`output:` paths does not change the
   key (see the README section on caching).
 - Software versions enter the key only through:
   - the content of `libepic.so` (`geometry_lib=find_epic_libraries()` input),
   - the `DD4HEP_HASH`, `NPSIM_HASH` and `EICRECON_HASH` params
     (`get_spack_package_hash()`: the Spack DAG hash of the first installation
     `spack find --json <package>` lists, which covers its version, commit,
     variants, compiler and dependencies),
   - upstream hashes of jobs whose key contains one of the above.
 - `libepic.so` embeds its Spack installation prefix in its RPATH. The prefix
   contains the Spack hash of the epic installation, which includes the epic
   commit (CI specs carry a `commit=<sha>` variant), e.g.
   `epic-git.<commit>_main-<hash>` or `epic-main-<hash>`. Its content, and
   therefore every key that includes it, changes with every epic commit.
   epic's compact XML files, field maps and calibration files are not hashed
   themselves.
 - Upgrading Snakemake can change every key (for example, newer versions hash
   inputs as `sha256:<hex>` instead of `<hex>`), which makes existing entries
   unreachable.


4. Cached rules
---------------

```ebnf
cached rule       = eligible rule - excluded rule ;
eligible rule     = ? a rule with "cache: True" (the 57 rules in the table below) ? ;
excluded rule     = ? an eligible rule whose params include a name listed in
                      $BENCHMARKS_CACHE_EXCLUDE_PARAMS ? ;
```

CI does not set `BENCHMARKS_CACHE_EXCLUDE_PARAMS` (section 6), so in CI every
eligible rule is cached.

The table lists every eligible rule in this repository.

 - **epic in key**: the rule has the `libepic.so` content in its own key, gets
   it through the named upstream rule, or does not depend on the epic build
   (—).
 - **Spack hash params**: software identifiers among the rule's `params`.

`fetch_epic` downloads input files (`EPIC/...`) from XRootD or S3. Its key
contains its shell template, which depends on `remote_provider` in
`snakemake.yml`, and the `PATH` param.

| Rule | Snakefile | Kind | epic in key | Spack hash params |
|---|---|---|---|---|
| `backgrounds_sim` | [benchmarks/backgrounds/Snakefile](../benchmarks/backgrounds/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `backwards_ecal_sim` | [benchmarks/backwards_ecal/Snakefile](../benchmarks/backwards_ecal/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `backwards_ecal_recon` | [benchmarks/backwards_ecal/Snakefile](../benchmarks/backwards_ecal/Snakefile) | reconstruction | via `backwards_ecal_sim` | `EICRECON_HASH` |
| `backwards_ecal_dis_recon` | [benchmarks/backwards_ecal/Snakefile](../benchmarks/backwards_ecal/Snakefile) | reconstruction | via `trk_dis_sim` | `EICRECON_HASH` |
| `emcal_barrel_particles` | [benchmarks/barrel_ecal/Snakefile](../benchmarks/barrel_ecal/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `beamline_steering_sim` | [benchmarks/beamline/Snakefile](../benchmarks/beamline/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `beamline_acceptance_sim` | [benchmarks/beamline/Snakefile](../benchmarks/beamline/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `calo_pid_sim` | [benchmarks/calo_pid/Snakefile](../benchmarks/calo_pid/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `calo_pid_recon` | [benchmarks/calo_pid/Snakefile](../benchmarks/calo_pid/Snakefile) | reconstruction | via `calo_pid_sim` | `EICRECON_HASH` |
| `ecal_gaps_sim` | [benchmarks/ecal_gaps/Snakefile](../benchmarks/ecal_gaps/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `ecal_gaps_recon` | [benchmarks/ecal_gaps/Snakefile](../benchmarks/ecal_gaps/Snakefile) | reconstruction | via `ecal_gaps_sim` | `EICRECON_HASH` |
| `far_forward_dvcs_sim` | [benchmarks/far_forward_dvcs/Snakefile](../benchmarks/far_forward_dvcs/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `far_forward_dvcs_reco` | [benchmarks/far_forward_dvcs/Snakefile](../benchmarks/far_forward_dvcs/Snakefile) | reconstruction | via `far_forward_dvcs_sim` | `EICRECON_HASH` |
| `lowq2_training_sim` | [benchmarks/lowq2_reconstruction/Snakefile](../benchmarks/lowq2_reconstruction/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `lowq2_reconstruction_particles_eicrecon` | [benchmarks/lowq2_reconstruction/Snakefile](../benchmarks/lowq2_reconstruction/Snakefile) | reconstruction | via `lowq2_training_sim` | `EICRECON_HASH` |
| `lowq2_reconstruction_tensors_eicrecon` | [benchmarks/lowq2_reconstruction/Snakefile](../benchmarks/lowq2_reconstruction/Snakefile) | reconstruction | via `lowq2_training_sim` | `EICRECON_HASH` |
| `lowq2_steering_reconstruction_particles_trained_eicrecon` | [benchmarks/lowq2_reconstruction/Snakefile](../benchmarks/lowq2_reconstruction/Snakefile) | reconstruction | via `lowq2_training_sim` | `EICRECON_HASH` |
| `tracking_performance_sim` | [benchmarks/tracking_performances/Snakefile](../benchmarks/tracking_performances/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `tracking_performance_recon` | [benchmarks/tracking_performances/Snakefile](../benchmarks/tracking_performances/Snakefile) | reconstruction | via `tracking_performance_sim` | `EICRECON_HASH` |
| `trk_dis_sim` | [benchmarks/tracking_performances_dis/Snakefile](../benchmarks/tracking_performances_dis/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `trk_dis_reco` | [benchmarks/tracking_performances_dis/Snakefile](../benchmarks/tracking_performances_dis/Snakefile) | reconstruction | via `trk_dis_sim` | `EICRECON_HASH` |
| `lfhcal_sim` | [benchmarks/lfhcal/Snakefile](../benchmarks/lfhcal/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `lfhcal_recon` | [benchmarks/lfhcal/Snakefile](../benchmarks/lfhcal/Snakefile) | reconstruction | via `lfhcal_sim` | `EICRECON_HASH` |
| `zdc_lyso_sim` | [benchmarks/zdc_lyso/Snakefile](../benchmarks/zdc_lyso/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `zdc_lyso_reco` | [benchmarks/zdc_lyso/Snakefile](../benchmarks/zdc_lyso/Snakefile) | reconstruction | via `zdc_lyso_sim` | `EICRECON_HASH` |
| `zdc_neutron_sim` | [benchmarks/zdc_neutron/Snakefile](../benchmarks/zdc_neutron/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `zdc_neutron_reco` | [benchmarks/zdc_neutron/Snakefile](../benchmarks/zdc_neutron/Snakefile) | reconstruction | via `zdc_neutron_sim` | `EICRECON_HASH` |
| `insert_muon_simulate` | [benchmarks/insert_muon/Snakefile](../benchmarks/insert_muon/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `insert_muon_recon` | [benchmarks/insert_muon/Snakefile](../benchmarks/insert_muon/Snakefile) | reconstruction | via `insert_muon_simulate` | `EICRECON_HASH` |
| `irt_sim` | [benchmarks/irt/Snakefile](../benchmarks/irt/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `irt_recon` | [benchmarks/irt/Snakefile](../benchmarks/irt/Snakefile) | reconstruction | via `irt_sim` | `EICRECON_HASH` |
| `zdc_lambda_simulate` | [benchmarks/zdc_lambda/Snakefile](../benchmarks/zdc_lambda/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `zdc_lambda_recon` | [benchmarks/zdc_lambda/Snakefile](../benchmarks/zdc_lambda/Snakefile) | reconstruction | via `zdc_lambda_simulate` | `EICRECON_HASH` |
| `zdc_photon_simulate` | [benchmarks/zdc_photon/Snakefile](../benchmarks/zdc_photon/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `zdc_photon_recon` | [benchmarks/zdc_photon/Snakefile](../benchmarks/zdc_photon/Snakefile) | reconstruction | via `zdc_photon_simulate` | `EICRECON_HASH` |
| `zdc_pi0_simulate` | [benchmarks/zdc_pi0/Snakefile](../benchmarks/zdc_pi0/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `zdc_pi0_recon` | [benchmarks/zdc_pi0/Snakefile](../benchmarks/zdc_pi0/Snakefile) | reconstruction | via `zdc_pi0_simulate` | `EICRECON_HASH` |
| `zdc_sigma_simulate` | [benchmarks/zdc_sigma/Snakefile](../benchmarks/zdc_sigma/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `zdc_sigma_recon` | [benchmarks/zdc_sigma/Snakefile](../benchmarks/zdc_sigma/Snakefile) | reconstruction | via `zdc_sigma_simulate` | `EICRECON_HASH` |
| `insert_neutron_simulate` | [benchmarks/insert_neutron/Snakefile](../benchmarks/insert_neutron/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `insert_neutron_recon` | [benchmarks/insert_neutron/Snakefile](../benchmarks/insert_neutron/Snakefile) | reconstruction | via `insert_neutron_simulate` | `EICRECON_HASH` |
| `insert_tau_simulate` | [benchmarks/insert_tau/Snakefile](../benchmarks/insert_tau/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `insert_tau_recon` | [benchmarks/insert_tau/Snakefile](../benchmarks/insert_tau/Snakefile) | reconstruction | via `insert_tau_simulate` | `EICRECON_HASH` |
| `femc_electron_simulate` | [benchmarks/femc_electron/Snakefile](../benchmarks/femc_electron/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `femc_electron_recon` | [benchmarks/femc_electron/Snakefile](../benchmarks/femc_electron/Snakefile) | reconstruction | via `femc_electron_simulate` | `EICRECON_HASH` |
| `femc_photon_simulate` | [benchmarks/femc_photon/Snakefile](../benchmarks/femc_photon/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `femc_photon_recon` | [benchmarks/femc_photon/Snakefile](../benchmarks/femc_photon/Snakefile) | reconstruction | via `femc_photon_simulate` | `EICRECON_HASH` |
| `femc_pi0_simulate` | [benchmarks/femc_pi0/Snakefile](../benchmarks/femc_pi0/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `femc_pi0_recon` | [benchmarks/femc_pi0/Snakefile](../benchmarks/femc_pi0/Snakefile) | reconstruction | via `femc_pi0_simulate` | `EICRECON_HASH` |
| `nhcal_acceptance_simulate` | [benchmarks/nhcal_acceptance/Snakefile](../benchmarks/nhcal_acceptance/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `nhcal_basic_distribution_simulate` | [benchmarks/nhcal_basic_distribution/Snakefile](../benchmarks/nhcal_basic_distribution/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `nhcal_sampling_fraction_simulate` | [benchmarks/nhcal_sampling_fraction/Snakefile](../benchmarks/nhcal_sampling_fraction/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `nhcal_dimuon_photoproduction_simulate` | [benchmarks/nhcal_dimuon_photoproduction/Snakefile](../benchmarks/nhcal_dimuon_photoproduction/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `nhcal_dimuon_photoproduction_recon` | [benchmarks/nhcal_dimuon_photoproduction/Snakefile](../benchmarks/nhcal_dimuon_photoproduction/Snakefile) | reconstruction | via `nhcal_dimuon_photoproduction_simulate` | `EICRECON_HASH` |
| `nhcal_pion_rejection_simulate` | [benchmarks/nhcal_pion_rejection/Snakefile](../benchmarks/nhcal_pion_rejection/Snakefile) | simulation | `libepic.so` content | `DD4HEP_HASH`, `NPSIM_HASH` |
| `nhcal_pion_rejection_recon` | [benchmarks/nhcal_pion_rejection/Snakefile](../benchmarks/nhcal_pion_rejection/Snakefile) | reconstruction | via `nhcal_pion_rejection_simulate` | `EICRECON_HASH` |
| `fetch_epic` | [Snakefile](../Snakefile) | download | — | — |

Rules without `cache: True` (event generators, `warmup_run`, analyses,
`compile_analysis`, `org2py`, `metadata`, …) are never cached. Their job hash
still enters the key of cached jobs that consume their outputs.


5. When the cache is consulted
------------------------------

Snakemake first decides which jobs have to run (output missing, older than an
input, or forced). Only for those jobs, and only for cached rules:

```ebnf
cached job run    = cache hit | cache miss ;
cache hit         = ? the entry exists: the output path becomes a symbolic link
                      to the entry and the job is not executed ? ;
cache miss        = ? the job is executed, its output is moved into the cache
                      directory and replaced by a symbolic link to the entry ? ;
```

Jobs whose output already exists and is up to date neither read nor write the
cache. A job that must run because an input was regenerated (for example
`warmup.edm4hep.root`, which each CI job recreates) does consult the cache, so
it normally gets the entry stored by an earlier job instead of re-running.

Before running Snakemake, every CI job removes symbolic links under
`sim_output/` and `EPIC/` whose cache entry was deleted, so that their outputs
are fetched or regenerated.


6. When the shared cache is skipped
-----------------------------------

The cache directory is chosen in the `before_script` of `.det_benchmark`
(inherited by every CI job that runs Snakemake), based on `GITHUB_REPOSITORY`
and `GITHUB_PR`, which the eic/epic and eic/EICrecon GitHub workflows pass
through eic/containers to the pipelines it triggers. `GITHUB_PR` is the pull
request number for `pull_request` events and empty for all other events
(merge queue, tags, branch pushes, manual runs):

```ebnf
policy            = epic policy | default policy ;

epic policy       = ( GITHUB_REPOSITORY = "eic/epic" ) , ( GITHUB_PR <> "" ) ,
                    ( cache directory = pipeline cache ) ;

default policy    = ? any other GITHUB_REPOSITORY (including empty), or GITHUB_PR = "" ? ,
                    ( cache directory = shared cache ) ;
```

Rationale: a pipeline that tests an epic or EICrecon pull request builds the
PR's head commit (`epic@git.<sha>=main`, `eicrecon@git.<sha>=main`), so its
keys contain that commit's `libepic.so` or `EICRECON_HASH`. The head commit
never becomes `main` as is (epic only allows squash merges through its merge
queue, and GitHub never fast-forwards a pull request), so no later pipeline
builds the same software and such entries could never be retrieved again.

 - For an epic PR, every cached output depends on the epic build, so the whole
   cache moves into the pipeline. Jobs of the same pipeline still share
   entries, as described in section 5.
 - An EICrecon PR uses the shared cache, including for reconstruction, even
   though its reconstruction entries are never retrieved by later pipelines
   (section 8). Its simulations do not depend on EICrecon and are shared with
   other pipelines. Reconstruction stays cached because jobs of the same
   pipeline share outputs under `sim_output/`: when a job reruns a cached
   reconstruction it only replaces the output with a link to the cache entry,
   but when it reruns an uncached one it rewrites the shared output in place
   while another job may be reading it.
 - An epic merge-queue run tests the squash commit that `main` is then
   fast-forwarded to, which the next nightly container build also builds, so
   it uses the shared cache. Tag, branch and manual runs also use the shared
   cache: they are rare, and their entries are evicted by garbage collection.

Complete list of conditions under which an output is neither read from nor
written to the shared cache:

| # | Condition | Effect |
|---|---|---|
| 1 | The rule has no `cache: True` | Never cached |
| 2 | Snakemake runs without `--cache`: any run outside CI unless `--cache` is passed, and the `RUN_EXTENDED_BENCHMARK` steps in `benchmarks/barrel_ecal/config.yml`, which do not pass `$SNAKEMAKE_FLAGS` | Nothing cached; Snakemake warns that the rules are "eligible for caching" |
| 3 | `--cache` is given but `SNAKEMAKE_OUTPUT_CACHE` is not set | Nothing cached; Snakemake warns that no cache location is set |
| 4 | The job's output exists and is up to date | Cache not consulted (section 5) |
| 5 | `GITHUB_REPOSITORY` is `eic/epic` and `GITHUB_PR` is not empty | All cached rules use the pipeline cache instead of the shared cache |
| 6 | The rule has a param listed in `BENCHMARKS_CACHE_EXCLUDE_PARAMS` (never set in CI) | The rule is not cached at all |


7. Removal of entries
---------------------

```ebnf
entry removal     = pipeline scratch removal | stale scratch removal | garbage collection ;

pipeline scratch removal = ? the cleanup:scratch job of a successful pipeline removes
                             /scratch/detector_benchmarks_<CI_PIPELINE_ID>, including
                             its pipeline cache and its symbolic links ? ;

stale scratch removal    = ? common:setup removes /scratch/detector_benchmarks_* and
                             /scratch/physics_benchmarks_* directories modified more than
                             72 hours ago, including their pipeline caches ? ;

garbage collection       = ? common:setup runs garbage-collect.py on the shared cache
                             with a 500G target size; see below ? ;
```

Stale scratch removal and garbage collection run under
`flock --nonblock /scratch/garbage-collect.lock`. They are skipped when another
pipeline holds the lock.

`garbage-collect.py` lists the shared cache, then:

 1. removes empty regular files (never valid entries) unconditionally,
 2. never removes entries modified less than 6 hours ago (grace period),
 3. never removes entries that a symbolic link anywhere under `/scratch`
    points to ("alive"),
 4. removes the remaining entries in decreasing order of
    `size × (now − access time)` until the total size of the cache is below
    the target.

The pipeline cache is not garbage collected; it lives exactly as long as the
pipeline's scratch directory.


8. Known gaps
-------------

 - A pipeline that tests an EICrecon PR stores its reconstruction outputs in
   the shared cache, where no later pipeline retrieves them, until garbage
   collection evicts them (section 6 explains why).
 - When two pipelines miss the cache for the same entry at the same time, the
   one that finishes last replaces the entry with its own, newer file. The
   other pipeline's outputs that depend on that entry then look outdated, so
   its later jobs rerun them, and through the regenerated
   `warmup.edm4hep.root` everything else that depends on `warmup_run`. For
   cached rules these reruns are cache hits.
 - If `ctypes.util.find_library("epic")` finds no library, `geometry_lib` is
   empty and keys silently lose the epic dependency.
 - If `spack` is missing or `spack find` fails, `get_spack_package_hash()`
   returns `""` (with a warning) and the key loses that package. If several
   installations of a package exist (e.g. the eic_xl container ships several
   epic versions), the first one listed is used, which need not be the one in
   `DETECTOR_PREFIX`.
 - physics_benchmarks writes to the same shared cache and does not apply the
   policy in section 6.
