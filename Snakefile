configfile: "config.yaml"

if config["remote"] == "S3":
    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    provider = S3RemoteProvider(
        endpoint_url="https://eics3.sdcc.bnl.gov:9000",
        access_key_id=os.environ["S3_ACCESS_KEY"],
        secret_access_key=os.environ["S3_SECRET_KEY"],
    )
    remote_path = lambda path: f"eictest/{path}"
elif config["remote"] == "XRootD":
    from snakemake.remote.XRootD import RemoteProvider as XRootDRemoteProvider
    provider = XRootDRemoteProvider(
        stay_on_remote=False,
    )
    remote_path = lambda path: f"root://dtn-eic.jlab.org//work/eic2/{path}"
else:
    raise ValueError(f"Unexpected config[\"remote\"] = {config['remote']}")

include: "benchmarks/backgrounds/Snakefile"
include: "benchmarks/barrel_ecal/Snakefile"
include: "benchmarks/ecal_gaps/Snakefile"
include: "benchmarks/material_scan/Snakefile"


rule warmup_run:
    output:
        "warmup/{DETECTOR_CONFIG}.edm4hep.root",
    message: "Ensuring that calibrations/fieldmaps are available for {wildcards.DETECTOR_CONFIG}"
    shell: """
ddsim \
  --runType batch \
  --numberOfEvents 1 \
  --compactFile "$DETECTOR_PATH/{wildcards.DETECTOR_CONFIG}.xml" \
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
