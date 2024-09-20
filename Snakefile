configfile: "snakemake.yml"

include: "benchmarks/backgrounds/Snakefile"
include: "benchmarks/backwards_ecal/Snakefile"
include: "benchmarks/barrel_ecal/Snakefile"
include: "benchmarks/ecal_gaps/Snakefile"
include: "benchmarks/material_scan/Snakefile"
include: "benchmarks/tracking_performances/Snakefile"
include: "benchmarks/tracking_performances_dis/Snakefile"
include: "benchmarks/zdc_lyso/Snakefile"
include: "benchmarks/insert_muon/Snakefile"


use_s3 = config["remote_provider"].lower() == "s3"
use_xrootd = config["remote_provider"].lower() == "xrootd"


def get_remote_path(path):
    if use_s3:
        return f"s3https://eics3.sdcc.bnl.gov:9000/eictest/{path}"
    elif use_xrootd:
        return f"root://dtn-eic.jlab.org//work/eic2/{path}"
    else:
        raise runtime_exception('Unexpected value for config["remote_provider"]: {config["remote_provider"]}')


rule fetch_epic:
    output:
        filepath="EPIC/{PATH}"
    cache: True
    retries: 3
    shell: """
xrdcp --debug 2 root://dtn-eic.jlab.org//work/eic2/{output.filepath} {output.filepath}
""" if use_xrootd else """
mc cp S3/eictest/{output.filepath} {output.filepath}
""" if use_s3 else f"""
echo 'Unexpected value for config["remote_provider"]: {config["remote_provider"]}'
exit 1
"""


rule warmup_run:
    output:
        "warmup/{DETECTOR_CONFIG}.edm4hep.root",
    message: "Ensuring that calibrations/fieldmaps are available for {wildcards.DETECTOR_CONFIG}"
    shell: """
set -m # monitor mode to prevent lingering processes
exec ddsim \
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
