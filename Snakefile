include: "benchmarks/backgrounds/Snakefile"
include: "benchmarks/barrel_ecal/Snakefile"
include: "benchmarks/ecal_gaps/Snakefile"
include: "benchmarks/material_scan/Snakefile"


rule fetch_epic:
    output:
        filepath="EPIC/{PATH}"
    shell: """
xrdcp root://dtn-eic.jlab.org//work/eic2/{output.filepath} {output.filepath}
"""


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
