# PNG plot outputs for benchmarks
# Benchmarks that produce PDFs will auto-convert to PNG

DETECTOR_CONFIGS = config.get("detector_configs", ["epic_craterlake"])
CAMPAIGNS = config.get("campaigns", ["local"])

# PDF to PNG conversion
rule pdf_to_png:
    input: "{path}.pdf"
    output: "{path}.png"
    shell: "convert -density 150 {input} -quality 90 {output}"

# Benchmark plot outputs (PNG only)
BENCHMARK_PLOTS = {
    "zdc_neutron": lambda cfg: [
        f"results/zdc_neutron/{cfg}/fwd_neutrons_geant.png",
        f"results/zdc_neutron/{cfg}/fwd_neutrons_recon.png",
    ],
    "beamline": lambda cfg: [
        f"benchmarks/beamline/analysis/beamspot_{CAMPAIGNS[0]}.png",
        f"benchmarks/beamline/analysis/x_px_{CAMPAIGNS[0]}.png",
        f"benchmarks/beamline/analysis/y_py_{CAMPAIGNS[0]}.png",
    ],
}

rule test_benchmark:
    """Run single benchmark: snakemake test_benchmark --config benchmark=zdc_neutron"""
    input:
        lambda wildcards: BENCHMARK_PLOTS.get(
            config.get("benchmark", "zdc_neutron"),
            lambda cfg: []
        )(config.get("detector_config", DETECTOR_CONFIGS[0]))

rule all_plots:
    """Generate all benchmark plots"""
    input:
        [plot for benchmark_func in BENCHMARK_PLOTS.values() 
         for cfg in DETECTOR_CONFIGS for plot in benchmark_func(cfg)]
         for plot in benchmark_func(cfg)]
