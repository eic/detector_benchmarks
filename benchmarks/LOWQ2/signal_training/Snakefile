# Snakefile

# Rule to run the simulation using Allpix2
rule run_simulation:
    input:
        "allpix2_config.cfg"  # Input file containing position and momentum
    output:
        "output/simulation_output.root"  # Output file from the simulation
    shell:
        "/home/simong/geant4/allpix-squared/bin/allpix -c {input} -o {output}"

# Rule to train the GAN using the simulation output
rule train_gan:
    input:
        "simulation_output.root"  # Output file from the simulation
    output:
        "trained_gan_model.h5"  # Trained GAN model file
    shell:
        "python train_gan_signal.py --input {input} --output {output}"

# Default rule to run the simulation and train the GAN
rule all:
    input:
        "trained_gan_model.h5"  # Trained GAN model file
    output:
        "final_results.root"  # Final results file
    shell:
        "python generate_results.py --model {input} --output {output}"