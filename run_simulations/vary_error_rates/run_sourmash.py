"""
Run sourmash on the simulated genomes

metagenomes are in: metagenomes/
name format: metagenome_<error_rate>_seed_<seed>.fastq

sourmash output files are in: sourmash_output/
name format: sourmash_output_<size>_seed_<seed_value>_k_<ksize>

error rates = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05]
num_runs = 30
seed = list(range(1, num_runs + 1))

sourmash command:
python ../../run_tools/run_sourmash/run_sourmash_wrapper.py 
--ksize <ksize> 
--threshold 1000 
--metagenome <metagenome_name> 
--kosig /scratch/mbr5797/KOs_sketched/KOs_sbt_scaled_1000_k_<ksize>.sbt.zip 
--gatherfile <gather_file_name> 
--outfile <ko_abund_file_name> 
--scaled 1000 
--metagenome_sketch_name <metagenome_sketch_name> 
--resource <resource_usage_filename>
"""

# we need to run sourmash for two k sizes: 11, and 15

def main():
    shell_script_name = "run_sourmash.sh"
    shell_script = open(shell_script_name, "w")
    shell_script.write("#!/bin/bash\n\n")
    shell_script.write("# This script was automatically generated by the script ")
    shell_script.write("run_sourmash.py\n\n")

    counter = 0
    for error_rate in [0.00, 0.01, 0.02, 0.03, 0.04, 0.05]:
        for i in range(30):
            seed = i + 1
            metagenome_filename = "metagenomes/metagenome_" + str(error_rate) + "_seed_" + str(seed) + ".fastq"
            for ksize in [11, 15]:
                gather_file = "sourmash_output/sourmash_gather_" + str(error_rate) + "_seed_" + str(seed) + "_k_" + str(ksize)
                ko_abund_file = "sourmash_output/ko_abund_" + str(error_rate) + "_seed_" + str(seed) + "_k_" + str(ksize)
                resource_file = "sourmash_output/resource_usage_" + str(error_rate) + "_seed_" + str(seed) + "_k_" + str(ksize)
                metagenome_sketch_name = "sourmash_output/metagenome_sketch_" + str(error_rate) + "_seed_" + str(seed) + "_k_" + str(ksize)
                shell_script.write("python ../../run_tools/run_sourmash/run_sourmash_wrapper.py ")
                shell_script.write("--ksize " + str(ksize) + " ")
                shell_script.write("--threshold 1000 ")
                shell_script.write("--metagenome " + metagenome_filename + " ")
                shell_script.write("--kosig /scratch/mbr5797/KOs_sketched/KOs_sbt_scaled_1000_k_" + str(ksize) + ".sbt.zip ")
                shell_script.write("--gatherfile " + gather_file + " ")
                shell_script.write("--outfile " + ko_abund_file + " ")
                shell_script.write("--scaled 1000 ")
                shell_script.write("--metagenome_sketch_name " + metagenome_sketch_name + " ")
                shell_script.write("--resource " + resource_file + " &\n")

                counter += 1
                if counter % 64 == 0:
                    shell_script.write("wait\n")

    shell_script.close()
    