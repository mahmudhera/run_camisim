"""
Run sourmash on the simulated genomes

metagenomes are in: metagenomes/
name format: metagenome_<num_genomes>_seed_<seed_value>.fastq

num_genomes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
seed_value = [1, 2, ..., 50]

sourmash command:
python ../../run_tools/run_sourmash/run_sourmash_wrapper.py 
--ksize <ksize> 
--threshold 1000 
--metagenome <metagenome_name> 
--kosig /scratch/mbr5797/KOs_sketched/KOs_sbt_scaled_1000_k_11.sbt.zip 
--gatherfile <gather_file_name> 
--outfile <ko_abund_file_name> 
--scaled 100
--metagenome_sketch_name <metagenome_sketch_name> 
--resource <resource_usage_filename>
"""

# we need to run sourmash for two k sizes: 11, and 15

def main():
    shell_script_name = "run_sourmash_sensitive.sh"
    shell_script = open(shell_script_name, "w")
    shell_script.write("#!/bin/bash\n\n")
    shell_script.write("# This script was automatically generated by the script ")
    shell_script.write("run_sourmash.py\n\n")

    counter = 0

    for num_genomes in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
        for i in range(50):
            seed = i + 1
            metagenome_filename = "metagenomes/metagenome_" + str(num_genomes) + "_seed_" + str(seed) + ".fastq"
            
            # ksize = 11
            gather_file = "sourmash_sensitive_output/gather_output_" + str(num_genomes) + "_seed_" + str(seed) + "_k_11"
            ko_abund_file = "sourmash_sensitive_output/ko_abund_output_" + str(num_genomes) + "_seed_" + str(seed) + "_k_11"
            resource_file = "sourmash_sensitive_output/resource_usage_" + str(num_genomes) + "_seed_" + str(seed) + "_k_11"
            metagenome_sketch_name = "sourmash_sensitive_output/metagenome_sketch_" + str(num_genomes) + "_seed_" + str(seed) + "_k_11"
            shell_script.write("python ../../run_tools/run_sourmash/run_sourmash_wrapper.py ")
            shell_script.write("--ksize 11 ")
            shell_script.write("--threshold 100 ")
            shell_script.write("--metagenome " + metagenome_filename + " ")
            shell_script.write("--kosig /scratch/mbr5797/KOs_sketched/KOs_sbt_scaled_1000_k_11.sbt.zip ")
            shell_script.write("--gatherfile " + gather_file + " ")
            shell_script.write("--outfile " + ko_abund_file + " ")
            shell_script.write("--scaled 1000 ")
            shell_script.write("--metagenome_sketch_name " + metagenome_sketch_name + " ")
            shell_script.write("--resource " + resource_file + " &\n")

            counter += 1

            # ksize = 15
            gather_file = "sourmash_sensitive_output/gather_output_" + str(num_genomes) + "_seed_" + str(seed) + "_k_15"
            ko_abund_file = "sourmash_sensitive_output/ko_abund_output_" + str(num_genomes) + "_seed_" + str(seed) + "_k_15"
            resource_file = "sourmash_sensitive_output/resource_usage_" + str(num_genomes) + "_seed_" + str(seed) + "_k_15"
            metagenome_sketch_name = "sourmash_sensitive_output/metagenome_sketch_" + str(num_genomes) + "_seed_" + str(seed) + "_k_15"
            shell_script.write("python ../../run_tools/run_sourmash/run_sourmash_wrapper.py ")
            shell_script.write("--ksize 15 ")
            shell_script.write("--threshold 100 ")
            shell_script.write("--metagenome " + metagenome_filename + " ")
            shell_script.write("--kosig /scratch/mbr5797/KOs_sketched/KOs_sbt_scaled_1000_k_15.sbt.zip ")
            shell_script.write("--gatherfile " + gather_file + " ")
            shell_script.write("--outfile " + ko_abund_file + " ")
            shell_script.write("--scaled 1000 ")
            shell_script.write("--metagenome_sketch_name " + metagenome_sketch_name + " ")
            shell_script.write("--resource " + resource_file + " &\n")

            counter += 1

            if counter % 32 == 0:
                shell_script.write("wait\n")

    shell_script.close()

if __name__:
    main()