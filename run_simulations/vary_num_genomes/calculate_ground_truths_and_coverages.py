# command to run simulation:
# python ../../run_camisim.py
#--outdir out_num_genomes_seed_xx 
#--gene_g_t gene_ground_truth_num_genomes_seed_xx
#--ko_g_t ko_ground_truth_num_genomes_seed_xx 
#--genome_cov_file genome_coverages_num_genomes_seed_xx

num_simulations = 50

def main():
    shell_script_name = "calculate_ground_truths.sh"
    shell_script = open(shell_script_name, "w")
    shell_script.write("#!/bin/bash\n\n")
    shell_script.write("# This script was automatically generated by the script ")
    shell_script.write("calculate_ground_truths_and_coverages.py\n\n")

    counter = 0

    for num_genomes in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
        for i in range(num_simulations):
            seed = i + 1
            shell_script.write("python ../../calculate_ground_truth.py ")
            shell_script.write("--outdir out_" + str(num_genomes) + "_seed_" + str(seed) + " ")
            shell_script.write("--gene_g_t gene_ground_truth_" + str(num_genomes) + "_seed_" + str(seed) + " ")
            shell_script.write("--ko_g_t ko_ground_truth_" + str(num_genomes) + "_seed_" + str(seed) + " ")
            shell_script.write("--genome_cov_file genome_coverages_" + str(num_genomes) + "_seed_" + str(seed) + " &\n")

            counter += 1
            if counter % 32 == 0:
                shell_script.write("wait\n\n")

    shell_script.close()

if __name__:
    main()