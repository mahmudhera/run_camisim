
"""
In this script, we will run simulations varying the num of novel genomes. 
We will run 30 simulations for each num of novel genomes. 
We will fix the size at 1.0 Gbp, and number of genomes at 64. We will use an error rate of 0.01.
"""

# command to run simulation:
# python ../../run_camisim.py <num_genomes> --seed xx --size 0.1 
#--config config_<error_rate>_seed_xx.ini 
#--outdir out_<error_rate>_seed_xx 
#--gene_g_t gene_ground_truth_<error_rate>_seed_xx
#--ko_g_t ko_ground_truth_<error_rate>_seed_xx 
#--metagenome_filename metagenome_<error_rate>_seed_xx.fastq 
#--genome_cov_file genome_coverages_<error_rate>_seed_xx
#-e <error_rate> -n wgsim
# -v <num_novel_genomes>

num_simulations = 20
num_genomes = 64

def main():
    shell_script_name = "run_sims_vary_novels.sh"
    shell_script = open(shell_script_name, "w")
    shell_script.write("#!/bin/bash\n\n")
    shell_script.write("# This script was automatically generated by the script ")
    shell_script.write("run_simulations.py\n\n")

    counter = 0

    num_novel_genomes = [0, 2, 4, 8, 16, 32]
    for num_novel_genome in num_novel_genomes:
        for i in range(num_simulations):
            seed = 1 + i
            outdir = f'out_{num_novel_genome}_seed_{seed}'
            gene_ground_truth = f'gene_ground_truth_{num_novel_genome}_seed_{seed}'
            ko_ground_truth = f'ko_ground_truth_{num_novel_genome}_seed_{seed}'
            metagenome_filename = f'metagenome_{num_novel_genome}_seed_{seed}.fastq'
            genome_coverages = f'genome_coverages_{num_novel_genome}_seed_{seed}'
            command = f'python ../../run_camisim.py {num_genomes} --seed {seed} --size 1.0 --config config_{num_novel_genome}_seed_{seed}.ini --outdir {outdir} --gene_g_t {gene_ground_truth} --ko_g_t {ko_ground_truth} --metagenome_filename {metagenome_filename} --genome_cov_file {genome_coverages} -e 0.01 -n wgsim -v {num_novel_genome} &\n'
            shell_script.write(command)

            counter += 1
            if counter % 4 == 0:
                shell_script.write("wait\n")

    shell_script.close()

if __name__:
    main()