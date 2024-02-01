"""
This script will run DIAMOND for each simulated metagenome, and record output and resource usages.
simulation files are in ./metagenomes
config files are in ./configs
output folders are in ./outputs
gene ground truth files are in ./gene_ground_truths
ko ground truth files are in ./ko_ground_truths
genome coverage files are in ./genome_coverages
"""
# command to run diamond:
# run_diamond.py <metagenome_name> <diamond_output_file> <gene_output_file> <ko_output_file> <resource_output_file> --verbose

def main():
    shell_script_name = "run_diamond.sh"
    shell_script = open(shell_script_name, "w")
    shell_script.write("#!/bin/bash\n\n")
    shell_script.write("# This script was automatically generated by the script ")
    shell_script.write("run_diamond.py\n\n")

    for num_genomes in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
        for i in range(50):
            seed = i + 1
            metagenome_filename = "metagenome_" + str(num_genomes) + "_seed_" + str(seed) + ".fastq"
            diamond_output = "diamond_output/diamond_output_" + str(num_genomes) + "_seed_" + str(seed)
            gene_output = "diamond_output/diamond_gene_results_" + str(num_genomes) + "_seed_" + str(seed)
            ko_output = "diamond_output/diamond_ko_results_" + str(num_genomes) + "_seed_" + str(seed)
            resource_output = "diamond_output/diamond_resource_usage_" + str(num_genomes) + "_seed_" + str(seed)
            shell_script.write("python ../../run_tools/run_diamond/run_diamond.py ")
            shell_script.write(metagenome_filename + " ")
            shell_script.write(diamond_output + " ")
            shell_script.write(gene_output + " ")
            shell_script.write(ko_output + " ")
            shell_script.write(resource_output + " ")
            shell_script.write("--verbose\n")

    shell_script.close()

if __name__:
    main()