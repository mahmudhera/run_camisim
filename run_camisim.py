import argparse
import datetime
import configparser
import subprocess
import os
import pysam
import pandas as pd

# hardcoded items
# designed for the GPU machine
#camisim_path = '/home/grads/mbr5797/camisim/CAMISIM-1.3/metagenome_from_profile.py'
camisim_path = '/home/grads/mbr5797/camisim/CAMISIM-1.3/metagenomesimulation.py'
kegg_genomes_directory = '/scratch/mbr5797/genomes_extracted_from_kegg'

def get_kegg_genome_mapping_filename(kegg_genome_name):
    return kegg_genomes_directory + '/' + kegg_genome_name + '/' + kegg_genome_name + '_mapping.csv'

def get_bam_filename(outdir, simulation_directory_name, genome_name):
    return outdir + '/' + simulation_directory_name + '/bam/' + genome_name + '.bam'

def main():
    # parse arguments
    parser = argparse.ArgumentParser(description='Process number of genomes.')
    parser.add_argument('number_of_genomes', type=int, help='Number of genomes')
    parser.add_argument('--config', type=str, help='Config file')
    parser.add_argument('--seed', type=int, help='Seed used in simulation', default=0)
    parser.add_argument('--outdir', type=str, help='Output directory', default='./out')
    parser.add_argument('--size', type=float, help='Size of the file in Gbp', default=0.1)
    parser.add_argument('--ground_truth', type=str, help='Ground truth file', default='./ground_truth.csv')
    args = parser.parse_args()

    # read the arguments
    number_of_genomes = args.number_of_genomes
    if args.config:
        config_filename = args.config
    else:
        timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        config_filename = 'config_' + timestamp + '.ini'
    seed = args.seed
    outdir = args.outdir
    size = args.size
    ground_truth_filename = args.ground_truth

    # read the sample config file
    config = configparser.ConfigParser()
    config.read('sample_config.ini')

    # update the num_genomes field with the command line argument
    # section: Main
    config.set('Main', 'seed', str(seed))
    config.set('Main', 'output_directory', outdir)

    # section: ReadSimulator
    config.set('ReadSimulator', 'size', str(size))

    # section: community0
    config.set('community0', 'genomes_total', str(number_of_genomes))
    config.set('community0', 'genomes_real', str(number_of_genomes))

    # write the updated config to the config filename
    with open(config_filename, 'w') as config_file:
        config.write(config_file)

    # run camisim
    #subprocess.run(['rm', '-rf', outdir])
    #subprocess.run(['mkdir', outdir])
    #cmd = f'python {camisim_path} -seed {seed} -id KL -p 64 {config_filename}'
    #print(cmd)
    #subprocess.run(cmd.split(' '))

    # find the simulation directory name
    for directory in os.listdir(outdir):
        if 'sample' in directory:
            simulation_directory_name = directory
            break

    #########################################################   
    # locate the bam files, then convert to a ground truth ##
    #########################################################
    # find the genomes that have been used in the simulation
    source_genomes_directory = outdir + '/source_genomes'
    genome_names_used_in_simulation = []
    filenames = os.listdir(source_genomes_directory)
    for filename in filenames:
        if filename.endswith('.fasta') or filename.endswith('.fna'):
            genome_names_used_in_simulation.append(filename.split('.')[0])
    
    # go into the kegg genomes directory and find the corresponding mapping files
    for used_genome_name in genome_names_used_in_simulation:
        # get the file names
        mapping_filename = get_kegg_genome_mapping_filename(used_genome_name)
        bam_filename = get_bam_filename(outdir, simulation_directory_name, used_genome_name)
        
        # check that the files exist
        assert os.path.exists(mapping_filename)
        assert os.path.exists(bam_filename)

        # read the bam file using pysam
        bamfile = pysam.AlignmentFile(bam_filename, "rb")

        # read the mapping file as a pandas dataframe
        mapping_df = pd.read_csv(mapping_filename)
        
        # iterate over gene_name, contig_id, start and end positions
        gene_name_list = mapping_df['gene_name'].tolist()
        contig_id_list = mapping_df['contig_id'].tolist()
        start_position_list = mapping_df['start_position'].tolist()
        end_position_list = mapping_df['end_position'].tolist()
        
        gene_name_to_num_reads_dict = {}
        for gene_name, contig_id, start_position, end_position in zip(gene_name_list, contig_id_list, start_position_list, end_position_list):
            # for each gene, query the bam file using these intervals and the contig id
            list_of_matches = list(bamfile.fetch(contig_id, start_position, end_position))
            
            # record the matches that have been found in the bam file
            num_reads_in_this_gene = len(list_of_matches)
            gene_name_to_num_reads_dict[gene_name] = num_reads_in_this_gene

        # use all the matches to create a ground truth file by studying what we did before
        # relative_abundance of genes by number of reads coming from that gene
        total_num_of_reads = sum(gene_name_to_num_reads_dict.values())
        gene_name_to_relative_abundance_by_num_reads_dict = {}
        for gene_name, num_reads in gene_name_to_num_reads_dict.items():
            gene_name_to_relative_abundance_by_num_reads_dict[gene_name] = 1.0 * num_reads / total_num_of_reads
    
    # sort the gene names based on decreasing order of the relative abundance
    gene_names_sorted_by_relative_abundance = sorted(gene_name_to_relative_abundance_by_num_reads_dict.keys(), key=lambda x: gene_name_to_relative_abundance_by_num_reads_dict[x], reverse=True)

    with open(ground_truth_filename, 'w') as ground_truth_file:
        for gene_name in gene_names_sorted_by_relative_abundance:
            relative_abundance = gene_name_to_relative_abundance_by_num_reads_dict[gene_name]
            ground_truth_file.write(f'{gene_name},{relative_abundance}\n')



if __name__ == '__main__':
    main()
