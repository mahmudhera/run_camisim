# in this script, we will recompute the ground truths for the different sizes of the datasets

import pandas as pd
import subprocess
import os
import pysam
import tqdm

kegg_genomes_directory = '/scratch/mbr5797/genomes_extracted_from_kegg'
gene_to_ko_mapping_filename = '/scratch/mbr5797/genomes_extracted_from_kegg/present_genes_and_koids.csv'


def get_kegg_genome_mapping_filename(kegg_genome_name):
    return kegg_genomes_directory + '/' + kegg_genome_name + '/' + kegg_genome_name + '_mapping.csv'

def get_bam_filename(outdir, simulation_directory_name, genome_name):
    return outdir + '/' + simulation_directory_name + '/bam/' + genome_name + '.bam'

sizes = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
seeds = range(1, 51)

command = 'mkdir new_ground_truths'
subprocess.run(command, shell=True)

# read the gene to ko mapping file
genes_to_kos_df = pd.read_csv(gene_to_ko_mapping_filename)
gene_name_to_ko_id_dict = {}
gene_name_list = genes_to_kos_df['gene_id'].tolist()
ko_id_list = genes_to_kos_df['ko_id'].tolist()
for gene_name, ko_id in zip(gene_name_list, ko_id_list):
    gene_name_to_ko_id_dict[gene_name] = ko_id

for variable in tqdm.tqdm(sizes, desc='Sizes'):
    for seed in tqdm.tqdm(seeds, desc='Seeds'):
        # locate the camisim output directory, format: outputs_of_camisim/out_<size>_seed_<seed>
        outdir = f'outputs_of_camisim/out_{variable}_seed_{seed}'

        # if outdir does not exist, skip
        if not os.path.exists(outdir):
            continue

        # new ground truth file name, format: new_ground_truths/ko_ground_truth_<size>_seed_<seed>
        ko_ground_truth_filename = f'new_ground_truths/ko_ground_truth_{variable}_seed_{seed}'

        # if this ground truth file already exists, skip
        if os.path.exists(ko_ground_truth_filename):
            continue

        # find the simulation directory name
        for directory in os.listdir(outdir):
            if 'sample' in directory:
                simulation_directory_name = directory
                break

        ######################################################################
        # locate the bam files, then convert to a ground truth at gene level #
        ######################################################################
        
        # find the genomes that have been used in the simulation
        source_genomes_directory = outdir + '/source_genomes'
        genome_names_used_in_simulation = []
        filenames = os.listdir(source_genomes_directory)
        for filename in filenames:
            if filename.endswith('.fasta') or filename.endswith('.fna'):
                genome_names_used_in_simulation.append(filename.split('.')[0])

        # records to keep for each gene
        # gene name, and num of reads that mapped to it
        gene_name_to_num_reads_dict = {}
        gene_name_mean_coverage_dict = {}
        gene_name_to_median_coverage_dict = {}
        gene_name_to_num_nucleotides_covered_dict = {}

        # go into the kegg genomes directory and ensure that mapping files exist for all the genomes used in the simulation
        for used_genome_name in genome_names_used_in_simulation:
            # get the file names
            mapping_filename = get_kegg_genome_mapping_filename(used_genome_name)
            bam_filename = get_bam_filename(outdir, simulation_directory_name, used_genome_name)
            
            # check that the files exist
            try:
                assert os.path.exists(mapping_filename)
                assert os.path.exists(bam_filename)
            except AssertionError:
                print(f'Error: {mapping_filename} or {bam_filename} does not exist')
                exit(-1)
        
        gene_name_to_start_end_dict = {}
        # go into the kegg genomes directory and find the corresponding mapping files
        for used_genome_name in genome_names_used_in_simulation:
            # get the file names
            mapping_filename = get_kegg_genome_mapping_filename(used_genome_name)
            bam_filename = get_bam_filename(outdir, simulation_directory_name, used_genome_name)

            # read the bam file using pysam
            bamfile = pysam.AlignmentFile(bam_filename, "rb")

            # read the mapping file as a pandas dataframe
            mapping_df = pd.read_csv(mapping_filename)
            
            # iterate over gene_name, contig_id, start and end positions
            gene_name_list = mapping_df['gene_name'].tolist()
            contig_id_list = mapping_df['contig_id'].tolist()
            start_position_list = mapping_df['start_position'].tolist()
            end_position_list = mapping_df['end_position'].tolist()

            # map gene name to start and end positions
            for gene_name, contig_id, start_position, end_position in zip(gene_name_list, contig_id_list, start_position_list, end_position_list):
                gene_name_to_start_end_dict[gene_name] = (contig_id, start_position, end_position)
            
            for gene_name, contig_id, start_position, end_position in zip(gene_name_list, contig_id_list, start_position_list, end_position_list):
                # for each gene, query the bam file using these intervals and the contig id
                list_of_matches = list(bamfile.fetch(contig_id, start_position, end_position))
                
                # record the matches that have been found in the bam file
                num_reads_in_this_gene = len(list_of_matches)
                if num_reads_in_this_gene == 0:
                    continue
                gene_name_to_num_reads_dict[gene_name] = num_reads_in_this_gene

                # record the mean and coverage of this gene
                coverage_list = [ pileupcolumn.n for pileupcolumn in bamfile.pileup(contig_id, start_position, end_position) ]
                num_nucleotides_covered = len(coverage_list)

                # calculate num zeros: gene length - num nucleotides covered
                gene_length = end_position - start_position + 1
                if gene_length < 0:
                    gene_length = -gene_length
                num_zeros = gene_length - num_nucleotides_covered
                
                # add zeros to the coverage list
                coverage_list += [0] * num_zeros
                
                gene_name_mean_coverage_dict[gene_name] = sum(coverage_list) / len(coverage_list)
                gene_name_to_median_coverage_dict[gene_name] = sorted(coverage_list)[len(coverage_list) // 2]
                gene_name_to_num_nucleotides_covered_dict[gene_name] = num_nucleotides_covered

        ###########################################################
        # use gene to ko mapping to get the ko level ground truth #
        ###########################################################

        ko_abundances_by_num_reads = {}
        ko_abundances_by_num_nts_in_reads = {}
        ko_abundances_by_mean_cov = {}
        ko_abundances_by_med_cov = {}

        koid_to_divide_by = {}

        # iterate over the genes
        for gene_name in gene_name_to_num_reads_dict.keys():
            koid = gene_name_to_ko_id_dict[gene_name]
            num_reads = gene_name_to_num_reads_dict[gene_name]
            num_nts_in_reads = gene_name_to_num_nucleotides_covered_dict[gene_name]
            mean_cov = gene_name_mean_coverage_dict[gene_name]
            med_cov = gene_name_to_median_coverage_dict[gene_name]

            contig_id, start_position, end_position = gene_name_to_start_end_dict[gene_name]
            gene_length = end_position - start_position + 1
            if gene_length < 0:
                gene_length = -gene_length

            # update the ko abundances
            if koid not in ko_abundances_by_num_reads:
                ko_abundances_by_num_reads[koid] = 0
            ko_abundances_by_num_reads[koid] += num_reads

            if koid not in ko_abundances_by_num_nts_in_reads:
                ko_abundances_by_num_nts_in_reads[koid] = 0
            ko_abundances_by_num_nts_in_reads[koid] += num_nts_in_reads

            if koid not in ko_abundances_by_mean_cov:
                ko_abundances_by_mean_cov[koid] = 0
            ko_abundances_by_mean_cov[koid] += mean_cov * gene_length
            
            if koid not in koid_to_divide_by:
                koid_to_divide_by[koid] = 0
            koid_to_divide_by[koid] += gene_length

            if koid not in ko_abundances_by_med_cov:
                ko_abundances_by_med_cov[koid] = 0
            ko_abundances_by_med_cov[koid] += med_cov

        # convert to abundance for mean cov
        for koid in ko_abundances_by_mean_cov.keys():
            ko_abundances_by_mean_cov[koid] /= (1.0*koid_to_divide_by[koid])

        # convert to relative abundance
        total_num_reads = sum(ko_abundances_by_num_reads.values())
        total_num_nts_in_reads = sum(ko_abundances_by_num_nts_in_reads.values())
        total_mean_cov = sum(ko_abundances_by_mean_cov.values())
        total_med_cov = sum(ko_abundances_by_med_cov.values())

        for koid in ko_abundances_by_num_reads.keys():
            ko_abundances_by_num_reads[koid] /= (1.0*total_num_reads)
            ko_abundances_by_num_nts_in_reads[koid] /= (1.0*total_num_nts_in_reads)
            ko_abundances_by_mean_cov[koid] /= (1.0*total_mean_cov)
            ko_abundances_by_med_cov[koid] /= (1.0*total_med_cov)

        # write the ko level ground truth
        with open(ko_ground_truth_filename, 'w') as ground_truth_file:
            ground_truth_file.write('ko_id,abund_by_num_reads,abund_by_num_nts,abund_by_mean_cov,abund_by_med_cov\n')
            for koid in ko_abundances_by_num_reads.keys():
                ground_truth_file.write(f'{koid},{ko_abundances_by_num_reads[koid]},{ko_abundances_by_num_nts_in_reads[koid]},{ko_abundances_by_mean_cov[koid]},{ko_abundances_by_med_cov[koid]}\n')