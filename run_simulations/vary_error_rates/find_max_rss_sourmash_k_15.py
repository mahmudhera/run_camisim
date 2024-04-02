# sourmash resource usage filename format: sourmash_batch_output/metagenome_<error_rate>_seed_<seed>.fastq_resource_usage_k_15
# resource usage file contains: walltime,cputime,max_rss,avg_rss

error_rates = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05]
seeds = range(2, 32)

# for every error rate, find the maximum RSS value
for error_rate in error_rates:
    max_rss_list = []
    for seed in seeds:
        filename = f'sourmash_batch_output/metagenome_{error_rate}_seed_{seed}.fastq_resource_usage_k_15'
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
        line = lines[0]
        walltime, cputime, max_rss_str, avg_rss = line.strip().split(',')
        max_rss = float(max_rss_str)
        max_rss_list.append(max_rss)
    print(f'Error rate: {error_rate}, max RSS: {max(max_rss_list)}')