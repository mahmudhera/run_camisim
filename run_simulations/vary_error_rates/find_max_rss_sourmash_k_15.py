# sourmash resource usage filename format: sourmash_batch_output/metagenome_<error_rate>_seed_<seed>.fastq_resource_usage_k_15
# resource usage file contains: walltime,cputime,max_rss,avg_rss

error_rates = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
seeds = range(2, 32)

# for every error rate, find the maximum RSS value
for error_rate in error_rates:
    max_rss = 0
    for seed in seeds:
        filename = f'sourmash_batch_output/metagenome_{error_rate}_seed_{seed}.fastq_resource_usage_k_15'
        with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                walltime, cputime, max_rss_str, avg_rss = line.strip().split(',')
                if float(max_rss_str) > max_rss:
                    max_rss = float(max_rss)
    print(f'Error rate: {error_rate}, max RSS: {max_rss}')