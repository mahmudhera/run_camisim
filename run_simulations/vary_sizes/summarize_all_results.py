"""
In this script, we will summarize the results of all the presence-absence metrics.
We have 50 runs for each number of genomes.
We will compute the mean and standard deviation of the following for each number of genomes.
- purity
- completeness
- kendall's tau
- pearsonr coefficient for common kos
- pearsonr coefficient for all kos
- weighted purity
- weighted completeness
- KL divergence of common kos in ground truth to predictions
- KL divergence of common kos in predictions to ground truth
- bray-curtis distance

sizes = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
seed = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ..., 50]

tools:
- diamond (fast) (performance metrics are in: diamond_fast_output)
- sourmash, k = 11 (performance metrics are in: sourmash_output)
- sourmash, k = 15 (performance metrics are in: sourmash_output)

file naming convention for resources used by each tool:
- diamond fast: diamond_fast_output/diamond_resource_usage_<size>_seed_<seed>
- sourmash: sourmash_output/resource_usage_<size>_seed_<seed>_k_<ksize>
these files are csv files. header: none, columns: walltime, cputime, max_rss, avg_rss

file naming convention for performance metrics:
- diamond fast: diamond_fast_output/diamond_performance_metrics_<size>_seed_<seed>
- sourmash: sourmash_output/sourmash_performance_metrics_<size>_seed_<seed>_k_<ksize>
these files are csv files, header: first line, column names: purity,completeness,kendalltau,pearsonr_common,pearsonr_all,wt_purity,wt_completeness,kl_div_common_gt_to_pred,kl_div_common_pred_to_gt,bray_curtis
"""

import os
import pandas as pd
import numpy as np

sizes = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
num_runs = 50
seed = list(range(1, num_runs + 1))

# create a dataframe to store the performance metrics
performance_metrics = pd.DataFrame(columns=['purity', 'completeness', 'kendalltau', 'pearsonr_common', 'pearsonr_all', 'wt_purity', 'wt_completeness', 'kl_div_common_gt_to_pred', 'kl_div_common_pred_to_gt', 'bray_curtis', 'size', 'tool'])

# read the performance metrics for diamond fast
for size in sizes:
    for s in seed:
        diamond_fast_file = "diamond_fast_output/diamond_performance_metrics_" + str(size) + "_seed_" + str(s)
        if os.path.exists(diamond_fast_file):
            df = pd.read_csv(diamond_fast_file)
            df['size'] = size
            df['tool'] = 'diamond_fast'
            performance_metrics = pd.concat( [performance_metrics, df], ignore_index=True)

# read the performance metrics for sourmash
for size in sizes:
    for s in seed:
        for ksize in [11, 15]:
            sourmash_file = "sourmash_output/sourmash_performance_metrics_" + str(size) + "_seed_" + str(s) + "_k_" + str(ksize)
            if os.path.exists(sourmash_file):
                df = pd.read_csv(sourmash_file)
                df['size'] = size
                df['tool'] = 'sourmash,k=' + str(ksize)
                performance_metrics = pd.concat( [performance_metrics, df], ignore_index=True)

# compute the mean and standard deviation of the performance metrics
mean_performance_metrics = performance_metrics.groupby(['size', 'tool']).mean().reset_index()
std_performance_metrics = performance_metrics.groupby(['size', 'tool']).std().reset_index()

# write the mean and standard deviation to a file
mean_performance_metrics.to_csv('mean_performance_metrics.csv', index=False)
std_performance_metrics.to_csv('std_performance_metrics.csv', index=False)

# create a dataframe to store the resources used
resources_used = pd.DataFrame(columns=['walltime', 'cputime', 'max_rss', 'avg_rss', 'size', 'tool'])

# read the resources used by diamond fast
for size in sizes:
    for s in seed:
        diamond_fast_file = "diamond_fast_output/diamond_resource_usage_" + str(size) + "_seed_" + str(s)
        if os.path.exists(diamond_fast_file):
            df = pd.read_csv(diamond_fast_file, header=None)
            df.columns = ['walltime', 'cputime', 'max_rss', 'avg_rss']
            df['size'] = size
            df['tool'] = 'diamond_fast'
            resources_used = pd.concat( [resources_used, df], ignore_index=True)

# read the resources used by sourmash for k = 11 and k = 15
for size in sizes:
    for s in seed:
        for ksize in [11, 15]:
            sourmash_file = "sourmash_output/resource_usage_" + str(size) + "_seed_" + str(s) + "_k_" + str(ksize)
            if os.path.exists(sourmash_file):
                df = pd.read_csv(sourmash_file, header=None)
                df.columns = ['walltime', 'cputime', 'max_rss', 'avg_rss']
                df['size'] = size
                df['tool'] = 'sourmash,k=' + str(ksize)
                resources_used = pd.concat( [resources_used, df], ignore_index=True)

# compute the mean and standard deviation of the resources used
mean_resources_used = resources_used.groupby(['size', 'tool']).mean().reset_index()
std_resources_used = resources_used.groupby(['size', 'tool']).std().reset_index()

# write the mean and standard deviation to a file
mean_resources_used.to_csv('mean_resources_used.csv', index=False)
std_resources_used.to_csv('std_resources_used.csv', index=False)