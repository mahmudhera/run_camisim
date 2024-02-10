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

Num genomes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
seed = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ..., 50]

tools:
- diamond (fast) (performance metrics are in: diamond_fast_output)
- diamond (sensitive) (performance metrics are in: diamond_sensitive_output)
- sourmash, k = 11 (performance metrics are in: sourmash_output)
- sourmash, k = 15 (performance metrics are in: sourmash_output)

file naming convention for resources used by each tool:
- diamond fast: diamond_fast_output/diamond_resource_usage_<num_genomes>_seed_<seed>
- diamond sensitive: diamond_sensitive_output/diamond_resource_usage_<num_genomes>_seed_<seed>
- sourmash: sourmash_output/resource_usage_<num_genomes>_seed_<seed>_k_<ksize>
these files are csv files. header: none, columns: walltime, cputime, max_rss, avg_rss

file naming convention for performance metrics:
- diamond fast: diamond_fast_output/diamond_performance_metrics_<num_genomes>_seed_<seed>
- diamond sensitive: diamond_sensitive_output/diamond_performance_metrics_<num_genomes>_seed_<seed>
- sourmash: sourmash_output/sourmash_performance_metrics_<num_genomes>_seed_<seed>_k_<ksize>
these files are csv files, header: first line, column names: purity,completeness,kendalltau,pearsonr_common,pearsonr_all,wt_purity,wt_completeness,kl_div_common_gt_to_pred,kl_div_common_pred_to_gt,bray_curtis
"""

import os
import pandas as pd
import numpy as np

num_genomes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
num_runs = 50
seed = list(range(1, num_runs + 1))

# create a dataframe to store the performance metrics
performance_metrics = pd.DataFrame(columns=['purity', 'completeness', 'kendalltau', 'pearsonr_common', 'pearsonr_all', 'wt_purity', 'wt_completeness', 'kl_div_common_gt_to_pred', 'kl_div_common_pred_to_gt', 'bray_curtis', 'num_genomes', 'tool'])

# read the performance metrics for diamond fast
for ng in num_genomes:
    for s in seed:
        diamond_fast_file = "diamond_fast_output/diamond_performance_metrics_" + str(ng) + "_seed_" + str(s)
        if os.path.exists(diamond_fast_file):
            df = pd.read_csv(diamond_fast_file)
            df['num_genomes'] = ng
            df['tool'] = 'diamond_fast'
            performance_metrics = pd.concat( [performance_metrics, df], ignore_index=True)

# read the performance metrics for diamond sensitive
for ng in num_genomes:
    for s in seed:
        diamond_sensitive_file = "diamond_sensitive_output/diamond_performance_metrics_" + str(ng) + "_seed_" + str(s)
        if os.path.exists(diamond_sensitive_file):
            df = pd.read_csv(diamond_sensitive_file)
            df['num_genomes'] = ng
            df['tool'] = 'diamond_sensitive'
            performance_metrics = pd.concat( [performance_metrics, df], ignore_index=True)

# read the performance metrics for sourmash
for ng in num_genomes:
    for s in seed:
        sourmash_file_11 = "sourmash_output/sourmash_performance_metrics_" + str(ng) + "_seed_" + str(s) + "_k_11"
        if os.path.exists(sourmash_file_11):
            df = pd.read_csv(sourmash_file_11)
            df['num_genomes'] = ng
            df['tool'] = 'sourmash,k=11'
            performance_metrics = pd.concat( [performance_metrics, df], ignore_index=True)
            
        sourmash_file_15 = "sourmash_output/sourmash_performance_metrics_" + str(ng) + "_seed_" + str(s) + "_k_15"
        if os.path.exists(sourmash_file_15):
            df = pd.read_csv(sourmash_file_15)
            df['num_genomes'] = ng
            df['tool'] = 'sourmash,k=15'
            performance_metrics = pd.concat( [performance_metrics, df], ignore_index=True)

# compute the mean and standard deviation of the performance metrics
mean_performance_metrics = performance_metrics.groupby(['num_genomes', 'tool']).mean().reset_index()
std_performance_metrics = performance_metrics.groupby(['num_genomes', 'tool']).std().reset_index()

# write the mean and standard deviation to a file
mean_performance_metrics.to_csv('mean_performance_metrics.csv', index=False)
std_performance_metrics.to_csv('std_performance_metrics.csv', index=False)

# create a dataframe to store the resources used
resources_used = pd.DataFrame(columns=['walltime', 'cputime', 'max_rss', 'avg_rss', 'num_genomes', 'tool'])

# read the resources used by diamond fast
for ng in num_genomes:
    for s in seed:
        diamond_fast_file = "diamond_fast_output/diamond_resource_usage_" + str(ng) + "_seed_" + str(s)
        if os.path.exists(diamond_fast_file):
            df = pd.read_csv(diamond_fast_file, header=None)
            df.columns = ['walltime', 'cputime', 'max_rss', 'avg_rss']
            df['num_genomes'] = ng
            df['tool'] = 'diamond_fast'
            resources_used = pd.concat( [resources_used, df], ignore_index=True)

# read the resources used by diamond sensitive
for ng in num_genomes:
    for s in seed:
        diamond_sensitive_file = "diamond_sensitive_output/diamond_resource_usage_" + str(ng) + "_seed_" + str(s)
        if os.path.exists(diamond_sensitive_file):
            df = pd.read_csv(diamond_sensitive_file, header=None)
            df.columns = ['walltime', 'cputime', 'max_rss', 'avg_rss']
            df['num_genomes'] = ng
            df['tool'] = 'diamond_sensitive'
            resources_used = pd.concat( [resources_used, df], ignore_index=True)

# read the resources used by sourmash for k = 11 and k = 15
for ng in num_genomes:
    for s in seed:
        sourmash_file_11 = "sourmash_output/resource_usage_" + str(ng) + "_seed_" + str(s) + "_k_11"
        if os.path.exists(sourmash_file_11):
            df = pd.read_csv(sourmash_file_11, header=None)
            df.columns = ['walltime', 'cputime', 'max_rss', 'avg_rss']
            df['num_genomes'] = ng
            df['tool'] = 'sourmash,k=11'
            resources_used = pd.concat( [resources_used, df], ignore_index=True)

        sourmash_file_15 = "sourmash_output/resource_usage_" + str(ng) + "_seed_" + str(s) + "_k_15"
        if os.path.exists(sourmash_file_15):
            df = pd.read_csv(sourmash_file_15, header=None)
            df.columns = ['walltime', 'cputime', 'max_rss', 'avg_rss']
            df['num_genomes'] = ng
            df['tool'] = 'sourmash,k=15'
            resources_used = pd.concat( [resources_used, df], ignore_index=True)

# compute the mean and standard deviation of the resources used
mean_resources_used = resources_used.groupby(['num_genomes', 'tool']).mean().reset_index()
std_resources_used = resources_used.groupby(['num_genomes', 'tool']).std().reset_index()

# write the mean and standard deviation to a file
mean_resources_used.to_csv('mean_resources_used.csv', index=False)
std_resources_used.to_csv('std_resources_used.csv', index=False)