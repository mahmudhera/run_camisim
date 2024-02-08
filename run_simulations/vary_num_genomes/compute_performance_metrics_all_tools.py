"""
In this script, we will compute the performance metrics for all tools.
Tool names:
- diamond (fast)
- diamond (sensitive)
- sourmash, k=11
- sourmash, k=15

num_genomes: 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024
num_runs: 50
seed: 1, 2, 3, ..., 50

diamond_fast ko files: diamond_fast_output/diamond_ko_results_<num_genomes>_seed_<seed>
diamond_sensitive ko files: diamond_sensitive_output/diamond_ko_results_<num_genomes>_seed_<seed>
sourmash ko files: ko_abund_output_/<num_genomes>_seed_<seed>_k_<ksize>
ground_truth_ko_files: ko_ground_truths/ko_ground_truth_<num_genomes>_seed_<seed>

Output file naming convention:
- diamond fast: diamond_fast_output/diamond_performance_metrics_<num_genomes>_seed_<seed>
- diamond sensitive: diamond_sensitive_output/diamond_performance_metrics_<num_genomes>_seed_<seed>
- sourmash: sourmash_output/sourmash_performance_metrics_<num_genomes>_seed_<seed>_k_<ksize>

command to compute perf metrics for one tool:
python ../../compute_presence_absence_metrics.py <ground_truth_ko_file> <predicted_ko_file> <output_file> --toolname <toolname>
the toolname argument can be: diamond, sourmash
"""

import os
import pandas as pd
import numpy as np

num_genomes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
num_runs = 50
seed = list(range(1, num_runs + 1))

# invoke runs for diamond fast
print("Invoking runs for diamond fast")
progress = 0
for ng in num_genomes:
    for s in seed:
        ground_truth_ko_file = "ko_ground_truths/ko_ground_truth_" + str(ng) + "_seed_" + str(s)
        predicted_ko_file = "diamond_fast_output/diamond_ko_results_" + str(ng) + "_seed_" + str(s)
        output_file = "diamond_fast_output/diamond_performance_metrics_" + str(ng) + "_seed_" + str(s)
        os.system("python ../../compute_presence_absence_metrics.py " + ground_truth_ko_file + " " + predicted_ko_file + " " + output_file + " --toolname diamond")
        progress += 1
        print("Progress: " + str(progress) + " / " + str(len(num_genomes) * len(seed)))

# invoke runs for diamond sensitive
print("Invoking runs for diamond sensitive")
progress = 0
for ng in num_genomes:
    for s in seed:
        ground_truth_ko_file = "ko_ground_truths/ko_ground_truth_" + str(ng) + "_seed_" + str(s)
        predicted_ko_file = "diamond_sensitive_output/diamond_ko_results_" + str(ng) + "_seed_" + str(s)
        output_file = "diamond_sensitive_output/diamond_performance_metrics_" + str(ng) + "_seed_" + str(s)
        os.system("python ../../compute_presence_absence_metrics.py " + ground_truth_ko_file + " " + predicted_ko_file + " " + output_file + " --toolname diamond")
        progress += 1
        print("Progress: " + str(progress) + " / " + str(len(num_genomes) * len(seed)))

# invoke runs for sourmash, k=11,15
print("Invoking runs for sourmash")
progress = 0
for ng in num_genomes:
    for s in seed:
        for ksize in [11, 15]:
            ground_truth_ko_file = "ko_ground_truths/ko_ground_truth_" + str(ng) + "_seed_" + str(s)
            predicted_ko_file = "sourmash_output/ko_abund_output_" + str(ng) + "_seed_" + str(s) + "_k_" + str(ksize)
            output_file = "sourmash_output/sourmash_performance_metrics_" + str(ng) + "_seed_" + str(s) + "_k_" + str(ksize)
            os.system("python ../../compute_presence_absence_metrics.py " + ground_truth_ko_file + " " + predicted_ko_file + " " + output_file + " --toolname sourmash")
            progress += 1
            print("Progress: " + str(progress) + " / " + str(len(num_genomes) * len(seed) * 2))
