"""
In this script, we will compute the performance metrics for all tools.
Tool names:
- diamond (fast)
- sourmash, k=11
- sourmash, k=15

num_novel_genomes = [0, 2, 4, 8, 16, 32]
num_runs: 20
seed: 1, 2, 3, ..., 20

diamond_fast ko files: diamond_fast_output/diamond_ko_results_<num_novel>_seed_<seed>
sourmash ko files: ko_abund_output_/ko_abund_<num_novel>_seed_<seed>_k_<ksize>
ground_truth_ko_files: ko_ground_truths/ko_ground_truth_<num_novel>_seed_<seed>

Output file naming convention:
- diamond fast: diamond_fast_output/diamond_performance_metrics_<num_novel>_seed_<seed>
- sourmash: sourmash_output/sourmash_performance_metrics_<num_novel>_seed_<seed>_k_<ksize>

command to compute perf metrics for one tool:
python ../../compute_presence_absence_metrics.py <ground_truth_ko_file> <predicted_ko_file> <output_file> --toolname <toolname>
the toolname argument can be: diamond, sourmash
"""

import os
import pandas as pd
import numpy as np

num_novel_genomes = [0, 2, 4, 8, 16, 32]
num_runs = 20
seed = list(range(1, num_runs + 1))

# invoke runs for diamond fast
print("Invoking runs for diamond fast")
progress = 0
for num_novel in num_novel_genomes:
    for s in seed:
        ground_truth_ko_file = "ko_ground_truths/ko_ground_truth_" + str(num_novel) + "_seed_" + str(s)
        predicted_ko_file = "diamond_fast_output/diamond_ko_results_" + str(num_novel) + "_seed_" + str(s)
        output_file = "diamond_fast_output/diamond_performance_metrics_" + str(num_novel) + "_seed_" + str(s)
        try:
            os.system("python ../../compute_presence_absence_metrics.py " + ground_truth_ko_file + " " + predicted_ko_file + " " + output_file + " --toolname diamond &")
        except:
            print("Error encountered for size: " + str(num_novel) + " seed: " + str(s))
        progress += 1
        print("Progress: " + str(progress) + " / " + str(len(num_novel) * len(seed)))

# invoke runs for sourmash, k=11,15
print("Invoking runs for sourmash")
progress = 0
for num_novel in num_novel_genomes:
    for s in seed:
        for ksize in [11, 15]:
            ground_truth_ko_file = "ko_ground_truths/ko_ground_truth_" + str(num_novel) + "_seed_" + str(s)
            predicted_ko_file = "sourmash_output/ko_abund_" + str(num_novel) + "_seed_" + str(s) + "_k_" + str(ksize)
            output_file = "sourmash_output/sourmash_performance_metrics_" + str(num_novel) + "_seed_" + str(s) + "_k_" + str(ksize)
            try:
                os.system("python ../../compute_presence_absence_metrics.py " + ground_truth_ko_file + " " + predicted_ko_file + " " + output_file + " --toolname sourmash &")
            except:
                print("Error encountered for size: " + str(num_novel) + " seed: " + str(s) + " ksize: " + str(ksize))
            progress += 1
            print("Progress: " + str(progress) + " / " + str(len(num_novel_genomes) * len(seed) * 2))