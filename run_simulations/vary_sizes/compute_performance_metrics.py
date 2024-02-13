"""
In this script, we will compute the performance metrics for all tools.
Tool names:
- diamond (fast)
- sourmash, k=11
- sourmash, k=15

sizes: 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4
num_runs: 50
seed: 1, 2, 3, ..., 50

diamond_fast ko files: diamond_fast_output/diamond_ko_results_<size>_seed_<seed>
sourmash ko files: ko_abund_output_/<size>_seed_<seed>_k_<ksize>
ground_truth_ko_files: ko_ground_truths/ko_ground_truth_<size>_seed_<seed>

Output file naming convention:
- diamond fast: diamond_fast_output/diamond_performance_metrics_<size>_seed_<seed>
- sourmash: sourmash_output/sourmash_performance_metrics_<size>_seed_<seed>_k_<ksize>

command to compute perf metrics for one tool:
python ../../compute_presence_absence_metrics.py <ground_truth_ko_file> <predicted_ko_file> <output_file> --toolname <toolname>
the toolname argument can be: diamond, sourmash
"""

import os
import pandas as pd
import numpy as np

sizes = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
num_runs = 50
seed = list(range(1, num_runs + 1))

# invoke runs for diamond fast
print("Invoking runs for diamond fast")
progress = 0
for size in sizes:
    for s in seed:
        ground_truth_ko_file = "ko_ground_truths/ko_ground_truth_" + str(size) + "_seed_" + str(s)
        predicted_ko_file = "diamond_fast_output/diamond_ko_results_" + str(size) + "_seed_" + str(s)
        output_file = "diamond_fast_output/diamond_performance_metrics_" + str(size) + "_seed_" + str(s)
        try:
            os.system("python ../../compute_presence_absence_metrics.py " + ground_truth_ko_file + " " + predicted_ko_file + " " + output_file + " --toolname diamond")
        except:
            print("Error encountered for size: " + str(size) + " seed: " + str(s))
        progress += 1
        print("Progress: " + str(progress) + " / " + str(len(sizes) * len(seed)))


# invoke runs for sourmash, k=11,15
print("Invoking runs for sourmash")
progress = 0
for size in sizes:
    for s in seed:
        for ksize in [11, 15]:
            ground_truth_ko_file = "ko_ground_truths/ko_ground_truth_" + str(size) + "_seed_" + str(s)
            predicted_ko_file = "sourmash_output/ko_abund_" + str(size) + "_seed_" + str(s) + "_k_" + str(ksize)
            output_file = "sourmash_output/sourmash_performance_metrics_" + str(size) + "_seed_" + str(s) + "_k_" + str(ksize)
            try:
                os.system("python ../../compute_presence_absence_metrics.py " + ground_truth_ko_file + " " + predicted_ko_file + " " + output_file + " --toolname sourmash")
            except:
                print("Error encountered for size: " + str(size) + " seed: " + str(s) + " ksize: " + str(ksize))
            progress += 1
            print("Progress: " + str(progress) + " / " + str(len(sizes) * len(seed) * 2))
