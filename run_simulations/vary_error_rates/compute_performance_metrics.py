"""
In this script, we will compute the performance metrics for all tools.
Tool names:
- diamond (fast)
- sourmash, k=11
- sourmash, k=15

error_rates = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05]
num_runs: 30
seed: 1, 2, 3, ..., 30

diamond_fast ko files: diamond_fast_output/diamond_ko_results_<error_rate>_seed_<seed>
sourmash ko files: ko_abund_output_/ko_abund_<error_rate>_seed_<seed>_k_<ksize>
ground_truth_ko_files: ko_ground_truths/ko_ground_truth_<error_rate>_seed_<seed>

Output file naming convention:
- diamond fast: diamond_fast_output/diamond_performance_metrics_<error_rate>_seed_<seed>
- sourmash: sourmash_output/sourmash_performance_metrics_<error_rate>_seed_<seed>_k_<ksize>

command to compute perf metrics for one tool:
python ../../compute_presence_absence_metrics.py <ground_truth_ko_file> <predicted_ko_file> <output_file> --toolname <toolname>
the toolname argument can be: diamond, sourmash
"""

import os
import pandas as pd
import numpy as np

error_rates = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05]
num_runs = 30
seed = list(range(2, num_runs + 2))

# invoke runs for diamond fast
print("Invoking runs for diamond fast")
progress = 0
for error_rate in error_rates:
    for s in seed:
        ground_truth_ko_file = "ko_ground_truths/ko_ground_truth_" + str(error_rate) + "_seed_" + str(s)
        predicted_ko_file = "diamond_fast_output/diamond_ko_results_" + str(error_rate) + "_seed_" + str(s)
        output_file = "diamond_fast_output/diamond_performance_metrics_" + str(error_rate) + "_seed_" + str(s)
        try:
            os.system("python ../../compute_presence_absence_metrics.py " + ground_truth_ko_file + " " + predicted_ko_file + " " + output_file + " --toolname diamond &")
        except:
            print("Error encountered for size: " + str(error_rate) + " seed: " + str(s))
        progress += 1
        print("Progress: " + str(progress) + " / " + str(len(error_rates) * len(seed)))

        if progress % 255 == 0:
            os.system("wait")


# invoke runs for sourmash, k=11,15
print("Invoking runs for sourmash")
progress = 0
for error_rate in error_rates:
    for s in seed:
        for ksize in [11, 15]:
            ground_truth_ko_file = "ko_ground_truths/ko_ground_truth_" + str(error_rate) + "_seed_" + str(s)
            predicted_ko_file = "sourmash_output/ko_abund_" + str(error_rate) + "_seed_" + str(s) + "_k_" + str(ksize)
            output_file = "sourmash_output/sourmash_performance_metrics_" + str(error_rate) + "_seed_" + str(s) + "_k_" + str(ksize)
            try:
                os.system("python ../../compute_presence_absence_metrics.py " + ground_truth_ko_file + " " + predicted_ko_file + " " + output_file + " --toolname sourmash")
            except:
                print("Error encountered for size: " + str(error_rate) + " seed: " + str(s) + " ksize: " + str(ksize))
            progress += 1
            print("Progress: " + str(progress) + " / " + str(len(error_rates) * len(seed) * 2))

            if progress % 255 == 0:
                os.system("wait")