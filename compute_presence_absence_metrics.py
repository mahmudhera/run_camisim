"""
In this script, we will compute the presence-absence metrics for a generated KO list.
Inputs:
    - KO ground truth
    - KO predictions by some method (DIAMOND, MMSeqs2, sourmash, etc.)
    - Output file name
Outputs:
    - Precision
    - Recall

Workflow:
    1. Read in the KO ground truth and KO predictions
    2. Filter KO predictions so that the filtered KOs make up top 90%, 95%, or 100% of predictions
    3. Compute precision and recall for each filtered set of KOs
    4. Write precision and recall to output file (three lines, one for each filter level, 100% first, 95% second, 90% last)
"""

import pandas as pd
import numpy as np
import argparse
from scipy.stats import kendalltau

def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Compute precision and recall for KO predictions')
    parser.add_argument('ground_truth', type=str, help='Path to KO ground truth')
    parser.add_argument('predictions', type=str, help='Path to KO predictions')
    parser.add_argument('output', type=str, help='Path to output file')
    parser.add_argument('--toolname', type=str, help='Name of tool used to generate predictions', default='diamond')
    args = parser.parse_args()

    # Read in the KO ground truth, headers: ko_id,abund_by_num_reads,abund_by_num_nts,abund_by_mean_cov,abund_by_med_cov
    ground_truth = pd.read_csv(args.ground_truth, sep=',')
    ground_truth_kos_set = set(ground_truth['ko_id'])
    ground_truth_kos_list = ground_truth['ko_id'].tolist()
    ground_truth_abundances = ground_truth['abund_by_num_reads'].tolist()

    # Sort the KO ground truth by abundance (highest to lowest)
    sorted_indices = np.argsort(ground_truth_abundances)[::-1]
    sorted_kos_ground_truth = [ground_truth_kos_list[i] for i in sorted_indices]

    # Read in the KO predictions: ko_id,num_reads,num_nucleotides_covered,relative_abundance_by_num_reads,relative_abundance_by_nucleotides_covered
    predictions = pd.read_csv(args.predictions, sep=',')
    predicted_kos = predictions['ko_id'].tolist()
    if args.toolname == 'sourmash':
        colname = 'abundance'
    if args.toolname == 'diamond':
        colname = 'relative_abundance_by_num_reads'
    relative_abundances = predictions[colname].tolist()

    # Sort the KO predictions by relative abundance (highest to lowest)
    sorted_indices = np.argsort(relative_abundances)[::-1]
    sorted_kos_predictions = [predicted_kos[i] for i in sorted_indices]
    predicted_kos_set = set(sorted_kos_predictions)

    # write headers to output file
    with open(args.output, 'w') as f:
        f.write('purity,completeness,kendalltau\n')

    # Compute precision and recall
    true_positives = len(ground_truth_kos_set.intersection(predicted_kos_set))
    false_positives = len(predicted_kos_set) - true_positives
    false_negatives = len(ground_truth_kos_set) - true_positives
    precision = true_positives / (true_positives + false_positives)
    recall = true_positives / (true_positives + false_negatives)

    # get the common kos between ground truth and predictions
    common_kos = ground_truth_kos_set.intersection(predicted_kos_set)

    # get the sorted order of these commons kos in the ground truth
    common_kos_sorted_indices = [sorted_kos_ground_truth.index(ko) for ko in common_kos]
    common_kos_sorted_gt = [sorted_kos_ground_truth[i] for i in common_kos_sorted_indices]

    # get the sorted order of these commons kos in the predictions
    common_kos_sorted_indices = [sorted_kos_predictions.index(ko) for ko in common_kos]
    common_kos_sorted_pred = [sorted_kos_predictions[i] for i in common_kos_sorted_indices]

    # compute kendall tau using sorted_kos_ground_truth and most_abundant_kos_list
    tau, p_value = kendalltau(common_kos_sorted_gt, common_kos_sorted_pred)

    # write as completeness and purity
    completeness = recall
    purity = precision

    # Write precision, recall, and tau to output file
    with open(args.output, 'a') as f:
        f.write(f'{purity},{completeness},{tau}\n')

if __name__ == '__main__':
    main()