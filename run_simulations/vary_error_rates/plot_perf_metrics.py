import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV files
#mean_performance_file = './mean_performance_metrics.csv'
#std_performance_file = './std_performance_metrics.csv'

mean_performance_file = './mean_performance_metrics_batch.csv'
std_performance_file = './std_performance_metrics_batch.csv'

# Reading the CSV files
mean_performance_data = pd.read_csv(mean_performance_file)
std_performance_data = pd.read_csv(std_performance_file)

# Merging the mean and standard deviation data
merged_data = pd.merge(mean_performance_data, std_performance_data, on=['size', 'tool'], suffixes=('_mean', '_std'))

# List of metrics to plot
metrics = ['purity', 'completeness', 'wt_purity', 'wt_completeness', 'pearsonr_common', 'pearsonr_all', 
            'kl_div_common_gt_to_pred', 
           'kl_div_common_pred_to_gt', 'kendalltau', 'bray_curtis']

metric_titles = ['Purity', 'Completeness', 'Weighted Purity', 'Weighted Completeness', 'Pearson Correlation (Common KOs)',
                 'Pearson Correlation (All KOs)', 'KL Divergence (Ground Truth to Predictions)', 
                 'KL Divergence (Predictions to Ground Truth)', 'Kendall Tau', 'Bray-Curtis Distance']

# Number of unique tools in the dataset
num_tools = merged_data['tool'].nunique()

# Setting up a colorblind-friendly palette and markers
palette = sns.color_palette("colorblind", num_tools)  # Colorblind-friendly palette with correct number of colors
markers = ["o", "s", "^", "D", "P", "*", "X", "<", ">"]  # Different markers for each tool

# Setting up the plotting
sns.set(style="whitegrid")
plt.figure(figsize=(12, 20))

# Plotting each metric in a separate subplot
for i, metric in enumerate(metrics, 1):
    # get title of the metric
    metric_title = metric_titles[i-1]

    plt.subplot(5, 2, i)

    # Creating a lineplot with colorblind-friendly colors and markers
    lineplot = sns.lineplot(
        x='size', y=f'{metric}_mean', hue='tool', style='tool',
        data=merged_data, palette=palette, markers=markers[:num_tools], markersize=8)

    # Adding error bars with matching colors
    lines = lineplot.get_lines()
    for j, tool in enumerate(merged_data['tool'].unique()):
        subset = merged_data[merged_data['tool'] == tool]
        plt.errorbar(subset['size'], subset[f'{metric}_mean'], yerr=subset[f'{metric}_std'], 
                     fmt='none', capsize=5, label='_nolegend_', color=lines[j].get_color())

    plt.xlabel('Error rate (%)')
    plt.ylabel(metric_title)
    plt.title(f'Mean {metric_title} with Error Bars')
    plt.xticks(subset['size'].unique(), labels=subset['size'].unique())
    plt.legend(title='Tool')

plt.tight_layout()
plt.savefig('perf_metrics_batch.pdf')
