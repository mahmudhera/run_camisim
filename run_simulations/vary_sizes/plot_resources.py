import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the CSV files
mean_resources_file = './mean_resources_used.csv'
std_resources_file = './std_resources_used.csv'

# Reading the CSV files
mean_resources_data = pd.read_csv(mean_resources_file)
std_resources_data = pd.read_csv(std_resources_file)

# Merging the mean and standard deviation data
merged_resources_data = pd.merge(mean_resources_data, std_resources_data, on=['size', 'tool'], suffixes=('_mean', '_std'))

# Apply logarithmic transformation to walltime and cputime (mean and std)
#for col in ['walltime', 'cputime']:
#    merged_resources_data[f'{col}_mean'] = np.log(merged_resources_data[f'{col}_mean'])
#    merged_resources_data[f'{col}_std'] = np.log(merged_resources_data[f'{col}_std'])

# Conversion factor from kbytes to Gbytes for max_rss
kbytes_to_gbytes = 1 / (1024**2)
merged_resources_data['max_rss_mean'] = merged_resources_data['max_rss_mean'] * kbytes_to_gbytes
merged_resources_data['max_rss_std'] = merged_resources_data['max_rss_std'] * kbytes_to_gbytes

# List of resource usage metrics to plot
resource_metrics = ['walltime', 'cputime', 'max_rss']

# Number of unique tools in the dataset
num_tools = merged_resources_data['tool'].nunique()

# Setting up a colorblind-friendly palette and markers
palette = sns.color_palette("colorblind", num_tools)  # Colorblind-friendly palette
markers = ["o", "s", "^", "D", "P", "*", "X", "<", ">"]  # Different markers for each tool

# Setting up the plotting
sns.set(style="whitegrid")
plt.figure(figsize=(20, 6))

# Plotting each resource usage metric in a separate subplot
for i, metric in enumerate(resource_metrics, 1):
    plt.subplot(1, 3, i)  # Adjusting layout for three columns and one row

    # Creating a lineplot for each metric
    lineplot = sns.lineplot(
        x='size', y=f'{metric}_mean', hue='tool', style='tool',
        data=merged_resources_data, palette=palette, markers=markers[:num_tools], markersize=8)

    # Adding error bars with matching colors
    lines = lineplot.get_lines()
    for j, tool in enumerate(merged_resources_data['tool'].unique()):
        subset = merged_resources_data[merged_resources_data['tool'] == tool]
        plt.errorbar(subset['size'], subset[f'{metric}_mean'], yerr=subset[f'{metric}_std'], 
                     fmt='none', capsize=5, label='_nolegend_', color=lines[j].get_color())

    plt.xlabel('Size of metagenome (Gbp)')
    #y_label = 'Log Seconds' if metric in ['walltime', 'cputime'] else 'Gbytes'
    y_label = 'Seconds' if metric in ['walltime', 'cputime'] else 'Gbytes'
    plt.ylabel(y_label)
    #plt.title(f'Mean Log {metric.title()} with Error Bars' if metric in ['walltime', 'cputime'] else f'Mean {metric.title()} with Error Bars')
    plt.title(f'Average {metric.title()} with Error Bars' if metric in ['walltime', 'cputime'] else f'Mean {metric.title()} with Error Bars')
    plt.xscale('log')
    plt.xticks(subset['size'].unique(), labels=subset['size'].unique())
    plt.legend(title='Tool')

plt.tight_layout()
plt.savefig('resources.pdf')
