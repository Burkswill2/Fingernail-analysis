import numpy as np
import pandas as pd
from scipy import stats

def process_data(file_path, num_bins, output_file):
    """
    Processes fluorescence data to analyze cortical and cytosolic expression in bins across a spatial gradient.

    Parameters:
    - file_path: The path to the CSV file containing the data.
    - num_bins: The number of bins to divide the xCortex range into for analysis.
    - output_file: The path for the output Excel file where the results will be saved.

    The function reads the data, divides it into bins based on xCortex values, calculates statistical measures for
    each bin, and saves the results to an Excel file.
    """
    # Read in the data from the CSV file.
    data = pd.read_csv(file_path)

    # Determine the range of xCortex values to establish bin edges.
    min_value = data['xCortex'].min()
    max_value = data['xCortex'].max()
    bin_range = max_value - min_value
    bin_size = bin_range / num_bins

    # Initialize lists to store results for each bin.
    bins, starts, ends = [], [], []
    cortical_bin_averages, cytosol_bin_averages = [], []
    cortical_bin_sds, cytosol_bin_sds = [], []
    ci_lowers, ci_uppers, mean_differences = [], [], []
    bin_ratios, plus_minus = [], []

    # Iterate over each bin to calculate statistical measures.
    for n in range(num_bins):
        bin_start = min_value + n * bin_size
        bin_end = min_value + (n + 1) * bin_size

        # Select data points that fall within the current bin's range of xCortex values.
        bin_data = data[(data['xCortex'] >= bin_start) & (data['xCortex'] < bin_end)]

        if not bin_data.empty:
            # Calculate statistical measures for the bin.
            cortical_mean, cytosol_mean, cortical_sd, cytosol_sd, ci_lower, ci_upper, mean_difference, bin_ratio = calculate_statistics(bin_data)

            # Store the calculated measures and bin information in lists.
            bins.append(f'Bin {n + 1}')
            starts.append(f'{bin_start:.2f}')
            ends.append(f'{bin_end:.2f}')
            cortical_bin_averages.append(cortical_mean)
            cytosol_bin_averages.append(cytosol_mean)
            cortical_bin_sds.append(cortical_sd)
            cytosol_bin_sds.append(cytosol_sd)
            ci_lowers.append(ci_lower)
            ci_uppers.append(ci_upper)
            mean_differences.append(mean_difference)
            bin_ratios.append(bin_ratio)
            plus_minus.append("+" if bin_ratio > 1 else "-")
        else:
            # Handle empty bins by appending None or an empty string for each measure.
            [lst.append(None) for lst in (cortical_bin_averages, cytosol_bin_averages, bin_ratios)]
            plus_minus.append('')

    # Combine the results into a DataFrame.
    results_df = pd.DataFrame({
        'Bin': bins,
        'Start': starts,
        'End': ends,
        'Cortical Average': cortical_bin_averages,
        'Cytosol Average': cytosol_bin_averages,
        'Cortical Std': cortical_bin_sds,
        'Cytosol Std': cytosol_bin_sds,
        'CI Lower': ci_lowers,
        'CI Upper': ci_uppers,
        'Mean Differences': mean_differences,
        'Bin Ratio': bin_ratios,
        '+/-': plus_minus
    })

    # Write the results DataFrame to an Excel file.
    results_df.to_excel(output_file, index=False)
    print(f'Results have been saved to {output_file}')



def calculate_statistics(data):
    """
    Calculates statistical measures for a given bin of data.

    Parameters:
    - data: A DataFrame containing the bin of data for which to calculate statistics.

    Returns:
    - A tuple containing the calculated measures: means, standard deviations, confidence interval bounds,
      mean difference, and the ratio of cortical to cytosol mean values.
    """
    cortical_mean = data['yCortex'].mean()
    cytosol_mean = data['yCytosol'].mean()
    bin_ratio = cortical_mean / cytosol_mean

    cortical_sd = data['yCortex'].std()
    cytosol_sd = data['yCytosol'].std()

    n_cortical = len(data['yCortex'])
    n_cytosol = len(data['yCytosol'])

    mean_difference = cortical_mean - cytosol_mean
    se_difference = np.sqrt((cortical_sd ** 2 / n_cortical) + (cytosol_sd ** 2 / n_cytosol))

    # Assuming 95% CI, for two-tailed t-test
    confidence_level = 0.95
    degrees_freedom = min(n_cortical - 1, n_cytosol - 1)
    t_critical = stats.t.ppf((1 + confidence_level) / 2, degrees_freedom)

    ci_lower = mean_difference - t_critical * se_difference
    ci_upper = mean_difference + t_critical * se_difference

    return cortical_mean, cytosol_mean, cortical_sd, cytosol_sd, ci_lower, ci_upper, mean_difference, bin_ratio


# Example call to process_data with output file path
file_path = '/Users/willburks/Desktop/data/FT1.csv'  # Replace with the actual file path
output_file = '/Users/willburks/Desktop/data/AnalysisResults.xlsx'  # Specify the desired output file path
process_data(file_path, 10, output_file)
