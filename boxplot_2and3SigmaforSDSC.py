import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def load_and_process_files(result_file, patient_info_file, apply_filter=False):
    # Load the CSV files
    result_df = pd.read_csv(result_file)
    patient_info_df = pd.read_csv(patient_info_file)
    
    # Extract patient ID from 'tstfname'
    result_df['patient_id'] = result_df['tstfname'].str.extract(r'(\d+)\.dcm')

    # Select relevant columns and filter out rows with Prostate_SurfaceDice1 < 0.1
    filtered_df = result_df[['patient_id', 'Prostate_SurfaceDice1', 'Prostate_SurfaceDice2', 
                             'Rectum_SurfaceDice1', 'Rectum_SurfaceDice2', 
                             'Bladder_SurfaceDice1', 'Bladder_SurfaceDice2',
                             'Femur_head_l_SurfaceDice1', 'Femur_head_l_SurfaceDice2',
                             'Femur_head_r_SurfaceDice1', 'Femur_head_r_SurfaceDice2']]
    
    filtered_df = filtered_df[filtered_df['Prostate_SurfaceDice1'] >= 0.1]

    # Rename the patient_id column in filtered_df to match patient_info_df
    filtered_df = filtered_df.rename(columns={'patient_id': 'PatientID'})

    # Convert PatientID columns to string
    filtered_df['PatientID'] = filtered_df['PatientID'].astype(str)
    patient_info_df['PatientID'] = patient_info_df['PatientID'].astype(str)

    # Merge the dataframes
    merged_df = pd.merge(filtered_df, patient_info_df, on='PatientID')

    if apply_filter:
        # Calculate mean and std for each organ's Surface Dice and filter low outliers
        for organ in ['Prostate', 'Rectum', 'Bladder', 'Femur_head_l', 'Femur_head_r']:
            surface_dice1_mean = merged_df[f'{organ}_SurfaceDice1'].mean()
            surface_dice1_std = merged_df[f'{organ}_SurfaceDice1'].std()
            surface_dice2_mean = merged_df[f'{organ}_SurfaceDice2'].mean()
            surface_dice2_std = merged_df[f'{organ}_SurfaceDice2'].std()

            merged_df = merged_df[
                (merged_df[f'{organ}_SurfaceDice1'] >= (surface_dice1_mean - 2 * surface_dice1_std)) & 
                (merged_df[f'{organ}_SurfaceDice2'] >= (surface_dice2_mean - 2 * surface_dice2_std))
            ]

    # Convert ScanDate to proper datetime format
    merged_df['ScanDate'] = pd.to_datetime(merged_df['ScanDate'], format='%Y%m%d')
    
    return merged_df

# File paths for the pairs of result and patient info files
file_pairs = [
    (
        r'J:\xxxx',
        r'J:\xxxx',
        True
    ),
    (
        r'J:\xxxx', 
        r'J://xxxx', 
        True
    )
]

# Load and process each pair of files
merged_dfs = [load_and_process_files(result, patient_info, apply_filter) 
              for result, patient_info, apply_filter in file_pairs]

# Combine all the processed dataframes into one
combined_df = pd.concat(merged_dfs, ignore_index=True)

# Exclude data points with Prostate_SurfaceDice1 < 0.1
filtered_combined_df = combined_df[combined_df['Prostate_SurfaceDice1'] >= 0.1].copy()

# Extract month and year for grouping
filtered_combined_df['YearMonth'] = filtered_combined_df['ScanDate'].dt.to_period('M')

# Filter out groups with less than 3 data points
grouped_df = filtered_combined_df.groupby('YearMonth').filter(lambda x: len(x) > 3)

# Get the subset of data from save1 folder only (unchanged downstream usage)
save1_df = merged_dfs[1]

# List of organs and their respective Surface Dice column names
organs = {
    'Prostate': {'surface_dice1': 'Prostate_SurfaceDice1', 'surface_dice2': 'Prostate_SurfaceDice2'},
    'Rectum': {'surface_dice1': 'Rectum_SurfaceDice1', 'surface_dice2': 'Rectum_SurfaceDice2'},
    'Bladder': {'surface_dice1': 'Bladder_SurfaceDice1', 'surface_dice2': 'Bladder_SurfaceDice2'},
    'Femur_head_l': {'surface_dice1': 'Femur_head_l_SurfaceDice1', 'surface_dice2': 'Femur_head_l_SurfaceDice2'},
    'Femur_head_r': {'surface_dice1': 'Femur_head_r_SurfaceDice1', 'surface_dice2': 'Femur_head_r_SurfaceDice2'}
}

# ===== NEW CODE TO EXCLUDE DATA AFTER 2023-08-20 FOR PLOTTING =====
plot_df = grouped_df[grouped_df['ScanDate'] <= pd.to_datetime('2023-08-20')].copy()

# Unique YearMonth values for plotting
year_months = plot_df['YearMonth'].dt.strftime('%Y-%m').unique()
year_months = np.sort(year_months)

# =======================
# CONTROL LIMITS: MANTEIA ONLY, FIRST 7 MONTHS BY MONTHLY MEAN
# =======================
# Build control limits from the Manteia cohort (merged_dfs[0]) using the earliest 7 distinct months
manteia_df = merged_dfs[0].copy()
manteia_df = manteia_df[manteia_df['Prostate_SurfaceDice1'] >= 0.1].copy()
manteia_df['YearMonth'] = manteia_df['ScanDate'].dt.to_period('M')

# Earliest 7 calendar months by ScanDate (from patient info)
control_months = np.array(sorted(manteia_df['YearMonth'].unique()))[:7]

# Precompute control means/stds (based on MONTHLY MEANS across the first 7 months)
control_stats = {}
for organ, cols in organs.items():
    monthly_means_1 = []
    monthly_means_2 = []
    for ym in control_months:
        m1 = manteia_df.loc[manteia_df['YearMonth'] == ym, cols['surface_dice1']].dropna()
        m2 = manteia_df.loc[manteia_df['YearMonth'] == ym, cols['surface_dice2']].dropna()
        if len(m1) > 0:
            monthly_means_1.append(m1.mean())
        if len(m2) > 0:
            monthly_means_2.append(m2.mean())
    monthly_means_1 = np.array(monthly_means_1, dtype=float)
    monthly_means_2 = np.array(monthly_means_2, dtype=float)

    ctrl_mean_1 = float(np.mean(monthly_means_1)) if monthly_means_1.size else np.nan
    ctrl_std_1  = float(np.std(monthly_means_1))  if monthly_means_1.size else np.nan
    ctrl_mean_2 = float(np.mean(monthly_means_2)) if monthly_means_2.size else np.nan
    ctrl_std_2  = float(np.std(monthly_means_2))  if monthly_means_2.size else np.nan

    control_stats[organ] = {
        'sd1': {'mean': ctrl_mean_1, 'std': ctrl_std_1},
        'sd2': {'mean': ctrl_mean_2, 'std': ctrl_std_2},
    }

# Prepare plotting for all organs
for organ, columns in organs.items():
    # Prepare data for SurfaceDice1 and SurfaceDice2 for the current organ, excluding NaN values
    surface_dice1_data = [
        plot_df[plot_df['YearMonth'].dt.strftime('%Y-%m') == ym][columns['surface_dice1']].dropna() 
        for ym in year_months
    ]
    surface_dice2_data = [
        plot_df[plot_df['YearMonth'].dt.strftime('%Y-%m') == ym][columns['surface_dice2']].dropna() 
        for ym in year_months
    ]

    # Plotting for the current organ
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)

    meanprops = {"marker": "x", "markerfacecolor": "black", "markeredgecolor": "black", "markersize": 5, "zorder": 3}

    # Boxplot for SurfaceDice1
    box1 = ax1.boxplot(surface_dice1_data, patch_artist=True, showmeans=True, meanprops=meanprops)
    for box in box1['boxes'][:13]:
        box.set_facecolor('orange')
    for i in range(13, len(box1['boxes'])):
        box1['boxes'][i].set_facecolor('lightblue')
    ax1.set_title(f'{organ} Surface DSC 1mm by Month')
    ax1.set_ylabel(f'{organ} SDSC 1')
    ax1.set_ylim([0.2, 1])
    ax1.grid(False)

    # === CONTROL LIMIT DISPLAY (SDSC1) ===
    ctrl1 = control_stats[organ]['sd1']
    if not np.isnan(ctrl1['mean']) and not np.isnan(ctrl1['std']):
        # Mean (green dashed)
        ax1.axhline(y=ctrl1['mean'], color='green', linestyle='--', label='Mean')
        # ±2σ (green shaded)
        ax1.fill_between(
            np.arange(0.5, len(year_months) + 1.5),
            ctrl1['mean'] - 2 * ctrl1['std'],
            ctrl1['mean'] + 2 * ctrl1['std'],
            color='lightgreen', alpha=0.5, label=r'±2$\sigma_{\bar{x}}$', zorder=2
        )
        # ±3σ (red dotted lines)
        ax1.axhline(y=ctrl1['mean'] - 3 * ctrl1['std'], color='red', linestyle=':', label=r'±3$\sigma_{\bar{x}}$')
        ax1.axhline(y=ctrl1['mean'] + 3 * ctrl1['std'], color='red', linestyle=':')

    #ax1.legend(loc='lower left', fontsize=8)
    ax1.legend()

    # Boxplot for SurfaceDice2
    box2 = ax2.boxplot(surface_dice2_data, patch_artist=True, showmeans=True, meanprops=meanprops)
    for box in box2['boxes'][:13]:
        box.set_facecolor('orange')
    for i in range(13, len(box2['boxes'])):
        box2['boxes'][i].set_facecolor('lightblue')
    ax2.set_title(f'{organ} Surface DSC 2mm by Month')
    ax2.set_ylabel(f'{organ} SDSC 2')
    ax2.set_ylim([0.2, 1])
    ax2.grid(False)

    # === CONTROL LIMIT DISPLAY (SDSC2) ===
    ctrl2 = control_stats[organ]['sd2']
    if not np.isnan(ctrl2['mean']) and not np.isnan(ctrl2['std']):
        # Mean (green dashed)
        ax2.axhline(y=ctrl2['mean'], color='green', linestyle='--', label='Mean')
        # ±2σ (green shaded)
        ax2.fill_between(
            np.arange(0.5, len(year_months) + 1.5),
            ctrl2['mean'] - 2 * ctrl2['std'],
            ctrl2['mean'] + 2 * ctrl2['std'],
            color='lightgreen', alpha=0.5, label=r'±2$\sigma_{\bar{x}}$', zorder=2
        )
        # ±3σ (red dotted lines)
        ax2.axhline(y=ctrl2['mean'] - 3 * ctrl2['std'], color='red', linestyle=':', label=r'±3$\sigma_{\bar{x}}$')
        ax2.axhline(y=ctrl2['mean'] + 3 * ctrl2['std'], color='red', linestyle=':')

    #ax2.legend(loc='lower left', fontsize=8)
    ax2.legend()

    # Set x-axis labels for Year-Month and rotate them
    plt.xticks(np.arange(1, len(year_months) + 1), year_months, rotation=45)
    plt.xlabel('Year-Month')

    # Adjust layout
    plt.tight_layout()

    # Show the plot for the current organ
    plt.show()

# -------------------------
# Monthly statistics export (unchanged)
# -------------------------
def extract_monthly_statistics(grouped_df, year_months, organs):
    monthly_statistics = []
    for year_month in year_months:
        month_data = {"YearMonth": year_month}
        for organ, columns in organs.items():
            dice_values = grouped_df[grouped_df['YearMonth'].dt.strftime('%Y-%m') == year_month][columns['surface_dice1']].dropna()
            hausdorff_values = grouped_df[grouped_df['YearMonth'].dt.strftime('%Y-%m') == year_month][columns['surface_dice2']].dropna()

            month_data[f"{organ}_surface_dice1_Mean"] = dice_values.mean()
            month_data[f"{organ}_surface_dice2_Mean"] = hausdorff_values.mean()
        monthly_statistics.append(month_data)

    stats_df = pd.DataFrame(monthly_statistics)
    return stats_df

stats_df = extract_monthly_statistics(grouped_df, year_months, organs)
print(stats_df)
stats_df.to_csv(r"C:\xxxx.csv", index=False)
