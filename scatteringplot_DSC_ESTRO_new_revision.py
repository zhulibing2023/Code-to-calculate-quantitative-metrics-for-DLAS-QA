import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.stats import mannwhitneyu

def load_and_process_files(result_file, patient_info_file, source_label):
    # Load the CSV files
    result_df = pd.read_csv(result_file)
    patient_info_df = pd.read_csv(patient_info_file)
    
    # Extract patient ID from 'tstfname'
    result_df['patient_id'] = result_df['tstfname'].str.extract(r'(\d+)\.dcm')

    # Select relevant columns and filter out rows with Prostate_dice < 0.1
    structures = ['Prostate', 'Bladder', 'Rectum', 'Femur_head_l', 'Femur_head_r']
    structure_cols = [f'{structure}_dice' for structure in structures] + [f'{structure}_hausdorff95' for structure in structures]
    filtered_df = result_df[['patient_id'] + structure_cols]
    
    # Apply filter to remove rows where any of the structure Dice scores are smaller than 0.1
    dice_cols = [f'{structure}_dice' for structure in structures]
    filtered_df = filtered_df[filtered_df[dice_cols].ge(0.1).all(axis=1)]

    # Add source label to differentiate between datasets
    filtered_df['Source'] = source_label

    # Rename the patient_id column in filtered_df to match patient_info_df
    filtered_df = filtered_df.rename(columns={'patient_id': 'PatientID'})

    # Convert PatientID columns to string
    filtered_df['PatientID'] = filtered_df['PatientID'].astype(str)
    patient_info_df['PatientID'] = patient_info_df['PatientID'].astype(str)

    # Merge the dataframes
    merged_df = pd.merge(filtered_df, patient_info_df, on='PatientID', how='inner')

    # Convert ScanDate to proper datetime format
    merged_df['ScanDate'] = pd.to_datetime(merged_df['ScanDate'], format='%Y%m%d', errors='coerce')
    
    return merged_df

# -----------------------
# Load data
# -----------------------
save1_file = r'J:\xxx.csv'
save1_patient_info_file = r'J:\xxx.csv'
save1_df = load_and_process_files(save1_file, save1_patient_info_file, 'save1')

manteia_file = r'J:\xxx.csv'
manteia_patient_info_file = r'J:\xxx.csv'
manteia_df = load_and_process_files(manteia_file, manteia_patient_info_file, 'manteia')

# Structures list for reuse
structures = ['Prostate', 'Bladder', 'Rectum', 'Femur_head_l', 'Femur_head_r']

# -----------------------
# CONTROL LIMITS from model A (earliest 150 scans by ScanDate)
# -----------------------
# Sort by ScanDate coming from patient-info merge, ensure unique PatientID if duplicates exist, then take first 150
manteia_sorted = (
    manteia_df
    .sort_values('ScanDate')
    .drop_duplicates(subset='PatientID', keep='first')
)

control_base = manteia_sorted.head(150).copy()
print(f"Base control cohort: {len(control_base)} "
      f"({control_base['ScanDate'].min().date()} → {control_base['ScanDate'].max().date()})")

structure_stats = {}
for structure in structures:
    dsc_col = f'{structure}_dice'
    hd_col  = f'{structure}_hausdorff95'

    # Exclude DSC<0.1 and HD95==0 for this structure's stats
    cohort_s = control_base[
        (control_base[dsc_col] >= 0.1) &
        (control_base[hd_col] < 5)
    ].dropna(subset=[dsc_col, hd_col])

    dice_mean = cohort_s[dsc_col].mean()
    dice_std  = cohort_s[dsc_col].std()
    hd_mean   = cohort_s[hd_col].mean()
    hd_std    = cohort_s[hd_col].std()

    structure_stats[structure] = {
        'dice_mean': dice_mean,
        'dice_std': dice_std,
        'hausdorff_mean': hd_mean,
        'hausdorff_std': hd_std,
        'dice_lower_bound': dice_mean - 2 * dice_std,
        'hausdorff_upper_bound': hd_mean + 2 * hd_std,
        'n_used': len(cohort_s)
    }
    
    print(f"{structure}: used {len(cohort_s)} of 150 for limits")


print(f"Control cohort size used: {len(cohort_s)} "
      f"({cohort_s['ScanDate'].min().date()} → {cohort_s['ScanDate'].max().date()})")



# -----------------------
# Combine & flag outliers using manteia-derived bounds
# -----------------------
combined_df = pd.concat([save1_df, manteia_df], ignore_index=True)

# Mark rows that are outliers/normal based on the manteia-derived bounds
for df in (manteia_df, save1_df, combined_df):
    for structure in structures:
        df[f'{structure}_Outlier'] = (
            (df[f'{structure}_dice'] < structure_stats[structure]['dice_lower_bound']) |
            (df[f'{structure}_hausdorff95'] > structure_stats[structure]['hausdorff_upper_bound'])
        )
        df[f'{structure}_Normal'] = (
            (df[f'{structure}_dice'] >= structure_stats[structure]['dice_lower_bound']) &
            (df[f'{structure}_hausdorff95'] <= structure_stats[structure]['hausdorff_upper_bound'])
        )
        print(structure_stats[structure])

# -----------------------
# Plotting helpers (unchanged from your structure, except they now rely on manteia-based bounds)
# -----------------------
# Get the min and max dates for shading the entire x-axis range
min_date = combined_df['ScanDate'].min()
max_date = combined_df['ScanDate'].max()
save1_part1 = save1_df[((save1_df['ScanDate'] >= '2023-02-01') & (save1_df['ScanDate'] <= '2023-08-05')) ]
max_date = save1_part1['ScanDate'].max()  # NOTE: this line narrows max_date; keep or remove as you prefer

def plot_structure_separate(structure_name):
    dice_lower_bound = structure_stats[structure_name]['dice_lower_bound']
    hausdorff_upper_bound = structure_stats[structure_name]['hausdorff_upper_bound']
    dice_mean = structure_stats[structure_name]['dice_mean']
    hausdorff_mean = structure_stats[structure_name]['hausdorff_mean']

    fig, axs = plt.subplots(2, 1, sharex=True)
    
    # Separate data
    manteia_data = manteia_df[manteia_df['Source'] == 'manteia']
    save1_part1 = save1_df[((save1_df['ScanDate'] >= '2023-02-01') & (save1_df['ScanDate'] <= '2023-08-05')) | 
                            ((save1_df['ScanDate'] >= '2024-05-02') & (save1_df['ScanDate'] <= '2024-11-01'))]
    save1_part2 = save1_df[(save1_df['ScanDate'] >= '2023-08-06') & (save1_df['ScanDate'] <= '2024-05-01')]
    
    # Exclude data after 2023-08-20 for plotting (kept as in your script)
    cutoff_date = pd.to_datetime("2023-08-20")
    manteia_data = manteia_data[manteia_data['ScanDate'] <= cutoff_date]
    save1_part1 = save1_part1[save1_part1['ScanDate'] <= cutoff_date]
    save1_part2 = save1_part2[save1_part2['ScanDate'] <= cutoff_date]
    
    # Plot Dice score
    axs[0].fill_between([min_date, max_date], dice_lower_bound, 1, color='lightgreen', alpha=0.3, label=f'{structure_name} Dice Bound ±2std')
    axs[0].axhline(y=dice_mean, color='green', linestyle='--', label=f'{structure_name} Mean Dice')
    
    axs[0].scatter(manteia_data[manteia_data[f'{structure_name}_Normal']]['ScanDate'], 
                   manteia_data[manteia_data[f'{structure_name}_Normal']][f'{structure_name}_dice'], 
                   c='orange', marker='^', s=20, label='Model A Normal')
    axs[0].scatter(manteia_data[manteia_data[f'{structure_name}_Outlier']]['ScanDate'], 
                   manteia_data[manteia_data[f'{structure_name}_Outlier']][f'{structure_name}_dice'], 
                   color='orange', marker='x', s=20, label='Model A Outlier')
    
    axs[0].scatter(save1_part1[save1_part1[f'{structure_name}_Normal']]['ScanDate'], 
                   save1_part1[save1_part1[f'{structure_name}_Normal']][f'{structure_name}_dice'], 
                   c='blue', marker='o', s=20, label='Model B Normal')
    axs[0].scatter(save1_part1[save1_part1[f'{structure_name}_Outlier']]['ScanDate'], 
                   save1_part1[save1_part1[f'{structure_name}_Outlier']][f'{structure_name}_dice'], 
                   color='blue', marker='x', s=20, label='Model B Outlier')
    
    axs[0].scatter(save1_part2[save1_part2[f'{structure_name}_Normal']]['ScanDate'], 
                   save1_part2[save1_part2[f'{structure_name}_Normal']][f'{structure_name}_dice'], 
                   c='orange', marker='^', s=20, label='Model A Normal')  # <- kept as in your code
    axs[0].scatter(save1_part2[save1_part2[f'{structure_name}_Outlier']]['ScanDate'], 
                   save1_part2[save1_part2[f'{structure_name}_Outlier']][f'{structure_name}_dice'], 
                   color='orange', marker='x', s=20, label='Model A Outlier')  # <- kept as in your code
    
    axs[0].set_ylim([0.5, 1])
    axs[0].set_ylabel('DSC', color='blue')
    axs[0].grid(False)
    axs[0].set_title(f'{structure_name} Monitoring')

    # Plot Hausdorff distance
    axs[1].fill_between([min_date, max_date], 0, hausdorff_upper_bound, color='lightgreen', alpha=0.3, label=f'{structure_name} Hausdorff Bound ±2std')
    axs[1].axhline(y=hausdorff_mean, color='green', linestyle='--', label=f'{structure_name} Mean Hausdorff')
    
    axs[1].scatter(manteia_data[manteia_data[f'{structure_name}_Normal']]['ScanDate'], 
                   manteia_data[manteia_data[f'{structure_name}_Normal']][f'{structure_name}_hausdorff95'], 
                   c='orange', marker='^', s=20, label='Model A Normal')
    axs[1].scatter(manteia_data[manteia_data[f'{structure_name}_Outlier']]['ScanDate'], 
                   manteia_data[manteia_data[f'{structure_name}_Outlier']][f'{structure_name}_hausdorff95'], 
                   color='orange', marker='x', s=20, label='Model A Outlier')
    
    axs[1].scatter(save1_part1[save1_part1[f'{structure_name}_Normal']]['ScanDate'], 
                   save1_part1[save1_part1[f'{structure_name}_Normal']][f'{structure_name}_hausdorff95'], 
                   c='blue', marker='o', s=20, label='Model B Normal')
    axs[1].scatter(save1_part1[save1_part1[f'{structure_name}_Outlier']]['ScanDate'], 
                   save1_part1[save1_part1[f'{structure_name}_Outlier']][f'{structure_name}_hausdorff95'], 
                   color='blue', marker='x', s=20, label='Model B Outlier')
    
    axs[1].scatter(save1_part2[save1_part2[f'{structure_name}_Normal']]['ScanDate'], 
                   save1_part2[save1_part2[f'{structure_name}_Normal']][f'{structure_name}_hausdorff95'], 
                   c='orange', marker='^', s=20, label='Model A Normal')  # <- kept as in your code
    axs[1].scatter(save1_part2[save1_part2[f'{structure_name}_Outlier']]['ScanDate'], 
                   save1_part2[save1_part2[f'{structure_name}_Outlier']][f'{structure_name}_hausdorff95'], 
                   color='orange', marker='x', s=20, label='Model A Outlier')  # <- kept as in your code
    
    axs[1].set_ylim([0, 4])
    axs[1].set_ylabel('HD 95 (cm)', color='purple')
    axs[1].grid(False)
    # Add legend with selected items
    handles, labels = axs[1].get_legend_handles_labels()
    selected_indices = [2, 3, 4, 5]  # Indices of the important legend items
    axs[1].legend([handles[i] for i in selected_indices], [labels[i] for i in selected_indices], loc='best')
    
    # Formatting x-axis
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    axs[1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    fig.autofmt_xdate()
    
    plt.tight_layout()
    plt.show()



    # -----------------------
# Outlier counts summary
# -----------------------

def outlier_counts_by_structure(df, structures):
    """Return dict: {structure: count_of_outliers} where outlier = (DSC < LB) OR (HD95 > UB)."""
    return {s: int(df[f'{s}_Outlier'].sum()) for s in structures}

# 150 selected control cohort (the exact 150 earliest manteia scans used for limits)
control_ids = set(control_base['PatientID'].astype(str))
manteia_control_150 = manteia_df[manteia_df['PatientID'].astype(str).isin(control_ids)]

counts_control_150 = outlier_counts_by_structure(manteia_control_150, structures)
counts_manteia_all = outlier_counts_by_structure(manteia_df, structures)
counts_save1_all   = outlier_counts_by_structure(save1_df, structures)

summary_counts = pd.DataFrame({
    '150_selected': [counts_control_150[s] for s in structures],
    'manteia_all':  [counts_manteia_all[s] for s in structures],
    'save1_all':    [counts_save1_all[s] for s in structures],
}, index=structures)

print("\nOutlier counts per organ (either DSC fails OR HD95 fails):")
print(summary_counts)


# Generate plots for each structure
for structure in structures:
    plot_structure_separate(structure)

# -----------------------
# Wilcoxon rank-sum (Mann–Whitney U) tests (unchanged)
# -----------------------
wilcoxon_results = []

# Model A = manteia_df, Model B = save1_df
for structure in structures:
    # DSC comparison
    dsc_model_a = manteia_df[f'{structure}_dice'].dropna()
    dsc_model_i = save1_df[f'{structure}_dice'].dropna()
    dsc_stat, dsc_p = mannwhitneyu(dsc_model_a, dsc_model_i, alternative='two-sided')
    
    # HD95 comparison
    hd_model_a = manteia_df[f'{structure}_hausdorff95'].dropna()
    hd_model_i = save1_df[f'{structure}_hausdorff95'].dropna()
    hd_stat, hd_p = mannwhitneyu(hd_model_a, hd_model_i, alternative='two-sided')
    
    wilcoxon_results.append({
        'Structure': structure,
        'DSC_stat': dsc_stat,
        'DSC_p_value': dsc_p,
        'HD95_stat': hd_stat,
        'HD95_p_value': hd_p
    })

# Convert to DataFrame for readability
wilcoxon_df = pd.DataFrame(wilcoxon_results)
print("\nWilcoxon Rank-Sum Test Results (Model A vs Model B):")
print(wilcoxon_df)
