import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------
# Data loading & preprocessing
# ----------------------------
def load_and_process_files(result_file, patient_info_file, apply_filter=False):
    result_df = pd.read_csv(result_file)
    patient_info_df = pd.read_csv(patient_info_file)

    # Extract patient ID from tstfname
    result_df['patient_id'] = result_df['tstfname'].str.extract(r'(\d+)\.dcm')

    # Keep needed columns
    filtered_df = result_df[['patient_id', 'Prostate_dice', 'Prostate_hausdorff95', 
                             'Rectum_dice', 'Rectum_hausdorff95', 
                             'Bladder_dice', 'Bladder_hausdorff95',
                             'Femur_head_l_dice', 'Femur_head_l_hausdorff95',
                             'Femur_head_r_dice', 'Femur_head_r_hausdorff95']].copy()

    # Minimal quality gate on prostate DSC
    filtered_df = filtered_df[filtered_df['Prostate_dice'] >= 0.1].copy()

    # Merge with patient info
    filtered_df = filtered_df.rename(columns={'patient_id': 'PatientID'})
    filtered_df['PatientID'] = filtered_df['PatientID'].astype(str)
    patient_info_df['PatientID'] = patient_info_df['PatientID'].astype(str)

    merged_df = pd.merge(filtered_df, patient_info_df, on='PatientID', how='inner').copy()

    # --- CHANGED: convert HD95 from cm to mm ---
    hd_cols = [
        'Prostate_hausdorff95', 'Rectum_hausdorff95', 'Bladder_hausdorff95',
        'Femur_head_l_hausdorff95', 'Femur_head_r_hausdorff95'
    ]
    for c in hd_cols:
        if c in merged_df.columns:
            merged_df[c] = merged_df[c].astype(float) * 10.0  # cm -> mm
    # -------------------------------------------

    # Optional 2σ/2σ filter per organ (applied to this file's rows only)
    if apply_filter:
        for organ in ['Prostate', 'Rectum', 'Bladder', 'Femur_head_l', 'Femur_head_r']:
            dice_col = f'{organ}_dice'
            hd_col   = f'{organ}_hausdorff95'
            if dice_col not in merged_df.columns or hd_col not in merged_df.columns:
                continue
            dice_mean = merged_df[dice_col].mean()
            dice_std  = merged_df[dice_col].std()
            hd_mean   = merged_df[hd_col].mean()
            hd_std    = merged_df[hd_col].std()
            merged_df = merged_df[
                (merged_df[dice_col] >= (dice_mean - 2 * dice_std)) &
                (merged_df[hd_col]   <= (hd_mean + 2 * hd_std))
            ].copy()

    # Parse date and drop invalid
    merged_df['ScanDate'] = pd.to_datetime(merged_df['ScanDate'], format='%Y%m%d', errors='coerce')
    merged_df = merged_df.dropna(subset=['ScanDate']).copy()
    return merged_df


# ----------------------------
# File pairs (edit as needed)
# ----------------------------
file_pairs = [
    (
        r'J:\xxx',
        r'J:\xxx',
        True
    ),
    (
        r'J:\xxx',
        r'J:\sxxx',
        True
    )
]

merged_dfs = [load_and_process_files(result, patient_info, apply_filter) 
              for result, patient_info, apply_filter in file_pairs]

combined_df = pd.concat(merged_dfs, ignore_index=True)

# Filter to prostate DSC >= 0.5
filtered_combined_df = combined_df[combined_df['Prostate_dice'] >= 0.5].copy()
filtered_combined_df['YearMonth'] = filtered_combined_df['ScanDate'].dt.to_period('M')

# Keep months with >3 cases
grouped_df = filtered_combined_df.groupby('YearMonth').filter(lambda x: len(x) > 3).copy()

# Limit to dates up to 2023-08-30
plotting_df = grouped_df[grouped_df['ScanDate'] <= pd.Timestamp('2023-08-30')].copy()

# Ordered month labels like ['2022-01', ...]
year_months = plotting_df['YearMonth'].dt.strftime('%Y-%m').unique()
year_months = np.sort(year_months)

# Organs & columns
organs = {
    'Prostate':      {'dice': 'Prostate_dice',      'hausdorff': 'Prostate_hausdorff95'},
    'Rectum':        {'dice': 'Rectum_dice',        'hausdorff': 'Rectum_hausdorff95'},
    'Bladder':       {'dice': 'Bladder_dice',       'hausdorff': 'Bladder_hausdorff95'},
    'Femur_head_l':  {'dice': 'Femur_head_l_dice',  'hausdorff': 'Femur_head_l_hausdorff95'},
    'Femur_head_r':  {'dice': 'Femur_head_r_dice',  'hausdorff': 'Femur_head_r_hausdorff95'}
}

# ----------------------------
# Control limits from first 7 months of Manteia (2022–2023) only
# ----------------------------
def compute_monthly_control_limits(df_manteia, organs):
    df = df_manteia.copy()
    df['YearMonth'] = df['ScanDate'].dt.to_period('M')
    baseline_months = sorted(df['YearMonth'].unique())[:7]
    ctrl = {}
    base_df = df[df['YearMonth'].isin(baseline_months)].copy()
    month_group = base_df.groupby('YearMonth')

    for organ, cols in organs.items():
        if cols['dice'] not in base_df.columns or cols['hausdorff'] not in base_df.columns:
            continue
        monthly_dice_means = month_group[cols['dice']].mean().values
        monthly_hd_means   = month_group[cols['hausdorff']].mean().values
        ctrl[organ] = {
            'dice_mean': float(np.mean(monthly_dice_means)),
            'dice_std':  float(np.std(monthly_dice_means)),
            'hd_mean':   float(np.mean(monthly_hd_means)),
            'hd_std':    float(np.std(monthly_hd_means))
        }
    return ctrl

manteia_df = merged_dfs[0].copy()
control_limits = compute_monthly_control_limits(manteia_df, organs)


# ----------------------------
# Helper: draw baseline divider + arrow + label
# ----------------------------
def add_baseline_marker(ax, year_months, y_frac=0.5, arrow_len_frac=0.10, label="Baseline"):
    """
    Draws a vertical SOLID line between 2022-07 and 2022-08 and a left-pointing arrow
    that STARTS at the line with a fixed axes-fraction length (arrow_len_frac).
    Places 'Baseline' above the arrow.

    ax:          matplotlib Axes
    year_months: array of '%Y-%m' strings in x order
    y_frac:      vertical position (axes fraction) for arrow
    arrow_len_frac: horizontal arrow length (axes fraction)
    label:       text to place above the arrow
    """
    ym_list = list(year_months)
    idx_jul = ym_list.index('2022-07') if '2022-07' in ym_list else None
    idx_aug = ym_list.index('2022-08') if '2022-08' in ym_list else None

    vline_x = None
    if idx_jul is not None:
        # boxes are centered at 1..N; boundary after July is (idx_jul + 1) + 0.5
        vline_x = (idx_jul + 1) + 0.5
    elif idx_aug is not None:
        # boundary before Aug is (idx_aug + 1) - 0.5
        vline_x = (idx_aug + 1) - 0.5
    else:
        # If neither month appears, do nothing
        return

    # Draw solid vertical line
    ax.axvline(vline_x, color='k', linestyle='--', linewidth=1.1)

    # Convert to axes fraction so arrow has constant length visually
    N = len(year_months)
    x_min, x_max = 0.5, N + 0.5
    start_frac = (vline_x - x_min) / (x_max - x_min)
    end_frac = max(0.0, start_frac - arrow_len_frac)

    # Arrow: starts at the vertical line, points left
    ax.annotate(
        "", xy=(end_frac, y_frac), xytext=(start_frac, y_frac),
        xycoords=ax.transAxes, textcoords=ax.transAxes,
        arrowprops=dict(arrowstyle="->", lw=1.5)
    )

    # Label above arrow, centered
    text_x = (start_frac + end_frac) / 2.0
    ax.text(text_x-0.01, y_frac + 0.05, label, transform=ax.transAxes,
            ha='center', va='bottom', fontsize=10)


# ----------------------------
# Plotting
# ----------------------------
for organ, columns in organs.items():
    if columns['dice'] not in plotting_df.columns:
        continue

    # Build monthly data lists in the same order as year_months
    dice_data = [
        plotting_df[plotting_df['YearMonth'].dt.strftime('%Y-%m') == ym][columns['dice']].dropna()
        for ym in year_months
    ]
    hd_data = [
        plotting_df[plotting_df['YearMonth'].dt.strftime('%Y-%m') == ym][columns['hausdorff']].dropna()
        for ym in year_months
    ]

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    meanprops = {"marker": "x", "markerfacecolor": "black", "markeredgecolor": "black",
                 "markersize": 5, "zorder": 3}

    # --- DSC (upper) ---
    box1 = ax1.boxplot(dice_data, patch_artist=True, showmeans=True, meanprops=meanprops)
    # Color first 13 months orange (model A), rest lightblue (model B)
    for box in box1['boxes'][:13]:
        box.set_facecolor('orange')
    for i in range(13, len(box1['boxes'])):
        box1['boxes'][i].set_facecolor('lightblue')

    ax1.set_title(f'{organ} DSC by Month')
    ax1.set_ylabel('DSC')
    ax1.set_ylim([0.7, 1])

    # Control limits
    ctrl = control_limits.get(organ)
    if ctrl:
        x_span = np.arange(0.5, len(dice_data) + 1.5)
        ax1.fill_between(
            x_span,
            ctrl['dice_mean'] - 2 * ctrl['dice_std'],
            ctrl['dice_mean'] + 2 * ctrl['dice_std'],
            color='lightgreen', alpha=0.5, label='±2σ region', zorder=2
        )
        ax1.axhline(ctrl['dice_mean'], color='green', linestyle='--', label='Mean')
        ax1.axhline(ctrl['dice_mean'] + 3 * ctrl['dice_std'], color='red', linestyle=':', label='±3σ')
        ax1.axhline(ctrl['dice_mean'] - 3 * ctrl['dice_std'], color='red', linestyle=':')

    # Add "(a)" to upper-left corner (upper subplot ONLY)
    ax1.text(-0.12, 0.99, '(a)', transform=ax1.transAxes,
             ha='left', va='top', fontsize=12, fontweight='bold')

    # Add baseline vertical divider + arrow + label on upper subplot
    add_baseline_marker(ax1, year_months, y_frac=0.8, arrow_len_frac=0.10, label="Baseline")

    # --- HD95 (lower) ---
    box2 = ax2.boxplot(hd_data, patch_artist=True, showmeans=True, meanprops=meanprops)
    for box in box2['boxes'][:13]:
        box.set_facecolor('orange')
    for i in range(13, len(box2['boxes'])):
        box2['boxes'][i].set_facecolor('lightblue')

    ax2.set_title(f'{organ} HD95 by Month')
    ax2.set_ylabel('HD95 (mm)')      # CHANGED: label to mm
    ax2.set_ylim([0, 20])            # CHANGED: was 0–2 cm; now 0–20 mm

    if ctrl:
        x_span = np.arange(0.5, len(hd_data) + 1.5)
        ax2.fill_between(
            x_span,
            ctrl['hd_mean'] - 2 * ctrl['hd_std'],
            ctrl['hd_mean'] + 2 * ctrl['hd_std'],
            color='lightgreen', alpha=0.5, label=r'±2$\sigma_{\bar{x}}$', zorder=2
        )
        ax2.axhline(ctrl['hd_mean'], color='green', linestyle='--', label='Mean')
        ax2.axhline(ctrl['hd_mean'] + 3 * ctrl['hd_std'], color='red', linestyle=':', label=r'±3$\sigma_{\bar{x}}$')
        ax2.axhline(ctrl['hd_mean'] - 3 * ctrl['hd_std'], color='red', linestyle=':')

    # Baseline marker on lower subplot too
    add_baseline_marker(ax2, year_months, y_frac=0.8, arrow_len_frac=0.10, label="Baseline")

    ax2.legend()

    plt.xticks(np.arange(1, len(year_months) + 1), year_months, rotation=45)
    plt.tight_layout()
    plt.show()
