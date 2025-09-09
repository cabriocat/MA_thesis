#!/usr/bin/env python3
"""
Signal-to-Noise Ratio (SNR) Analysis
=====================================

This script performs comprehensive SNR analysis including cumulative SNR
computation, statistical testing, and data export for R analysis.

Tasks performed:
- Compute cumulative SNR across repetitions
- Perform paired t-tests and TOST equivalence testing
- Apply Holm-Bonferroni correction
- Generate SNR and F-value plots
- Export results for R/ggplot analysis
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
import mne

# Suppress MNE info messages
mne.set_log_level('WARNING')


def split_by_rep(epo, max_rep=6):
    """
    Split epochs by repetition number.
    
    Parameters
    ----------
    epo : mne.Epochs
        Epochs object
    max_rep : int
        Maximum repetition number to consider
        
    Returns
    -------
    dict
        Dictionary mapping repetition number to Epochs object
    """
    out = {}
    
    # Use metadata if available
    if epo.metadata is not None and 'rep' in epo.metadata.columns:
        for rep in range(1, max_rep + 1):
            try:
                sel = epo[f'rep == {rep}']
            except Exception:
                sel = epo.copy()[[]]  # Empty epochs
            if len(sel) > 0:
                out[rep] = sel
        return out
    
    # Fall back to event_id name pattern
    if epo.event_id is None:
        return out
    
    # Map event_id -> name for quick lookup
    inv = {v: k for k, v in epo.event_id.items()}
    
    # For each rep, pick events whose name contains "/rep/"
    for rep in range(1, max_rep + 1):
        wanted_codes = [code for code, name in inv.items() if f"/{rep}/" in name]
        if not wanted_codes:
            continue
        
        # Select epochs whose event code is in wanted_codes
        picks = np.isin(epo.events[:, 2], wanted_codes)
        if picks.any():
            out[rep] = epo[picks]
    
    return out


def compute_cumulative_snr(epochs_dict, subjects, roi, baseline, signal_win, max_rep=6):
    """
    Compute cumulative SNR across repetitions for all subjects.
    
    Parameters
    ----------
    epochs_dict : dict
        Dictionary of subject epochs
    subjects : list
        List of subject IDs
    roi : list
        List of ROI channel names
    baseline : tuple
        Baseline time window (start, end) in seconds
    signal_win : tuple
        Signal time window (start, end) in seconds
    max_rep : int
        Maximum number of repetitions
        
    Returns
    -------
    pd.DataFrame
        DataFrame with SNR results
    """
    
    print("Computing cumulative SNR...")
    
    # Pre-compute masks (shared timing across all subjects)
    _times = next(iter(epochs_dict.values())).times
    bl_mask = (_times >= baseline[0]) & (_times <= baseline[1])
    sg_mask = (_times >= signal_win[0]) & (_times <= signal_win[1])
    
    records = []
    
    for subj in subjects:
        if subj not in epochs_dict:
            continue
            
        epo = epochs_dict[subj].copy().pick_channels(roi, ordered=True)
        by_rep = split_by_rep(epo, max_rep=max_rep)
        
        cumul = None
        n_added = 0  # Track actual cumulative count
        
        for rep in range(1, max_rep + 1):
            rep_epo = by_rep.get(rep, None)
            if rep_epo is None or len(rep_epo) == 0:
                continue
            
            # Accumulate epochs
            cumul = rep_epo if cumul is None else mne.concatenate_epochs([cumul, rep_epo])
            n_added += 1
            
            # Compute evoked response
            evoked = cumul.copy().apply_baseline(baseline).average()
            
            # Convert to µV
            data_uv = evoked.data * 1e6
            bl_data = data_uv[:, bl_mask]
            sg_data = data_uv[:, sg_mask]
            
            # Per-channel SNR, then average across ROI channels
            rms_noise = np.sqrt(np.mean(bl_data**2, axis=1))
            mean_sig = np.abs(sg_data).mean(axis=1)
            snr = np.nanmean(mean_sig / rms_noise)
            
            records.append({
                'subject': subj,
                'cum_reps': n_added,
                'snr': snr
            })
    
    df = pd.DataFrame(records)
    print(f"Computed SNR for {len(df['subject'].unique())} subjects")
    
    return df


def holm_decisions(pvals_dict, alpha=0.05):
    """
    Apply Holm-Bonferroni step-down correction.
    
    Parameters
    ----------
    pvals_dict : dict
        Dictionary mapping test names to p-values
    alpha : float
        Significance level
        
    Returns
    -------
    dict
        Dictionary mapping test names to rejection decisions
    """
    items = sorted(pvals_dict.items(), key=lambda kv: kv[1])  # Ascending by p-value
    m = len(items)
    decisions = {lvl: False for lvl, _ in items}
    
    # Step-down: stop at first failure
    for i, (lvl, p) in enumerate(items, start=1):
        thresh = alpha / (m - i + 1)
        if p <= thresh:
            decisions[lvl] = True
        else:
            break
    
    return decisions


def paired_t_one_sided(x, y, alternative):
    """
    Perform one-sided paired t-test.
    
    Parameters
    ----------
    x, y : array-like
        Paired observations
    alternative : str
        'greater' or 'less'
        
    Returns
    -------
    float
        One-sided p-value
    """
    d = np.asarray(x) - np.asarray(y)
    d = d[~np.isnan(d)]
    n = d.size
    
    if n < 2:
        return np.nan
    
    mean_d = d.mean()
    sd_d = d.std(ddof=1)
    
    if sd_d == 0:
        if ((alternative == 'greater' and mean_d > 0) or 
            (alternative == 'less' and mean_d < 0)):
            return 0.0
        else:
            return 1.0
    
    t = mean_d / (sd_d / np.sqrt(n))
    dfree = n - 1
    
    if alternative == 'greater':
        return stats.t.sf(t, dfree)  # P(T > t)
    else:
        return stats.t.cdf(t, dfree)  # P(T < t)


def perform_statistical_tests(df, equivalence_pct=0.05, alpha=0.05):
    """
    Perform statistical tests on SNR data.
    
    Parameters
    ----------
    df : pd.DataFrame
        SNR dataframe
    equivalence_pct : float
        Equivalence margin as percentage
    alpha : float
        Significance level
        
    Returns
    -------
    tuple
        (summary_df, tests_df)
    """
    
    print("Performing statistical tests...")
    
    # Summarize across subjects
    summary = (
        df.groupby('cum_reps', as_index=False)['snr']
        .agg(mean='mean', std='std', n='count')
    )
    summary['sem'] = summary['std'] / np.sqrt(summary['n'])
    summary['delta_abs'] = summary['mean'].diff()
    summary['delta_pct'] = 100 * summary['mean'].pct_change()
    
    # Prepare wide table for paired tests
    wide = df.pivot_table(index='subject', columns='cum_reps', values='snr', aggfunc='mean')
    available_reps = sorted([c for c in wide.columns if isinstance(c, (int, np.integer))])
    
    # Determine equivalence margin based on N=3 group mean
    if 3 in summary['cum_reps'].values:
        snr3_mean = float(summary.loc[summary['cum_reps'] == 3, 'mean'])
        eq_margin_plateau = equivalence_pct * snr3_mean
    else:
        snr3_mean = np.nan
        eq_margin_plateau = np.nan
    
    # Paired t-tests and TOST
    ttest_pvals = {}
    rows = []
    
    for i, N in enumerate(available_reps):
        if i == 0:  # Skip first repetition
            rows.append({
                'cum_reps': N,
                'paired_t_p': np.nan,
                'holm_sig': np.nan,
                'tost_p_lower': np.nan,
                'tost_p_upper': np.nan,
                'tost_equivalent': np.nan,
                'compared_to': np.nan,
                'f_value': np.nan
            })
            continue
        
        N_prev = available_reps[i - 1]  # Previous cumulative repetition count
        
        # Use only subjects with both N and N_prev
        submask = wide[[N_prev, N]].dropna()
        if submask.shape[0] < 2:
            rows.append({
                'cum_reps': N,
                'paired_t_p': np.nan,
                'holm_sig': np.nan,
                'tost_p_lower': np.nan,
                'tost_p_upper': np.nan,
                'tost_equivalent': np.nan,
                'compared_to': N_prev,
                'f_value': np.nan
            })
            continue
        
        x = submask[N].values      # Current N
        y = submask[N_prev].values # Previous N
        
        # Two-sided paired t-test
        t_stat, p_two_sided = stats.ttest_rel(x, y, nan_policy='omit')
        ttest_pvals[N] = p_two_sided
        f_value = t_stat**2  # F-value is t-squared for paired t-test
        
        # TOST equivalence testing
        prev_mean = submask[N_prev].mean()
        eq_margin_adaptive = equivalence_pct * prev_mean
        
        if np.isnan(eq_margin_adaptive) or eq_margin_adaptive == 0:
            p_lower = np.nan
            p_upper = np.nan
            equivalent = np.nan
        else:
            diff = x - y
            p_lower = paired_t_one_sided(
                diff, 
                np.zeros_like(diff) - eq_margin_adaptive, 
                alternative='greater'
            )
            p_upper = paired_t_one_sided(
                diff, 
                np.zeros_like(diff) + eq_margin_adaptive, 
                alternative='less'
            )
            equivalent = (p_lower < alpha) and (p_upper < alpha)
        
        rows.append({
            'cum_reps': N,
            'paired_t_p': p_two_sided,
            'tost_p_lower': p_lower,
            'tost_p_upper': p_upper,
            'tost_equivalent': bool(equivalent) if equivalent == equivalent else np.nan,
            'compared_to': N_prev,
            'f_value': f_value
        })
    
    # Holm-Bonferroni correction
    holm = holm_decisions({k: v for k, v in ttest_pvals.items() if np.isfinite(v)}, alpha=alpha)
    for r in rows:
        N = r['cum_reps']
        r['holm_sig'] = holm.get(N, np.nan)
    
    tests = pd.DataFrame(rows).sort_values('cum_reps').reset_index(drop=True)
    
    # Merge tests into summary
    summary = summary.merge(tests, how='left', on='cum_reps')
    
    # Plateau rule
    if 3 in summary['cum_reps'].values and np.isfinite(snr3_mean):
        thresh = equivalence_pct * snr3_mean
        summary['plateau_rule'] = summary['delta_abs'].abs() < thresh
    else:
        summary['plateau_rule'] = np.nan
    
    return summary, tests


def create_snr_plots(summary):
    """
    Create SNR and F-value plots.
    
    Parameters
    ----------
    summary : pd.DataFrame
        Summary statistics dataframe
        
    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    
    # Create subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
    
    # Top plot: SNR
    ax1.errorbar(
        summary['cum_reps'], 
        summary['mean'], 
        yerr=summary['sem'], 
        marker='o', 
        color='blue'
    )
    
    for x, y in zip(summary['cum_reps'], summary['mean']):
        ax1.text(x, y, f'{y:.2f}', ha='center', va='bottom', fontsize=8)
    
    ax1.set_xlabel('Cumulative repetitions (actual count)')
    ax1.set_ylabel('Mean SNR (µV/µV)')
    ax1.set_title('Cumulative SNR across subjects (ROI)')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.grid(True)
    
    # Bottom plot: F-values
    f_vals = summary['f_value'].dropna()
    cum_reps_f = summary.loc[summary['f_value'].notna(), 'cum_reps']
    
    if len(f_vals) > 0:
        ax2.bar(cum_reps_f, f_vals, alpha=0.7, color='red', width=0.6)
        for x, y in zip(cum_reps_f, f_vals):
            ax2.text(x, y, f'{y:.2f}', ha='center', va='bottom', fontsize=8)
    
    ax2.set_xlabel('Cumulative repetitions (actual count)')
    ax2.set_ylabel('F-value (t²)')
    ax2.set_title('F-values from paired t-tests (next-neighbor comparisons)')
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    return fig


def export_results(summary, project_root):
    """
    Export analysis results for R/ggplot.
    
    Parameters
    ----------
    summary : pd.DataFrame
        Summary statistics
    project_root : str
        Project root directory
    """
    
    print("Exporting results...")
    
    # Create output directory
    export_path = os.path.join(project_root, "data", "3_analysis")
    os.makedirs(export_path, exist_ok=True)
    
    # Export SNR data
    snr_data = summary[['cum_reps', 'mean', 'sem', 'delta_abs', 'delta_pct']].copy()
    snr_data.rename(columns={
        'cum_reps': 'cumulative_repetitions',
        'mean': 'mean_snr',
        'sem': 'sem_snr',
        'delta_abs': 'delta_absolute',
        'delta_pct': 'delta_percent'
    }, inplace=True)
    
    snr_export_file = os.path.join(export_path, "cumulative_snr_data.csv")
    snr_data.to_csv(snr_export_file, index=False)
    print(f"SNR data exported to: {snr_export_file}")
    
    # Export F-value data
    f_data = summary[summary['f_value'].notna()][
        ['cum_reps', 'f_value', 'paired_t_p', 'compared_to']
    ].copy()
    f_data.rename(columns={
        'cum_reps': 'cumulative_repetitions',
        'f_value': 'f_statistic',
        'paired_t_p': 'p_value',
        'compared_to': 'comparison_baseline'
    }, inplace=True)
    
    f_export_file = os.path.join(export_path, "f_values_data.csv")
    f_data.to_csv(f_export_file, index=False)
    print(f"F-value data exported to: {f_export_file}")
    
    # Export complete summary
    complete_export_file = os.path.join(export_path, "complete_snr_analysis.csv")
    summary.to_csv(complete_export_file, index=False)
    print(f"Complete analysis exported to: {complete_export_file}")


def main():
    """Main execution function."""
    
    print("=" * 60)
    print("SNR ANALYSIS")
    print("=" * 60)
    
    # Setup parameters (inline to avoid import issues)
    project_root = os.path.abspath(os.path.join(os.getcwd(), ".."))
    data_path = os.path.join(project_root, "data", "1_raw")
    
    subjects = [
        "sub-bravoc001", "sub-bravoc002", "sub-bravoc003", "sub-bravoc004",
        "sub-bravoc006", "sub-bravoc007", "sub-bravoc008", "sub-bravoc009",
        "sub-bravoc011", "sub-bravoc012", "sub-bravoc013", "sub-bravoc014",
        "sub-bravoc016", "sub-bravoc017", "sub-bravoc018", "sub-bravoc019",
        "sub-bravoc020", "sub-bravoc022", "sub-bravoc025", "sub-bravoc026",
        "sub-bravoc027", "sub-bravoc028", "sub-bravoc029", "sub-bravoc030",
        "sub-bravoc032", "sub-bravoc033", "sub-bravoc035", "sub-bravoc036",
        "sub-bravoc037", "sub-bravoc040", "sub-bravoc041", "sub-bravoc042",
        "sub-bravoc043", "sub-bravoc044", "sub-bravoc045", "sub-bravoc046",
        "sub-bravoc049", "sub-bravoc050", "sub-bravoc051", "sub-bravoc052"
    ]
    
    baseline = (-0.100, 0.0)
    signal_win = (0.285, 0.345)
    
    roi = [
        "FC1", "FCz", "FC2", "FCC1h", "FCC2h",
        "C1", "Cz", "C2", "CCP1h", "CCP2h", 
        "CP1", "CPz", "CP2", "CPP1h", "CPP2h"
    ]
    
    print(f"Analyzing {len(subjects)} subjects")
    print(f"ROI: {len(roi)} channels")
    print(f"Signal window: {signal_win[0]*1000:.0f}-{signal_win[1]*1000:.0f} ms")
    
    # Load epochs data
    epochs_dict = {}
    included_subjects = []
    
    for subject in subjects:
        try:
            epochs_file = os.path.join(
                data_path,
                f"{subject}/preprocessed-epo/{subject}_task-nouns_eeg-epo_100_70_cutoff.fif"
            )
            
            if os.path.exists(epochs_file):
                epo = mne.read_epochs(epochs_file)
                epochs_dict[subject] = epo
                included_subjects.append(subject)
            else:
                print(f"Warning: File not found for {subject}")
                
        except Exception as e:
            print(f"Error loading {subject}: {str(e)}")
            continue
    
    print(f"Loaded {len(included_subjects)} subjects")
    
    # Compute cumulative SNR
    df = compute_cumulative_snr(
        epochs_dict=epochs_dict,
        subjects=included_subjects,
        roi=roi,
        baseline=baseline,
        signal_win=signal_win
    )
    
    # Perform statistical tests
    summary, tests = perform_statistical_tests(df)
    
    # Create plots
    fig = create_snr_plots(summary)
    
    # Save plot
    plots_dir = os.path.join(project_root, "plots", "snr")
    os.makedirs(plots_dir, exist_ok=True)
    plot_path = os.path.join(plots_dir, "snr_analysis.pdf")
    fig.savefig(plot_path, format='pdf', dpi=300, bbox_inches='tight')
    print(f"SNR plot saved: {plot_path}")
    
    plt.show()
    
    # Export results
    export_results(summary, project_root)
    
    # Display summary table
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    
    cols = ['cum_reps', 'mean', 'sem', 'delta_abs', 'delta_pct',
            'paired_t_p', 'holm_sig', 'tost_equivalent', 'plateau_rule']
    
    display_cols = [col for col in cols if col in summary.columns]
    print(summary[display_cols].to_string(index=False, float_format='%.4f'))
    
    print("\n" + "=" * 60)
    print("SNR ANALYSIS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
