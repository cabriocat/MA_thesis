#!/usr/bin/env python3
"""
SNR (Signal-to-Noise Ratio) Analysis Script

This script handles SNR analysis:
- Cumulative SNR computation across repetitions
- Statistical testing (paired t-tests, TOST equivalence testing)
- Holm-Bonferroni correction for multiple comparisons
- Data export for R visualization

Author: [Your Name]
"""

import os
import numpy as np
import pandas as pd
import mne
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats

def split_epochs_by_repetition(epo, max_rep=6):
    """
    Split epochs by repetition number.
    
    Parameters:
    -----------
    epo : mne.Epochs
        MNE Epochs object
    max_rep : int
        Maximum number of repetitions
        
    Returns:
    --------
    dict
        Dictionary mapping repetition number to Epochs
    """
    out = {}
    
    # Try metadata approach first
    if epo.metadata is not None and 'rep' in epo.metadata.columns:
        for rep in range(1, max_rep+1):
            try:
                sel = epo['rep == @rep']
            except Exception:
                sel = epo.copy()[[]]  # empty
            if len(sel) > 0:
                out[rep] = sel
        return out

    # Fall back to event_id name pattern
    if epo.event_id is None:
        return out
    
    # Map event_id->name for quick lookup
    inv = {v: k for k, v in epo.event_id.items()}
    
    # For each rep, pick events whose name contains "/rep/"
    for rep in range(1, max_rep+1):
        wanted_codes = [code for code, name in inv.items() if f"/{rep}/" in name]
        if not wanted_codes:
            continue
        picks = np.isin(epo.events[:, 2], wanted_codes)
        if picks.any():
            out[rep] = epo[picks]
    
    return out

def compute_cumulative_snr(epochs_dict, subjects, roi, baseline, signal_win, max_rep=6):
    """
    Compute cumulative SNR for all subjects.
    
    Parameters:
    -----------
    epochs_dict : dict
        Dictionary mapping subject to Epochs
    subjects : list
        List of subject identifiers
    roi : list
        List of ROI channel names
    baseline : tuple
        Baseline time window
    signal_win : tuple
        Signal time window
    max_rep : int
        Maximum number of repetitions
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with SNR data
    """
    
    # Pre-compute masks (shared timing across all subjects)
    _times = next(iter(epochs_dict.values())).times
    bl_mask = (_times >= baseline[0]) & (_times <= baseline[1])
    sg_mask = (_times >= signal_win[0]) & (_times <= signal_win[1])

    records = []
    
    for subj in subjects:
        if subj not in epochs_dict:
            continue
            
        epo = epochs_dict[subj].copy().pick_channels(roi)
        by_rep = split_epochs_by_repetition(epo, max_rep=max_rep)

        cumul = None
        n_added = 0
        
        for rep in range(1, max_rep+1):
            rep_epo = by_rep.get(rep, None)
            if rep_epo is None or len(rep_epo) == 0:
                continue

            cumul = rep_epo if cumul is None else mne.concatenate_epochs([cumul, rep_epo])
            n_added += 1

            evoked = cumul.copy().apply_baseline(baseline).average()

            data_uv = evoked.data * 1e6  # µV
            bl_data = data_uv[:, bl_mask]
            sg_data = data_uv[:, sg_mask]

            # Per-channel SNR, then average across ROI channels
            rms_noise = np.sqrt(np.mean(bl_data**2, axis=1))
            mean_sig = np.abs(sg_data).mean(axis=1)
            snr = np.nanmean(mean_sig / rms_noise)

            records.append({'subject': subj, 'cum_reps': n_added, 'snr': snr})

    return pd.DataFrame(records)

def compute_summary_statistics(df):
    """Compute summary statistics across subjects."""
    summary = (
        df.groupby('cum_reps', as_index=False)['snr']
          .agg(mean='mean', std='std', n='count')
    )
    summary['sem'] = summary['std'] / np.sqrt(summary['n'])
    summary['delta_abs'] = summary['mean'].diff()
    summary['delta_pct'] = 100 * summary['mean'].pct_change()
    
    return summary

def holm_bonferroni_correction(pvals_dict, alpha=0.05):
    """
    Apply Holm-Bonferroni step-down correction.
    
    Parameters:
    -----------
    pvals_dict : dict
        Dictionary mapping level to p-value
    alpha : float
        Significance level
        
    Returns:
    --------
    dict
        Dictionary mapping level to rejection decision
    """
    items = sorted(pvals_dict.items(), key=lambda kv: kv[1])  # ascending by p
    m = len(items)
    decisions = {lvl: False for lvl, _ in items}
    
    for i, (lvl, p) in enumerate(items, start=1):
        thresh = alpha / (m - i + 1)
        if p <= thresh:
            decisions[lvl] = True
        else:
            break
    
    return decisions

def paired_t_one_sided(x, y, alternative):
    """
    One-sided paired t-test.
    
    Parameters:
    -----------
    x, y : array-like
        Paired samples
    alternative : str
        'greater' or 'less'
        
    Returns:
    --------
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
        if (alternative == 'greater' and mean_d > 0) or (alternative == 'less' and mean_d < 0):
            return 0.0
        else:
            return 1.0
    
    t = mean_d / (sd_d / np.sqrt(n))
    dfree = n - 1
    
    if alternative == 'greater':
        return stats.t.sf(t, dfree)     # P(T > t)
    else:
        return stats.t.cdf(t, dfree)    # P(T < t)

def perform_statistical_tests(df, summary, alpha=0.05, equivalence_pct=0.05):
    """
    Perform statistical tests on SNR data.
    
    Parameters:
    -----------
    df : pd.DataFrame
        SNR data per subject
    summary : pd.DataFrame
        Summary statistics
    alpha : float
        Significance level
    equivalence_pct : float
        Equivalence margin as percentage
        
    Returns:
    --------
    pd.DataFrame
        Summary with test results
    """
    
    # Prepare wide table: rows=subject, cols=cum_reps, values=snr
    wide = df.pivot_table(index='subject', columns='cum_reps', values='snr', aggfunc='mean')
    available_reps = sorted([c for c in wide.columns if isinstance(c, (int, np.integer))])

    # Determine equivalence margin based on N=3 group mean
    if 3 in summary['cum_reps'].values:
        snr3_mean = float(summary.loc[summary['cum_reps']==3, 'mean'])
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
                'cum_reps': N, 'paired_t_p': np.nan, 'holm_sig': np.nan,
                'tost_p_lower': np.nan, 'tost_p_upper': np.nan,
                'tost_equivalent': np.nan, 'compared_to': np.nan, 'f_value': np.nan
            })
            continue
        
        N_prev = available_reps[i-1]
        
        # Use only subjects with both N and N_prev
        submask = wide[[N_prev, N]].dropna()
        if submask.shape[0] < 2:
            rows.append({
                'cum_reps': N, 'paired_t_p': np.nan, 'holm_sig': np.nan,
                'tost_p_lower': np.nan, 'tost_p_upper': np.nan,
                'tost_equivalent': np.nan, 'compared_to': N_prev, 'f_value': np.nan
            })
            continue

        x = submask[N].values
        y = submask[N_prev].values

        # Two-sided paired t-test
        t_stat, p_two_sided = stats.ttest_rel(x, y, nan_policy='omit')
        ttest_pvals[N] = p_two_sided
        f_value = t_stat**2

        # TOST equivalence testing
        prev_mean = submask[N_prev].mean()
        eq_margin_adaptive = equivalence_pct * prev_mean
        
        if np.isnan(eq_margin_adaptive) or eq_margin_adaptive == 0:
            p_lower = np.nan
            p_upper = np.nan
            equivalent = np.nan
        else:
            diff = x - y
            p_lower = paired_t_one_sided(diff, np.zeros_like(diff) - eq_margin_adaptive, 
                                       alternative='greater')
            p_upper = paired_t_one_sided(diff, np.zeros_like(diff) + eq_margin_adaptive, 
                                       alternative='less')
            equivalent = (p_lower < alpha) and (p_upper < alpha)

        rows.append({
            'cum_reps': N, 'paired_t_p': p_two_sided, 'tost_p_lower': p_lower,
            'tost_p_upper': p_upper, 'tost_equivalent': bool(equivalent) if equivalent==equivalent else np.nan,
            'compared_to': N_prev, 'f_value': f_value
        })

    # Holm-Bonferroni correction
    holm = holm_bonferroni_correction({k: v for k, v in ttest_pvals.items() if np.isfinite(v)}, 
                                     alpha=alpha)
    for r in rows:
        N = r['cum_reps']
        r['holm_sig'] = holm.get(N, np.nan)

    tests = pd.DataFrame(rows).sort_values('cum_reps').reset_index(drop=True)
    
    # Merge tests into summary
    summary = summary.merge(tests, how='left', on='cum_reps')

    # Optional plateau flag
    if 3 in summary['cum_reps'].values and np.isfinite(snr3_mean):
        thresh = equivalence_pct * snr3_mean
        summary['plateau_rule'] = summary['delta_abs'].abs() < thresh
    else:
        summary['plateau_rule'] = np.nan

    return summary

def create_snr_plots(summary):
    """Create SNR and F-value plots."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

    # Top plot: SNR
    ax1.errorbar(summary['cum_reps'], summary['mean'], yerr=summary['sem'], 
                marker='o', color='blue')
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
    plt.show()
    
    return fig

def export_data_for_r(summary, project_root):
    """Export data for R visualization."""
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
    f_data = summary[summary['f_value'].notna()][['cum_reps', 'f_value', 'paired_t_p', 'compared_to']].copy()
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

def perform_snr_analysis(data_objects, max_rep=6, equivalence_pct=0.05, alpha=0.05):
    """
    Main SNR analysis function.
    
    Parameters:
    -----------
    data_objects : dict
        Dictionary containing data objects from 01_data_setup.py
    max_rep : int
        Maximum number of repetitions
    equivalence_pct : float
        Equivalence margin as percentage
    alpha : float
        Significance level
    """
    
    print("=== SNR Analysis ===")
    
    # Extract data objects
    epochs_dict = data_objects['epochs_dict']
    included_subjects = data_objects['included_subjects']
    roi = data_objects['roi']
    baseline = data_objects['baseline']
    signal_win = data_objects['signal_win']
    project_root = data_objects['project_root']
    
    print(f"Number of subjects: {len(included_subjects)}")
    print(f"ROI channels: {len(roi)}")
    print(f"Baseline: {baseline}")
    print(f"Signal window: {signal_win}")
    print(f"Max repetitions: {max_rep}")
    
    # Compute cumulative SNR
    df = compute_cumulative_snr(epochs_dict, included_subjects, roi, baseline, signal_win, max_rep)
    print(f"SNR data computed: {df.shape}")
    
    # Summary statistics
    summary = compute_summary_statistics(df)
    print("Summary statistics computed")
    
    # Statistical tests
    summary = perform_statistical_tests(df, summary, alpha, equivalence_pct)
    print("Statistical tests completed")
    
    # Create plots
    fig = create_snr_plots(summary)
    
    # Export data for R
    export_data_for_r(summary, project_root)
    
    # Display results
    cols = ['cum_reps','mean','sem','delta_abs','delta_pct',
            'paired_t_p','holm_sig','tost_p_lower','tost_p_upper','tost_equivalent',
            'plateau_rule','compared_to','f_value']
    print("\n=== SNR Analysis Results ===")
    print(summary[cols].to_string(index=False))
    
    return summary, fig

def main():
    """Main execution function - requires data from 01_data_setup.py"""
    print("=== SNR Analysis Script ===")
    print("This script requires data objects from 01_data_setup.py")
    print("Run this script by importing and calling functions with your data objects.")
    print("\nExample usage:")
    print("from data_setup import main as load_data")
    print("data = load_data()")
    print("summary, fig = perform_snr_analysis(data)")

if __name__ == "__main__":
    main()
