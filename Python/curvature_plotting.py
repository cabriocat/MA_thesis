#!/usr/bin/env python3
"""
ERP Curvature Analysis
======================

This script computes and visualizes the curvature of ERP waveforms across conditions.
Curvature is calculated as the second derivative of the signal, showing how rapidly
the ERP is changing direction at each time point.

This is an exploratory analysis to identify time points where neural processing
shows rapid transitions or inflection points.

Tasks performed:
- Load raw epochs data for all subjects
- Compute grand-average ERPs per condition
- Calculate curvature (second derivative) for each condition
- Visualize curvature across time for all conditions
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import mne
from mne import grand_average

# Suppress MNE info messages
mne.set_log_level('WARNING')


# ============================================================================
# CONFIGURATION
# ============================================================================

def setup_parameters():
    """Define experimental parameters and paths."""
    
    # Project setup
    project_root = os.path.abspath(os.getcwd())
    
    # Data path
    data_path = "/Users/johannberger/Documents/Archive/MA_thesis/data/1_raw"
    
    # Validate path
    if not os.path.exists(data_path):
        raise FileNotFoundError(
            f"Data path not found: {data_path}\n"
            f"Please update the data_path variable."
        )
    
    # Subjects
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
    
    # Baseline correction window
    baseline = (-0.100, 0.0)
    
    # Condition codes
    conditions = ['a', 'f', 't', 'c', 'e', 's']
    
    # Label mapping
    label_mapping = {
        'a': 'Animal',
        'f': 'Food',
        't': 'Tool',
        'c': 'Communication',
        'e': 'Emotion',
        's': 'Social'
    }
    
    # ROI channels (frontocentral region)
    roi = [
        "FC1", "FCz", "FC2", "FCC1h", "FCC2h",
        "C1", "Cz", "C2", "CCP1h", "CCP2h",
        "CP1", "CPz", "CP2", "CPP1h", "CPP2h"
    ]
    
    # Curvature computation parameter
    # sample_distance controls the window size for averaging before computing derivatives
    # sample_distance = 1: no averaging, compute on raw samples (finest resolution)
    # sample_distance = 10: average samples 1-10, 11-20, etc., then compute curvature
    #                       (~10ms windows for 1000Hz sampling, smoother result)
    sample_distance = 1
    
    return {
        'project_root': project_root,
        'data_path': data_path,
        'subjects': subjects,
        'baseline': baseline,
        'conditions': conditions,
        'label_mapping': label_mapping,
        'roi': roi,
        'sample_distance': sample_distance
    }


# ============================================================================
# DATA LOADING
# ============================================================================

def load_epochs_data(data_path, subjects, baseline):
    """
    Load epochs data for all subjects and apply baseline correction.
    
    Parameters
    ----------
    data_path : str
        Path to raw data directory
    subjects : list
        List of subject identifiers
    baseline : tuple
        Baseline time window
        
    Returns
    -------
    epochs_dict : dict
        Dictionary mapping subjects to Epochs objects
    included_subjects : list
        Successfully loaded subjects
    """
    
    epochs_dict = {}
    included_subjects = []
    
    print("Loading epochs data...")
    
    for subject in subjects:
        try:
            epochs_file = os.path.join(
                data_path,
                f"{subject}/preprocessed-epo/{subject}_task-nouns_eeg-epo_100_70_cutoff.fif"
            )
            
            if not os.path.exists(epochs_file):
                print(f"Warning: File not found for {subject}")
                continue
            
            # Load and apply baseline
            epo = mne.read_epochs(epochs_file, preload=True)
            epo.apply_baseline(baseline)
            
            epochs_dict[subject] = epo
            included_subjects.append(subject)
            print(f"✓ Loaded: {subject}")
            
        except Exception as e:
            print(f"✗ Error loading {subject}: {str(e)}")
            continue
    
    print(f"\nSuccessfully loaded {len(included_subjects)} subjects")
    return epochs_dict, included_subjects


def build_condition_evoked(epochs_dict, conditions, included_subjects):
    """
    Build grand-average evoked responses for each condition.
    
    Parameters
    ----------
    epochs_dict : dict
        Dictionary of subject epochs
    conditions : list
        List of condition codes
    included_subjects : list
        List of included subjects
        
    Returns
    -------
    grand_avg_dict : dict
        Dictionary mapping conditions to grand-average evoked responses
    """
    
    # Collect evoked responses per condition
    evoked_by_condition = {cond: [] for cond in conditions}
    
    print("\nBuilding condition-specific evoked responses...")
    
    for subject in included_subjects:
        epo = epochs_dict[subject]
        all_event_names = list(epo.event_id.keys())
        
        for cond in conditions:
            # Get all trials for this condition (across all repetitions)
            condition_events = [
                name for name in all_event_names 
                if f"/{cond}/n" in name
            ]
            
            if not condition_events:
                continue
            
            # Average across all trials for this subject and condition
            evoked = epo[condition_events].average()
            evoked_by_condition[cond].append(evoked)
    
    # Create grand averages
    grand_avg_dict = {}
    for cond in conditions:
        if evoked_by_condition[cond]:
            grand_avg_dict[cond] = grand_average(evoked_by_condition[cond])
            n_subjects = len(evoked_by_condition[cond])
            print(f"  Condition {cond}: {n_subjects} subjects included")
    
    return grand_avg_dict


# ============================================================================
# CURVATURE COMPUTATION
# ============================================================================

def compute_curvature(signal, sampling_rate, sample_distance=1):
    """
    Compute curvature as the second derivative of the signal.
    
    The curvature measures how quickly the signal changes direction.
    A positive curvature indicates convex bending (upward curve),
    while negative curvature indicates concave bending (downward curve).
    
    Parameters
    ----------
    signal : np.ndarray
        1D time-series signal
    sampling_rate : float
        Sampling rate in Hz
    sample_distance : int, optional
        Window size for averaging. If sample_distance=1, uses neighboring samples
        (finest resolution). If sample_distance=10, averages samples in windows of
        10 (e.g., 1-10, 11-20, 21-30, etc.) before computing derivatives.
        Default is 1.
        
    Returns
    -------
    curvature : np.ndarray
        Curvature values (second derivative), averaged over windows
    """
    
    if sample_distance == 1:
        # No averaging needed, compute directly
        dt = 1.0 / sampling_rate
        first_derivative = np.gradient(signal, dt)
        second_derivative = np.gradient(first_derivative, dt)
        return second_derivative
    
    else:
        # Average signal in windows, then compute curvature
        n_samples = len(signal)
        n_windows = n_samples // sample_distance
        
        # Truncate to fit complete windows
        signal_truncated = signal[:n_windows * sample_distance]
        
        # Reshape into windows and average
        signal_windowed = signal_truncated.reshape(n_windows, sample_distance)
        signal_averaged = signal_windowed.mean(axis=1)
        
        # Compute derivatives on averaged signal
        # Effective time step is the window size
        dt = sample_distance / sampling_rate
        first_derivative = np.gradient(signal_averaged, dt)
        second_derivative = np.gradient(first_derivative, dt)
        
        return second_derivative


def extract_roi_curvature(grand_avg_dict, roi, conditions, sample_distance=1):
    """
    Extract ROI-averaged curvature for each condition.
    
    Parameters
    ----------
    grand_avg_dict : dict
        Grand-average evoked responses per condition
    roi : list
        ROI channel names
    conditions : list
        Condition codes
    sample_distance : int, optional
        Window size for averaging curvature. Default is 1 (no averaging).
        
    Returns
    -------
    curvature_dict : dict
        Dictionary mapping conditions to curvature arrays
    times : np.ndarray
        Time points (downsampled to window centers)
    """
    
    print(f"\nComputing curvature for each condition (sample_distance={sample_distance})...")
    
    # Get first condition to extract time info
    first_cond = conditions[0]
    times_full = grand_avg_dict[first_cond].times
    sampling_rate = grand_avg_dict[first_cond].info['sfreq']
    
    # Compute times for window centers
    if sample_distance == 1:
        times = times_full
    else:
        n_samples = len(times_full)
        n_windows = n_samples // sample_distance
        # Use center time of each window
        window_indices = np.arange(n_windows) * sample_distance + sample_distance // 2
        times = times_full[window_indices]
    
    curvature_dict = {}
    
    for cond in conditions:
        evoked = grand_avg_dict[cond]
        
        # Pick ROI channels
        picks = mne.pick_channels(evoked.ch_names, include=roi)
        
        # Average across ROI channels to get single waveform
        roi_signal = evoked.data[picks, :].mean(axis=0)
        
        # Convert to µV for interpretability
        roi_signal_uv = roi_signal * 1e6
        
        # Compute curvature with specified window size
        curvature = compute_curvature(roi_signal_uv, sampling_rate, sample_distance)
        
        curvature_dict[cond] = curvature
        print(f"  ✓ Computed curvature for condition: {cond}")
    
    return curvature_dict, times


# ============================================================================
# VISUALIZATION
# ============================================================================

def plot_curvature(curvature_dict, times, label_mapping, sample_distance=1, save_path=None):
    """
    Create curvature plot for all conditions.
    
    Parameters
    ----------
    curvature_dict : dict
        Curvature values per condition
    times : np.ndarray
        Time points in seconds
    label_mapping : dict
        Condition labels
    sample_distance : int, optional
        Sample distance used in computation (for title). Default is 1.
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        Created figure
    """
    
    # Okabe-Ito colorblind-friendly palette
    colors = {
        'a': '#388E3C',  # Green - Animal
        'f': '#1976D2',  # Blue - Food
        't': '#D32F2F',  # Red - Tool
        'c': '#d9a00f',  # Orange - Communication
        'e': '#616161',  # Gray - Emotion
        's': '#ff00ff'   # Magenta - Social
    }
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot curvature for each condition
    for cond in curvature_dict.keys():
        ax.plot(
            times * 1000,  # Convert to milliseconds
            curvature_dict[cond],
            label=label_mapping[cond],
            color=colors[cond],
            linewidth=2,
            alpha=0.85
        )
    
    # Add zero reference line
    ax.axhline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # Add stimulus onset line
    ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.7)
    
    # Styling
    ax.set_xlabel('Time (ms)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Curvature (µV/s²)', fontsize=13, fontweight='bold')
    
    # Add sample distance info to title
    if sample_distance == 1:
        title_suffix = '(no averaging)'
    else:
        title_suffix = f'(averaged over {sample_distance}-sample windows)'
    
    ax.set_title(f'ERP Curvature Analysis by Semantic Category\n(ROI: Frontocentral) {title_suffix}', 
                 fontsize=14, fontweight='bold', pad=15)
    
    # Legend
    ax.legend(
        loc='best',
        frameon=True,
        framealpha=0.95,
        fontsize=11,
        title='Condition'
    )
    
    # Grid
    ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.5)
    
    # Tight layout
    plt.tight_layout()
    
    # Save if path provided
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Saved plot to: {save_path}")
    
    return fig


def plot_erp_and_curvature(grand_avg_dict, curvature_dict, times, roi, 
                           label_mapping, sample_distance=1, save_path=None):
    """
    Create two-panel plot showing both ERP waveforms and curvature.
    
    Parameters
    ----------
    grand_avg_dict : dict
        Grand-average evoked responses
    curvature_dict : dict
        Curvature values per condition
    times : np.ndarray
        Time points
    roi : list
        ROI channel names
    label_mapping : dict
        Condition labels
    sample_distance : int, optional
        Sample distance used in computation. Default is 1.
    save_path : str, optional
        Save path
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        Created figure
    """
    
    # Colors
    colors = {
        'a': '#388E3C',
        'f': '#1976D2',
        't': '#D32F2F',
        'c': '#d9a00f',
        'e': '#616161',
        's': '#ff00ff'
    }
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    # Get full time axis for ERP plot
    first_cond = list(grand_avg_dict.keys())[0]
    times_full = grand_avg_dict[first_cond].times
    
    # Panel 1: ERP waveforms
    for cond in grand_avg_dict.keys():
        evoked = grand_avg_dict[cond]
        picks = mne.pick_channels(evoked.ch_names, include=roi)
        roi_signal = evoked.data[picks, :].mean(axis=0) * 1e6  # µV
        
        ax1.plot(
            times_full * 1000,
            roi_signal,
            label=label_mapping[cond],
            color=colors[cond],
            linewidth=2,
            alpha=0.85
        )
    
    ax1.axhline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax1.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.7)
    ax1.set_ylabel('Amplitude (µV)', fontsize=13, fontweight='bold')
    ax1.set_title('ERP Waveforms by Semantic Category', 
                  fontsize=14, fontweight='bold', pad=15)
    ax1.legend(loc='best', frameon=True, framealpha=0.95, fontsize=11)
    ax1.grid(True, alpha=0.25, linestyle=':', linewidth=0.5)
    
    # Panel 2: Curvature
    for cond in curvature_dict.keys():
        ax2.plot(
            times * 1000,
            curvature_dict[cond],
            label=label_mapping[cond],
            color=colors[cond],
            linewidth=2,
            alpha=0.85
        )
    
    ax2.axhline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax2.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.7)
    ax2.set_xlabel('Time (ms)', fontsize=13, fontweight='bold')
    ax2.set_ylabel('Curvature (µV/s²)', fontsize=13, fontweight='bold')
    
    # Add sample distance info to title
    if sample_distance == 1:
        title_suffix = '(no averaging)'
    else:
        title_suffix = f'(averaged over {sample_distance}-sample windows)'
    
    ax2.set_title(f'ERP Curvature (Second Derivative) {title_suffix}', 
                  fontsize=14, fontweight='bold', pad=15)
    ax2.grid(True, alpha=0.25, linestyle=':', linewidth=0.5)
    
    plt.tight_layout()
    
        # Save if path provided
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved combined plot to: {save_path}")
    
    return fig


def plot_curvature_heatmap(curvature_dict, times, label_mapping, 
                           sample_distance=1, save_path=None):
    """
    Create heatmap visualization of curvature across conditions and time.
    
    This shows the intensity/magnitude of curvature changes as a color-coded
    heatmap, making it easier to identify temporal patterns and differences
    between conditions.
    
    Parameters
    ----------
    curvature_dict : dict
        Curvature values per condition
    times : np.ndarray
        Time points in seconds
    label_mapping : dict
        Condition labels
    sample_distance : int, optional
        Sample distance used. Default is 1.
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        Created figure
    """
    
    # Prepare data matrix (conditions x time)
    conditions = ['a', 'f', 't', 'c', 'e', 's']
    condition_labels = [label_mapping[c] for c in conditions]
    
    # Stack curvature data
    curvature_matrix = np.array([curvature_dict[c] for c in conditions])
    
    # Create figure
    fig = plt.figure(figsize=(14, 8))
    gs = GridSpec(2, 1, height_ratios=[1, 0.3], hspace=0.3)
    
    # Main heatmap
    ax_heat = fig.add_subplot(gs[0])
    
    # Use diverging colormap (blue-white-red) centered at zero
    vmax = np.abs(curvature_matrix).max()
    vmin = -vmax
    
    im = ax_heat.imshow(
        curvature_matrix,
        aspect='auto',
        cmap='RdBu_r',  # Red for positive, Blue for negative
        vmin=vmin,
        vmax=vmax,
        extent=[times[0]*1000, times[-1]*1000, len(conditions)-0.5, -0.5],
        interpolation='bilinear'
    )
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax_heat, orientation='vertical', pad=0.02)
    cbar.set_label('Curvature (µV/s²)', fontsize=12, fontweight='bold')
    
    # Set y-axis labels
    ax_heat.set_yticks(range(len(conditions)))
    ax_heat.set_yticklabels(condition_labels, fontsize=11)
    ax_heat.set_ylabel('Semantic Category', fontsize=12, fontweight='bold')
    
    # Add stimulus onset line
    ax_heat.axvline(0, color='black', linestyle='-', linewidth=2, alpha=0.8)
    
    # Title
    if sample_distance == 1:
        title_suffix = '(no averaging)'
    else:
        title_suffix = f'(averaged over {sample_distance}-sample windows)'
    
    ax_heat.set_title(
        f'ERP Curvature Heatmap Across Semantic Categories\n{title_suffix}',
        fontsize=14, fontweight='bold', pad=15
    )
    
    # Bottom panel: Mean absolute curvature over time
    ax_mean = fig.add_subplot(gs[1])
    
    # Okabe-Ito colors
    colors = {
        'a': '#388E3C',
        'f': '#1976D2',
        't': '#D32F2F',
        'c': '#d9a00f',
        'e': '#616161',
        's': '#ff00ff'
    }
    
    for cond in conditions:
        abs_curv = np.abs(curvature_dict[cond])
        ax_mean.plot(
            times * 1000,
            abs_curv,
            label=label_mapping[cond],
            color=colors[cond],
            linewidth=2,
            alpha=0.8
        )
    
    ax_mean.axvline(0, color='black', linestyle='-', linewidth=1.5, alpha=0.7)
    ax_mean.set_xlabel('Time (ms)', fontsize=12, fontweight='bold')
    ax_mean.set_ylabel('Abs. Curvature', fontsize=11, fontweight='bold')
    ax_mean.set_title('Absolute Curvature Magnitude', fontsize=12, fontweight='bold')
    ax_mean.legend(loc='upper right', ncol=6, fontsize=9, framealpha=0.9)
    ax_mean.grid(True, alpha=0.25, linestyle=':', linewidth=0.5)
    
    # Save if path provided
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved heatmap to: {save_path}")
    
    return fig


def plot_curvature_statistics(curvature_dict, times, label_mapping,
                              sample_distance=1, save_path=None):
    """
    Create statistical summary plots of curvature distributions.
    
    Shows time windows with highest curvature variability and
    overall curvature distributions per condition.
    
    Parameters
    ----------
    curvature_dict : dict
        Curvature values per condition
    times : np.ndarray
        Time points in seconds
    label_mapping : dict
        Condition labels
    sample_distance : int, optional
        Sample distance used. Default is 1.
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        Created figure
    """
    
    conditions = ['a', 'f', 't', 'c', 'e', 's']
    
    # Okabe-Ito colors
    colors = {
        'a': '#388E3C',
        'f': '#1976D2',
        't': '#D32F2F',
        'c': '#d9a00f',
        'e': '#616161',
        's': '#ff00ff'
    }
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel 1: Standard deviation of curvature over time (variability)
    ax1 = axes[0, 0]
    curvature_matrix = np.array([curvature_dict[c] for c in conditions])
    std_over_conditions = np.std(curvature_matrix, axis=0)
    
    ax1.fill_between(times * 1000, 0, std_over_conditions, alpha=0.3, color='gray')
    ax1.plot(times * 1000, std_over_conditions, linewidth=2, color='black')
    ax1.axvline(0, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    ax1.set_xlabel('Time (ms)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('SD of Curvature', fontsize=11, fontweight='bold')
    ax1.set_title('Curvature Variability Across Conditions', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.25)
    
    # Panel 2: Distribution of curvature values per condition (violin plot)
    ax2 = axes[0, 1]
    violin_data = [curvature_dict[c] for c in conditions]
    violin_labels = [label_mapping[c] for c in conditions]
    
    parts = ax2.violinplot(violin_data, positions=range(len(conditions)),
                           widths=0.7, showmeans=True, showmedians=True)
    
    # Color the violins
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[conditions[i]])
        pc.set_alpha(0.7)
    
    ax2.set_xticks(range(len(conditions)))
    ax2.set_xticklabels(violin_labels, rotation=45, ha='right')
    ax2.set_ylabel('Curvature (µV/s²)', fontsize=11, fontweight='bold')
    ax2.set_title('Curvature Distribution per Category', fontsize=12, fontweight='bold')
    ax2.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax2.grid(True, alpha=0.25, axis='y')
    
    # Panel 3: Peak positive curvature timing
    ax3 = axes[1, 0]
    
    for cond in conditions:
        # Find peaks in positive curvature
        pos_curv = np.copy(curvature_dict[cond])
        pos_curv[pos_curv < 0] = 0
        
        ax3.plot(times * 1000, pos_curv, label=label_mapping[cond],
                color=colors[cond], linewidth=2, alpha=0.8)
    
    ax3.axvline(0, color='black', linestyle='-', linewidth=1.5, alpha=0.7)
    ax3.set_xlabel('Time (ms)', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Positive Curvature', fontsize=11, fontweight='bold')
    ax3.set_title('Positive Curvature (Upward Bending)', fontsize=12, fontweight='bold')
    ax3.legend(loc='best', fontsize=9)
    ax3.grid(True, alpha=0.25)
    
    # Panel 4: Peak negative curvature timing
    ax4 = axes[1, 1]
    
    for cond in conditions:
        # Find peaks in negative curvature (make positive for display)
        neg_curv = np.copy(curvature_dict[cond])
        neg_curv[neg_curv > 0] = 0
        neg_curv = -neg_curv  # Flip for visualization
        
        ax4.plot(times * 1000, neg_curv, label=label_mapping[cond],
                color=colors[cond], linewidth=2, alpha=0.8)
    
    ax4.axvline(0, color='black', linestyle='-', linewidth=1.5, alpha=0.7)
    ax4.set_xlabel('Time (ms)', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Negative Curvature (abs)', fontsize=11, fontweight='bold')
    ax4.set_title('Negative Curvature (Downward Bending)', fontsize=12, fontweight='bold')
    ax4.legend(loc='best', fontsize=9)
    ax4.grid(True, alpha=0.25)
    
    # Overall title
    if sample_distance == 1:
        title_suffix = '(no averaging)'
    else:
        title_suffix = f'(averaged over {sample_distance}-sample windows)'
    
    fig.suptitle(f'Curvature Statistical Analysis {title_suffix}',
                fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    # Save if path provided
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved statistics plot to: {save_path}")
    
    return fig


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function."""
    
    print("=" * 70)
    print("ERP CURVATURE ANALYSIS")
    print("=" * 70)
    
    # Load parameters
    params = setup_parameters()
    
    # Load epochs data
    epochs_dict, included_subjects = load_epochs_data(
        params['data_path'],
        params['subjects'],
        params['baseline']
    )
    
    if not epochs_dict:
        print("No data loaded. Exiting.")
        return
    
    # Build grand averages per condition
    grand_avg_dict = build_condition_evoked(
        epochs_dict,
        params['conditions'],
        included_subjects
    )
    
    # Create plots directory
    plots_dir = os.path.join(params['project_root'], 'plots', 'curvature')
    os.makedirs(plots_dir, exist_ok=True)
    
    # Generate exploratory curvature visualizations for selected sample distances
    print("\nGenerating curvature heatmaps and statistical plots...")
    
    # Use a subset of sample distances for detailed visualization
    explore_distances = [1, 10, 25, 50]
    
    for sd in explore_distances:
        print(f"\n  Processing sample_distance = {sd}...")
        
        # Compute curvature for this sample distance
        curvature_dict, times = extract_roi_curvature(
            grand_avg_dict,
            params['roi'],
            params['conditions'],
            sd
        )
        
        # Determine filename suffix
        if sd == 1:
            filename_suffix = ''
        else:
            filename_suffix = f'_sd{sd}'
        
        # Generate heatmap visualization
        print(f"    Creating heatmap...")
        fig_heat = plot_curvature_heatmap(
            curvature_dict,
            times,
            params['label_mapping'],
            sd,
            save_path=os.path.join(plots_dir, f'curvature_heatmap{filename_suffix}.pdf')
        )
        plt.close(fig_heat)
        
        # Generate statistical analysis plots
        print(f"    Creating statistical analysis...")
        fig_stats = plot_curvature_statistics(
            curvature_dict,
            times,
            params['label_mapping'],
            sd,
            save_path=os.path.join(plots_dir, f'curvature_statistics{filename_suffix}.pdf')
        )
        plt.close(fig_stats)
    
    print("\n" + "=" * 70)
    print("CURVATURE ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"Plots saved to: {plots_dir}")
    print(f"Generated plots for window sizes: {explore_distances}")
    print("\nVisualization types created:")
    print("  - Heatmaps: Color-coded curvature intensity across time and conditions")
    print("  - Statistics: Variability, distributions, and positive/negative curvature")
    print("\nInterpretation:")
    print("  - Positive curvature (red) = Signal curving upward (acceleration)")
    print("  - Negative curvature (blue) = Signal curving downward (deceleration)")
    print("  - Bright colors = Strong curvature (rapid direction changes)")
    print("  - Dark/white = Weak/zero curvature (linear or stable trend)")
    print("\n  Note: Larger window sizes smooth the curvature by averaging")
    print("        signal values within each window before computing derivatives.")


if __name__ == "__main__":
    main()
