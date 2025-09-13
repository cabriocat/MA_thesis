#!/usr/bin/env python3
"""
ERP Visualization
=================

This script creates grand-average ERP plots for different numbers of repetitions.
It uses the ROI-averaged data and applies the Okabe-Ito color palette for
accessibility.

Tasks performed:
- Generate grand-average evoked responses
- Create ROI-averaged ERP plots
- Export plots in PDF format
- Support for variable number of repetitions (1-6)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mne import grand_average
import mne

# Suppress MNE info messages
mne.set_log_level('WARNING')


# ============================================================================
# DATA LOADING FUNCTIONS (Embedded for independence)
# ============================================================================

def setup_parameters():
    """Define experimental parameters and paths."""
    
    # Project setup - use current directory (MA_thesis) as project root
    project_root = os.path.abspath(os.getcwd())
    
    # Data path - Update this to match your actual data location
    data_path = "/Users/johannberger/Documents/thesis/data/1_raw"
    
    # Validate that the data path exists
    if not os.path.exists(data_path):
        raise FileNotFoundError(
            f"Data path not found: {data_path}\n"
            f"Please update the data_path variable in the setup_parameters function "
            f"to point to your actual data directory."
        )
    
    # Subjects list
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
    
    # Time windows (seconds)
    baseline = (-0.100, 0.0)
    signal_win = (0.285, 0.345)
    
    # Tasks and conditions
    tasks = ["nouns"]
    conditions = {"nouns": ['a', 'f', 't', 'c', 'e', 's']}
    
    # ROI channels (frontocentral region)
    roi = [
        "FC1", "FCz", "FC2", "FCC1h", "FCC2h",
        "C1", "Cz", "C2", "CCP1h", "CCP2h",
        "CP1", "CPz", "CP2", "CPP1h", "CPP2h"
    ]
    
    # Label mapping for conditions
    label_mapping_nouns = {
        "a": "Animal",
        "f": "Food", 
        "t": "Tool",
        "c": "Communication",
        "e": "Emotion",
        "s": "Social"
    }
    
    return {
        'project_root': project_root,
        'data_path': data_path,
        'subjects': subjects,
        'baseline': baseline,
        'signal_win': signal_win,
        'tasks': tasks,
        'conditions': conditions,
        'roi': roi,
        'label_mapping_nouns': label_mapping_nouns
    }


def load_epochs_data(data_path, subjects):
    """
    Load epochs data for all subjects.
    
    Parameters
    ----------
    data_path : str
        Path to the raw data directory
    subjects : list
        List of subject identifiers
        
    Returns
    -------
    epochs_dict : dict
        Dictionary mapping subject IDs to Epochs objects
    included_subjects : list
        List of successfully loaded subjects
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
                
            epo = mne.read_epochs(epochs_file)
            epochs_dict[subject] = epo
            included_subjects.append(subject)
            print(f"✓ Loaded: {subject}")
            
        except Exception as e:
            print(f"✗ Error loading {subject}: {str(e)}")
            continue
    
    print(f"\nSuccessfully loaded {len(included_subjects)} subjects")
    return epochs_dict, included_subjects


def build_evoked_dict(epochs_dict, conditions, included_subjects):
    """
    Build evoked response dictionary organized by repetition and condition.
    
    Parameters
    ----------
    epochs_dict : dict
        Dictionary of subject epochs
    conditions : dict
        Experimental conditions
    included_subjects : list
        List of subjects to include
        
    Returns
    -------
    evoked_dict : dict
        Nested dictionary: repetition -> condition -> list of evoked responses
    """
    
    # Initialize evoked dictionary
    evoked_dict = {
        str(rep): {cond: [] for cond in conditions['nouns']} 
        for rep in range(1, 7)
    }
    
    print("Building evoked response dictionary...")
    
    for subject in included_subjects:
        epo = epochs_dict[subject]
        all_names = list(epo.event_id.keys())
        
        for rep in range(1, 7):
            for cond in conditions['nouns']:
                # Find matching event names for this repetition and condition
                names = [n for n in all_names if f"/{rep}/{cond}/n" in n]
                
                if not names:
                    continue
                    
                # Average across trials for this subject/rep/condition
                ev = epo[names].average()
                evoked_dict[str(rep)][cond].append(ev)
    
    # Report statistics
    print("\nEvoked dictionary statistics:")
    for rep in range(1, 7):
        for cond in conditions['nouns']:
            n_subjects = len(evoked_dict[str(rep)][cond])
            print(f"  Rep {rep}, Condition {cond}: {n_subjects} subjects")
    
    return evoked_dict


# ============================================================================
# ERP PLOTTING FUNCTIONS
# ============================================================================


def create_erp_plot(grand_evokeds, roi, signal_win, n_reps_to_use, 
                   label_mapping_nouns, conditions, save_path=None):
    """
    Create ERP plot with ROI averaging.
    
    Parameters
    ----------
    grand_evokeds : dict
        Grand-averaged evoked responses by condition
    roi : list
        List of ROI channel names
    signal_win : tuple
        Signal window in seconds
    n_reps_to_use : int
        Number of repetitions used
    label_mapping_nouns : dict
        Mapping from condition codes to labels
    conditions : list
        List of condition codes
    save_path : str, optional
        Path to save the plot
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure
    """
    
    # Okabe-Ito color palette (colorblind friendly)
    colors = ['#388E3C', '#1976D2', '#D32F2F', '#d9a00f', '#616161', '#ff00ff']
    
    # ROI averaging
    pick_roi = mne.pick_channels(grand_evokeds['a'].ch_names, include=roi)
    roi_data = {
        cond: grand_evokeds[cond].data[pick_roi, :].mean(axis=0) * 1e6  # Convert to µV
        for cond in conditions
    }
    
    # Create plot
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Plot ERP traces
    for color, cond in zip(colors, conditions):
        ax.plot(
            grand_evokeds[cond].times * 1e3,  # Convert to ms
            roi_data[cond],
            label=label_mapping_nouns[cond], 
            color=color, 
            linewidth=1.7
        )
    
    # Shade the signal window
    ax.axvspan(
        signal_win[0] * 1e3, 
        signal_win[1] * 1e3,
        color='grey', 
        alpha=0.25, 
        zorder=0
    )
    
    # Formatting
    ax.invert_yaxis()  # Negative up convention
    ax.set_ylim(4, -4)
    
    # Zero lines
    ax.axhline(0, color='k', linewidth=0.8)
    ax.axvline(0, color='k', linewidth=0.8)
    
    # Labels and styling
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Amplitude (µV)')
    ax.set_title("")
    ax.legend(loc='upper right', fontsize=8)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
    
    fig.tight_layout()
    
    # Save if path provided
    if save_path:
        fig.savefig(save_path, format='pdf', dpi=300, bbox_inches='tight')
        print(f"Saved plot: {save_path}")
    
    return fig


def generate_multiple_erp_plots(evoked_dict, params, output_dir):
    """
    Generate ERP plots for different numbers of repetitions.
    
    Parameters
    ----------
    evoked_dict : dict
        Evoked responses dictionary
    params : dict
        Experimental parameters
    output_dir : str
        Directory to save plots
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("Generating ERP plots for different repetition counts...")
    
    # Generate plots for 1, 3, and 6 repetitions
    rep_counts = [1, 2, 3, 4, 5, 6]
    
    for n_reps_to_use in rep_counts:
        print(f"\nCreating plot for {n_reps_to_use} repetition(s)...")
        
        # Build evokeds for this repetition count
        evokeds_for_plot = {cond: [] for cond in params['conditions']['nouns']}
        
        for rep in range(1, n_reps_to_use + 1):
            for cond in params['conditions']['nouns']:
                evokeds_for_plot[cond].extend(evoked_dict[str(rep)][cond])
        
        # Check if we have data
        total_evokeds = sum(len(evokeds_for_plot[cond]) for cond in params['conditions']['nouns'])
        if total_evokeds == 0:
            print(f"Warning: No data available for {n_reps_to_use} repetitions")
            continue
        
        # Grand-average per category
        grand_evokeds = {}
        for cond in params['conditions']['nouns']:
            if len(evokeds_for_plot[cond]) > 0:
                grand_evokeds[cond] = grand_average(evokeds_for_plot[cond])
                print(f"  {cond}: {len(evokeds_for_plot[cond])} evoked responses")
        
        # Create plot
        save_path = os.path.join(output_dir, f'ERP_plot_rep_{n_reps_to_use}.pdf')
        
        fig = create_erp_plot(
            grand_evokeds=grand_evokeds,
            roi=params['roi'],
            signal_win=params['signal_win'],
            n_reps_to_use=n_reps_to_use,
            label_mapping_nouns=params['label_mapping_nouns'],
            conditions=params['conditions']['nouns'],
            save_path=save_path
        )
        
        # Display plot
        plt.show()
        plt.close(fig)


def main():
    """Main execution function."""
    
    print("=" * 60)
    print("ERP VISUALIZATION")
    print("=" * 60)
    
    # Load data and setup
    print("Loading data...")
    params = setup_parameters()
    epochs_dict, included_subjects = load_epochs_data(
        params['data_path'], 
        params['subjects']
    )
    evoked_dict = build_evoked_dict(
        epochs_dict, 
        params['conditions'], 
        included_subjects
    )
    
    # Create output directory
    output_dir = os.path.join(params['project_root'], "plots", "erp")
    
    # Generate plots
    generate_multiple_erp_plots(evoked_dict, params, output_dir)
    
    print("\n" + "=" * 60)
    print("ERP VISUALIZATION COMPLETE") 
    print("=" * 60)
    print(f"Plots saved to: {output_dir}")


if __name__ == "__main__":
    main()
