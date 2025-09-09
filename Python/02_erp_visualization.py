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

# Import from data loading script
try:
    from data_loading_setup import setup_parameters, load_epochs_data, build_evoked_dict
except ImportError:
    # Alternative: run data loading script first or define functions locally
    print("Warning: Could not import from data_loading_setup. Run 01_data_loading_setup.py first.")
    from importlib import import_module
    import sys
    sys.path.append('.')
    data_loading = import_module('01_data_loading_setup')
    setup_parameters = data_loading.setup_parameters
    load_epochs_data = data_loading.load_epochs_data
    build_evoked_dict = data_loading.build_evoked_dict


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
    ax.set_ylim(3, -3)
    
    # Zero lines
    ax.axhline(0, color='k', linewidth=0.8)
    ax.axvline(0, color='k', linewidth=0.8)
    
    # Labels and styling
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Amplitude (µV)')
    ax.set_title(f'Grand-average ERP (ROI) – {n_reps_to_use} repetitions')
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
    rep_counts = [1, 3, 6]
    
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
