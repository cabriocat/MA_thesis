#!/usr/bin/env python3
"""
ERP Plotting and Visualization Script

This script handles ERP plotting and visualization:
- Creating grand-average ERPs across repetitions
- ROI-based ERP plotting with customizable repetition numbers
- Generating publication-ready plots

Author: [Your Name]
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import mne
from mne import grand_average

def create_erp_plot(evoked_dict, conditions, roi, label_mapping_nouns, n_reps_to_use=6, save_path=None):
    """
    Create ERP plot with specified number of repetitions.
    
    Parameters:
    -----------
    evoked_dict : dict
        Dictionary mapping repetition -> condition -> list of evokeds
    conditions : dict
        Dictionary mapping task -> list of conditions
    roi : list
        List of ROI channel names
    label_mapping_nouns : dict
        Mapping from condition codes to readable labels
    n_reps_to_use : int
        Number of repetitions to include (1-6)
    save_path : str, optional
        Path to save the plot
    """
    
    assert 1 <= n_reps_to_use <= 6, "n_reps_to_use must be between 1 and 6"
    
    print(f"Creating ERP plot with {n_reps_to_use} repetitions")
    
    # Build a dict with the evokeds we actually want
    evokeds_for_plot = {cond: [] for cond in conditions['nouns']}
    for rep in range(1, n_reps_to_use + 1):
        for cond in conditions['nouns']:
            evokeds_for_plot[cond].extend(evoked_dict[str(rep)][cond])

    # Grand-average per category
    grand_evokeds = {cond: grand_average(evokeds_for_plot[cond])
                     for cond in conditions['nouns']}

    # Okabe–Ito palette
    colors = ['#388E3C', '#1976D2', '#D32F2F', '#d9a00f',
              '#616161', '#ff00ff']

    # ROI AVERAGE
    pick_roi = mne.pick_channels(grand_evokeds['a'].ch_names, include=roi)
    roi_data = {
        cond: grand_evokeds[cond].data[pick_roi, :].mean(axis=0) * 1e6   # µV
        for cond in conditions['nouns']
    }

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 5))
    for color, cond in zip(colors, conditions['nouns']):
        ax.plot(grand_evokeds[cond].times * 1e3, roi_data[cond],
                label=label_mapping_nouns[cond], color=color, linewidth=1.7)

    # Shade the signal window (285-345ms)
    signal_win = (0.285, 0.345)
    ax.axvspan(signal_win[0]*1e3, signal_win[1]*1e3,
               color='grey', alpha=0.25, zorder=0)

    # Flip y-axis (negative up) + manual limits
    ax.invert_yaxis()
    ax.set_ylim(3, -3)

    # Zero lines
    ax.axhline(0, color='k', linewidth=0.8)
    ax.axvline(0, color='k', linewidth=0.8)

    # Labels & finish
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Amplitude (µV)')
    ax.set_title(f'Grand-average ERP (ROI) – {n_reps_to_use} repetitions')
    ax.legend(loc='upper right', fontsize=8)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
    fig.tight_layout()
    
    # Save plot if path provided
    if save_path:
        fig.savefig(save_path, format='pdf')
        print(f"Plot saved to: {save_path}")
    
    plt.show()
    return fig

def create_multiple_erp_plots(evoked_dict, conditions, roi, label_mapping_nouns, output_dir=None):
    """
    Create ERP plots for different repetition numbers.
    
    Parameters:
    -----------
    evoked_dict : dict
        Dictionary mapping repetition -> condition -> list of evokeds
    conditions : dict
        Dictionary mapping task -> list of conditions
    roi : list
        List of ROI channel names
    label_mapping_nouns : dict
        Mapping from condition codes to readable labels
    output_dir : str, optional
        Directory to save plots
    """
    
    repetition_numbers = [1, 3, 6]  # Key repetition numbers to plot
    
    for n_reps in repetition_numbers:
        save_path = None
        if output_dir:
            import os
            os.makedirs(output_dir, exist_ok=True)
            save_path = os.path.join(output_dir, f'ERP_plot_{n_reps}_reps.pdf')
        
        create_erp_plot(
            evoked_dict=evoked_dict,
            conditions=conditions,
            roi=roi,
            label_mapping_nouns=label_mapping_nouns,
            n_reps_to_use=n_reps,
            save_path=save_path
        )

def main():
    """Main execution function - requires data from 01_data_setup.py"""
    print("=== ERP Plotting Script ===")
    print("This script requires data objects from 01_data_setup.py")
    print("Run this script by importing and calling functions with your data objects.")
    print("\nExample usage:")
    print("from data_setup import main as load_data")
    print("data = load_data()")
    print("create_erp_plot(data['evoked_dict'], data['conditions'], data['roi'], data['label_mapping_nouns'])")

if __name__ == "__main__":
    main()
