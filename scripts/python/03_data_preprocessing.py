#!/usr/bin/env python3
"""
Data Preprocessing and CSV Export Script

This script handles data preprocessing and CSV export:
- Converting epochs to DataFrame format
- Applying time window filtering
- Exporting per-subject and combined CSV files
- Data transformation for statistical analysis

Author: [Your Name]
"""

import os
import numpy as np
import pandas as pd
import mne

def setup_csv_export_path(project_root, time_window):
    """Setup CSV export directory based on time window."""
    time_start_ms = int(time_window[0] * 1000)
    time_end_ms = int(time_window[1] * 1000)
    
    save_to_csv_path = os.path.join(
        project_root, "data", "2_preprocessed", 
        f"{time_start_ms}-{time_end_ms}ms"
    )
    os.makedirs(save_to_csv_path, exist_ok=True)
    
    return save_to_csv_path, time_start_ms, time_end_ms

def process_subject_epochs(subject, data_path, experiment, baseline_window, time_window):
    """
    Process epochs for a single subject.
    
    Parameters:
    -----------
    subject : str
        Subject identifier
    data_path : str
        Path to raw data
    experiment : str
        Experiment name (e.g., 'nouns')
    baseline_window : tuple
        Baseline time window in seconds
    time_window : tuple
        Analysis time window in seconds
        
    Returns:
    --------
    pd.DataFrame or None
        Processed dataframe for the subject
    """
    
    epochs_file = os.path.join(
        data_path,
        subject,
        f"preprocessed-epo/{subject}_task-{experiment}_eeg-epo_100_70_cutoff.fif"
    )
    
    if not os.path.exists(epochs_file):
        print(f"File not found: {epochs_file}")
        return None

    print(f"Processing: {subject}")
    epochs = mne.read_epochs(epochs_file, preload=True)
    epochs.apply_baseline(baseline_window)

    # To long-format DataFrame
    df = epochs.to_data_frame()
    df.insert(0, 'subject', subject)
    df['repetition'] = df['condition'].apply(lambda x: x.split("/")[1])
    
    # Normalize time column to seconds (handles both ms and s representations)
    df['time'] = df['time'].astype(float)
    if df['time'].abs().max() > 10:  # heuristic: large values imply milliseconds
        df['time'] = df['time'] / 1000.0

    # Filter for window and compute mean
    df_w = df[(df['time'] >= time_window[0]) & (df['time'] <= time_window[1])]
    df_w_mean = df_w.groupby(
        ['epoch', 'subject', 'condition'],
        as_index=False
    ).mean(numeric_only=True)

    return df_w_mean

def export_subject_csvs(included_subjects, data_path, experiment, baseline_window, 
                       time_window, save_to_csv_path, time_start_ms, time_end_ms):
    """
    Export per-subject CSV files and collect data.
    
    Returns:
    --------
    list
        List of DataFrames for all subjects
    """
    
    df_list = []
    
    for subject in included_subjects:
        df_w_mean = process_subject_epochs(
            subject, data_path, experiment, baseline_window, time_window
        )
        
        if df_w_mean is not None:
            # Save individual CSV
            out_file = os.path.join(
                save_to_csv_path,
                f"{subject}_task-{experiment}_{time_start_ms}-{time_end_ms}ms.csv"
            )
            df_w_mean.to_csv(out_file, index=False)
            print(f"Saved: {out_file}")
            df_list.append(df_w_mean)
    
    return df_list

def combine_and_transform_data(df_list, experiment, save_to_csv_path, time_start_ms, time_end_ms):
    """
    Combine all subjects and transform data format.
    
    Parameters:
    -----------
    df_list : list
        List of subject DataFrames
    experiment : str
        Experiment name
    save_to_csv_path : str
        Directory to save combined CSV
    time_start_ms : int
        Start time in milliseconds
    time_end_ms : int
        End time in milliseconds
    """
    
    if not df_list:
        print("No dataframes collected; check your subject list.")
        return None
    
    # Concatenate
    df_all = pd.concat(df_list, ignore_index=True)
    
    # Check for bad conditions
    bad_conditions = df_all[~df_all['condition'].str.contains("/", regex=False)]
    if not bad_conditions.empty:
        print("⚠️ Bad condition rows:")
        print(bad_conditions)

    # Parse condition into separate columns
    parts = df_all['condition'].str.split('/', expand=True)
    df_all['item'] = parts.iloc[:, 0]
    df_all['repetition'] = parts.iloc[:, 1].astype(int)
    df_all['category'] = parts.iloc[:, 2].map({
        'a': 'animal',
        'c': 'commun',
        'e': 'emotion',
        'f': 'food',
        's': 'social',
        't': 'tool'
    })

    # Drop unwanted columns
    df_all.drop(columns=['condition', 'epoch', 'time'], inplace=True)

    # Pivot to long format (channels->voltage)
    id_vars = ['subject', 'item', 'repetition', 'category']
    channel_cols = [c for c in df_all.columns if c not in id_vars]
    df_long = df_all.melt(
        id_vars=id_vars,
        value_vars=channel_cols,
        var_name='channel',
        value_name='voltage'
    )

    # Save combined CSV
    combined_file = os.path.join(
        save_to_csv_path,
        f"combined_task-{experiment}_{time_start_ms}-{time_end_ms}ms.csv"
    )
    df_long.to_csv(combined_file, index=False)
    print(f"Saved combined file: {combined_file}")
    
    return df_long

def preprocess_and_export_data(data_objects):
    """
    Main preprocessing function using data objects from setup script.
    
    Parameters:
    -----------
    data_objects : dict
        Dictionary containing data objects from 01_data_setup.py
    """
    
    # Extract data objects
    included_subjects = data_objects['included_subjects']
    project_root = data_objects['project_root']
    data_path = data_objects['data_path']
    baseline = data_objects['baseline']
    signal_win = data_objects['signal_win']
    tasks = data_objects['tasks']
    
    # Setup
    experiment = tasks[0]
    time_window = signal_win
    baseline_window = baseline
    
    print(f"=== Data Preprocessing and Export ===")
    print(f"Experiment: {experiment}")
    print(f"Time window: {time_window}")
    print(f"Baseline: {baseline_window}")
    print(f"Number of subjects: {len(included_subjects)}")
    
    # Setup CSV export path
    save_to_csv_path, time_start_ms, time_end_ms = setup_csv_export_path(
        project_root, time_window
    )
    print(f"Export directory: {save_to_csv_path}")
    
    # Export per-subject CSVs
    df_list = export_subject_csvs(
        included_subjects, data_path, experiment, baseline_window,
        time_window, save_to_csv_path, time_start_ms, time_end_ms
    )
    
    # Combine and transform data
    df_combined = combine_and_transform_data(
        df_list, experiment, save_to_csv_path, time_start_ms, time_end_ms
    )
    
    print(f"\nPreprocessing complete!")
    print(f"Processed {len(df_list)} subjects")
    if df_combined is not None:
        print(f"Combined dataset shape: {df_combined.shape}")
    
    return df_combined

def main():
    """Main execution function - requires data from 01_data_setup.py"""
    print("=== Data Preprocessing Script ===")
    print("This script requires data objects from 01_data_setup.py")
    print("Run this script by importing and calling functions with your data objects.")
    print("\nExample usage:")
    print("from data_setup import main as load_data")
    print("data = load_data()")
    print("df_combined = preprocess_and_export_data(data)")

if __name__ == "__main__":
    main()
