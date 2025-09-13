#!/usr/bin/env python3
"""
Data Preparation and Export
============================

This script processes the loaded epochs data and exports it in CSV format
for further analysis. It handles time window extraction, baseline correction,
and data transformation to long format.

Tasks performed:
- Apply baseline correction to epochs
- Extract data for specified time window
- Transform to long format with proper condition parsing
- Export individual subject CSVs and combined dataset
"""

import os
import numpy as np
import pandas as pd
import mne

# Suppress MNE info messages
mne.set_log_level('WARNING')


def process_subject_data(subject, epochs_file, baseline_window, time_window, 
                        experiment, project_root):
    """
    Process epochs data for a single subject.
    
    Parameters
    ----------
    subject : str
        Subject identifier
    epochs_file : str
        Path to subject's epochs file
    baseline_window : tuple
        Baseline time window in seconds
    time_window : tuple
        Signal time window in seconds
    experiment : str
        Experiment name
    project_root : str
        Project root directory
        
    Returns
    -------
    df_w_mean : pd.DataFrame
        Processed dataframe for this subject
    """
    
    if not os.path.exists(epochs_file):
        print(f"File not found: {epochs_file}")
        return None
    
    print(f"Processing: {subject}")
    
    # Load and preprocess epochs
    epochs = mne.read_epochs(epochs_file, preload=True)
    epochs.apply_baseline(baseline_window)
    
    # Convert to DataFrame
    df = epochs.to_data_frame()
    df.insert(0, 'subject', subject)
    df['repetition'] = df['condition'].apply(lambda x: x.split("/")[1])
    
    # Normalize time column to seconds
    df['time'] = df['time'].astype(float)
    if df['time'].abs().max() > 10:  # Heuristic: large values imply milliseconds
        df['time'] = df['time'] / 1000.0
    
    # Filter for time window and compute mean
    df_w = df[(df['time'] >= time_window[0]) & (df['time'] <= time_window[1])]
    df_w_mean = df_w.groupby(
        ['epoch', 'subject', 'condition'],
        as_index=False
    ).mean(numeric_only=True)
    
    # Save individual CSV
    time_start_ms = int(time_window[0] * 1000)
    time_end_ms = int(time_window[1] * 1000)
    save_to_csv_path = os.path.join(
        project_root, "data", "2_preprocessed", f"{time_start_ms}-{time_end_ms}ms"
    )
    os.makedirs(save_to_csv_path, exist_ok=True)
    
    out_file = os.path.join(
        save_to_csv_path,
        f"{subject}_task-{experiment}_{time_start_ms}-{time_end_ms}ms.csv"
    )
    df_w_mean.to_csv(out_file, index=False)
    print(f"Saved: {out_file}")
    
    return df_w_mean


def combine_and_transform_data(df_list, experiment, time_window, project_root):
    """
    Combine all subject dataframes and transform to long format.
    
    Parameters
    ----------
    df_list : list
        List of subject dataframes
    experiment : str
        Experiment name
    time_window : tuple
        Time window in seconds
    project_root : str
        Project root directory
        
    Returns
    -------
    df_long : pd.DataFrame
        Combined long-format dataframe
    """
    
    if not df_list:
        print("No dataframes collected; check your subject list.")
        return None
    
    print("Combining and transforming data...")
    
    # Concatenate all subjects
    df_all = pd.concat(df_list, ignore_index=True)
    
    # Check for bad conditions
    bad_conditions = df_all[~df_all['condition'].str.contains("/", regex=False)]
    if not bad_conditions.empty:
        print("⚠️  Warning: Found bad condition rows:")
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
    
    # Pivot to long format (channels -> voltage)
    id_vars = ['subject', 'item', 'repetition', 'category']
    channel_cols = [c for c in df_all.columns if c not in id_vars]
    df_long = df_all.melt(
        id_vars=id_vars,
        value_vars=channel_cols,
        var_name='channel',
        value_name='voltage'
    )
    
    # Save combined CSV
    time_start_ms = int(time_window[0] * 1000)
    time_end_ms = int(time_window[1] * 1000)
    save_to_csv_path = os.path.join(
        project_root, "data", "2_preprocessed", f"{time_start_ms}-{time_end_ms}ms"
    )
    
    combined_file = os.path.join(
        save_to_csv_path,
        f"combined_task-{experiment}_{time_start_ms}-{time_end_ms}ms.csv"
    )
    df_long.to_csv(combined_file, index=False)
    print(f"Saved combined file: {combined_file}")
    
    return df_long


def main():
    """Main execution function."""
    
    print("=" * 60)
    print("DATA PREPARATION AND EXPORT")
    print("=" * 60)
    
    # Setup parameters (inline to avoid import issues)
    project_root = os.path.abspath(os.getcwd())
    
    # Data path - Update this to match your actual data location
    data_path = "/Users/johannberger/Documents/thesis/data/1_raw"
    
    # Validate path exists
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data path not found: {data_path}")
    
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
    
    # Analysis parameters
    experiment = "nouns"
    baseline_window = (-0.100, 0.0)
    time_window = (0.285, 0.345)
    
    print(f"Processing {len(subjects)} subjects")
    print(f"Time window: {time_window[0]*1000:.0f}-{time_window[1]*1000:.0f} ms")
    print(f"Baseline: {baseline_window[0]*1000:.0f}-{baseline_window[1]*1000:.0f} ms")
    
    # Process each subject
    df_list = []
    
    for subject in subjects:
        epochs_file = os.path.join(
            data_path,
            subject,
            f"preprocessed-epo/{subject}_task-{experiment}_eeg-epo_100_70_cutoff.fif"
        )
        
        df_subject = process_subject_data(
            subject=subject,
            epochs_file=epochs_file,
            baseline_window=baseline_window,
            time_window=time_window,
            experiment=experiment,
            project_root=project_root
        )
        
        if df_subject is not None:
            df_list.append(df_subject)
    
    # Combine and transform data
    df_long = combine_and_transform_data(
        df_list=df_list,
        experiment=experiment,
        time_window=time_window,
        project_root=project_root
    )
    
    print("\n" + "=" * 60)
    print("DATA PREPARATION COMPLETE")
    print("=" * 60)
    print(f"Processed {len(df_list)} subjects successfully")
    if df_long is not None:
        print(f"Combined dataset shape: {df_long.shape}")
        print(f"Categories: {sorted(df_long['category'].unique())}")
        print(f"Repetitions: {sorted(df_long['repetition'].unique())}")
        print(f"Subjects: {len(df_long['subject'].unique())}")


if __name__ == "__main__":
    main()
