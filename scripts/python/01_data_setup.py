#!/usr/bin/env python3
"""
Data Setup and Epoch Loading Script

This script handles the initial setup of the EEG analysis pipeline:
- Loading MNE epochs from preprocessed data
- Building evoked dictionaries by repetition and condition
- Setting up project paths and parameters

Author: [Your Name]
"""

import os
import numpy as np
import mne
from copy import deepcopy

# Set MNE logging level
mne.set_log_level('WARNING')

def setup_project_paths():
    """Setup and return project paths."""
    project_root = os.path.abspath(os.path.join(os.getcwd(), "..", ".."))
    data_path = os.path.join(project_root, "data", "1_raw")
    return project_root, data_path

def define_parameters():
    """Define experimental parameters."""
    # Time windows (seconds)
    baseline = (-0.100, 0.0)
    signal_win = (0.285, 0.345)
    
    # Tasks and conditions
    tasks = ["nouns"]
    conditions = {"nouns": ['a', 'f', 't', 'c', 'e', 's']}
    
    # ROI channels
    roi = [
        "FC1", "FCz", "FC2", "FCC1h", "FCC2h",
        "C1", "Cz", "C2", "CCP1h", "CCP2h",
        "CP1", "CPz", "CP2", "CPP1h", "CPP2h"
    ]
    
    # Label mapping
    label_mapping_nouns = {
        "a": "Animal",
        "f": "Food",
        "t": "Tool",
        "c": "Communication",
        "e": "Emotion",
        "s": "Social"
    }
    
    return baseline, signal_win, tasks, conditions, roi, label_mapping_nouns

def define_subjects():
    """Define subject list."""
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
    return subjects

def load_epochs_and_build_evoked_dict(subjects, data_path, conditions):
    """
    Load epochs for all subjects and build evoked dictionary.
    
    Returns:
        epochs_dict: Dictionary mapping subject -> Epochs
        evoked_dict: Dictionary mapping repetition -> condition -> list of evokeds
        included_subjects: List of successfully processed subjects
    """
    # Storage
    epochs_dict = {}
    evoked_dict = {str(rep): {cond: [] for cond in conditions['nouns']} for rep in range(1, 7)}
    included_subjects = []

    for subject in subjects:
        try:
            epo_path = os.path.join(
                data_path, 
                f"{subject}/preprocessed-epo/{subject}_task-nouns_eeg-epo_100_70_cutoff.fif"
            )
            epo = mne.read_epochs(epo_path)
            epochs_dict[subject] = epo

            all_names = list(epo.event_id.keys())
            for rep in range(1, 7):
                for cond in conditions['nouns']:
                    names = [n for n in all_names if f"/{rep}/{cond}/n" in n]
                    if not names:
                        continue
                    ev = epo[names].average()
                    evoked_dict[str(rep)][cond].append(ev)

            included_subjects.append(subject)
            print(f"Loaded: {subject}")
            
        except Exception as e:
            print(f"Failed to load {subject}: {e}")
            continue

    print(f"\nSuccessfully loaded {len(included_subjects)} subjects")
    return epochs_dict, evoked_dict, included_subjects

def main():
    """Main execution function."""
    print("=== EEG Data Setup ===")
    
    # Setup paths and parameters
    project_root, data_path = setup_project_paths()
    baseline, signal_win, tasks, conditions, roi, label_mapping_nouns = define_parameters()
    subjects = define_subjects()
    
    print(f"Project root: {project_root}")
    print(f"Data path: {data_path}")
    print(f"Number of subjects: {len(subjects)}")
    print(f"Conditions: {conditions}")
    print(f"Signal window: {signal_win}")
    print(f"Baseline: {baseline}")
    print(f"ROI channels: {len(roi)} channels")
    
    # Load data
    epochs_dict, evoked_dict, included_subjects = load_epochs_and_build_evoked_dict(
        subjects, data_path, conditions
    )
    
    # Save processed data (optional - could use pickle or other format)
    print("\nData loading complete!")
    print("Data structures created:")
    print(f"  - epochs_dict: {len(epochs_dict)} subjects")
    print(f"  - evoked_dict: {len(evoked_dict)} repetitions")
    print(f"  - included_subjects: {len(included_subjects)} subjects")
    
    return {
        'epochs_dict': epochs_dict,
        'evoked_dict': evoked_dict,
        'included_subjects': included_subjects,
        'project_root': project_root,
        'data_path': data_path,
        'baseline': baseline,
        'signal_win': signal_win,
        'tasks': tasks,
        'conditions': conditions,
        'roi': roi,
        'label_mapping_nouns': label_mapping_nouns
    }

if __name__ == "__main__":
    data_objects = main()
