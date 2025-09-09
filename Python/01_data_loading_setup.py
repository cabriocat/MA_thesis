#!/usr/bin/env python3
"""
Data Loading and Setup
======================

This script handles the initial data loading, subject configuration, and 
creation of evoked response dictionaries for EEG analysis.

Tasks performed:
- Load epochs data for all subjects
- Build evoked response dictionaries organized by repetition and condition
- Set up experimental parameters and ROI definitions
"""

import os
import numpy as np
import mne
from mne import grand_average
from copy import deepcopy

# Suppress MNE info messages
mne.set_log_level('WARNING')

# ============================================================================
# CONFIGURATION SECTION - MODIFY THESE PATHS FOR YOUR DATA
# ============================================================================

# Path to your raw data directory containing subject folders
# Each subject folder should contain a 'preprocessed-epo' subdirectory
# with .fif files named like: {subject}_task-nouns_eeg-epo_100_70_cutoff.fif
DATA_PATH = "/Users/johannberger/Documents/thesis/data/1_raw"

# Alternative: Use relative path if data is in your repository
# DATA_PATH = os.path.join(os.getcwd(), "..", "data", "1_raw")

# ============================================================================


def setup_parameters():
    """Define experimental parameters and paths."""
    
    # Project setup
    project_root = os.path.abspath(os.path.join(os.getcwd(), ".."))
    
    # Use the configured data path
    data_path = DATA_PATH
    
    # Validate that the data path exists
    if not os.path.exists(data_path):
        raise FileNotFoundError(
            f"Data path not found: {data_path}\n"
            f"Please update the DATA_PATH variable at the top of this script "
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


def main():
    """Main execution function."""
    
    print("=" * 60)
    print("DATA LOADING AND SETUP")
    print("=" * 60)
    
    # Setup parameters
    params = setup_parameters()
    print(f"Project root: {params['project_root']}")
    print(f"Data path: {params['data_path']}")
    print(f"Number of subjects: {len(params['subjects'])}")
    print(f"ROI channels: {len(params['roi'])} channels")
    
    # Load epochs data
    epochs_dict, included_subjects = load_epochs_data(
        params['data_path'], 
        params['subjects']
    )
    
    # Build evoked dictionary
    evoked_dict = build_evoked_dict(
        epochs_dict, 
        params['conditions'], 
        included_subjects
    )
    
    print("\n" + "=" * 60)
    print("SETUP COMPLETE")
    print("=" * 60)
    print(f"Loaded data for {len(included_subjects)} subjects")
    print(f"Created evoked responses for {len(params['conditions']['nouns'])} conditions")
    print("Ready for analysis and visualization")
    
    # Return all necessary variables for subsequent scripts
    return {
        'params': params,
        'epochs_dict': epochs_dict,
        'evoked_dict': evoked_dict,
        'included_subjects': included_subjects
    }


if __name__ == "__main__":
    results = main()
