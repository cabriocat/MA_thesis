# Configuration Template
# =====================
# 
# This file shows the paths that need to be updated in each script
# to run the analysis on your own data.

# ==============================================================================
# PYTHON SCRIPTS CONFIGURATION
# ==============================================================================

# Update these paths in ALL Python scripts (01, 02, 03, 04):

# Path to raw EEG data directory containing subject folders
# Each subject folder should have: preprocessed-epo/{subject}_task-nouns_eeg-epo_100_70_cutoff.fif
DATA_PATH = "/path/to/your/data/1_raw"

# Expected directory structure:
# /path/to/your/data/1_raw/
# ├── sub-bravoc001/
# │   └── preprocessed-epo/
# │       └── sub-bravoc001_task-nouns_eeg-epo_100_70_cutoff.fif
# ├── sub-bravoc002/
# │   └── preprocessed-epo/
# │       └── sub-bravoc002_task-nouns_eeg-epo_100_70_cutoff.fif
# └── ...

# ==============================================================================
# R SCRIPTS CONFIGURATION  
# ==============================================================================

# Update this path in R/01_data_preparation.R:

# Path to the combined CSV file created by Python script 03_data_preparation.py
DATA_CSV_PATH <- "/path/to/your/data/2_preprocessed/285-345ms/combined_task-nouns_285-345ms.csv"

# Update these paths in R/04_data_visualization.R:

# Paths to SNR analysis results created by Python script 04_snr_analysis.py
SNR_DATA_PATH <- "/path/to/your/data/3_analysis/cumulative_snr_data.csv"
F_DATA_PATH <- "/path/to/your/data/3_analysis/f_values_data.csv"

# Output directory for plots
OUTPUT_DIR <- "/path/to/your/data/4_plots"

# ==============================================================================
# SUBJECT LIST CONFIGURATION
# ==============================================================================

# If you have different subjects, update the subjects list in Python scripts:

SUBJECTS = [
    "sub-yourprefix001", "sub-yourprefix002", "sub-yourprefix003",
    # ... add your subject IDs here
]

# ==============================================================================
# FILE NAMING CONVENTIONS
# ==============================================================================

# The scripts expect specific file naming conventions:

# Raw EEG files:
# {subject}_task-nouns_eeg-epo_100_70_cutoff.fif

# If your files have different names, update the filename pattern in:
# - Python/01_data_loading_setup.py (line ~140)
# - Python/03_data_preparation.py (line ~190)  
# - Python/04_snr_analysis.py (line ~520)

# ==============================================================================
# QUICK SETUP CHECKLIST
# ==============================================================================

# 1. Update DATA_PATH in all Python scripts
# 2. Update DATA_CSV_PATH in R/01_data_preparation.R
# 3. Update SNR data paths in R/04_data_visualization.R
# 4. If needed, update subject list and file naming patterns
# 5. Create output directories if they don't exist:
#    - data/2_preprocessed/285-345ms/
#    - data/3_analysis/
#    - data/4_plots/
