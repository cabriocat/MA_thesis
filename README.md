# MA Thesis Analysis Scripts

This repository contains reorganized analysis scripts for EEG/ERP data analysis, converted from Jupyter notebooks to structured Python and R scripts.

## Structure

### Python Scripts (in `Python/` folder)

1. **`01_data_loading_setup.py`** - Data Loading and Setup
   - Loads epochs data for all subjects
   - Creates evoked response dictionaries
   - Sets up experimental parameters and ROI definitions

2. **`02_erp_visualization.py`** - ERP Visualization
   - Creates grand-average ERP plots
   - Supports different repetition counts (1, 3, 6)
   - Uses colorblind-friendly Okabe-Ito palette

3. **`03_data_preparation.py`** - Data Preparation and Export
   - Processes epochs data and applies baseline correction
   - Exports data in CSV format for R analysis
   - Transforms data to long format

4. **`04_snr_analysis.py`** - Signal-to-Noise Ratio Analysis
   - Computes cumulative SNR across repetitions
   - Performs statistical testing (paired t-tests, TOST)
   - Exports results for R visualization

### R Scripts (in `R/` folder)

1. **`01_data_preparation.R`** - Data Preparation and Setup
   - Loads preprocessed CSV data from Python
   - Creates cumulative datasets (1, 1-2, 1-3, ..., 1-6 repetitions)
   - Sets up contrast coding for statistical analysis

2. **`02_model_selection.R`** - Model Selection
   - Fits multiple candidate models for each dataset
   - Compares models using likelihood ratio tests
   - Selects best model based on AIC/BIC and significance

3. **`03_statistical_analysis.R`** - Statistical Analysis
   - Conducts Type III ANOVA tests
   - Performs pairwise comparisons with multiple comparison corrections
   - Summarizes statistical findings across repetition counts

4. **`04_data_visualization.R`** - Data Visualization
   - Creates ERP bar plots with significance brackets
   - Generates line plots showing trends across repetitions
   - Imports and visualizes SNR analysis results from Python
   - Exports plots in PDF and SVG formats

## Usage

### Python Workflow

Run the Python scripts in order:

```bash
cd Python/
python 01_data_loading_setup.py
python 02_erp_visualization.py
python 03_data_preparation.py
python 04_snr_analysis.py
```

### R Workflow

Run the R scripts in order:

```bash
cd R/
Rscript 01_data_preparation.R
Rscript 02_model_selection.R
Rscript 03_statistical_analysis.R
Rscript 04_data_visualization.R
```

### Dependencies

**Python:**
- mne
- numpy
- pandas
- matplotlib
- scipy

**R:**
- tidyverse
- lme4
- lmerTest
- emmeans
- ggplot2
- broom.mixed
- glue
- svglite (optional, for SVG export)

## Key Features

### Modular Design
- Each script handles one specific task
- Clear separation of concerns
- Easy to run individual components

### Error Handling
- Robust error handling and validation
- Informative progress messages
- Graceful handling of missing data

### Documentation
- Comprehensive docstrings/comments
- Clear function signatures
- Usage examples

### Output Management
- Organized output directory structure
- Multiple export formats (CSV, PDF, SVG)
- Consistent file naming conventions

## Data Flow

1. **Python scripts** process raw EEG data:
   - Load and validate epochs
   - Create visualizations
   - Export processed data to `data/2_preprocessed/`
   - Perform SNR analysis and export to `data/3_analysis/`

2. **R scripts** perform statistical analysis:
   - Load processed data from Python
   - Fit statistical models
   - Conduct hypothesis tests
   - Create publication-ready plots in `data/4_plots/`

## Directory Structure

```
MA_thesis/
├── Python/
│   ├── 01_data_loading_setup.py
│   ├── 02_erp_visualization.py
│   ├── 03_data_preparation.py
│   └── 04_snr_analysis.py
├── R/
│   ├── 01_data_preparation.R
│   ├── 02_model_selection.R
│   ├── 03_statistical_analysis.R
│   └── 04_data_visualization.R
└── data/
    ├── 1_raw/              # Raw EEG data
    ├── 2_preprocessed/     # Processed data (CSV)
    ├── 3_analysis/         # Analysis results
    └── 4_plots/            # Generated plots
```

## Notes

- Update file paths in scripts to match your directory structure
- The lint warnings in R scripts are expected (tidyverse syntax)
- Scripts can be run individually or as part of the complete workflow
- Results from Python SNR analysis are imported into R for visualization

## Original Research Question

The analysis investigates the optimal number of stimulus repetitions for ERP experiments, comparing statistical power and effect sizes across 1-6 repetitions for semantic category detection tasks.
