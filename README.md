# EEG Analysis Pipeline - Restructured Scripts

This repository contains a reorganized and structured analysis pipeline for EEG data processing and statistical analysis. The previous Jupyter notebooks have been converted into modular Python and R scripts for better reproducibility and maintainability.

## 📁 Directory Structure

```
MA_thesis/
├── scripts/
│   ├── python/           # Python analysis pipeline
│   │   ├── data_setup.py         # Data loading and epoch processing
│   │   ├── erp_plotting.py       # ERP visualization and plotting
│   │   ├── data_preprocessing.py # Data export and CSV generation
│   │   ├── snr_analysis.py       # Signal-to-noise ratio analysis
│   │   └── main_pipeline.py      # Main orchestration script
│   └── R/                # R statistical analysis pipeline
│       ├── 01_data_preparation.R     # Data loading and setup
│       ├── 02_model_selection.R      # Statistical model comparison
│       ├── 03_statistical_analysis.R # ANOVA and pairwise tests
│       ├── 04_visualization.R        # Publication-ready plots
│       └── main_r_pipeline.R         # Main orchestration script
├── data/
│   ├── 1_raw/            # Raw EEG data (epochs)
│   ├── 2_preprocessed/   # Preprocessed CSV files
│   ├── 3_analysis/       # Analysis results and exports
│   └── 4_plots/          # Generated plots and figures
├── Python/               # Original Python notebook (archived)
└── R/                    # Original R notebook (archived)
```

## 🐍 Python Pipeline

The Python pipeline handles EEG data processing, visualization, and SNR analysis:

### Core Scripts

1. **`data_setup.py`** - Data Loading and Setup
   - Loads MNE epochs from preprocessed data
   - Builds evoked dictionaries by repetition and condition
   - Sets up project paths and experimental parameters

2. **`erp_plotting.py`** - ERP Visualization
   - Creates grand-average ERPs across repetitions
   - Generates ROI-based ERP plots
   - Produces publication-ready figures

3. **`data_preprocessing.py`** - Data Export and Preprocessing
   - Converts epochs to DataFrame format
   - Applies time window filtering
   - Exports per-subject and combined CSV files

4. **`snr_analysis.py`** - Signal-to-Noise Ratio Analysis
   - Computes cumulative SNR across repetitions
   - Performs statistical testing (paired t-tests, TOST)
   - Applies Holm-Bonferroni correction
   - Exports data for R visualization

5. **`main_pipeline.py`** - Main Orchestration Script
   - Coordinates the entire Python pipeline
   - Provides interactive and batch execution modes

### Usage

#### Full Pipeline
```bash
cd scripts/python
python main_pipeline.py --full
```

#### Interactive Mode
```bash
cd scripts/python
python main_pipeline.py --interactive
```

#### Individual Scripts
```python
# In Python
from data_setup import main as setup_data
from erp_plotting import create_multiple_erp_plots
from data_preprocessing import preprocess_and_export_data
from snr_analysis import perform_snr_analysis

# Load data
data_objects = setup_data()

# Create plots
create_multiple_erp_plots(
    evoked_dict=data_objects['evoked_dict'],
    conditions=data_objects['conditions'],
    roi=data_objects['roi'],
    label_mapping_nouns=data_objects['label_mapping_nouns']
)

# Preprocess and export
df_combined = preprocess_and_export_data(data_objects)

# Perform SNR analysis
summary, fig = perform_snr_analysis(data_objects)
```

## 📊 R Pipeline

The R pipeline handles statistical modeling, analysis, and publication-ready visualization:

### Core Scripts

1. **`01_data_preparation.R`** - Data Setup and Model Preparation
   - Loads required libraries (lme4, tidyverse, etc.)
   - Reads preprocessed EEG data from Python pipeline
   - Creates cumulative datasets for different repetition ranges
   - Sets up experimental parameters and validation

2. **`02_model_selection.R`** - Statistical Model Selection
   - Fits mixed-effects models for each repetition subset
   - Compares models using likelihood ratio tests
   - Selects best-fitting models based on AIC/BIC and significance

3. **`03_statistical_analysis.R`** - Statistical Analysis
   - ANOVA tests for main effect of category across repetitions
   - Pairwise comparisons with multiple comparison corrections
   - Interaction analysis for repetition × category effects
   - Comprehensive summary tables and interpretation

4. **`04_visualization.R`** - Data Visualization
   - ERP bar plots with significance brackets
   - Line plots showing trends across repetitions
   - SNR analysis plots from Python data export
   - Interaction plots for repetition × category effects
   - Publication-ready themes and styling

5. **`main_r_pipeline.R`** - Main Orchestration Script
   - Coordinates the entire R analysis pipeline
   - Provides interactive and batch execution modes

### Usage

#### Full Pipeline
```bash
cd scripts/R
Rscript main_r_pipeline.R --full
```

#### Interactive Mode
```r
# In R
setwd("scripts/R")
source("main_r_pipeline.R")
run_interactive_r_pipeline()
```

#### Individual Scripts
```r
# Source all scripts
source("01_data_preparation.R")
source("02_model_selection.R")
source("03_statistical_analysis.R")
source("04_visualization.R")

# Run pipeline steps
setup_results <- setup_data_and_models()
model_results <- run_model_selection(setup_results$datasets)
analysis_results <- run_statistical_analysis(model_results$final_models)
plots <- create_all_plots(setup_results$datasets, analysis_results$pairwise_results, setup_results$params)
```

## 🔄 Complete Workflow

### Prerequisites
- **Python**: MNE-Python, pandas, numpy, matplotlib, scipy
- **R**: tidyverse, lme4, lmerTest, emmeans, ggplot2

### Step-by-Step Execution

1. **Run Python Pipeline** (generates preprocessed data)
   ```bash
   cd scripts/python
   python main_pipeline.py --full
   ```

2. **Run R Pipeline** (performs statistical analysis)
   ```bash
   cd scripts/R
   Rscript main_r_pipeline.R --full
   ```

3. **Check Results**
   - Preprocessed data: `data/2_preprocessed/`
   - Analysis results: `data/3_analysis/`
   - Generated plots: `data/4_plots/`

## 📈 Key Outputs

### Python Pipeline Outputs
- **ERP Plots**: Grand-average ERPs for different repetition numbers
- **SNR Analysis**: Cumulative SNR plots and statistical tests
- **Preprocessed Data**: CSV files ready for R analysis
- **Export Files**: SNR data and F-values for R visualization

### R Pipeline Outputs
- **Model Comparison**: Best-fitting models for each repetition subset
- **Statistical Analysis**: ANOVA results and pairwise comparisons
- **Visualization**: Publication-ready ERP plots with significance brackets
- **Summary Tables**: Comprehensive analysis summaries

## 🔍 Key Findings

The analysis pipeline supports the following key findings:

1. **Category Effect**: Semantic category shows significant effects from the first repetition
2. **Optimal Repetitions**: Statistical power stabilizes after 3 repetitions
3. **Diminishing Returns**: No substantial gain in SNR or statistical power beyond 3 repetitions
4. **Recommendation**: Use 3 repetitions for optimal efficiency in future experiments

## 🛠️ Advantages of Structured Approach

### Over Notebooks
- **Modularity**: Each script handles a specific task
- **Reusability**: Functions can be easily imported and reused
- **Version Control**: Better tracking of changes in individual components
- **Testing**: Easier to test individual components
- **Collaboration**: Multiple people can work on different scripts simultaneously

### Over Monolithic Scripts
- **Maintainability**: Easier to debug and modify specific functionalities
- **Scalability**: New analyses can be added as separate scripts
- **Documentation**: Each script is self-documented with clear purposes
- **Flexibility**: Can run individual steps or full pipeline as needed

## 📝 Notes

- The original notebooks are preserved in the `Python/` and `R/` directories
- All file paths are designed to work from the script directories
- The pipeline creates necessary directories automatically
- Error handling and validation are built into each script
- Progress indicators and logging help track execution

## 🆘 Troubleshooting

### Common Issues
1. **Missing Dependencies**: Install required Python/R packages
2. **Path Issues**: Ensure you're running scripts from the correct directories
3. **Data Not Found**: Run Python pipeline before R pipeline
4. **Memory Issues**: Large datasets may require more RAM

### Getting Help
- Check individual script documentation
- Use `--help` flags for command-line usage
- Review error messages for specific issues
- Ensure all prerequisites are installed
