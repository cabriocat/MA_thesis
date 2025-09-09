#!/usr/bin/env python3
"""
Main Pipeline Script

This script orchestrates the entire Python analysis pipeline:
- Data setup and loading
- ERP plotting
- Data preprocessing and export
- SNR analysis

Author: [Your Name]
"""

import sys
import os

# Add the current directory to the Python path to import local modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import our modules
from data_setup import main as setup_data
from erp_plotting import create_erp_plot, create_multiple_erp_plots
from data_preprocessing import preprocess_and_export_data
from snr_analysis import perform_snr_analysis

def run_full_pipeline():
    """Run the complete analysis pipeline."""
    
    print("========================================")
    print("    EEG ANALYSIS PIPELINE - PYTHON     ")
    print("========================================")
    print()
    
    # Step 1: Data Setup
    print("STEP 1: Data Setup and Loading")
    print("-" * 30)
    data_objects = setup_data()
    print("✓ Data setup complete\n")
    
    # Step 2: ERP Plotting
    print("STEP 2: ERP Plotting")
    print("-" * 20)
    output_dir = os.path.join(data_objects['project_root'], "data", "4_plots")
    create_multiple_erp_plots(
        evoked_dict=data_objects['evoked_dict'],
        conditions=data_objects['conditions'],
        roi=data_objects['roi'],
        label_mapping_nouns=data_objects['label_mapping_nouns'],
        output_dir=output_dir
    )
    print("✓ ERP plotting complete\n")
    
    # Step 3: Data Preprocessing
    print("STEP 3: Data Preprocessing and Export")
    print("-" * 37)
    df_combined = preprocess_and_export_data(data_objects)
    print("✓ Data preprocessing complete\n")
    
    # Step 4: SNR Analysis
    print("STEP 4: SNR Analysis")
    print("-" * 20)
    summary, fig = perform_snr_analysis(data_objects)
    print("✓ SNR analysis complete\n")
    
    print("========================================")
    print("    PIPELINE EXECUTION COMPLETE!       ")
    print("========================================")
    print()
    print("Output files have been created in:")
    print(f"  - ERP plots: {output_dir}")
    print(f"  - Preprocessed data: {os.path.join(data_objects['project_root'], 'data', '2_preprocessed')}")
    print(f"  - Analysis results: {os.path.join(data_objects['project_root'], 'data', '3_analysis')}")
    
    return {
        'data_objects': data_objects,
        'df_combined': df_combined,
        'snr_summary': summary,
        'snr_figure': fig
    }

def run_individual_steps():
    """Run individual pipeline steps interactively."""
    
    print("========================================")
    print("  EEG ANALYSIS PIPELINE - INTERACTIVE  ")
    print("========================================")
    print()
    print("Available steps:")
    print("1. Data Setup")
    print("2. ERP Plotting")
    print("3. Data Preprocessing")
    print("4. SNR Analysis")
    print("5. Run Full Pipeline")
    print()
    
    while True:
        choice = input("Enter step number (1-5) or 'q' to quit: ").strip()
        
        if choice.lower() == 'q':
            break
        elif choice == '1':
            data_objects = setup_data()
            print("Data objects stored in 'data_objects' variable")
        elif choice == '2':
            if 'data_objects' in locals():
                create_multiple_erp_plots(
                    evoked_dict=data_objects['evoked_dict'],
                    conditions=data_objects['conditions'],
                    roi=data_objects['roi'],
                    label_mapping_nouns=data_objects['label_mapping_nouns']
                )
            else:
                print("Please run Step 1 (Data Setup) first!")
        elif choice == '3':
            if 'data_objects' in locals():
                df_combined = preprocess_and_export_data(data_objects)
                print("Combined dataframe stored in 'df_combined' variable")
            else:
                print("Please run Step 1 (Data Setup) first!")
        elif choice == '4':
            if 'data_objects' in locals():
                summary, fig = perform_snr_analysis(data_objects)
                print("SNR results stored in 'summary' and 'fig' variables")
            else:
                print("Please run Step 1 (Data Setup) first!")
        elif choice == '5':
            return run_full_pipeline()
        else:
            print("Invalid choice. Please enter 1-5 or 'q'.")

def main():
    """Main function with command line interface."""
    
    if len(sys.argv) > 1:
        if sys.argv[1] == '--full':
            return run_full_pipeline()
        elif sys.argv[1] == '--interactive':
            return run_individual_steps()
        elif sys.argv[1] == '--help':
            print("Usage:")
            print("  python main_pipeline.py --full        # Run full pipeline")
            print("  python main_pipeline.py --interactive # Interactive mode")
            print("  python main_pipeline.py --help        # Show this help")
            return
    else:
        print("Usage:")
        print("  python main_pipeline.py --full        # Run full pipeline")
        print("  python main_pipeline.py --interactive # Interactive mode")
        print("  python main_pipeline.py --help        # Show this help")
        print()
        print("Running in interactive mode by default...")
        return run_individual_steps()

if __name__ == "__main__":
    results = main()
