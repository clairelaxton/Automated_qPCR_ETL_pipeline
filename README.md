# Automated_qPCR_ETL_pipeline
A custom Python3 extract-transform-load pipeline for qPCR data processing. This pipeline allows for rapid compilation and QC of multiple raw qPCR data files into one clean data sheet, ready for downstream analysis.

This pipeline is designed to automate the processing of raw qpcr export files, perform QC and collate into a clean sheet. We have designed it to be used by people with minimal coding experience.
Instructions
1.	Install Python and required libraries

a. Go to https://realpython.com/installing-python/ and install

b. Open your terminal and run: python3 -m pip install pandas openpyxl.

This installs pandas for spreadsheet and table operations and openpyxl for saving Excel (.xlsx) files


2.	Create a new directory, [your_project_folder] and add your files using the following structure (don’t edit file names):
<img width="468" height="185" alt="image" src="https://github.com/user-attachments/assets/07f74432-9c97-4723-9073-5252edfffb8b" />

Include a new directory within [your_project_folder] called ‘raw’ and drop in all your raw exported .csv qPCR files. Raw files must be in .csv format and names much match the names in the plate maps exactly.

<img width="489" height="324" alt="image" src="https://github.com/user-attachments/assets/d7c73606-d930-4ef5-a781-9d3e51cfe028" />
Figure 1: Instructions for preparing plate_maps.csv file.

<img width="409" height="153" alt="image" src="https://github.com/user-attachments/assets/93c6dd6f-e609-4e50-a7b6-6992cfd2339e" />
Figure 2: Instructions for preparing the target_mapping.csv file.


3.	Run the script: double click on the shortcut


4.	Output files
Files saved into [your_project_name] directory:
- your_project_name _qpcr_combined_results.csv. Contains all qpcr results collated together according to plate and target maps)
- your_project_name _qpcr_clean_summary.xlsx. Contains all cleaned qPCR data that passed QC. Results are collapsed by Sample name to give one row per Sample, per target. If enabled, correction to an internal control is completed at this stage. All inconclusive results flagged and excluded from this sheet.
- your_project_name _qpcr_exclusions.xlsx. Contains all samples which were marked as ‘Inconclusive’ and excluded, as well as the POS and NEG controls, as these are not required for analysis.
