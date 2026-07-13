CyTOF Debarcoding Pipeline
This repository contains the R code used to debarcode mass cytometry (CyTOF) data generated in:
"BCG vaccination attenuates Schistosoma mansoni-associated rural–urban immune variation" 

Overview
PBMC samples were multiplexed using a six-choose-three cadmium barcode scheme (Cd106, Cd110, Cd111, Cd112, Cd114 and Cd116), acquired on a Helios™ mass cytometer, normalized using the Fluidigm bead-normalization workflow, and exported as concatenated FCS files.
This script performs debarcoding using the CATALYST R package by:
1. Importing concatenated FCS files
2. Reading barcode definitions from an Excel barcode key
3. Assigning events to barcode populations
4. Estimating barcode separation cutoffs
5. Applying barcode-specific cutoffs
6. Exporting debarcoded FCS files for downstream analysis
Important: This pipeline assumes that all input FCS files have already undergone manual preprocessing in FlowJo™. EQ beads, non-Gaussian events, dead cells, doublets, and CD45-negative events were removed prior to export. Consequently, the input files contain only live, singlet, CD45+ cells and are provided as concatenated batch-level FCS files for debarcoding.
Requirements
The script was developed using:
• R (v4.5.1)
• CATALYST (v1.32.1)
• flowCore
• readr
• openxlsx

Input Files
Concatenated FCS file
Example:
batch1_panel1_cd45pos_concat_1.fcs

Barcode key
Example:
batch1_BarcodeOverview.xlsx

Rows correspond to samples and columns correspond to barcode channels. Barcode assignments are encoded using:
x = barcode present
blank = barcode absent

Once the script is ran, it will generate quality-control summaries, estimate barcode separation thresholds, and export individual debarcoded FCS files.

Output files are written as:
sample_1.fcs
sample_2.fcs
sample_3.fcs
...

Each file contains events assigned to a single barcode population following application of barcode-specific separation cutoffs.

Quality Control
The script provides:
table(re0$bc_id)

to summarize barcode assignments and to visualize barcode separation and debarcoding yield.
plotYields(re)
Users should inspect these outputs before downstream analyses.
