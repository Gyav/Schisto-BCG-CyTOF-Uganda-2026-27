# CyTOF Debarcoding Pipeline
# Purpose: Debarcode concatenated CyTOF FCS files using CATALYST and export individual FCS files for each barcode.
# Inputs:
#   - batch1_panel1_cd45pos_concat_1.fcs
#   - batch1_BarcodeOverview.xlsx
# Outputs:
#   - sample_<barcodeID>.fcs
# Associated manuscript: BCG vaccination attenuates Schistosoma mansoni-associated rural–urban immune variation 
# Author: Gyaviira Nkurunungi

# ---------------------------------------------------------
# Load required packages
library(CATALYST)
library(flowCore)
library(readr)

# Define input files
fcs_file <- "batch1_panel1_cd45pos_concat_1.fcs"
barcode_file <- "batch1_BarcodeOverview.xlsx"
# Basic input checks
stopifnot(file.exists(fcs_file))
stopifnot(file.exists(barcode_file))

# Read concatenated FCS file (Transformation is disabled because barcode intensities are used directly during debarcoding)
fcs <- read.FCS(
  fcs_file,
  transformation = FALSE,
  truncate_max_range = FALSE
)
fcs
pData(parameters(fcs)) # View marker and feature names

# ------------------------------------------------------------------
# Prepare barcode key. Barcode channels in the barcode key must correspond to the barcode channel indices in the FCS file.
key <- openxlsx::read.xlsx(
  barcode_file,
  rowNames = TRUE
)
key[is.na(key)] <- 0 # Replace missing values with 0
barcode_ids <- rownames(key)
# Convert "x" entries to binary barcode matrix
key <- apply(
  key,
  2,
  function(x) ifelse(x == "x", 1, 0)
)
key <- as.data.frame(key)
key[] <- lapply(key, as.numeric)
rownames(key) <- barcode_ids
key

# ------------------------------------------------------------------
# Debarcoding
sce <- prepData(fcs)
# Convert channel names to numeric isotope masses
rowData(sce)$channel_name <- readr::parse_number(
  rowData(sce)$channel_name
)
rowData(sce)
# Preliminary barcode assignment
re0 <- assignPrelim(
  x = sce,
  bc_key = key,
  verbose = TRUE
)
re0
table(re0$bc_id) # Number of events assigned to each barcode population

# ------------------------------------------------------------------
# Estimate and apply barcode separation cutoffs
re <- estCutoffs(x = re0)
metadata(re)$sep_cutoffs
# Apply population-specific cutoffs
re <- applyCutoffs(x = re)
# Visualise barcode separation and debarcoding yield
plotYields(re, which = 0)
metadata(re)

# ------------------------------------------------------------------
# Retain successfully assigned events
sce <- re[, re$bc_id != 0]
(fs <- sce2fcs(
  sce,
  split_by = "bc_id"
))
dim(sce)

# ------------------------------------------------------------------
# Export debarcoded FCS files
(ids <- fsApply(fs, identifier))
for (id in ids) {
  ff <- fs[[id]]
  out_file <- sprintf(
    "sample_%s.fcs",
    id
  )
  write.FCS(
    ff,
    out_file
  )
}