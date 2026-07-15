CytoNorm Batch Normalization of CyTOF Data

This repository contains the R code used for batch correction of CyTOF data using the CytoNorm package.

Overview

Exported FCS files were initially analysed in FlowJo™ (v10). Events were gated to exclude EQ beads (140Ce), followed by gating on Gaussian-derived parameters (residual, offset, centre, and width), live cells (195Pt), singlets (191Ir), and finally CD45⁺ cells.

Batched live, singlet, CD45⁺ cells were subsequently debarcoded and compensated using the CATALYST package in R. Batch effects arising from samples stained and acquired at different times were corrected using CytoNorm, employing a common reference sample included in every batch.
