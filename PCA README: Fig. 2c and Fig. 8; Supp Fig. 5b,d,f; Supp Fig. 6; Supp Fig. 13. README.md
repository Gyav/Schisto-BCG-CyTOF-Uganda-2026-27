**PCA, PERMANOVA and Principal Curve Analysis**

This repository contains the R code used to perform principal component analysis (PCA), PERMANOVA testing, and principal curve analysis of PBMC immune phenotype frequencies in:

_BCG vaccination attenuates Schistosoma mansoni-associated rural–urban immune variation_

**Overview**

The script was used to explore global differences in immune phenotypes between:

• Rural S. mansoni infected (Sm+)

• Rural S. mansoni uninfected (Sm−)

• Urban participants
**
The workflow performs:**

1. Principal component analysis (PCA) using prcomp
   
2. Visualisation of PCA embeddings
   
3. PERMANOVA testing using Euclidean distances from PCA coordinates
   
4. Principal curve fitting across the PCA space
   
5. Assignment of lambda values representing position along the inferred trajectory
   
6. Statistical comparison of lambda values between groups using generalized least squares models and estimated marginal means
