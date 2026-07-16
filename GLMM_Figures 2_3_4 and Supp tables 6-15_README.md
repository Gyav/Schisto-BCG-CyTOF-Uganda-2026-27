**Binomial GLMM Analysis of Immune Cell Cluster Frequencies**

This repository contains R code used to analyse immune cell cluster frequencies in the study:


**BCG vaccination attenuates Schistosoma mansoni-associated rural–urban immune variation.**

**Overview**

Immune cell cluster frequencies (% of total CD45+ PBMCs) were analysed using binomial generalized linear mixed models (GLMMs).

**Models included:**

• Fixed effects: study group, age, and sex

• Random effect: participant ID

• Link function: logit

• Optimizer: BOBYQA


**Analyses Performed**

1. Rural vs Urban comparisons
   
    • Differences in immune cell cluster frequencies between all rural and urban participants.
   
2. Linear trend tests
   
    • Assessment of ordered immune variation across:
   
    • Rural Sm+ → Rural Sm− → Urban
   
    • Performed using orthogonal polynomial contrasts via emmeans.
   
3. Pairwise comparisons
   
    • Rural Sm− vs Rural Sm+
   
    • Urban vs Rural Sm+
   
    • Urban vs Rural Sm−
   
4. Multiple-testing correction
   
    • Benjamini–Hochberg false discovery rate (FDR).
   
