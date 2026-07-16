# Volcano Plot Pipeline
# Purpose: Visualise differential immune phenotype frequencies and highlight significantly enriched or depleted cell clusters.
# Input:
# Excel file containing columns with:
#     - cluster identifiers
#     - model estimates
#     - age- and sex-adjusted p-values
#     - FDR-adjusted p-values
#     - major cell lineage names
# Associated manuscript: BCG vaccination attenuates Schistosoma mansoni-associated rural–urban immune variation
# Author: Gyaviira Nkurunungi

# ---------------------------------------------------------
# Load required packages

library(ggplot2)
library(readxl)
library(writexl)
library(ggrepel)

# ---------------------------------------------------------
# Input and output files
input_file <- "B_ex vivo_wk0_setting_age_sex.xlsx"
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)
stopifnot(file.exists(input_file))
# ---------------------------------------------------------
# Import data
dataVolcano <- read_excel(input_file)
dataVolcano$estimate <- as.numeric(dataVolcano$estimate)
dataVolcano$estimateC <- round(dataVolcano$estimate, 2)
# ---------------------------------------------------------
# Define differential abundance status
dataVolcano$diffexpressed <- "NO"
dataVolcano$diffexpressed[dataVolcano$estimate > 0 & dataVolcano$FDRpval < 0.05] <- "UP"
dataVolcano$diffexpressed[dataVolcano$estimate < 0 & dataVolcano$FDRpval < 0.05] <- "DOWN"
# ---------------------------------------------------------
# Cell-lineage colour scheme
category_colors <- c(
  "B cells" = "red", "CD4+ T cells" = "#ebc713", "CD8+ T cells" = "turquoise", 
  "gd T cells" = "#3E7B27", "Myeloid cells" = "blue", "Innate lymphoid cells (NK, ILCs)" = "magenta",
  "NO" = "black" 
)
dataVolcano$color_group <- ifelse(dataVolcano$diffexpressed == "NO", "NO", as.character(dataVolcano$`Major cell lineages`))
# ---------------------------------------------------------
# Volcano plot
p <- ggplot(data=dataVolcano, aes(x=estimate, y=-log10(FDRpval), color=color_group)) + 
  geom_point(size=5) + 
  geom_vline(xintercept=c(0), linetype="dashed", color="black") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") + 
  scale_color_manual(values=category_colors) + 
  geom_text_repel(data=dataVolcano[dataVolcano$diffexpressed != "NO", ], 
                  aes(label=cluster2), size=7, max.overlaps=35, 
                  #segment.color = NA, 
                  color="black") + 
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    legend.position = "none", #remove legend if necessary
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(hjust = 0.5, size = 18, color = "black"),
    axis.title.x = element_text(size = 25),
    axis.text.y = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 25),
  ) +
  labs( 
    x="Model estimate (log odds ratio)", y="-Log10(FDR)", 
    color="Cell lineages")+
  ylim(0,8)
p

# ---------------------------------------------------------
# Export figures
ggsave(file.path(output_dir, "B_ex vivo_wk0_setting_age_sex.pdf"), p, width = 7, height = 8, dpi = 600)
ggsave(file.path(output_dir, "B_ex vivo_wk0_setting_age_sex.png"), p, width = 7, height = 8, dpi = 600)
