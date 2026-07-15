# Principal Component Analysis (PCA) Pipeline
# Purpose:
#   Explore global immune phenotype variation between rural Schistosoma mansoni-infected (Sm+), rural Schistosoma mansoni-uninfected (Sm−), and urban participants.
# Analyses:
#   1. PCA
#   2. Scree plots
#   3. PCA visualisation
#   4. PERMANOVA
#   5. Principal curve analysis
#   6. Lambda trajectory analysis
#
# Associated manuscript: BCG vaccination attenuates Schistosoma mansoni-associated rural–urban immune variation 
# Author: Gyaviira Nkurunungi

#Load the required packages
library(stats)
library(ggplot2)
library(ggfortify)
library(readxl)
library(writexl)
library(vegan)
library(tidyverse)
library(here)
library(glue)
library(patchwork)
library(Hmisc)
library(nlme)
library(emmeans)
library(PMCMRplus)
library(princurve)
library(ggpubr)
library(ggbeeswarm)
library(ggnewscale)

# Input and output locations
input_file <- "masterdata_cd45.xlsx"
input_sheet <- "week0"
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)
stopifnot(file.exists(input_file))

# Import data
# Columns 11:135 contain immune phenotype frequencies used for PCA.
data <- read_excel(input_file, sheet = input_sheet)
dfdata <- as.data.frame(data)

#PCA
pca <- prcomp(dfdata[11:135], center = TRUE, scale = TRUE)
summary(pca)

#eigenvalues and eigenvectors:
(eigenvalues <- pca$sdev^2)
(eigenvectors <- pca$rotation)

#scree plots
screeplot(pca, npcs=10, type= "barplot", main = "Scree Plot")
screeplot(pca, npcs=10, type= "lines", main = "Scree Plot")

#PCA visualisation: Rural vs Urban comparison:
plot <- autoplot(
  pca,
  data = dfdata,
  colour = "setting",
  frame = TRUE,
  frame.type = "norm",  
  frame.level = 0.95,   
  size = 5             
) +
  scale_color_manual(values = c("darkgreen", "magenta")) +  
  scale_fill_manual(values = c("transparent", "transparent")) +  
  theme_classic() +     
  theme(
    axis.title.x = element_text(color="black", size=25),
    axis.text.x = element_text(color="black", size=18),
    axis.title.y = element_text(color="black", size=25), 
    axis.text.y = element_text(color="black", size=18), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 30)
  ) 

print(plot)
ggsave(file.path(output_dir, "allRuralvsallUrban_pca_exvivo.png"), plot = plot, width = 8, height = 5, dpi = 600)

##PERMANOVA: rural vs urban comparison
pca_scores <- as.data.frame(pca$x[, 1:2])
pca_scores$setting <- dfdata$setting
dist_matrix <- dist(pca_scores, method = "euclidean")
set.seed(1234)
permanova_result <- adonis2(dist_matrix ~ setting, data = pca_scores, permutations = 1000)
print(permanova_result)

####PCA visualisation: 3 group comparison
# arm3: rural_sm_inf   = Rural Sm+ ; rural_sm_uninf = Rural Sm− ; # urban = Urban

plot <- autoplot(
  pca,
  data = dfdata,
  colour = "arm3",
  frame = TRUE,
  frame.type = "norm", 
  frame.level = 0.95,   
  size = 5              
) +
  scale_color_manual(values = c("darkgreen", "blue", "magenta")) +  
  scale_fill_manual(values = c("transparent", "transparent", "transparent")) +  
  theme_classic() +    
  theme(
    axis.title.x = element_text(color="black", size=25),
    axis.text.x = element_text(color="black", size=18),
    axis.title.y = element_text(color="black", size=25), 
    axis.text.y = element_text(color="black", size=18),  # y-axis text size
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 30)
  ) 

plot$layers[[2]]$aes_params$linewidth <- 0.2 
plot$layers[[2]]$aes_params$linetype <- "dashed" 

print(plot)
ggsave(file.path(output_dir, "allRuralvsallUrban_pca_exvivo.png"), plot = plot, width = 8, height = 5, dpi = 600)


#PERMANOVA: Three-group comparison
pca_scores <- as.data.frame(pca$x[, 1:2])
pca_scores$arm3 <- dfdata$arm3
dist_matrix <- dist(pca_scores, method = "euclidean")
set.seed(1234)
permanova_result <- adonis2(dist_matrix ~ arm3, data = pca_scores, permutations = 1000)
print(permanova_result)

# Principal curve analysis ################

dfdata$pc1 <- pca$x[,1] # get the PC coordinate 1
dfdata$pc2 <- pca$x[,2] # now for coordinate 2

set.seed(45)
pcurve <- principal_curve(as.matrix(dfdata[,c("pc1", "pc2")]), maxit=100, thresh = 0.0005)
dfdata$lambda <- pcurve$lambda

#reorder curve coordinates for plotting
lf.reord <- dfdata
ix <- order(lf.reord$lambda)
pcurve.coords <- pcurve$s[ix,]
lambda <- lf.reord$lambda[ix]

dfdata$arm3 <- factor(dfdata$arm3, levels = c("rural_sm_inf", "rural_sm_uninf", "urban"))
color_groups <- c(
  "rural_sm_inf"   = "darkgreen", 
  "rural_sm_uninf"   = "blue",
  "urban" = "magenta" 
)

p2 <- lf.reord |>
  ggplot(aes(pc1, pc2, color=arm3)) +
  geom_point(cex=5, aes(shape=arm3, color=arm3, fill=arm3), show.legend=FALSE) +
  geom_segment(x=lf.reord$pc1[ix], xend=pcurve.coords[,1],
               y=lf.reord$pc2[ix], yend=pcurve.coords[,2],
               inherit.aes=FALSE, color="gray80", show.legend=FALSE) +
  scale_fill_manual(values=color_groups) +
  scale_color_manual(values=color_groups) +
  scale_shape_manual(values=c(21,22,23)) +
  
  new_scale_color() +
  geom_path(data=pcurve.coords, aes(pc1, pc2, color=lambda), linewidth=3, inherit.aes=FALSE) +
  scale_color_viridis_c(name="Lambda", limits=range(lambda, na.rm=TRUE), breaks=range(lambda, na.rm=TRUE), labels=c("Min", "Max")) +
  tidydr::theme_dr() +
  theme(panel.grid = element_blank(),
        legend.position = c("inside"),
        legend.position.inside = c(0.3,0.3),
        legend.justification = c(1,1),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25)) +
  labs(x="PC1", y="PC2") +
  guides(color=guide_colorbar(direction="horizontal", legend.ticks=element_blank(), title.position="top")) +
  
  coord_equal()

p2
ggsave(file.path(output_dir, "plot_pcurve_PCA_WK0.png"), p2, width = 8, height = 7, dpi = 600)


fit <- gls(rank(lambda) ~ arm3, data=dfdata)
emm.trend <- emmeans(fit, poly ~ arm3, adjust="none")
(trend.pval <- as.data.frame(emm.trend$contrasts)[1,"p.value"]) # lambda trend-test P-value

emm.pairs <- emmeans(fit, pairwise ~ arm3)
pvals <- as.data.frame(emm.pairs$contrasts) |>
  separate(contrast, " - ", into=c("xmin", "xmax")) |>
  mutate(y.position = max(dfdata$lambda)) |>
  mutate(p.value = rstatix::p_round(p.value))

fit <- gls(lambda ~ arm3, data=dfdata)
means <- data.frame(emmeans(fit, pairwise ~ arm3)$emmeans)

p3 <- dfdata |>
  ggplot(aes(arm3, lambda)) +
  geom_quasirandom(aes(color=arm3), size=4, shape=1, stroke=1.2) +
  geom_errorbar(data=means, aes(x=arm3, ymin=lower.CL, ymax=upper.CL), inherit.aes=FALSE, width=.5, linewidth=1.5) +
  geom_point(data=means, aes(x=arm3, y=emmean, shape=arm3, fill=arm3), inherit.aes=FALSE, cex=5, stroke=1.5) +
  ggpubr::geom_bracket(data=pvals, step.increase=.08, tip.length=0.01, label.size = 6) +
  labs(x = NULL, y = "Lambda", subtitle = paste0("Trend p = ", signif(trend.pval, 3))) +
  
  # coord_cartesian(clip="off") +
  scale_fill_manual(values=color_groups) +
  scale_color_manual(values=color_groups) +
  scale_shape_manual(values=c(21,22,23)) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 25),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), 
        axis.text.x = element_blank(),
        plot.subtitle = element_text(size = 18), 
        panel.background = element_blank(),
        axis.line = element_line(linewidth = 0.5, colour = "black"))

p3
ggsave(file.path(output_dir, "plot_lambda_plot_PCA.png"), p3, width = 4, height = 8, dpi = 600)
