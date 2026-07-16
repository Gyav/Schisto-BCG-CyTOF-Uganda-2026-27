# ---------------------------------------------------------
# Binomial GLMM Analysis of Immune Cell Cluster Frequencies
#
# Purpose:
#   Identify differences in immune cell cluster frequencies between study groups using binomial generalized linear mixed models (GLMMs).
#
# Analyses:
#   1. Rural versus urban comparisons
#   2. Linear trend tests across: Rural Sm+ → Rural Sm− → Urban
#   3. Pairwise group comparisons
#   4. FDR correction across clusters
#
# Associated manuscript: BCG vaccination attenuates Schistosoma mansoni-associated rural–urban immune variation
# Author: Gyaviira Nkurunungi

# ---------------------------------------------------------
# Load required packages

library(ggplot2)
library(ggbreak)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(dplyr)
library(readxl)
library(writexl)
library(emmeans)

getwd()
dir<- setwd("~/Documents/POPVAC/CyToF and Aurora/CYTOF ANALYSIS 2021_22/CYTOVAC OMIQ ANALYSES/vaccine responses/ex vivo/glmm results")

# ---------------------------------------------------------
# Input and output files
input_file <- "masterdata_cd45.xlsx"
input_sheet <- "week0"
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)
stopifnot(file.exists(input_file))

# ---------------------------------------------------------
# Import data
meta <- read_excel(input_file, sheet = input_sheet)
dfdata <- data.frame(meta)

######ASSOCIATIONS BETWEEN WEEK 0 IMMUNE PHENOTYPE AND SETTING [all rural vs all urban] #########

#Example glmm model using cluster 7
dfdata$cluster7_count <- round((dfdata$cluster7 / 100)*dfdata$CD45pos) 
glmm <- glmer(cbind(cluster7_count, CD45pos - cluster7_count)  ~ setting + age + factor(sex)+ (1 | file),
              data = dfdata,
              family = binomial(link = "logit"),
              control = glmerControl(optimizer = "bobyqa"))

summary(glmm)
results <- tidy(glmm, effects = "fixed") %>%
  mutate(p.value = 2 * (1 - pnorm(abs(statistic)))) # Calculate p-value from z-statistic
print(results)
(pvalue<-results$p.value[2])
(logOR<-results$estimate[2])

# loop over clusters 1 to 125. 
all_results <- data.frame()
for (i in 1:125) {
  cluster_var <- paste0("cluster", i)           
  count_var   <- paste0(cluster_var, "_count")  
  dfdata[[count_var]] <- round(dfdata[[cluster_var]] / 100 * dfdata$CD45pos)
  response_part <- paste0("cbind(", count_var, ", CD45pos - ", count_var, ")")
  form <- as.formula(paste0(response_part, " ~ setting + age + factor(sex) + (1 | file)"))
 
  glmm <- glmer(form,
                data = dfdata,
                family = binomial(link = "logit"),
                control = glmerControl(optimizer = "bobyqa"))
  
  results <- tidy(glmm, effects = "fixed") %>%
    mutate(p.value = 2 * (1 - pnorm(abs(statistic))))  
  
  setting_row <- results %>% filter(grepl("^setting", term)) 
  logOR  <- setting_row$estimate[1]
  p.value <- setting_row$p.value[1]
  
  all_results <- rbind(all_results, data.frame(cluster = cluster_var, logOR = logOR, p.value = p.value))
  
  print(paste("Cluster:", i, "logOR:", logOR, "P-value:", p.value))
}

print(all_results)

#fdr calculation 
(pvalues <- c(all_results$p.value))
FDRpval<- p.adjust(pvalues,method="BH")
(all_results2<- cbind(all_results,FDRpval))
write_xlsx(all_results2, file.path(output_dir, "ex_vivo_wk4_setting.xlsx"))


######trend test p values using GLMMs#########


##run the glmm model, using one cluster to check the code
dfdata$cluster7_count <- round((dfdata$cluster7 / 100)*dfdata$CD45pos) 
dfdata$arm3 <- factor(dfdata$arm3, ordered = TRUE) 
glmm <- glmer(cbind(cluster7_count, CD45pos - cluster7_count)  ~ arm3 + age + factor(sex) +  (1 | file),
              data = dfdata,
              family = binomial(link = "logit"),
              control = glmerControl(optimizer = "bobyqa"))

summary(glmm)
results <- tidy(glmm, effects = "fixed") %>%
  mutate(p.value = 2 * (1 - pnorm(abs(statistic))))
print(results)

# Get estimated marginal means
emm <- emmeans(glmm, ~ arm3)
# Test for linear trend (orthogonal polynomial contrast)
trend_test <- contrast(emm, method = "poly")  
print(trend_test)
confint(trend_test) #confidence intervals
# Extract just the linear trend p-value
p_trend <- summary(trend_test)$p.value[1]
cat("P for linear trend:", p_trend, "\n")

# Define pairwise contrasts with lower group as reference
pairwise_ref <- contrast(emm,
                         method = list(
                           "rural_sm_uninf vs rural_sm_inf" = c(-1, 1, 0),
                           "urban vs rural_sm_inf" = c(-1, 0, 1),
                           "urban vs rural_sm_uninf" = c(0, -1, 1)
                         ),
                         adjust = "none")

# Summarize pairwise contrasts
(pairwise_summary <- summary(pairwise_ref))
confint(pairwise_ref)

#####Now for all 125 clusters#################

all_results <- data.frame()
# Loop through clusters 1 to 125
for (i in 1:125) {
  cluster_var <- paste0("cluster", i)           
  count_var   <- paste0(cluster_var, "_count")  
  dfdata[[count_var]] <- round(dfdata[[cluster_var]] / 100 * dfdata$CD45pos)
  response_part <- paste0("cbind(", count_var, ", CD45pos - ", count_var, ")")
  form <- as.formula(paste0(response_part, " ~ arm3 + age + factor(sex) + (1 | file)"))
  
  glmm <- glmer(form,
                data = dfdata,
                family = binomial(link = "logit"),
                control = glmerControl(optimizer = "bobyqa"))
 
  emm        <- emmeans(glmm, ~ arm3)
  trend_test <- contrast(emm, method = "poly", max.degree = 1)  
  trend_sum  <- summary(trend_test, infer = TRUE)             
  
  logOR     <- trend_sum$estimate[1]
  ci_lower  <- trend_sum$asymp.LCL[1]
  ci_upper  <- trend_sum$asymp.UCL[1]
  p_trend   <- trend_sum$p.value[1]
  
  # Store results
  all_results <- rbind(
    all_results,
    data.frame(
      cluster   = cluster_var,
      logOR     = logOR,           
      OR        = exp(logOR),      
      ci_lower  = ci_lower,
      ci_upper  = ci_upper,
      p_trend   = p_trend,
      converged = !any(grepl("failed to converge|singular", 
                             summary(glmm)$optinfo$message, ignore.case = TRUE))
    )
  )
  
  cat(sprintf(
    "Cluster %3d | logOR: %6.3f  [%5.3f, %5.3f]  OR: %5.2f  p-trend: %8.5f\n",
    i, logOR, ci_lower, ci_upper, exp(logOR), p_trend
  ))
}

print(all_results)

#fdr calculation 
(ptrends <- c(all_results$p_trend))
FDRpval_trend<- p.adjust(ptrends,method="BH")
(all_results2<- cbind(all_results,FDRpval_trend))
write_xlsx(all_results2, file.path(output_dir, "ptrends_age_sex.xlsx"))


############Now loop for pairwise comparisons###############

# Initialize empty data frames for the three pairwise tables
table1 <- data.frame(Cluster = character(), estimate = numeric(), p.value = numeric())
table2 <- data.frame(Cluster = character(), estimate = numeric(), p.value = numeric())
table3 <- data.frame(Cluster = character(), estimate = numeric(), p.value = numeric())

# Loop through all clusters
for (i in 1:125) {
  cluster_name <- paste0("cluster", i)           
  count_var   <- paste0(cluster_name, "_count")  
  dfdata[[count_var]] <- round(dfdata[[cluster_name]] / 100 * dfdata$CD45pos)
  response_part <- paste0("cbind(", count_var, ", CD45pos - ", count_var, ")")
  form <- as.formula(paste0(response_part, " ~ arm3 + age + factor(sex) + (1 | file)"))
 
  glmm <- glmer(form,
                data = dfdata,
                family = binomial(link = "logit"),
                control = glmerControl(optimizer = "bobyqa"))
  
  emm <- emmeans(glmm, ~ arm3)
  # Define pairwise contrasts with lower group as reference
  pairwise_ref <- contrast(emm,
                           method = list(
                             "rural_sm_uninf vs rural_sm_inf" = c(-1, 1, 0),
                             "urban vs rural_sm_inf" = c(-1, 0, 1),
                             "urban vs rural_sm_uninf" = c(0, -1, 1)
                           ),
                           adjust = "none")
  
  # Summarize pairwise contrasts
  pairwise_summary <- summary(pairwise_ref)
  ci <-confint(pairwise_ref)
  logOR <-pairwise_summary$estimate[1]
  
  # Append to the appropriate tables
  table1 <- rbind(table1, data.frame(Cluster = cluster_name,
                                     logOR = pairwise_summary$estimate[1],
                                     ci_lower = summary(ci)$asymp.LCL[1],
                                     ci_upper = summary(ci)$asymp.UCL[1],
                                     p.value = pairwise_summary$p.value[1]))
  
  table2 <- rbind(table2, data.frame(Cluster = cluster_name,
                                     logOR = pairwise_summary$estimate[2],
                                     ci_lower = summary(ci)$asymp.LCL[2],
                                     ci_upper = summary(ci)$asymp.UCL[2],
                                     p.value = pairwise_summary$p.value[2]))
  
  table3 <- rbind(table3, data.frame(Cluster = cluster_name,
                                     logOR = pairwise_summary$estimate[3],
                                     ci_lower = summary(ci)$asymp.LCL[3],
                                     ci_upper = summary(ci)$asymp.UCL[3],
                                     p.value = pairwise_summary$p.value[3]))
}

# View results
head(table1)
head(table2)
head(table3)

print(table1)

#fdr calculation 
(pvals <- c(table1$p.value))
FDRpval<- p.adjust(pvals, method="BH")
(table1b<- cbind(table1,FDRpval))
write_xlsx(all_results2, file.path(output_dir, "ruralSMneg_vs_ruralSMpos_age_sex.xlsx"))


print(table2)
#fdr calculation using the p.adjust() command
(pvals <- c(table2$p.value))
FDRpval<- p.adjust(pvals, method="BH")
(table2b<- cbind(table2,FDRpval))
write_xlsx(all_results2, file.path(output_dir, "ruralSMpos_vs_URB_age_sex.xlsx"))


# Display full results
print(table3)
#fdr calculation using the p.adjust() command
(pvals <- c(table3$p.value))
FDRpval<- p.adjust(pvals, method="BH")
(table3b<- cbind(table3,FDRpval))
write_xlsx(all_results2, file.path(output_dir, "ruralSMneg_vs_URB_age_sex.xlsx"))



##################box plots##################
plot <- ggplot(dfdata, aes(y=cluster7, x=setting)) +
  geom_boxplot(width = 0.5, aes(color = setting), outlier.shape = NA) +
  geom_jitter(aes(color = setting), alpha = 0.5, size = 2.5, position = position_jitterdodge(jitter.width = 0.3)) +
  scale_color_manual(values = c("rural" = "blue", "urban" = "magenta")) +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 0, face = "bold"), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 15),  
        axis.title.y = element_text(size = 23, face = "bold"), 
        panel.background = element_rect(fill = "transparent"),
        legend.position = "none")+ 
  scale_y_continuous(limits = c(0, 13)) +  
  labs(x = "",
       y = "CD56+CD16+ (% of total CD45+)")

print(plot)
# Save plot
ggsave("cd56cd16_setting.png", plot, width = 4, height = 7, dpi = 600)

#split into 3 groups
dfdata$arm <- factor(dfdata$arm, levels = c("STD", "INT", "URB")) # Reorder the factor levels of 'arm'
plot <- ggplot(dfdata, aes(y=cluster7, x=arm)) +
  geom_boxplot(width = 0.5, aes(color = arm), outlier.shape = NA) +
  geom_jitter(aes(color = arm), alpha = 0.5, size = 2.5, position = position_jitterdodge(jitter.width = 0.9)) +
  scale_color_manual(
    values = c("INT" = "blue", "STD" = "darkgreen", "URB" = "magenta"),
    labels = c("STD" = "Rural, Sm+", "INT" = "Rural, Sm-", "URB" = "Urban")  
  ) +
  scale_x_discrete(labels = c(
    "STD" = "Rural, Sm+",
    "INT" = "Rural, Sm-",
    "URB" = "Urban"
  )) +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.line.y.right = element_blank(),  
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(), 
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 0, face = "bold"),  
        axis.ticks.x = element_blank(),  
        axis.text.y = element_text(size = 15),  
        axis.title.y = element_text(size = 23, face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        legend.position = "none")+ 
  labs(x = "",
       y = "CD56+CD16+ (% of total CD45+)")

print(plot)
# Save plot 
ggsave("cd56cd16_arm.png", plot, width = 4, height = 7, dpi = 600)

