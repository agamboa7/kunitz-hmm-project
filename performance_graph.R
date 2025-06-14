#----LABORATORY OF BIOINFORMATICS----
#----KUNITZ DOMAIN PROJECT----
# Programmer: Andrea Gamboa
# Date: 14/06/2025
# R Version: 4.5.0

# This script reads performance metrics in .tsv format for a profile HMM
# applied to two datasets using two scoring modes: full sequence and best domain.
# It visualizes the effect of different e-value thresholds on MCC, comparing results
# across datasets and modes using line plots with facets.
# Script was made to analyze the performance data obtained from the Kunitz domain
# HMM developed for the Laboratory of Bioinformatics I project as part of the  
# Bioinformatics Masters program at the University of Bologna.

sessionInfo()

setwd("C:/Users/andyg/Desktop/Unibo/LAB1") # Change to your wd 
getwd()

# Install the required package
install.packages("ggplot2")
install.packages("tidyverse")

# Load required libraries
library(tidyverse)
library(readr)

# Step 1: Read data
df1 <- read_tsv("set_1_performance.tsv")
df2 <- read_tsv("set_2_performance.tsv")

# Step 2: Preprocess both dataframes
df1 <- df1 %>%
  mutate(
    mode = ifelse(fullseq, "fullseq", "domain"),
    threshold_numeric = as.numeric(threshold),
    minus_log10_threshold = -log10(threshold_numeric),
    dataset = "set_1"
  )

df2 <- df2 %>%
  mutate(
    mode = ifelse(fullseq, "fullseq", "domain"),
    threshold_numeric = as.numeric(threshold),
    minus_log10_threshold = -log10(threshold_numeric),
    dataset = "set_2"
  )

# Step 3: Combine datasets
df_all <- bind_rows(df1, df2)

# Step 4: Plot MCC vs -log10(E-value), colored by dataset, faceted by mode
ggplot(df_all, aes(x = minus_log10_threshold, y = mcc, color = dataset)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = c("set_1" = "purple", "set_2" = "hotpink")) +
  scale_x_continuous(
    breaks = 1:12,
    labels = paste0("1e-", 1:12),
    name = "-log10(e-value)"
  ) +
  ylab("MCC") +
  facet_wrap(~mode, ncol = 1) +
  labs(
    title = "MCC vs e-value Threshold",
    color = "Dataset:"
  ) +
  theme_minimal(base_size = 14)
