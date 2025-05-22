# some_stat
some_stat

# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)

# Load the data
df <- read_csv("synthetic_rhabdomyosarcoma_data.csv")

# Check the structure
glimpse(df)

# Step 1: Calculate per-patient statistics
patient_stats <- df %>%
  group_by(patient_id, subtype) %>%
  summarise(
    total_cells = n(),
    ki67_pos = sum(Ki67 == 1),
    lactb_pos = sum(LACTB == 1),
    lactb_in_ki67_pos = sum(LACTB == 1 & Ki67 == 1),
    lactb_in_ki67_neg = sum(LACTB == 1 & Ki67 == 0)
  ) %>%
  mutate(
    pct_ki67_pos = ki67_pos / total_cells * 100,
    pct_lactb_pos = lactb_pos / total_cells * 100,
    pct_lactb_in_ki67_pos = lactb_in_ki67_pos / ki67_pos * 100,
    pct_lactb_in_ki67_neg = lactb_in_ki67_neg / (total_cells - ki67_pos) * 100
  )

# Step 2: Boxplots per tumor subtype
ggplot(patient_stats, aes(x = subtype, y = pct_ki67_pos)) +
  geom_boxplot(fill = "skyblue") +
  labs(title = "% of Ki67+ cells by Tumor Subtype", y = "% Ki67+", x = "Tumor Subtype")

ggplot(patient_stats, aes(x = subtype, y = pct_lactb_in_ki67_pos)) +
  geom_boxplot(fill = "salmon") +
  labs(title = "% of LACTB+ in Ki67+ cells", y = "% LACTB+ | Ki67+", x = "Tumor Subtype")

ggplot(patient_stats, aes(x = subtype, y = pct_lactb_in_ki67_neg)) +
  geom_boxplot(fill = "palegreen") +
  labs(title = "% of LACTB+ in Ki67− cells", y = "% LACTB+ | Ki67−", x = "Tumor Subtype")

# Step 3: Kruskal-Wallis test (non-parametric ANOVA)
kruskal.test(pct_lactb_in_ki67_pos ~ subtype, data = patient_stats)
kruskal.test(pct_lactb_in_ki67_neg ~ subtype, data = patient_stats)
kruskal.test(pct_ki67_pos ~ subtype, data = patient_stats)
