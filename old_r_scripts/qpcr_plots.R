############################################
## 16S qPCR: bacterial load on skin swabs
## Gnotobiotic Xenopus project
############################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# ---- paths ----
dir <- "~/Documents/GitHub/Gnotobiotic_Xenopus/"
qpcr_file <- file.path(dir, "qpcr_data.txt")
map_file  <- file.path(dir, "molly_xenopus_map.txt")
out_dir   <- dir

# ---- data ----
qpcr <- read.delim(qpcr_file)
meta <- read.delim(map_file)

# Keep sample wells only (NTC rows have Content == "NTC" and Cq == NA)
ntc  <- filter(qpcr, Content == "NTC")
qpcr <- filter(qpcr, Content == "Unkn", SampleID %in% meta$SampleID)
cat("NTC wells:", nrow(ntc), "| amplified:", sum(!is.na(ntc$Cq)), "\n")

# Treatment groups (antibiotic x Fusarium)
qpcr <- qpcr %>%
  select(-Antibiotic) %>%
  left_join(select(meta, SampleID, Antibiotic, Fusarium, FusariumDose), by = "SampleID") %>%
  mutate(Group = case_when(
           Fusarium == "Yes" & Antibiotic == "AllAB" ~ "AllAB + Fk",
           Fusarium == "Yes"                         ~ "Fk alone",
           Antibiotic == "Gentamycin"                ~ "Gentamicin",
           TRUE                                      ~ Antibiotic),
         Group = factor(Group, levels = c("Control", "Ampicillin", "AllAB", "PSN",
                                          "Sulfamethazine", "Gentamicin", "Trimethoprim",
                                          "Fk alone", "AllAB + Fk")))

# ---- technical replicate agreement ----
reps <- qpcr %>%
  group_by(SampleID) %>% arrange(Plate, Well, .by_group = TRUE) %>%
  mutate(rep = paste0("rep", row_number())) %>% ungroup() %>%
  select(SampleID, rep, LogCopy, Plate) %>%
  pivot_wider(names_from = rep, values_from = c(LogCopy, Plate)) %>%
  mutate(layout = ifelse(Plate_rep1 == Plate_rep2, "same plate", "split across plates"))

fit <- anova(lm(LogCopy ~ factor(SampleID), data = qpcr))
icc <- (fit$`Mean Sq`[1] - fit$`Mean Sq`[2]) / (fit$`Mean Sq`[1] + fit$`Mean Sq`[2])
mad <- mean(abs(reps$LogCopy_rep1 - reps$LogCopy_rep2))
r_by <- reps %>% group_by(layout) %>%
  summarise(n = n(), r = cor(LogCopy_rep1, LogCopy_rep2)) %>%
  mutate(lab = sprintf("%s (n = %d, r = %.2f)", layout, n, r))

ggplot(reps, aes(LogCopy_rep1, LogCopy_rep2, colour = layout)) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(size = 2) +
  coord_equal() +
  scale_colour_manual(values = c("same plate" = "#332288", "split across plates" = "#CC6677"),
                      labels = setNames(r_by$lab, r_by$layout), name = NULL) +
  labs(x = "Replicate 1 (log10 16S copies)", y = "Replicate 2 (log10 16S copies)",
       title = "Technical replicates",
       subtitle = sprintf("ICC = %.2f; mean |difference| = %.2f log", icc, mad)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "vertical")

# ---- per-sample load (mean of technical replicates) ----
load <- qpcr %>%
  group_by(SampleID, Group, Antibiotic, Fusarium, FusariumDose) %>%
  summarise(LogCopy = mean(LogCopy), .groups = "drop")
write.csv(load, file.path(out_dir, "qpcr_per_sample.csv"), row.names = FALSE)

# ---- statistics ----
ab_only <- filter(load, Fusarium == "No")
kw <- kruskal.test(LogCopy ~ Group, data = ab_only)

contrast <- function(g1, g2) {
  x <- load$LogCopy[load$Group == g1]; y <- load$LogCopy[load$Group == g2]
  data.frame(contrast = paste(g1, "vs", g2), n1 = length(x), n2 = length(y),
             delta_log10 = round(median(x) - median(y), 3),
             fold = round(10^(median(x) - median(y)), 2),
             p = signif(wilcox.test(x, y, exact = FALSE)$p.value, 3))
}
stats <- bind_rows(
  lapply(setdiff(levels(droplevels(ab_only$Group)), "Control"), contrast, g2 = "Control"),
  contrast("Fk alone", "Control"),
  contrast("AllAB + Fk", "AllAB"),
  contrast("AllAB + Fk", "Fk alone"))
fk <- filter(load, Group == "Fk alone")
dose <- cor.test(log10(as.numeric(fk$FusariumDose)), fk$LogCopy, method = "spearman", exact = FALSE)
print(kw); print(stats)

#Kruskal-Wallis chi-squared = 7.9254, df = 6, p-value = 0.2436

#contrast n1 n2 delta_log10 fold      p
#1     Ampicillin vs Control  4  5      -0.129 0.74 0.2700
#2          AllAB vs Control  3  5       0.220 1.66 0.3710
#3            PSN vs Control  3  5       0.023 1.05 1.0000
#4 Sulfamethazine vs Control  3  5      -0.028 0.94 0.5510
#5     Gentamicin vs Control  5  5      -0.223 0.60 0.2100
#6   Trimethoprim vs Control  8  5      -0.043 0.90 0.5100
#7       Fk alone vs Control 17  5      -0.373 0.42 0.0344
#8       AllAB + Fk vs AllAB  5  3      -0.752 0.18 0.0369
#9    AllAB + Fk vs Fk alone  5 17      -0.159 0.69 0.3470

cat(sprintf("Fk-alone load vs log10 dose: rho = %.2f, p = %.2f\n", dose$estimate, dose$p.value))
#Fk-alone load vs log10 dose: rho = -0.03, p = 0.90

write.csv(stats, file.path(out_dir, "qpcr_stats.csv"), row.names = FALSE)

# ---- load by treatment group ----
ctrl_med <- median(load$LogCopy[load$Group == "Control"])
cmp <- list(c("Control", "Fk alone"), c("AllAB", "AllAB + Fk"))

p_load <- ggplot(load, aes(Group, LogCopy)) +
  geom_hline(yintercept = ctrl_med, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey40") +
  geom_jitter(width = 0.12, size = 1.8, colour = "#332288", seed = 1) +
  stat_compare_means(comparisons = cmp, method = "wilcox.test",
                     method.args = list(exact = FALSE), label = "p.format", size = 3.2) +
  labs(x = NULL, y = "log10 16S copies per reaction\n(mean of technical replicates)",
       title = "Bacterial load by treatment",
       subtitle = sprintf("Antibiotics vs Control (no Fusarium): Kruskal-Wallis p = %.2f", kw$p.value)) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# ---- load vs Fusarium dose ----
dose_df <- filter(load, Fusarium == "Yes") %>%
  mutate(Dose = as.numeric(FusariumDose),
         Pretreatment = ifelse(Antibiotic == "AllAB", "AllAB", "None"))

p_dose <- ggplot(dose_df, aes(Dose, LogCopy, colour = Pretreatment)) +
  geom_hline(yintercept = ctrl_med, linetype = 2, colour = "grey50") +
  geom_jitter(width = 0.05, height = 0, size = 2.2, seed = 1) +
  scale_x_log10(breaks = c(1e2, 1e3, 1e4, 1e5),
                labels = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5))) +
  scale_colour_manual(values = c(None = "#88CCEE", AllAB = "#332288")) +
  labs(x = "Fusarium dose (spores/mL)", y = "log10 16S copies per reaction",
       title = "Load vs Fusarium dose",
       subtitle = sprintf("Fk alone: Spearman rho = %.2f, p = %.2f\nDashed line = Control median",
                          dose$estimate, dose$p.value)) +
  theme_bw(base_size = 11)

