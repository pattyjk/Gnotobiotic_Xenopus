###################################
## XL Size (snout-vent length, SVL)
###################################
dir <- "~/Documents/GitHub/Gnotobiotic_Xenopus/"

library(ggplot2)
library(dplyr)
library(patchwork)
source(file.path(dir, "pub_theme.R"))   # theme_pub(), recode_abx(), fk_dose_factor(), palettes, save_pub()

# load data
meta <- read.delim(file.path(dir, "molly_xenopus_map.txt"), header = TRUE,
                   na.strings = c("NA", "N/A"))

# Drop PCR control (by name, Control == "Yes") and any rows missing SVL
meta <- meta %>%
  filter(Control == "No", !is.na(SVL_mm)) %>%
  mutate(
    SVL_mm       = as.numeric(SVL_mm),
    FusariumDose = as.numeric(FusariumDose),
    # BUG FIX: old paste0("Fk ", 100000) gave "Fk 1e+05", which was not a factor
    # level -> NA, silently dropping all 10^5 animals from the dose ANOVA.
    DoseLabel    = fk_dose_factor(FusariumDose),
    Antibiotic   = factor(recode_abx(Antibiotic), levels = c(abx_levels, "None")),
    # antibiotic pretreatment of Fusarium-exposed animals (FuskNoAB -> "None")
    Pretreatment = factor(ifelse(Antibiotic == "None", "None", as.character(Antibiotic)),
                          levels = c("None", "AllAB"))
  ) %>%
  droplevels()
stopifnot(!anyNA(meta$DoseLabel[meta$Fusarium == "Yes"]))

meta_noFusk <- meta %>% filter(Fusarium == "No") %>% droplevels()
meta_Fusk   <- meta %>% filter(Fusarium == "Yes") %>% droplevels()

#---------------------------------------
# Statistics
#---------------------------------------
# ANOVA: SVL ~ Antibiotic (no Fusarium)  [unchanged data: n = 31, 7 groups]
print(bartlett.test(SVL_mm ~ Antibiotic, data = meta_noFusk))
# Bartlett's K-squared = 4.3924, df = 6, p-value = 0.6237
aov_abx <- aov(SVL_mm ~ Antibiotic, data = meta_noFusk)
print(summary(aov_abx))
#            Df Sum Sq Mean Sq F value Pr(>F)
# Antibiotic  6  2.728  0.4546   0.701  0.651
# Residuals  24 15.560  0.6483
print(TukeyHSD(aov_abx))
kw_abx <- kruskal.test(SVL_mm ~ Antibiotic, data = meta_noFusk)
print(kw_abx)

# ANOVA: SVL ~ FusariumDose (Fusarium samples only)
# NOTE: the old comment (df = 2, 11 residual df, n = 14) excluded the eight 10^5
# animals (NA bug). With all 22 Fk-exposed animals there are 4 dose levels.
print(bartlett.test(SVL_mm ~ DoseLabel, data = meta_Fusk))
aov_fusk <- aov(SVL_mm ~ DoseLabel, data = meta_Fusk)
print(summary(aov_fusk))
print(TukeyHSD(aov_fusk))
kw_fusk <- kruskal.test(SVL_mm ~ DoseLabel, data = meta_Fusk)
print(kw_fusk)
sp_fusk <- cor.test(log10(meta_Fusk$FusariumDose), meta_Fusk$SVL_mm,
                    method = "spearman", exact = FALSE)
print(sp_fusk)

# Two-way ANOVA: SVL ~ Antibiotic * Fusarium (full dataset)
# (only AllAB has both Fusarium levels, so the interaction is not estimable)
aov_full <- aov(SVL_mm ~ Antibiotic * Fusarium, data = meta)
print(summary(aov_full))
#            Df Sum Sq Mean Sq F value Pr(>F)
# Antibiotic  7   4.72  0.6750   0.901  0.514
# Fusarium    1   0.07  0.0701   0.094  0.761
# Residuals  44  32.95  0.7489

#---------------------------------------
# Figure: (A) SVL by antibiotic, no Fusarium; (B) SVL vs Fusarium dose
#---------------------------------------
y_lim <- range(meta$SVL_mm) + c(-0.3, 0.3)
n_lab <- meta_noFusk %>% count(Antibiotic)

pA <- ggplot(meta_noFusk, aes(Antibiotic, SVL_mm, colour = Antibiotic)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.35, colour = "grey30",
               fill = NA) +
  geom_jitter(width = 0.12, height = 0, size = 1.2, alpha = 0.9) +
  geom_text(data = n_lab, aes(Antibiotic, y_lim[1], label = paste0("n=", n)),
            colour = "grey30", size = 6 / .pt, vjust = 0) +
  annotate("text", x = 0.6, y = y_lim[2], hjust = 0, vjust = 1, size = 6.5 / .pt,
           label = paste0("Kruskal-Wallis ", fmt_p(kw_abx$p.value))) +
  scale_colour_manual(values = abx_cols, guide = "none") +
  scale_y_continuous(limits = y_lim) +
  labs(x = NULL, y = "SVL (mm)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ctrl_med <- median(meta_noFusk$SVL_mm[meta_noFusk$Antibiotic == "Control"])
pB <- ggplot(meta_Fusk, aes(FusariumDose, SVL_mm, shape = Pretreatment, colour = Pretreatment)) +
  geom_hline(yintercept = ctrl_med, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  annotate("text", x = 10^3.35, y = ctrl_med, label = "Control median", hjust = 0,
           vjust = -0.4, size = 6 / .pt, colour = "grey40") +
  geom_point(position = position_jitter(width = 0.04, height = 0, seed = 1), size = 1.5) +
  annotate("text", x = 10^1.8, y = y_lim[2], hjust = 0, vjust = 1, size = 6.5 / .pt, parse = TRUE,
           label = paste0("Spearman~rho == '", sprintf("%.2f", sp_fusk$estimate), ", ",
                          fmt_p(sp_fusk$p.value), "'")) +
  annotate("text", x = 10^1.8, y = y_lim[2] - 0.28, hjust = 0, vjust = 1, size = 6.5 / .pt,
           label = paste0("Kruskal-Wallis ", fmt_p(kw_fusk$p.value))) +
  scale_x_log10(breaks = 10^(2:5), labels = function(x) parse(text = paste0("10^", log10(x))),
                limits = c(10^1.8, 10^5.25)) +
  scale_shape_manual(values = c(None = 16, AllAB = 17), name = "Antibiotic pretreatment") +
  scale_colour_manual(values = c(None = "#E34A33", AllAB = unname(abx_cols["AllAB"])),
                      name = "Antibiotic pretreatment") +
  scale_y_continuous(limits = y_lim) +
  labs(x = "Fusarium dose (spores/mL)", y = "SVL (mm)") +
  theme_pub() +
  theme(legend.position = "top", legend.justification = "right",
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -4, 0),
        legend.title = element_text(size = 6.5), legend.text = element_text(size = 6.5))

p_svl <- pA + pB + plot_layout(widths = c(1.15, 1))
print(p_svl)
save_pub(p_svl, file.path(dir, "SVL_antibiotic_fusarium.pdf"), width_mm = 180, height_mm = 75)
