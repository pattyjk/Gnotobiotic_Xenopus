############################################
## Microbiome analysis: rarefaction, Bray-Curtis PCoA, PERMANOVA,
## dispersion, alpha diversity (ASV richness)
############################################
library(vegan)
library(ggplot2)
library(dplyr)

dir <- "~/Documents/GitHub/Gnotobiotic_Xenopus/"
source(file.path(dir, "pub_theme.R"))   # theme_pub(), abx_cols, fk_dose_factor(), save_pub(), ...
fig_dir <- "~/Documents/GitHub/Gnotobiotic_Xenopus/figures"                         # figures are written next to the data

# ---- data ----
meta <- read.delim(file.path(dir, "molly_xenopus_map.txt"), header = TRUE,
                   stringsAsFactors = FALSE)
otu_no_contam <- read.delim(file.path(dir, "no_contam_asv_table.txt"),
                            row.names = 1, header = TRUE, check.names = FALSE)

# Remove PCR control
pcr_ctrl <- meta$SampleID[meta$Control == "Yes"]          # "PCRCtrl"
otu_no_contam <- otu_no_contam[, !colnames(otu_no_contam) %in% c(pcr_ctrl, "PCRCtrl")]
meta <- meta[meta$Control != "Yes", ]
meta$FusariumDose <- as.numeric(meta$FusariumDose)   

# Canonical names (Gentamycin -> Gentamicin, FuskNoAB -> None) and plotting variables.
# Renaming factor levels does not change any test statistic.
meta <- meta %>%
  mutate(Antibiotic_raw = Antibiotic,
         Antibiotic     = recode_abx(Antibiotic),           # "None" = Fusarium, no antibiotic
         Pretreatment   = factor(ifelse(Antibiotic == "AllAB", "AllAB", "None"),
                                 levels = c("None", "AllAB")),
         FkExposure     = factor(ifelse(Fusarium == "Yes", "Exposed", "Not exposed"),
                                 levels = c("Not exposed", "Exposed")),
         Dose           = fk_dose_factor(FusariumDose),     # keeps 10^5 (was NA)
         Group          = ifelse(Fusarium == "Yes", paste(Antibiotic, "+ Fk"), Antibiotic))

# Rarefy 
# diversity 
set.seed(515)
cort_rare <- rrarefy(t(otu_no_contam), sample = 2832)

# ---- PCoA on Bray-Curtis ----
ko_pcoa <- capscale(cort_rare ~ 1, distance = "bray")
ko.coords <- as.data.frame(scores(ko_pcoa, display = "sites"))
ko.coords$SampleID <- row.names(ko.coords)
ko.coords <- merge(ko.coords, meta, by = "SampleID")

pct1 <- 100 * round(ko_pcoa$CA$eig[1] / sum(ko_pcoa$CA$eig), 3)
pct2 <- 100 * round(ko_pcoa$CA$eig[2] / sum(ko_pcoa$CA$eig), 3)
lab_x <- paste0("PCoA 1 (", pct1, "% of variation)")
lab_y <- paste0("PCoA 2 (", pct2, "% of variation)")

# ---- PERMANOVA (adonis2) ----
meta_ordered <- meta[match(rownames(cort_rare), meta$SampleID), ]

# Full model: Antibiotic + Fusarium + interaction
# (no interaction row is estimable: Fusarium exposure is nested in the Antibiotic levels None / AllAB)
adonis_full <- adonis2(cort_rare ~ Antibiotic * Fusarium, data = meta_ordered,
                       distance = "bray", permutations = 999, by = "terms")
print(adonis_full)
# Previous run:
#            Df SumOfSqs      R2      F Pr(>F)
#Antibiotic  7   6.2377 0.47615 5.9506  0.001 ***
#Fusarium    1   0.2736 0.02089 1.8272  0.055 .
#Residual   44   6.5890 0.50296
#Total      52  13.1004 1.00000

adonis_abx <- adonis2(cort_rare ~ Antibiotic, data = meta_ordered,
                      distance = "bray", permutations = 999)
print(adonis_abx)
#Model     7   6.2377 0.47615 5.8432  0.001 ***

adonis_fusk <- adonis2(cort_rare ~ Fusarium, data = meta_ordered,
                       distance = "bray", permutations = 999)
print(adonis_fusk)
#Model     1   1.9754 0.15079 9.056  0.001 ***
#Residual 45   6.8627 0.52385                  
#Total    52  13.1004 1.00000

# ---- Betadisper ----
bray_dist <- vegdist(cort_rare, method = "bray")
bd_abx  <- betadisper(bray_dist, meta_ordered$Antibiotic)
bd_fusk <- betadisper(bray_dist, meta_ordered$Fusarium)
print(permutest(bd_abx))
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
#Groups     7 0.46578 0.066540 8.4975    999  0.001 ***
#Residuals 45 0.35238 0.007831 

print(permutest(bd_fusk))
#Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.01168 0.0116786 1.4172    999  0.233
#Residuals 51 0.42027 0.0082406 

# Sub-models used to annotate the split PCoA panels (added)
idx_ab <- meta_ordered$Fusarium == "No"
adonis_ab_only <- adonis2(cort_rare[idx_ab, ] ~ Antibiotic, data = meta_ordered[idx_ab, ],
                          distance = "bray", permutations = 999)
print(adonis_ab_only)
#.         Df SumOfSqs      R2      F Pr(>F)    
#Model     6   4.3218 0.63443 6.9417  0.001 ***
#Residual 24   2.4903 0.36557                  
#Total    30   6.8121 1.00000 

idx_fk <- meta_ordered$Fusarium == "Yes"
adonis_fk_only <- adonis2(cort_rare[idx_fk, ] ~ Pretreatment + Dose,
                          data = droplevels(meta_ordered[idx_fk, ]),
                          distance = "bray", permutations = 999, by = "terms")
print(adonis_fk_only)
#Df SumOfSqs      R2     F Pr(>F)
#Pretreatment  1   0.2142 0.04966 1.034  0.413
#Dose          3   0.5773 0.13385 0.929  0.569
#Residual     17   3.5215 0.81649             
#Total        21   4.3129 1.00000 

dp <- function(e) paste(deparse(e, width.cutoff = 500L), collapse = "")
r2p <- function(tab, term, label = term) {
  bquote(.(label) ~ italic(R)^2 == .(sprintf("%.2f", tab[term, "R2"])) * "," ~
           italic(p) == .(sprintf("%.3g", tab[term, "Pr(>F)"])))
}

# ---- shared PCoA helpers ----
ell_min_n <- 4                               
ellipse_data <- function(d, grp) {
  d %>% group_by(across(all_of(grp))) %>% filter(n() >= ell_min_n) %>% ungroup()
}
abx_present <- abx_levels[abx_levels %in% ko.coords$Antibiotic]
abx_legend_labels <- function(x) ifelse(x == "Control", "None", x)
shape_fk <- c("Not exposed" = 21, "Exposed" = 24)

# ---- Fig: PCoA, all samples ----
ko.coords$Antibiotic <- factor(ko.coords$Antibiotic, levels = c(abx_levels, "None"))
ko.coords$EllGroup <- paste(ko.coords$Antibiotic, ko.coords$FkExposure)

ann_all <- c(
  dp(r2p(adonis_full, "Antibiotic", "Antibiotic:")),
  dp(r2p(adonis_full, "Fusarium", "Fusarium:")))

p_all <- ggplot(ko.coords, aes(MDS1, MDS2)) +
  stat_ellipse(data = ellipse_data(ko.coords, "EllGroup"),
               aes(group = EllGroup, colour = Antibiotic, linetype = FkExposure),
               type = "t", level = 0.95, linewidth = 0.4, show.legend = FALSE) +
  geom_point(aes(fill = Antibiotic, shape = FkExposure), size = 2, stroke = 0.3,
             colour = "black") +
  scale_fill_manual(values = abx_cols, breaks = abx_present, labels = abx_legend_labels,
                    name = "Antibiotic", drop = TRUE) +
  scale_colour_manual(values = abx_cols, guide = "none") +
  scale_shape_manual(values = shape_fk, name = "Fusarium") +
  scale_linetype_manual(values = c("Not exposed" = "solid", "Exposed" = "dashed"), guide = "none") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2.2), order = 1),
         shape = guide_legend(override.aes = list(fill = "grey60", size = 2.2), order = 2)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = c(1.3, 2.7),
           label = ann_all, parse = TRUE, size = 6.5 / .pt) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  labs(x = lab_x, y = lab_y) +
  theme_pub()

save_pub(p_all, file.path(fig_dir, "pcoa_all_ellipses"), width_mm = 120, height_mm = 85)

# ---- Fig: PCoA, antibiotic treatments only (no Fusarium) ----
ko_noFusk <- droplevels(ko.coords[ko.coords$Fusarium == "No", ])

p_noFusk <- ggplot(ko_noFusk, aes(MDS1, MDS2)) +
  stat_ellipse(data = ellipse_data(ko_noFusk, "Antibiotic"),
               aes(group = Antibiotic, colour = Antibiotic),
               type = "t", level = 0.95, linewidth = 0.4, show.legend = FALSE) +
  geom_point(aes(fill = Antibiotic), shape = 21, size = 2, stroke = 0.3, colour = "black") +
  scale_fill_manual(values = abx_cols, breaks = abx_levels[abx_levels %in% ko_noFusk$Antibiotic],
                    name = "Antibiotic") +
  scale_colour_manual(values = abx_cols, guide = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 1.3,
           label = dp(r2p(adonis_ab_only, "Model", "PERMANOVA:")),
           parse = TRUE, size = 6.5 / .pt) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  labs(x = lab_x, y = lab_y) +
  theme_pub()

save_pub(p_noFusk, file.path(fig_dir, "pcoa_antibiotics_ellipses"), width_mm = 110, height_mm = 80)

# ---- Fig: PCoA, Fusarium-exposed animals (colour = dose, shape = pretreatment) ----
ko_Fusk <- ko.coords[ko.coords$Fusarium == "Yes", ]
stopifnot(!any(is.na(ko_Fusk$Dose)))           # 10^5 must not be NA
dose_present <- levels(droplevels(ko_Fusk$Dose))

ann_fk <- c(dp(r2p(adonis_fk_only, "Pretreatment", "Pretreatment:")),
            dp(r2p(adonis_fk_only, "Dose", "Dose:")))

p_Fusk <- ggplot(ko_Fusk, aes(MDS1, MDS2)) +
  # ellipses by pretreatment (n = 17 None, 5 AllAB); dose groups mix pretreatments
  stat_ellipse(data = ellipse_data(ko_Fusk, "Pretreatment"),
               aes(group = Pretreatment, linetype = Pretreatment),
               colour = "grey40", type = "t", level = 0.95, linewidth = 0.4,
               show.legend = FALSE) +
  geom_point(aes(fill = Dose, shape = Pretreatment), size = 2.2, stroke = 0.3,
             colour = "black") +
  scale_fill_manual(values = fk_dose_cols, breaks = dose_present, labels = fk_labeller,
                    name = "Fusarium dose\n(spores/mL)", drop = TRUE) +
  scale_shape_manual(values = c("None" = 21, "AllAB" = 24),
                     labels = c("None" = "None", "AllAB" = "AllAB"),
                     name = "Antibiotic\npretreatment") +
  scale_linetype_manual(values = c("None" = "solid", "AllAB" = "dashed"), guide = "none") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2.4), order = 1),
         shape = guide_legend(override.aes = list(fill = "grey60", size = 2.2), order = 2)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = c(1.3, 2.7),
           label = ann_fk, parse = TRUE, size = 6.5 / .pt) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  labs(x = lab_x, y = lab_y) +
  theme_pub()

save_pub(p_Fusk, file.path(fig_dir, "pcoa_fusarium_ellipses"), width_mm = 110, height_mm = 80)

# ---- Alpha diversity: ASV richness on the same rarefied table ----
cort.alph <- data.frame(SampleID = rownames(cort_rare),
                        Richness = specnumber(cort_rare))
cort.alph <- merge(cort.alph, meta, by = "SampleID")

# Richness by antibiotic, no Fusarium (points over boxes; n = 3-8 per group)
cort.alph2 <- cort.alph %>% filter(Fusarium == "No") %>%
  mutate(Antibiotic = factor(Antibiotic, levels = rev(abx_levels[abx_levels %in% Antibiotic])))

kw_abx <- kruskal.test(Richness ~ Antibiotic, data = cort.alph2)
print(kw_abx)

#Kruskal-Wallis chi-squared = 19.893, df = 6, p-value =0.002894

wil <- cort.alph2 %>% filter(Antibiotic != "Control") %>% group_by(Antibiotic) %>%
  summarise(p = suppressWarnings(wilcox.test(
    Richness, cort.alph2$Richness[cort.alph2$Antibiotic == "Control"], exact = FALSE)$p.value),
    .groups = "drop") %>%
  mutate(p_adj = p.adjust(p, method = "BH"),
         label = ifelse(p_adj < 0.001, "<0.001", sprintf("%.2g", p_adj)))
print(wil)
#  Antibiotic          p  p_adj label
#1 Trimethoprim   0.374  0.374  0.37 
#2 Gentamicin     0.0570 0.0855 0.086
#3 Sulfamethazine 0.0806 0.0968 0.097
#4 PSN            0.0314 0.0649 0.065
#5 AllAB          0.0325 0.0649 0.065
#6 Ampicillin     0.0179 0.0649 0.065

y_lab_pos <- max(cort.alph2$Richness) * 1.08

p_alpha_abx <- ggplot(cort.alph2, aes(Richness, Antibiotic)) +
  geom_boxplot(aes(fill = Antibiotic), outlier.shape = NA, width = 0.6,
               linewidth = 0.3, alpha = 0.5) +
  geom_point(aes(fill = Antibiotic), shape = 21, size = 1.6, stroke = 0.3,
             position = position_jitter(height = 0.12, width = 0, seed = 1)) +
  geom_text(data = wil, aes(x = y_lab_pos, y = Antibiotic, label = label),
            size = 6.5 / .pt, hjust = 0) +
  annotate("text", x = y_lab_pos, y = length(levels(cort.alph2$Antibiotic)) + 0.75,
           label = "italic(p)[adj]", parse = TRUE, size = 6.5 / .pt, hjust = 0) +
  scale_fill_manual(values = abx_cols, guide = "none") +
  scale_x_continuous(limits = c(0, y_lab_pos * 1.2), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(clip = "off") +
  labs(x = "ASV richness (rarefied to 2,832 reads)", y = NULL,
       caption = paste0("Kruskal-Wallis ", fmt_p(kw_abx$p.value),
                        "; Wilcoxon vs Control, BH-adjusted")) +
  theme_pub() +
  theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey25"),
        plot.caption.position = "plot",
        plot.margin = margin(10, 6, 4, 4))

save_pub(p_alpha_abx, file.path(fig_dir, "alpha_antibiotics"), width_mm = 85, height_mm = 70)

# Richness vs Fusarium dose (continuous log10 dose; all Fk-exposed animals)
alphy <- cort.alph %>% filter(Fusarium == "Yes")
stopifnot(!any(is.na(alphy$Dose)))
sp_all  <- cor.test(alphy$FusariumDose, alphy$Richness, method = "spearman", exact = FALSE)
sp_none <- with(subset(alphy, Pretreatment == "None"),
                cor.test(FusariumDose, Richness, method = "spearman", exact = FALSE))
print(sp_all); print(sp_none)
sp_lab <- function(ct, who) {
  dp(bquote(.(who) * ":" ~ rho == .(sprintf("%.2f", ct$estimate)) * "," ~
                   italic(p) == .(sprintf("%.2g", ct$p.value))))
}

ctrl_rich <- median(cort.alph$Richness[cort.alph$Antibiotic == "Control" & cort.alph$Fusarium == "No"])
# Styled to match the other dose panels (qPCR load, SVL): linear y on the same
# 0-60 scale as the antibiotic richness panel, colour/shape = pretreatment,
# dashed line = Control median.
p_alpha_fusk_pt <- ggplot(alphy, aes(FusariumDose, Richness, fill = Pretreatment, shape = Pretreatment)) +
  geom_hline(yintercept = ctrl_rich, linetype = 2, colour = "grey60", linewidth = 0.4) +
  annotate("text", x = 10^1.95, y = ctrl_rich, label = "Control median", hjust = 0,
           vjust = -0.5, size = 6 / .pt, colour = "grey40") +
  geom_point(size = 1.9, stroke = 0.3, colour = "black",
             position = position_jitter(width = 0.06, height = 0, seed = 2)) +
  scale_x_log10(breaks = 10^(2:5), labels = scales::label_math(10^.x, format = log10),
                limits = c(10^1.85, 10^5.2)) +
  scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = c("None" = "#E34A33", "AllAB" = unname(abx_cols["AllAB"])),
                    name = "Antibiotic\npretreatment") +
  scale_shape_manual(values = c("None" = 21, "AllAB" = 24), name = "Antibiotic\npretreatment") +
  annotate("text", x = 10^1.95, y = 60, hjust = 0, vjust = c(1, 2.3),
           label = c(sp_lab(sp_all, "All"), sp_lab(sp_none, "No antibiotic")),
           parse = TRUE, size = 6.5 / .pt) +
  labs(x = "Fusarium dose (spores/mL)", y = "ASV richness (rarefied to 2,832 reads)") +
  theme_pub() +
  theme(legend.position = "top", legend.justification = "right",
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -4, 0),
        legend.title = element_text(size = 6.5), legend.text = element_text(size = 6.5))

save_pub(p_alpha_fusk_pt, file.path(fig_dir, "alpha_fusarium_dose"), width_mm = 85, height_mm = 70)
# NOTE: the former box-plot version (alpha_fusarium_boxplot.pdf) duplicated this figure and
# showed 10^5 as "NA"; it has been removed in favour of the continuous-dose plot above.

# ---- Stats: homogeneity of variance + ANOVA + Tukey (all samples, unchanged model) ----
print(bartlett.test(Richness ~ Antibiotic, data = cort.alph))
#Bartlett's K-squared = 4.7918, df = 7, p-value = 0.6853   (previous run, second rarefaction draw)
#Bartlett's K-squared = 6.4029, df = 7, p-value = 0.4936   (this script, single rarefied table)
print(bartlett.test(Richness ~ Fusarium, data = cort.alph))
#Bartlett's K-squared = 0.06422, df = 1, p-value = 0.7999  (previous run)
#Bartlett's K-squared = 0.0762,  df = 1, p-value = 0.7825  (this script)

# ANOVA: Richness ~ Antibiotic + Fusarium   (this line was an un-commented header -> syntax error)
aov_abx <- aov(Richness ~ Antibiotic + Fusarium, data = cort.alph)
print(summary(aov_abx))
#Antibiotic  7 1722.0  246.00   5.525 0.000132 ***   (previous run)
#Fusarium    1   10.2   10.21   0.229 0.634436
#Antibiotic  7 1748.1  249.73   5.530 0.000131 ***   (this script, single rarefied table)
#Fusarium    1    9.1    9.07   0.201 0.656136
print(TukeyHSD(aov_abx))
