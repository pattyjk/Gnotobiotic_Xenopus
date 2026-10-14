###################################
##XL Size
###################################
library(ggplot2)
library(dplyr)
library(ggpubr)

#laod data
meta<-read.delim('~/Documents/GitHub/Gnotobiotic_Xenopus/molly_xenopus_map.txt', header=T)

# Drop PCR control and any rows missing SVL
meta <- meta %>%
  filter(Control == "No", !is.na(SVL_mm)) %>%
  mutate(
    SVL_mm       = as.numeric(SVL_mm),
    FusariumDose = as.integer(FusariumDose),
    DoseLabel    = factor(paste0("Fk ", FusariumDose),
                          levels = paste0("Fk ", c(0, 100, 1000, 10000, 100000))),
    Antibiotic   = factor(Antibiotic,
                          levels = c("Control","Trimethoprim","Gentamycin",
                                     "Ampicillin","PSN","AllAB",
                                     "Sulfamethazine","FuskNoAB"))
  )

cb_palette <- c(
  "#332288","#88CCEE","#44AA99","#117733",
  "#DDCC77","#CC6677","#AA4499","#882255"
)

#SVL by antibiotic — no Fusarium samples
meta_noFusk <- meta %>% filter(Fusarium == "No")

p1 <- ggplot(meta_noFusk, aes(x = Antibiotic, y = SVL_mm, fill = Antibiotic)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = cb_palette) +
  coord_flip() +
  theme_bw(base_size = 13) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("SVL (mm)") +
  labs(title = "Snout-to-vent length by antibiotic treatment",
       subtitle = "No Fusarium exposure")

print(p1)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/SVL_antibiotics.pdf", p1, width = 7, height = 5, dpi = 150)

#SVL by Fusarium dose 
meta_Fusk <- meta %>% filter(Fusarium == "Yes")

p2 <- ggplot(meta_Fusk, aes(x = DoseLabel, y = SVL_mm, fill = DoseLabel)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none") +
  xlab("Fusarium dose (spores/mL)") +
  ylab("SVL (mm)") +
  labs(title = "Snout-to-vent length by Fusarium dose",
       subtitle = "Fusarium-exposed samples only")

print(p2)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/SVL_fusarium_dose.pdf", p2, width = 6, height = 5, dpi = 150)

#SVL vs Fusarium dose continuous 
p3 <- ggplot(meta_Fusk, aes(x = as.numeric(FusariumDose), y = SVL_mm,
                            color = Antibiotic)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = 0.05, inherit.aes = FALSE,
           aes(x = as.numeric(FusariumDose), y = SVL_mm)) +
  scale_x_log10() +
  scale_color_manual(values = cb_palette) +
  theme_bw(base_size = 13) +
  xlab("Fusarium dose (spores/mL, log scale)") +
  ylab("SVL (mm)") +
  labs(title = "SVL vs Fusarium dose",
       color = "Antibiotic")

print(p3)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/SVL_fusarium_continuous.pdf", p3, width = 7, height = 5, dpi = 150)

#SVL by Fusarium (Yes/No) + Antibiotic — full dataset ─────────────
p4 <- ggplot(meta, aes(x = Antibiotic, y = SVL_mm, fill = Fusarium)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,
               position = position_dodge(width = 0.75), width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
             size = 1.8, alpha = 0.7) +
  scale_fill_manual(values = c("No" = "#88CCEE", "Yes" = "#CC6677")) +
  coord_flip() +
  theme_bw(base_size = 13) +
  xlab("") + ylab("SVL (mm)") +
  labs(title = "SVL by antibiotic and Fusarium exposure",
       fill = "Fusarium")

print(p4)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/SVL_antibiotic_fusarium.pdf", p4, width = 8, height = 5, dpi = 150)

#ANOVA: SVL ~ Antibiotic (no Fusarium)
print(bartlett.test(SVL_mm ~ Antibiotic, data = meta_noFusk))
#Bartlett's K-squared = 4.3924, df = 6, p-value = 0.6237

aov_abx <- aov(SVL_mm ~ Antibiotic, data = meta_noFusk)
print(summary(aov_abx))
#            Df Sum Sq Mean Sq F value Pr(>F)
#Antibiotic   6  2.728  0.4546   0.701  0.651
#Residuals   24 15.560  0.6483 

print(TukeyHSD(aov_abx))

#ANOVA: SVL ~ FusariumDose (Fusarium samples only)
print(bartlett.test(SVL_mm ~ DoseLabel, data = meta_Fusk))
#Bartlett's K-squared = 4.0939, df = 2, p-value = 0.1291

aov_fusk <- aov(SVL_mm ~ DoseLabel, data = meta_Fusk)
print(summary(aov_fusk))
#            Df Sum Sq Mean Sq F value Pr(>F)
#DoseLabel    2  1.863  0.9315   0.998    0.4
#Residuals   11 10.272  0.9338    

print(TukeyHSD(aov_fusk))

#Two-way ANOVA: SVL ~ Antibiotic * Fusarium
aov_full <- aov(SVL_mm ~ Antibiotic * Fusarium, data = meta)
print(summary(aov_full))
#            Df Sum Sq Mean Sq F value Pr(>F)
#Antibiotic   7   4.72  0.6750   0.901  0.514
#Fusarium     1   0.07  0.0701   0.094  0.761
#Residuals   44  32.95  0.7489 