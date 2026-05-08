############################################
##Microbiome analysis
############################################
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)

#data
meta <- read.delim("~/Documents/GitHub/Gnotobiotic_Xenopus/molly_xenopus_map.txt", header = TRUE)
otu_no_contam <- read.delim("~/Documents/GitHub/Gnotobiotic_Xenopus/no_contam_asv_table.txt",
                            row.names = 1, header = TRUE)
# Remove PCR control
otu_no_contam <- otu_no_contam[, -54]

# Rarefy
set.seed(515)
cort_rare <- rrarefy(t(otu_no_contam), sample = 2832)

#PCoA on Bray-Curtis
ko_pcoa <- capscale(cort_rare ~ 1, distance = "bray")

ko.scores <- scores(ko_pcoa)
ko.coords <- as.data.frame(ko.scores$sites)
ko.coords$SampleID <- row.names(ko.coords)
ko.coords <- merge(ko.coords, meta, by = "SampleID")

pct1 <- 100 * round(ko_pcoa$CA$eig[1] / sum(ko_pcoa$CA$eig), 3)
pct2 <- 100 * round(ko_pcoa$CA$eig[2] / sum(ko_pcoa$CA$eig), 3)

cb_palette <- c(
  "#332288", "#88CCEE", "#44AA99", "#117733",
  "#DDCC77", "#CC6677", "#AA4499", "#882255"
)

#adonis
#subset metadata to match rarefied sample order
meta_ordered <- meta[match(rownames(cort_rare), meta$SampleID), ]

# Full model: Antibiotic + Fusarium + interaction
adonis_full <- adonis2(
  cort_rare ~ Antibiotic * Fusarium,
  data = meta_ordered,
  distance = "bray",
  permutations = 999,
  by = "terms"
)
print(adonis_full)
#            Df SumOfSqs      R2      F Pr(>F)    
#Antibiotic  7   6.2377 0.47615 5.9506  0.001 ***
#Fusarium    1   0.2736 0.02089 1.8272  0.055 .  
#Residual   44   6.5890 0.50296                  
#Total      52  13.1004 1.00000 

# Separate models for each factor
adonis_abx <- adonis2(cort_rare ~ Antibiotic, data = meta_ordered,
                      distance = "bray", permutations = 999)
print(adonis_abx)
#         Df SumOfSqs      R2      F Pr(>F)    
#Model     7   6.2377 0.47615 5.8432  0.001 ***
#Residual 45   6.8627 0.52385                  
#Total    52  13.1004 1.00000 

adonis_fusk <- adonis2(cort_rare ~ Fusarium, data = meta_ordered,
                       distance = "bray", permutations = 999)
print(adonis_fusk)
#         Df SumOfSqs      R2     F Pr(>F)    
#Model     1   1.9754 0.15079 9.056  0.001 ***
#Residual 51  11.1250 0.84921                 
#Total    52  13.1004 1.00000 

#Betadisper
bray_dist <- vegdist(cort_rare, method = "bray")
bd_abx  <- betadisper(bray_dist, meta_ordered$Antibiotic)
bd_fusk <- betadisper(bray_dist, meta_ordered$Fusarium)
print(permutest(bd_abx))
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#Groups     7 0.46578 0.066540 8.4975    999  0.001 ***
#Residuals 45 0.35238 0.007831 
print(permutest(bd_fusk))
#          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.01168 0.0116786 1.4172    999  0.243
#Residuals 51 0.42027 0.0082406

# PCoA plot

p_all <- ggplot(ko.coords, aes(MDS1, MDS2, color = Antibiotic, shape = Fusarium)) +
  stat_ellipse(aes(group = Antibiotic), type = "t", level = 0.95,
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = cb_palette) +
  xlab(paste0("PC1 – ", pct1, "%")) +
  ylab(paste0("PC2 – ", pct2, "%"))

print(p_all)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/pcoa_all_ellipses.pdf", p_all, width = 8, height = 6, dpi = 150)

#PCoA split: No Fusarium samples only
ko_noFusk <- ko.coords[ko.coords$Fusarium == "No", ]

p_noFusk <- ggplot(ko_noFusk, aes(MDS1, MDS2, color = Antibiotic)) +
  stat_ellipse(aes(group = Antibiotic), type = "t", level = 0.95,
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = cb_palette) +
  xlab(paste0("PC1 – ", pct1, "%")) +
  ylab(paste0("PC2 – ", pct2, "%")) +
  labs(title = "PCoA – Antibiotic treatments (no Fusarium)")

print(p_noFusk)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/pcoa_antibiotics_ellipses.pdf", p_noFusk, width = 7, height = 5, dpi = 150)

#PCoA split: Fusarium samples only
ko_Fusk <- ko.coords[ko.coords$Fusarium == "Yes", ]
ko_Fusk$DoseLabel <- factor(paste0("Fk ", as.integer(ko_Fusk$FusariumDose)),
                            levels = paste0("Fk ", c(100, 1000, 10000, 100000)))

p_Fusk <- ggplot(ko_Fusk, aes(MDS1, MDS2, color = DoseLabel, shape = Antibiotic)) +
  stat_ellipse(aes(group = DoseLabel), type = "t", level = 0.95,
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_brewer(palette = "RdYlBu", direction = -1) +
  xlab(paste0("PC1 – ", pct1, "%")) +
  ylab(paste0("PC2 – ", pct2, "%")) +
  labs(title    = "PCoA – Fusarium treatments",
       color    = "Fusarium dose",
       shape    = "Antibiotic")

print(p_Fusk)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/pcoa_fusarium_ellipses.pdf", p_Fusk, width = 7, height = 5, dpi = 150)

#Alpha diversity
asv_table <- read.delim("~/Documents/GitHub/Gnotobiotic_Xenopus/no_contam_asv_table.txt",
                        row.names = 1, header = TRUE)
asv_table <- asv_table[, -54]
cort_rare2 <- rrarefy(t(asv_table), sample = 2832)

cort.alph <- as.data.frame(specnumber(cort_rare2))
cort.alph$SampleID <- row.names(cort.alph)
cort.alph <- merge(cort.alph, meta, by = "SampleID")
names(cort.alph)[2] <- "Richness"

# No-Fusarium: richness by antibiotic
cort.alph2 <- cort.alph[cort.alph$Fusarium == "No", ]
p_alpha_abx <- ggplot(cort.alph2, aes(Antibiotic, Richness)) +
  geom_boxplot() +
  theme_bw() + xlab("") + coord_flip() + ylab("sOTU Richness") +
  labs(title = "Alpha diversity – antibiotic treatments")
print(p_alpha_abx)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/alpha_antibiotics.pdf", p_alpha_abx, width = 7, height = 5, dpi = 150)

# Fusarium: richness vs dose
alphy <- cort.alph[cort.alph$Fusarium == "Yes", ]
alphy$DoseLabel <- factor(paste0("Fk ", as.integer(alphy$FusariumDose)),
                          levels = paste0("Fk ", c(100, 1000, 10000, 100000)))

p_alpha_fusk_pt <- ggplot(alphy, aes(as.numeric(FusariumDose), Richness)) +
  geom_point() + scale_y_log10() + scale_x_log10() +
  stat_cor(method = "spearman") + geom_smooth(method = "lm") +
  theme_bw() + ylab("sOTU Richness") + xlab("Fusarium dose (spores/mL)") +
  labs(title = "Alpha diversity vs Fusarium dose")
print(p_alpha_fusk_pt)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/alpha_fusarium_dose.pdf", p_alpha_fusk_pt, width = 6, height = 5, dpi = 150)

p_alpha_fusk_box <- ggplot(alphy, aes(DoseLabel, Richness)) +
  geom_boxplot() + scale_y_log10() +
  theme_bw() + ylab("sOTU Richness") + xlab("Fusarium dose (spores/mL)") +
  labs(title = "Alpha diversity by Fusarium dose")
print(p_alpha_fusk_box)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/alpha_fusarium_boxplot.pdf", p_alpha_fusk_box, width = 6, height = 5, dpi = 150)

#Stats: homogeneity of variance + ANOVA + Tukey
#Bartlett test: Antibiotic
print(bartlett.test(cort.alph$Richness ~ cort.alph$Antibiotic))
#Bartlett's K-squared = 4.7918, df = 7, p-value = 0.6853

#Bartlett test: Fusarium
print(bartlett.test(cort.alph$Richness ~ cort.alph$Fusarium))
#Bartlett's K-squared = 0.06422, df = 1, p-value = 0.7999

ANOVA: Richness ~ Antibiotic
aov_abx <- aov(cort.alph$Richness ~ cort.alph$Antibiotic + cort.alph$Fusarium)
print(summary(aov_abx))
#Df Sum Sq Mean Sq F value   Pr(>F)    
#cort.alph$Antibiotic  7 1722.0  246.00   5.525 0.000132 ***
#cort.alph$Fusarium    1   10.2   10.21   0.229 0.634436    
#Residuals            44 1959.1   44.53

print(TukeyHSD(aov_abx))

