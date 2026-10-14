# ==============================================================================
# Differential abundance of Xenopus skin ASVs ~ Antibiotic / Fusarium
# ==============================================================================
#load libraries
#BiocManager::install(c("phyloseq", "ANCOMBC", "ALDEx2"))
library(phyloseq)
library(ANCOMBC)
library(ALDEx2)
library(tidyverse)
library(microbiome)

#read in data
# The ASV header has no label for the ID column, so its 54 fields are one fewer
# than the 55 in each data row -> R takes column 1 as row names automatically.
asv <- t(read.delim('~/Documents/GitHub/Gnotobiotic_Xenopus/genus_taxonomy.txt', row.names = 1, check.names = FALSE))
asv <- as.matrix(asv)
storage.mode(asv) <- "integer"

map <- read.delim('~/Documents/GitHub/Gnotobiotic_Xenopus/molly_xenopus_map.txt', row.names = 1, check.names = FALSE,
                  stringsAsFactors = FALSE)
dim(asv)                               
#[1] 82 54

#Drop the PCR control, light QC
asv_bio <- asv[, colnames(asv) != "PCRCtrl"]
map_bio <- map[colnames(asv_bio), ]

# drop all-zero ASVs
asv_bio <- asv_bio[rowSums(asv_bio) > 0, ]        

# library sizes
libsize <- colSums(asv_bio)
summary(libsize)                                 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2832    4370    5104    7744    8088   25121

# lowest-depth samples
sort(libsize)[1:5]      
#S37  S36  S43  S19  S38 
#2832 3070 3350 3375 3533

# typed factors with Control as the reference level
map_bio$Antibiotic   <- relevel(factor(map_bio$Antibiotic), ref = "Control")
map_bio$Fusarium     <- factor(map_bio$Fusarium, levels = c("No", "Yes"))
map_bio$FusariumDose <- as.numeric(map_bio$FusariumDose)

#Build phyloseq object
ps <- phyloseq(
  otu_table(asv_bio, taxa_are_rows = TRUE),
  sample_data(map_bio)
)

# ==============================================================================
# 5. ANTIBIOTIC EFFECT  (Fusarium = No subset; each antibiotic vs Control)
# ==============================================================================
ps_ab <- subset_samples(ps, Fusarium == "No")
ps_ab <- prune_taxa(taxa_sums(ps_ab) > 0, ps_ab)
sample_data(ps_ab)$Antibiotic <- droplevels(sample_data(ps_ab)$Antibiotic)


set.seed(515)
out_ab <- ancombc2(
  data        = ps_ab,
  fix_formula = "Antibiotic",
  group       = "Antibiotic",
  prv_cut     = 0.10,        # keep ASVs in >=10% of samples
  lib_cut     = 0,
  struc_zero  = TRUE,
  neg_lb      = TRUE,
  alpha       = 0.05,
  global      = TRUE,        # omnibus: differs across ANY antibiotic
  dunnet      = TRUE,        # each antibiotic vs Control reference
  pairwise    = FALSE,
  n_cl        = 1,
  verbose     = TRUE
)

#ASVs differing across the antibiotic panel
ab_global <- out_ab$res_global
subset(ab_global, diff_abn)          

# Per-antibiotic vs Control
ab_dunn <- out_ab$res_dunn

#write to file for safe keeping
OUT_DIR<-'~/Documents/GitHub/Gnotobiotic_Xenopus/diff_abun'
write_tsv(out_ab$res,        file.path(OUT_DIR, "antibiotic_res.tsv"))
write_tsv(ab_global,         file.path(OUT_DIR, "antibiotic_global.tsv"))
write_tsv(ab_dunn,           file.path(OUT_DIR, "antibiotic_dunnett.tsv"))


# ==============================================================================
#FUSARIUM EFFECT, UNDER ALL-ANTIBIOTICS  (AllAB: Yes vs No)
# ==============================================================================
ps_allab <- subset_samples(ps, Antibiotic == "AllAB")
ps_allab <- prune_taxa(taxa_sums(ps_allab) > 0, ps_allab)
table(sample_data(ps_allab)$Fusarium)             # 3 No, 5 Yes (small!)

out_fus_allab <- ancombc2(
  data        = ps_allab,
  fix_formula = "Fusarium",
  group       = "Fusarium",
  prv_cut     = 0.10, lib_cut = 0,
  struc_zero  = TRUE, neg_lb = TRUE, alpha = 0.05,
  n_cl = 1, verbose = TRUE
)
subset(out_fus_allab$res, get(grep("diff_FusariumYes", names(out_fus_allab$res),
                                   value = TRUE)[1]))

write_tsv(out_fus_allab$res, file.path(OUT_DIR, "fusarium_within_allab.tsv"))


# ==============================================================================
#FUSARIUM EFFECT, NO ANTIBIOTICS  (FuskNoAB[Yes] vs Control[No])
# ==============================================================================
ps_noab <- subset_samples(ps, Antibiotic %in% c("Control", "FuskNoAB"))
ps_noab <- prune_taxa(taxa_sums(ps_noab) > 0, ps_noab)
# Antibiotic perfectly aliases Fusarium here, so model Fusarium only:
table(sample_data(ps_noab)$Fusarium)              # 5 No, 17 Yes

set.seed(1)
out_fus_noab <- ancombc2(
  data        = ps_noab,
  fix_formula = "Fusarium",
  group       = "Fusarium",
  prv_cut     = 0.10, lib_cut = 0,
  struc_zero  = TRUE, neg_lb = TRUE, alpha = 0.05,
  n_cl = 1, verbose = TRUE
)
write_tsv(out_fus_noab$res, file.path(OUT_DIR, "fusarium_noab.tsv"))


# ==============================================================================
#Fusarium DOSE-RESPONSE within FuskNoAB
# ==============================================================================
ps_dose <- subset_samples(ps, Antibiotic == "FuskNoAB")
ps_dose <- prune_taxa(taxa_sums(ps_dose) > 0, ps_dose)
sample_data(ps_dose)$logDose <- log10(sample_data(ps_dose)$FusariumDose)
table(sample_data(ps_dose)$FusariumDose)

out_dose <- ancombc2(
  data        = ps_dose,
  fix_formula = "logDose",
  prv_cut     = 0.10, lib_cut = 0,
  struc_zero  = FALSE, alpha = 0.05,
  n_cl = 1, verbose = TRUE
)
write_tsv(out_dose$res, file.path(OUT_DIR, "fusarium_dose_trend.tsv"))

# ==============================================================================
# Figures for differential abundance output - Fusarium forest plot of LFC
# ==============================================================================

#libraries
library(tidyverse)
library(patchwork)
library(scales)

df <- read_tsv('~/Documents/GitHub/Gnotobiotic_Xenopus/diff_abun/fusarium_noab.tsv', show_col_types = FALSE)
names(df)[1] <- "taxon"

d <- df %>%
  pivot_longer(matches("^(lfc|se|q)_"),
               names_to = c("stat", "term"),
               names_pattern = "^(lfc|se|q)_(.*)$") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(term != "(Intercept)") %>%       # binary contrast -> one term remains
  mutate(sig   = !is.na(q) & q < ALPHA,
         label = str_sub(taxon, 1, 8))

#filter out non-sig
ALPHA <- 0.05
TOP_N <- 25
sig_d <- filter(d, sig)
if (nrow(sig_d) > 0) {
  plot_d  <- sig_d
  subtxt  <- paste0("ASVs at q < ", ALPHA)
} else {
  plot_d  <- slice_max(d, abs(lfc), n = TOP_N)
  subtxt  <- paste0("no ASV reached q < ", ALPHA,
                    "; showing top ", nrow(plot_d), " by |lfc|")
}
plot_d <- mutate(plot_d, label = fct_reorder(label, lfc))

p <- ggplot(plot_d, aes(lfc, taxon, colour = lfc > 0)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_errorbarh(aes(xmin = lfc - se, xmax = lfc + se), height = 0) +
  geom_point(aes(alpha = sig), size = 2.6) +
  scale_colour_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"),
                      labels = c("lower with Fusarium", "higher with Fusarium"),
                      name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35),
                     guide = "none") +
  ggtitle("Differential taxa- Fusarium effect") +
  xlab('Log-fold change')+
  ylab('')+
  theme_bw(base_size = 12)

ggsave('~/Documents/GitHub/Gnotobiotic_Xenopus/Figures/diff_taxa_fusarium.pdf', p, width =12, height = 7)

# ==============================================================================
# Figures for differential abundance output - heatmap for antibiotics
# ==============================================================================

library(tidyverse)
ALPHA <- 0.05
TOP_N <- 40   

#read and QC data
df <- read_tsv('~/Documents/GitHub/Gnotobiotic_Xenopus/diff_abun/antibiotic_dunnett.tsv', show_col_types = FALSE)
names(df)[1] <- "taxon"

d <- df %>%
  pivot_longer(matches("^(lfc|q)_"),
               names_to = c("stat", "term"),
               names_pattern = "^(lfc|q)_(.*)$") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(term != "(Intercept)") %>%
  mutate(term  = str_remove(term, "^Antibiotic"),   # tidy contrast labels
         sig   = !is.na(q) & q < ALPHA,
         label = str_sub(taxon, 1, 8))


#choose ASVs: sig in >=1 contrast
keep <- d %>% group_by(taxon) %>%
  summarise(any_sig = any(sig), max_abs = max(abs(lfc), na.rm = TRUE),
            .groups = "drop")

if (any(keep$any_sig)) {
  keep_taxa <- filter(keep, any_sig)$taxon
  subtxt    <- paste0(length(keep_taxa), " ASVs significant in >=1 antibiotic (q < ", ALPHA, ")")
} else {
  keep_taxa <- slice_max(keep, max_abs, n = TOP_N)$taxon
  subtxt    <- paste0("no ASV reached q < ", ALPHA, "; showing top ", length(keep_taxa), " by |lfc|")
}
d <- filter(d, taxon %in% keep_taxa)


#heatmap
lim <- max(abs(d$lfc), na.rm = TRUE)

p <- ggplot(d, aes(term, taxon, fill = lfc)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(data = filter(d, sig), aes(label = "*"), vjust = 0.75, size = 5) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim, lim),
                       name = "lfc vs\nControl", na.value = "grey92") +
  ggtitle("Antibiotic effects on Xenopus skin ASVs")+
  
  theme_minimal(base_size = 12) +
  ylab('')+
  xlab('')+
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid  = element_blank())

ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/Figures/diff_taxa_antibiotics.pdf", p, width = 12, height = 6)
