library(phyloseq)
library(ANCOMBC)
library(ALDEx2)
library(tidyverse)

KO_FILE  <- "~/Documents/GitHub/Gnotobiotic_Xenopus/picrust_out/KO_metagenome_out/pred_metagenome_unstrat.tsv"
MAP_FILE <- "~/Documents/GitHub/Gnotobiotic_Xenopus/molly_xenopus_map.txt"
OUT_DIR  <- "~/Documents/GitHub/Gnotobiotic_Xenopus/ko_results"
dir.create(OUT_DIR, showWarnings = FALSE)


# ---- 1. Read + round + align -------------------------------------------------
ko <- read.delim(KO_FILE, row.names = 1, check.names = FALSE)   # 'function' -> rownames
ko <- round(as.matrix(ko))                                      # predictions -> integer
storage.mode(ko) <- "integer"

map <- read.delim(MAP_FILE, row.names = 1, check.names = FALSE,
                  stringsAsFactors = FALSE)

common  <- intersect(colnames(ko), rownames(map))
dropped <- setdiff(rownames(map), colnames(ko))
message("dropped (not in KO table): ", paste(dropped, collapse = ", "))
ko  <- ko[, common]
map <- map[common, ]
ko  <- ko[rowSums(ko) > 0, ]                                    # drop empty KOs
dim(ko)                                                         # ~6.2k KOs x 51

map$Antibiotic   <- relevel(factor(map$Antibiotic), ref = "Control")
map$Fusarium     <- factor(map$Fusarium, levels = c("No", "Yes"))
map$FusariumDose <- as.numeric(map$FusariumDose)
table(map$Antibiotic, map$Fusarium)                             # confirm cell n's


#make phyloseq object
ps <- phyloseq(otu_table(ko, taxa_are_rows = TRUE), sample_data(map))

# ==============================================================================
#ANTIBIOTIC EFFECT (Fusarium = No; each antibiotic vs Control)
# ==============================================================================
ps_ab <- subset_samples(ps, Fusarium == "No")
ps_ab <- prune_taxa(taxa_sums(ps_ab) > 0, ps_ab)
sample_data(ps_ab)$Antibiotic <- droplevels(sample_data(ps_ab)$Antibiotic)

set.seed(515)
out_ab <- ancombc2(
  data = ps_ab, fix_formula = "Antibiotic", group = "Antibiotic",
  prv_cut = 0.10,        # raise (e.g. 0.25) to cut KOs / speed up
  lib_cut = 0, struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05,
  global = TRUE, dunnet = TRUE, pairwise = FALSE, n_cl = 1, verbose = TRUE)
names(out_ab)
write_tsv(out_ab$res,        file.path(OUT_DIR, "antibiotic_res.tsv"))
write_tsv(out_ab$res_global, file.path(OUT_DIR, "antibiotic_global.tsv"))
write_tsv(out_ab$res_dunn,   file.path(OUT_DIR, "antibiotic_dunnett.tsv"))


# ==============================================================================
# 4. FUSARIUM EFFECT under all-antibiotics (AllAB: Yes vs No)
# ==============================================================================
ps_allab <- subset_samples(ps, Antibiotic == "AllAB")
ps_allab <- prune_taxa(taxa_sums(ps_allab) > 0, ps_allab)
table(sample_data(ps_allab)$Fusarium)                           # 3 No, 5 Yes
set.seed(1)
out_fus_allab <- ancombc2(
  data = ps_allab, fix_formula = "Fusarium", group = "Fusarium",
  prv_cut = 0.10, lib_cut = 0, struc_zero = TRUE, neg_lb = TRUE,
  alpha = 0.05, n_cl = 1, verbose = TRUE)
write_tsv(out_fus_allab$res, file.path(OUT_DIR, "fusarium_within_allab.tsv"))


# ==============================================================================
# 5. FUSARIUM EFFECT no antibiotics (FuskNoAB[Yes] vs Control[No])
# ==============================================================================
ps_noab <- subset_samples(ps, Antibiotic %in% c("Control", "FuskNoAB"))
ps_noab <- prune_taxa(taxa_sums(ps_noab) > 0, ps_noab)
table(sample_data(ps_noab)$Fusarium)                            # 5 No, 16 Yes
set.seed(1)
out_fus_noab <- ancombc2(
  data = ps_noab, fix_formula = "Fusarium", group = "Fusarium",
  prv_cut = 0.10, lib_cut = 0, struc_zero = TRUE, neg_lb = TRUE,
  alpha = 0.05, n_cl = 1, verbose = TRUE)
write_tsv(out_fus_noab$res, file.path(OUT_DIR, "fusarium_noab.tsv"))