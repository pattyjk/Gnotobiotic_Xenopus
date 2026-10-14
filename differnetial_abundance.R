# ==============================================================================
# Differential abundance of Xenopus skin bacterial genera ~ Antibiotic / Fusarium
# (ANCOM-BC2 on the genus-collapsed, decontaminated count table)
# ==============================================================================
#BiocManager::install(c("phyloseq", "ANCOMBC", "ALDEx2", "microbiome"))
library(phyloseq)
library(ANCOMBC)
library(tidyverse)
library(patchwork)
library(scales)

# ---- paths: edit `dir` only --------------------------------------------------
dir     <- "~/Documents/GitHub/Gnotobiotic_Xenopus/"
OUT_DIR <- file.path(dir, "diff_abun")
FIG_DIR <- file.path(dir, "Figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(dir, "pub_theme.R"))   # theme_pub(), abx_levels/abx_cols, recode_abx(), save_pub()

# ---- read data ---------------------------------------------------------------
# genus_taxonomy.txt is samples x genera (first column "index" = SampleID);
# transpose to genera x samples for phyloseq.
asv <- t(read.delim(file.path(dir, "genus_taxonomy.txt"), row.names = 1, check.names = FALSE))
asv <- as.matrix(asv)
storage.mode(asv) <- "integer"

map <- read.delim(file.path(dir, "molly_xenopus_map.txt"), row.names = 1, check.names = FALSE,
                  stringsAsFactors = FALSE)
dim(asv)
#[1] 66 53   (66 genus-level taxa x 52 samples + PCR control)

# ---- drop the PCR control BY NAME (metadata Control == "Yes"), light QC -------
pcr_ctrl <- rownames(map)[map$Control == "Yes"]                 # "PCRCtrl"
asv_bio  <- asv[, !(colnames(asv) %in% pcr_ctrl)]
map_bio  <- map[colnames(asv_bio), ]
stopifnot(identical(rownames(map_bio), colnames(asv_bio)))

# drop all-zero taxa
asv_bio <- asv_bio[rowSums(asv_bio) > 0, ]

# library sizes
libsize <- colSums(asv_bio)
summary(libsize)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   2875    4334    5097    7875    8233   25221
# (genus-collapsed table; the original script comment, 2832/4370/5104/7744/8088/25121, came from the ASV table)
sort(libsize)[1:5]
#  S37  S36  S43  S19  S38
# 2875 3105 3371 3392 3667

# ---- canonical treatment names (Gentamycin -> Gentamicin, FuskNoAB -> None) ---
map_bio$Antibiotic   <- recode_abx(map_bio$Antibiotic)
map_bio$Antibiotic   <- relevel(factor(map_bio$Antibiotic), ref = "Control")
map_bio$Fusarium     <- factor(map_bio$Fusarium, levels = c("No", "Yes"))
map_bio$FusariumDose <- as.numeric(map_bio$FusariumDose)
table(map_bio$Antibiotic, map_bio$Fusarium)

# ---- phyloseq object ---------------------------------------------------------
ps <- phyloseq(otu_table(asv_bio, taxa_are_rows = TRUE), sample_data(map_bio))

# ==============================================================================
# 1. ANTIBIOTIC EFFECT  (Fusarium = No subset; each antibiotic vs Control)
# ==============================================================================
ps_ab <- subset_samples(ps, Fusarium == "No")
ps_ab <- prune_taxa(taxa_sums(ps_ab) > 0, ps_ab)
sample_data(ps_ab)$Antibiotic <- droplevels(sample_data(ps_ab)$Antibiotic)

set.seed(515)
out_ab <- ancombc2(
  data        = ps_ab,
  fix_formula = "Antibiotic",
  group       = "Antibiotic",
  prv_cut     = 0.10,        # keep taxa in >=10% of samples
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

# taxa differing across the antibiotic panel
ab_global <- out_ab$res_global
subset(ab_global, diff_abn)

# per-antibiotic vs Control
ab_dunn <- out_ab$res_dunn

write_tsv(out_ab$res, file.path(OUT_DIR, "antibiotic_res.tsv"))
write_tsv(ab_global,  file.path(OUT_DIR, "antibiotic_global.tsv"))
write_tsv(ab_dunn,    file.path(OUT_DIR, "antibiotic_dunnett.tsv"))
# taxa with structural zeros (absent from >=1 group) are dropped from the tests above
# (struc_zero = TRUE); keep the presence/absence table, it is itself informative
write_tsv(out_ab$zero_ind, file.path(OUT_DIR, "antibiotic_structural_zeros.tsv"))

# ==============================================================================
# 2. FUSARIUM EFFECT UNDER ALL-ANTIBIOTICS  (AllAB: Fusarium Yes vs No)
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
# 3. FUSARIUM EFFECT, NO ANTIBIOTICS  (Fk-only [Yes] vs Control [No])
# ==============================================================================
ps_noab <- subset_samples(ps, Antibiotic %in% c("Control", "None"))   # "None" = FuskNoAB
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
# 4. FUSARIUM DOSE-RESPONSE within Fk-only animals (log10 dose as continuous)
# ==============================================================================
ps_dose <- subset_samples(ps, Antibiotic == "None")
ps_dose <- prune_taxa(taxa_sums(ps_dose) > 0, ps_dose)
sample_data(ps_dose)$logDose <- log10(sample_data(ps_dose)$FusariumDose)
table(fk_dose_factor(sample_data(ps_dose)$FusariumDose))   # all 5 levels incl. 10^5

out_dose <- ancombc2(
  data        = ps_dose,
  fix_formula = "logDose",
  prv_cut     = 0.10, lib_cut = 0,
  struc_zero  = FALSE, alpha = 0.05,
  n_cl = 1, verbose = TRUE
)
write_tsv(out_dose$res, file.path(OUT_DIR, "fusarium_dose_trend.tsv"))


# ==============================================================================
# FIGURES
# ==============================================================================
ALPHA <- 0.05

# ---- short taxon labels from SILVA strings ------------------------------------
# "d__Bacteria;p__..;f__Comamonadaceae;g__Comamonas;__"  -> italic Comamonas
# "...;f__Comamonadaceae;__;__"                          -> Comamonadaceae (unclassified genus)
# "...;o__Burkholderiales;__;__;__"                      -> Burkholderiales (unclassified family)
# genus + informative species epithet -> "Genus species" (keeps duplicates distinct)
short_taxon <- function(x) {
  get_rank <- function(s, r) {
    v <- sub(paste0("^", r, "__"), "", grep(paste0("^", r, "__"), s, value = TRUE))
    if (length(v) == 0 || v == "" || grepl("^uncultured", v)) NA_character_ else v
  }
  vapply(strsplit(x, ";"), function(s) {
    g <- get_rank(s, "g"); f <- get_rank(s, "f"); o <- get_rank(s, "o"); sp <- get_rank(s, "s")
    if (!is.na(sp) && grepl("_", sp) && !grepl("bacterium|metagenome", sp))
      sp <- sub("^[^_]+_", "", sp) else sp <- NA_character_
    if (!is.na(g)) {
      g <- sub("^Allorhizobium-Neorhizobium-Pararhizobium-", "", g)     # SILVA composite -> Rhizobium
      g <- sub("^Burkholderia-Caballeronia-Paraburkholderia$", "Burkholderia s.l.", g)
      if (!is.na(sp)) paste(g, sp) else g
    } else if (!is.na(f)) paste0(f, " (unclassified genus)")
    else if (!is.na(o)) paste0(o, " (unclassified family)")
    else x
  }, character(1))
}
# plotmath: italicise genus/species part, keep the qualifier upright
taxon_expr <- function(lab) {
  is_unc <- grepl(" \\(unclassified", lab)
  nm  <- sub(" \\(unclassified.*$", "", lab)
  qual <- sub("^.* \\((unclassified[^)]*)\\)$", "\\1", lab)
  ifelse(is_unc,
         sprintf('italic("%s")~"(%s)"', nm, qual),
         sprintf('italic("%s")', lab))
}
lab_parse <- function(x) parse(text = taxon_expr(x))

# ------------------------------------------------------------------------------
# Figure: Fusarium forest plot (LFC +/- SE), two contrasts side by side
#   left : Fk-only vs Control (no antibiotic; n = 17 vs 5)
#   right: AllAB + Fk vs AllAB alone (n = 5 vs 3)
# rows = taxa significant (q < 0.05) in either contrast; else top N by |LFC|
# ------------------------------------------------------------------------------
TOP_N <- 25
tidy_res <- function(file, contrast) {
  df <- read_tsv(file, show_col_types = FALSE); names(df)[1] <- "taxon"
  df %>%
    pivot_longer(matches("^(lfc|se|q)_"), names_to = c("stat", "term"),
                 names_pattern = "^(lfc|se|q)_(.*)$") %>%
    pivot_wider(id_cols = c(taxon, term), names_from = stat, values_from = value) %>%
    filter(term != "(Intercept)") %>%                # binary contrast -> one term remains
    mutate(contrast = contrast, sig = !is.na(q) & q < ALPHA)
}
fd <- bind_rows(
  tidy_res(file.path(OUT_DIR, "fusarium_noab.tsv"),         "No antibiotic:\nFusarium vs Control"),
  tidy_res(file.path(OUT_DIR, "fusarium_within_allab.tsv"), "AllAB background:\nAllAB + Fusarium vs AllAB")) %>%
  mutate(contrast = fct_inorder(contrast), label = short_taxon(taxon))
stopifnot(!anyDuplicated(unique(fd[c("taxon", "label")])$label))

keep_f <- fd %>% group_by(taxon) %>% summarise(any_sig = any(sig), mx = max(abs(lfc)), .groups = "drop")
fus_caption <- if (any(keep_f$any_sig)) {
  paste0("Taxa with q < ", ALPHA, " in either contrast (ANCOM-BC2). No point = not tested in that contrast.")
} else {
  paste0("No taxon reached q < ", ALPHA, " in either contrast; all ", min(nrow(keep_f), TOP_N),
         " taxa passing the prevalence filter are shown (ANCOM-BC2).\nNo point = taxon not tested in that contrast (prevalence filter / structural zero).")
}
keep_f <- if (any(keep_f$any_sig)) filter(keep_f, any_sig) else slice_max(keep_f, mx, n = TOP_N)
fd_plot <- fd %>% filter(taxon %in% keep_f$taxon)
# order rows by the no-antibiotic effect
row_ord <- fd_plot %>% filter(contrast == levels(contrast)[1]) %>% arrange(lfc) %>% pull(label)
row_ord <- c(setdiff(unique(fd_plot$label), row_ord), row_ord)
fd_plot <- fd_plot %>%
  mutate(label = factor(label, levels = row_ord),
         dir   = factor(ifelse(lfc > 0, "Higher with Fusarium", "Lower with Fusarium"),
                        levels = c("Lower with Fusarium", "Higher with Fusarium")),
         sigf  = factor(ifelse(sig, paste0("q < ", ALPHA), "n.s."), levels = c(paste0("q < ", ALPHA), "n.s.")))

p_fus <- ggplot(fd_plot, aes(lfc, label, colour = dir)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lfc - se, xmax = lfc + se), height = 0, linewidth = 0.4) +
  geom_point(aes(shape = sigf), size = 1.6, fill = "white", stroke = 0.6, show.legend = TRUE) +
  scale_colour_manual(values = c("Lower with Fusarium" = "#2166AC", "Higher with Fusarium" = "#B2182B"),
                      name = NULL) +
  scale_shape_manual(values = setNames(c(16, 21), c(paste0("q < ", ALPHA), "n.s.")), name = NULL, drop = FALSE) +
  scale_y_discrete(labels = lab_parse) +
  facet_wrap(~ contrast, nrow = 1) +
  guides(shape = guide_legend(override.aes = list(colour = "black"), order = 2),
         colour = guide_legend(order = 1)) +
  labs(x = "Log fold change (ln, \u00B1 SE)", y = NULL, caption = fus_caption) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0), panel.spacing.x = unit(4, "mm"),
        strip.text = element_text(face = "plain", size = 7.5),
        plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey25"))
h_fus <- 25 + 4.2 * nlevels(fd_plot$label)
save_pub(p_fus, file.path(FIG_DIR, "diff_taxa_fusarium.pdf"), width_mm = 180, height_mm = h_fus)

# ------------------------------------------------------------------------------
# Figure: antibiotic heatmap, LFC vs Control (Dunnett), significant taxa only
# ------------------------------------------------------------------------------
TOP_N <- 40
df <- read_tsv(file.path(OUT_DIR, "antibiotic_dunnett.tsv"), show_col_types = FALSE)
names(df)[1] <- "taxon"

d <- df %>%
  pivot_longer(matches("^(lfc|q)_"), names_to = c("stat", "term"),
               names_pattern = "^(lfc|q)_(.*)$") %>%
  pivot_wider(id_cols = c(taxon, term), names_from = stat, values_from = value) %>%
  filter(term != "(Intercept)") %>%
  mutate(term  = recode_abx(str_remove(term, "^Antibiotic")),   # tidy contrast labels
         sig   = !is.na(q) & q < ALPHA,
         label = short_taxon(taxon))

# choose taxa: significant in >= 1 antibiotic contrast
keep <- d %>% group_by(taxon, label) %>%
  summarise(any_sig = any(sig), n_sig = sum(sig), max_abs = max(abs(lfc), na.rm = TRUE),
            .groups = "drop")
if (any(keep$any_sig)) {
  keep_taxa <- filter(keep, any_sig)$taxon
} else {
  keep_taxa <- slice_max(keep, max_abs, n = TOP_N)$taxon
}
d <- filter(d, taxon %in% keep_taxa)

# rows clustered on LFC profile; columns in canonical antibiotic order
m <- d %>% dplyr::select(label, term, lfc) %>%
  pivot_wider(names_from = term, values_from = lfc) %>%
  column_to_rownames("label") %>% as.matrix()
m[is.na(m)] <- 0
row_ord <- if (nrow(m) > 2) rownames(m)[hclust(dist(m), method = "average")$order] else rownames(m)
d <- d %>% mutate(term  = factor(term, levels = intersect(abx_levels, unique(term))),
                  label = factor(label, levels = rev(row_ord)))

# omnibus (global) test across the antibiotic panel -> extra annotation column
gl <- read_tsv(file.path(OUT_DIR, "antibiotic_global.tsv"), show_col_types = FALSE)
names(gl)[1] <- "taxon"
gl <- gl %>% filter(taxon %in% keep_taxa) %>%
  mutate(label = factor(short_taxon(taxon), levels = levels(d$label)),
         txt   = paste0(ifelse(q_val < 0.001, "<0.001", formatC(q_val, format = "fg", digits = 2)),
                        ifelse(diff_abn & !passed_ss, "\u00A7", "")),
         face  = ifelse(diff_abn, "bold", "plain"))
n_ab <- nlevels(d$term)

lim <- ceiling(max(abs(d$lfc), na.rm = TRUE))
any_cell_sig <- any(d$sig)
ab_caption <- paste0(
  "Cells: ANCOM-BC2 LFC vs Control; * Dunnett-adjusted q < ", ALPHA,
  if (!any_cell_sig) " (none reached this threshold)" else "", ".\n",
  "Global q: omnibus test across antibiotics (bold = q < ", ALPHA,
  "; \u00A7 = failed pseudo-count sensitivity analysis).\n",
  "Taxa absent from any group (structural zeros) are excluded from testing.")

p_ab <- ggplot(d, aes(term, label, fill = lfc)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = filter(d, sig), shape = 8, size = 0.9, stroke = 0.4, colour = "black") +
  geom_text(data = gl, aes(x = n_ab + 0.9, y = label, label = txt, fontface = face),
            inherit.aes = FALSE, hjust = 0, size = 7 / .pt) +
  annotate("text", x = n_ab + 0.9, y = nlevels(d$label) + 0.5, label = "Global q",
           hjust = 0, vjust = -0.6, size = 7 / .pt, fontface = "italic") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-lim, lim), oob = squish,
                       name = "LFC vs\nControl (ln)", na.value = "grey92") +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.6)), position = "top") +
  scale_y_discrete(labels = lab_parse, expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, caption = ab_caption) +
  theme_pub() +
  theme(axis.text.x.top = element_text(angle = 40, hjust = 0, vjust = 0),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key.height = unit(5, "mm"), legend.key.width = unit(2.5, "mm"),
        plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey25"),
        plot.caption.position = "plot")
h_ab <- 45 + 5.0 * nlevels(d$label)
save_pub(p_ab, file.path(FIG_DIR, "diff_taxa_antibiotics.pdf"), width_mm = 130, height_mm = h_ab)
