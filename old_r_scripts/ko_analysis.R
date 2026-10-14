library(phyloseq)
library(ANCOMBC)
library(ALDEx2)
library(tidyverse)

KO_FILE  <- "~/Documents/GitHub/Gnotobiotic_Xenopus/picrust_out/KO_metagenome_out/pred_metagenome_unstrat.tsv"
MAP_FILE <- "~/Documents/GitHub/Gnotobiotic_Xenopus/molly_xenopus_map.txt"
OUT_DIR  <- "~/Documents/GitHub/Gnotobiotic_Xenopus/ko_results"
dir.create(OUT_DIR, showWarnings = FALSE)


#Read + round + align
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
#FUSARIUM EFFECT under all-antibiotics (AllAB: Yes vs No)
# ==============================================================================
ps_allab <- subset_samples(ps, Antibiotic == "AllAB")
ps_allab <- prune_taxa(taxa_sums(ps_allab) > 0, ps_allab)

out_fus_allab <- ancombc2(
  data = ps_allab, fix_formula = "Fusarium", group = "Fusarium",
  prv_cut = 0.10, lib_cut = 0, struc_zero = TRUE, neg_lb = TRUE,
  alpha = 0.05, n_cl = 1, verbose = TRUE)
write_tsv(out_fus_allab$res, file.path(OUT_DIR, "fusarium_within_allab.tsv"))


# ==============================================================================
#FUSARIUM EFFECT no antibiotics (FuskNoAB[Yes] vs Control[No])
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

######################
#pull KOs
######################

kos <- c(
  "K00098","K00869","K00949","K00980","K01034","K01035","K01597","K01728",
  "K01771","K01815","K01820","K01880","K02022","K02351","K02526","K03921",
  "K03931","K03932","K04516","K04518","K05593","K06442","K06881","K06974",
  "K06998","K07086","K07456","K07710","K07714","K08678","K09749","K09992",
  "K10118","K10119","K10548","K10924","K10926","K10927","K12267","K12276",
  "K12278","K12279","K12282","K12283","K12284","K12285","K12286","K13018",
  "K13650","K15923","K16138","K16898","K16899","K18384","K18552","K18682",
  "K18785","K19419","K19572","K19701"
)

# --- 1. definitions: KEGG /list in chunks of 10 (API limit) -------------------
fetch_chunk <- function(ids) {
  url <- paste0("https://rest.kegg.jp/list/", paste(ids, collapse = "+"))
  txt <- tryCatch(readLines(url, warn = FALSE), error = function(e) character(0))
  if (length(txt) == 0) return(NULL)
  p <- strsplit(txt, "\t", fixed = TRUE)
  data.frame(taxon = vapply(p, `[`, "", 1),
             defn  = vapply(p, function(x) if (length(x) > 1) x[2] else NA_character_, ""),
             stringsAsFactors = FALSE)
}
chunks <- split(kos, ceiling(seq_along(kos) / 10))
res <- do.call(rbind, lapply(chunks, function(ch) { Sys.sleep(0.34); fetch_chunk(ch) }))
res$name <- sub(";.*$", "", res$defn)
res$desc <- trimws(sub("^[^;]*;\\s*", "", res$defn))
res$desc[res$desc == res$name] <- res$defn[res$desc == res$name]

# --- 2. KO hierarchy levels: parse the br:ko00001 BRITE tree ------------------
parse_brite <- function(kos) {
  hx <- tryCatch(readLines(url("https://rest.kegg.jp/get/br:ko00001"), warn = FALSE),
                 error = function(e) character(0))
  if (!length(hx)) { warning("could not fetch br:ko00001"); return(NULL) }
  strip <- function(s) trimws(gsub("<[^>]+>", "", s))     # drop any HTML tags
  A <- B <- C <- Acode <- NA_character_
  out <- list(); want <- kos
  for (ln in hx) {
    if (!nzchar(ln)) next
    tag <- substr(ln, 1, 1)
    if (!tag %in% c("A", "B", "C", "D")) next             # skip #, !, + comment lines
    body <- strip(substring(ln, 2))
    if (tag == "A") {
      Acode <- if (grepl("^[0-9]{5}", body)) sub("^([0-9]{5}).*", "\\1", body) else NA_character_
      A <- sub("^[0-9]{5}\\s*", "", body); B <- C <- NA_character_
    } else if (tag == "B") {
      if (nzchar(body)) { B <- sub("^[0-9]{5}\\s*", "", body); C <- NA_character_ }
    } else if (tag == "C") {
      C <- sub("^[0-9]{5}\\s*", "", body)
    } else if (tag == "D") {
      ko <- sub("^(K[0-9]{5}).*", "\\1", body)
      if (grepl("^K[0-9]{5}$", ko) && ko %in% want) {
        pid   <- if (grepl("\\[PATH:ko[0-9]+\\]", C)) sub(".*\\[PATH:(ko[0-9]+)\\].*", "\\1", C) else NA_character_
        pname <- trimws(sub("\\s*\\[(PATH|BR):[^]]+\\]", "", C))
        out[[length(out) + 1]] <- data.frame(
          taxon = ko, A_code = Acode, level_A = A, level_B = B,
          pathway_id = pid, pathway_name = pname, stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(out)) return(NULL)
  unique(do.call(rbind, out))
}
brite <- parse_brite(kos)

# --- 3. collapse to one row per KO (pathway maps only: drop 09180/09190) ------
collapse_by <- function(df, col, key)
  tapply(df[[col]], df[[key]], function(x) paste(unique(x[!is.na(x) & nzchar(x)]), collapse = "; "))

ko_desc <- res[, c("taxon", "name", "desc")]
if (!is.null(brite)) {
  pw <- brite[is.na(brite$A_code) | !(brite$A_code %in% c("09180", "09190")), ]
  ut <- unique(pw$taxon)
  summ <- data.frame(
    taxon    = ut,
    level_A  = collapse_by(pw, "level_A", "taxon")[ut],
    level_B  = collapse_by(pw, "level_B", "taxon")[ut],
    pathways = collapse_by(pw, "pathway_name", "taxon")[ut],
    stringsAsFactors = FALSE)
  ko_desc <- merge(ko_desc, summ, by = "taxon", all.x = TRUE)
  write.table(brite, "ko_brite_long.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  message("wrote ko_brite_long.tsv (", nrow(brite), " KO x pathway rows)")
}

# --- 4. write + report --------------------------------------------------------
missing <- setdiff(kos, res$taxon)
if (length(missing)) message("No KEGG record for: ", paste(missing, collapse = ", "))
no_level <- setdiff(kos, ko_desc$taxon[!is.na(ko_desc$level_A)])
if (length(no_level)) message("No pathway-map level (BRITE-only/unclassified): ",
                              paste(no_level, collapse = ", "))

#write.table(ko_desc, "~/Documents/GitHub/Gnotobiotic_Xenopus/ko_desc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# Heatmap (KOs): antibiotic effects, lfc vs Control (Dunnett) — labelled
# ==============================================================================

# ==============================================================================
# Heatmap (KOs): antibiotic effects vs Control, FACETED by KEGG level-2 category
# Standalone. Reads ko_results/antibiotic_dunnett.tsv + ko_desc.tsv -> ko_figures/.
# ==============================================================================
library(tidyverse)

IN <- "~/Documents/GitHub/Gnotobiotic_Xenopus/ko_results/antibiotic_dunnett.tsv" 
OUT <- "~/Documents/GitHub/Gnotobiotic_Xenopus/figures/heatmap_antibiotic_lfc"

# ---- config ----
ALPHA <- 0.05; LFC_MIN <- 1; MAX_ROWS <- 60
FACET_COL <- "level_B"          # "level_B" (KEGG level 2) or "pathways" (level 3)

df <- read_tsv(IN, show_col_types = FALSE); names(df)[1] <- "taxon"
d <- df %>%
  pivot_longer(matches("^(lfc|q)_"), names_to = c("stat","term"),
               names_pattern = "^(lfc|q)_(.*)$") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = str_remove(term, "^Antibiotic"),
         sig = !is.na(q) & q < ALPHA, label = taxon, cat = "Unclassified")

# ---- labels + level-2 category from ko_desc.tsv ------------------------------
if (file.exists("~/Documents/GitHub/Gnotobiotic_Xenopus/ko_desc.tsv")) {
  kd <- read_tsv("ko_desc.tsv", show_col_types = FALSE)
  d <- left_join(d, kd, by = "taxon") %>%
    mutate(sym   = ifelse(is.na(name) | name == taxon, NA, sub(",.*$", "", name)),
           label = ifelse(is.na(sym), taxon, paste0(sym, " (", taxon, ")")))
  if (FACET_COL %in% names(d)) {
    d <- mutate(d, cat = trimws(sub(";.*$", "", .data[[FACET_COL]])),  # primary category
                cat = ifelse(is.na(cat) | cat == "", "Unclassified", cat))
  }
}

# ---- keep KOs with a real signal --------------------------------------------
rank_tbl <- d %>% group_by(taxon, label, cat) %>%
  summarise(score = max(ifelse(sig & abs(lfc) >= LFC_MIN, abs(lfc), 0), na.rm = TRUE),
            .groups = "drop") %>%
  filter(score > 0) %>% slice_max(score, n = MAX_ROWS)
if (nrow(rank_tbl) == 0) stop("No KO passes sig & |lfc| >= ", LFC_MIN, " — relax filters.")
d <- semi_join(d, rank_tbl, by = "taxon")

# ---- column order by clustering; facet & row order by effect -----------------
m <- d %>% dplyr::select(label, term, lfc) %>%      # dplyr:: avoids AnnotationDbi mask
  pivot_wider(names_from = term, values_from = lfc) %>%
  column_to_rownames("label") %>% as.matrix()
m0 <- m; m0[is.na(m0)] <- 0
col_ord <- if (ncol(m0) > 2) colnames(m0)[hclust(dist(t(m0)))$order] else colnames(m0)

ord <- d %>% group_by(cat) %>% mutate(gmax = max(abs(lfc), na.rm = TRUE)) %>% ungroup()
cat_levels <- ord %>% distinct(cat, gmax) %>% arrange(desc(gmax)) %>% pull(cat)
row_levels <- ord %>% group_by(cat, label) %>% summarise(rm = mean(lfc, na.rm = TRUE), .groups = "drop") %>%
  mutate(cat = factor(cat, levels = cat_levels)) %>% arrange(cat, desc(rm)) %>% pull(label)
d <- mutate(d, term  = factor(term, levels = col_ord),
            cat   = factor(cat, levels = cat_levels),
            label = factor(label, levels = row_levels))

lim <- max(abs(d$lfc), na.rm = TRUE)
p <- ggplot(d, aes(term, label, fill = lfc)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(data = filter(d, sig), aes(label = "*"), vjust = 0.75, size = 5) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       limits = c(-lim, lim), name = "lfc vs\nControl", na.value = "grey92") +
  facet_grid(cat ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = label_wrap_gen(18)) +
  labs(title = "KO: antibiotic effects by KEGG level-2 category",
       subtitle = paste0("|lfc| >= ", LFC_MIN, ", q < ", ALPHA, ", top ", nrow(rank_tbl), " KOs"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
        strip.background.y = element_rect(fill = "grey92", colour = NA),
        panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.4),
        panel.spacing.y = unit(8, "pt"))

h <- max(5, 0.22 * length(unique(d$label)) + 0.25 * length(unique(d$cat)) + 1.5)
ggsave(paste0(OUT, ".pdf"), p, width = 10, height = h)
p

# ==============================================================================
# Forest plot (KOs): KO: Fusarium effect (no antibiotics), FACETED by KEGG level-2 category
# ==============================================================================
library(tidyverse)

IN <- "~/Documents/GitHub/Gnotobiotic_Xenopus/ko_results/fusarium_noab.tsv"; 
OUT <- "~/Documents/GitHub/Gnotobiotic_Xenopus/figures/forest_fusarium_noab"
ALPHA <- 0.05; TOP_N <- 30
FACET_COL <- "level_B"          # "level_B" (KEGG level 2) or "pathways" (level 3)

df <- read_tsv(IN, show_col_types = FALSE); names(df)[1] <- "taxon"
d <- df %>%
  pivot_longer(matches("^(lfc|se|q)_"), names_to = c("stat","term"),
               names_pattern = "^(lfc|se|q)_(.*)$") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = !is.na(q) & q < ALPHA, label = taxon, cat = "Unclassified")

# labels + level-2 category from ko_desc.tsv
if (file.exists("~/Documents/GitHub/Gnotobiotic_Xenopus/ko_desc.tsv")) {
  kd <- read_tsv("ko_desc.tsv", show_col_types = FALSE)
  d <- left_join(d, kd, by = "taxon") %>%
    mutate(sym   = ifelse(is.na(name) | name == taxon, NA, sub(",.*$", "", name)),
           label = ifelse(is.na(sym), taxon, paste0(sym, " (", taxon, ")")))
  if (FACET_COL %in% names(d)) {
    d <- mutate(d, cat = trimws(sub(";.*$", "", .data[[FACET_COL]])),
                cat = ifelse(is.na(cat) | cat == "", "Unclassified", cat))
  }
}

sig_d <- filter(d, sig)
if (nrow(sig_d) > 0) { plot_d <- sig_d; subtxt <- paste0("KOs at q < ", ALPHA)
} else { plot_d <- slice_max(d, abs(lfc), n = TOP_N)
subtxt <- paste0("no KO reached q < ", ALPHA, "; top ", nrow(plot_d), " by |lfc|") }

# facet order by group mean lfc; rows by lfc
cat_levels <- plot_d %>% group_by(cat) %>% summarise(m = mean(lfc), .groups = "drop") %>%
  arrange(m) %>% pull(cat)
plot_d <- plot_d %>% mutate(cat = factor(cat, levels = cat_levels),
                            label = fct_reorder(label, lfc))

p <- ggplot(plot_d, aes(lfc, label, colour = lfc > 0)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_errorbarh(aes(xmin = lfc - se, xmax = lfc + se), height = 0) +
  geom_point(aes(alpha = sig), size = 2.6) +
  scale_colour_manual(values = c(`TRUE`="#b2182b", `FALSE`="#2166ac"),
                      labels = c("lower with Fusarium","higher with Fusarium"), name = NULL) +
  scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.35), guide = "none") +
  facet_grid(cat ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = label_wrap_gen(18)) +
  labs(title = "KO: Fusarium effect (no antibiotics)", subtitle = subtxt, x = "log fold change", y = NULL) +
  theme_bw(base_size = 12) +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
        panel.spacing.y = unit(2, "pt"))

h <- max(4, 0.24 * nrow(plot_d) + 0.25 * length(unique(plot_d$cat)) + 1.5)
ggsave(paste0(OUT, ".pdf"), p, width = 8.5, height = h)
