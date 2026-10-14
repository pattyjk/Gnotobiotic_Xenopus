# Xenopus microbiome – Order-level relative abundance
# Main figure : group-mean composition per treatment
#               (A) antibiotic regimens, no Fusarium; (B) Fusarium challenge
# Supplement  : per-sample bars for the same two sets
# One order palette (top 10 + Other) is shared by every panel.

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(ggh4x)       # facet_nested(): merged "AllAB + Fk" / "Fk" strips above dose strips

# ---- paths: edit this one line to run elsewhere ----
dir <- "~/Documents/GitHub/Gnotobiotic_Xenopus/"
source(file.path(dir, "pub_theme.R"))   # theme_pub(), recode_abx(), fk_dose_factor(), save_pub()
fig_dir <- dir

# ---- Load data ----
tax  <- read.table(file.path(dir, "order_taxonomy.txt"), header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
meta <- read.table(file.path(dir, "molly_xenopus_map.txt"), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)

# Drop PCR control(s) by name (metadata Control == "Yes"), never by position
pcr_ctrl <- meta$SampleID[meta$Control == "Yes"]
tax  <- tax[!rownames(tax) %in% c(pcr_ctrl, "PCRCtrl"), , drop = FALSE]
meta <- meta[meta$Control != "Yes", ]
meta$FusariumDose <- as.numeric(meta$FusariumDose)

missing_in_tax <- setdiff(meta$SampleID, rownames(tax))
if (length(missing_in_tax))
  message("Samples in metadata but not in order_taxonomy.txt: ",
          paste(missing_in_tax, collapse = ", "))

# ---- Clean up column names: extract just the order label ----
extract_order <- function(x) {
  parts <- strsplit(x, ";")[[1]]
  o <- grep("^o__", parts, value = TRUE)
  if (length(o) == 0) {                      # unresolved at order level
    last <- sub("^[a-z]__", "", trimws(tail(parts, 1)))
    return(paste0("Unclassified ", last))
  }
  sub("^o__", "", trimws(o))
}
colnames(tax) <- sapply(colnames(tax), extract_order)
tax <- as.data.frame(t(rowsum(t(as.matrix(tax)), colnames(tax))))  # merge duplicate labels

# ---- Top 10 orders by total read count across all samples (unchanged rule) ----
order_totals <- colSums(tax)
top10 <- names(sort(order_totals, decreasing = TRUE))[1:10]
orders_all <- c(top10, "Other")

tax_top <- tax[, top10]
tax_top$Other <- rowSums(tax[, !colnames(tax) %in% top10, drop = FALSE])
tax_top <- tax_top / rowSums(tax_top) * 100           # relative abundance (%)
tax_top$SampleID <- rownames(tax_top)

# ---- Merge with metadata, define treatment groups ----
df <- inner_join(tax_top, meta[, c("SampleID", "Antibiotic", "Fusarium", "FusariumDose")],
                 by = "SampleID") %>%
  mutate(Antibiotic = recode_abx(Antibiotic),          # Gentamycin -> Gentamicin, FuskNoAB -> None
         Dose       = fk_dose_factor(FusariumDose))    # keeps 10^5 (was NA)

# Fusarium set: Control / AllAB / AllAB + Fk (by dose) / Fk (by dose)
fk_set <- df %>%
  filter(Fusarium == "Yes" | Antibiotic %in% c("Control", "AllAB")) %>%
  mutate(GroupKey = case_when(
    Fusarium == "No"      ~ Antibiotic,
    Antibiotic == "AllAB" ~ paste0("AllAB + Fk ", Dose),
    TRUE                  ~ paste0("Fk ", Dose)))
fk_dose_present <- levels(droplevels(fk_set$Dose))
fk_keys <- c("Control", "AllAB",
             paste0("AllAB + Fk ", fk_dose_present), paste0("Fk ", fk_dose_present))
fk_keys <- fk_keys[fk_keys %in% fk_set$GroupKey]
fk_set <- fk_set %>% mutate(GroupKey = factor(GroupKey, levels = fk_keys))
stopifnot(!any(is.na(fk_set$GroupKey)))

# "AllAB + Fk 10^3" -> plotmath with superscript dose
key_expr <- function(k) {
  e <- sub("^(.*) (10\\^[0-9])$", "\\1~\\2", k)
  gsub("AllAB \\+ Fk", "AllAB~'+'~Fk", e)
}

ab_set <- df %>% filter(Fusarium == "No") %>%
  mutate(GroupKey = factor(Antibiotic, levels = abx_levels[abx_levels %in% Antibiotic]))

to_long <- function(d) {
  d %>% pivot_longer(all_of(orders_all), names_to = "Order", values_to = "RelAbund") %>%
    mutate(Order = factor(Order, levels = rev(orders_all)))
}
long_abx  <- to_long(ab_set)
long_fusk <- to_long(fk_set)

# ---- Shared order palette (top 10 + Other), used by every panel ----
order_colors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
                  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#79706E", "#D3D3D3")
names(order_colors) <- orders_all
fill_scale <- scale_fill_manual(values = order_colors, breaks = orders_all,
                                name = "Order", drop = FALSE)
y_scale <- scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25))
y_coord <- coord_cartesian(ylim = c(0, 100), clip = "off")

# ---- Main figure: group means ----
group_mean <- function(lng) {
  lng %>% group_by(GroupKey, Order) %>%
    summarise(RelAbund = mean(RelAbund), n = n(), .groups = "drop")
}
mean_abx  <- group_mean(long_abx)
mean_fusk <- group_mean(long_fusk)
n_lab <- function(m) m %>% distinct(GroupKey, n)

mean_plot <- function(m, xlabels = waiver()) {
  ggplot(m, aes(GroupKey, RelAbund, fill = Order)) +
    geom_col(width = 0.8, colour = "white", linewidth = 0.1) +
    geom_text(data = n_lab(m), aes(GroupKey, 101.5, label = paste0("n = ", n)),
              inherit.aes = FALSE, size = 6 / .pt, vjust = 0) +
    fill_scale + y_scale + y_coord +
    scale_x_discrete(labels = xlabels) +
    labs(x = NULL, y = "Mean relative abundance (%)") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.margin = margin(10, 4, 4, 4))
}
pA <- mean_plot(mean_abx)
pB <- mean_plot(mean_fusk, xlabels = function(k) parse(text = key_expr(k)))

p_main <- (pA + pB + plot_layout(widths = c(nlevels(mean_abx$GroupKey),
                                            nlevels(mean_fusk$GroupKey)),
                                 guides = "collect")) &
  theme(legend.position = "right")
save_pub(p_main, file.path(fig_dir, "xenopus_orders_group_means"), width_mm = 180, height_mm = 85)

# ---- Supplementary: per-sample bars ----
sample_order <- function(lng) {        # within a group, order samples by the top order
  lng %>% filter(Order == top10[1]) %>% arrange(GroupKey, desc(RelAbund)) %>% pull(SampleID)
}
per_sample_plot <- function(lng, facet) {
  lng <- lng %>% mutate(SampleID = factor(SampleID, levels = sample_order(lng)))
  ggplot(lng, aes(SampleID, RelAbund, fill = Order)) +
    geom_col(width = 0.9, colour = "white", linewidth = 0.1) +
    facet +
    fill_scale + y_scale + y_coord +
    labs(x = "Sample", y = "Relative abundance (%)") +
    theme_pub() +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          panel.spacing.x = unit(1.2, "mm"),
          strip.text   = element_text(size = 7, face = "bold", margin = margin(1.5, 0, 1.5, 0)))
}

# Antibiotics: one strip per regimen; the long name in the 2-animal facet is hyphenated
abx_strip <- function(x) ifelse(x == "Sulfamethazine", "Sulfa-\nmethazine", x)
p1 <- per_sample_plot(long_abx,
        facet_grid(~ GroupKey, scales = "free_x", space = "free_x",
                   labeller = as_labeller(abx_strip))) +
  theme(strip.clip = "off")            # let "Sulfa-/methazine" extend into the panel gap
save_pub(p1, file.path(fig_dir, "xenopus_orders_antibiotics"), width_mm = 180, height_mm = 80)

# Fusarium: nested strips -> treatment on top, dose (10^2 .. 10^5) underneath
long_fusk_nest <- long_fusk %>%
  mutate(key = as.character(GroupKey),
         Treatment = case_when(grepl("^AllAB \\+ Fk", key) ~ "AllAB~'+'~Fk",
                               grepl("^Fk", key)            ~ "Fk",
                               TRUE                         ~ key),
         Treatment = factor(Treatment, levels = c("Control", "AllAB", "AllAB~'+'~Fk", "Fk")),
         DoseLab   = ifelse(grepl("Fk", key), sub("^.* (10\\^[0-9])$", "\\1", key), "''"),
         DoseLab   = factor(DoseLab, levels = c("''", paste0("10^", 2:6))))
p2 <- per_sample_plot(long_fusk_nest,
        facet_nested(~ Treatment + DoseLab, scales = "free_x", space = "free_x",
                     labeller = label_parsed,
                     nest_line = element_line(colour = "black", linewidth = 0.3)))
save_pub(p2, file.path(fig_dir, "xenopus_orders_fusarium"), width_mm = 180, height_mm = 85)

# ---- Top 10 orders and group sizes for reference ----
print(data.frame(Order = top10, TotalReads = order_totals[top10]))
print(n_lab(mean_abx)); print(n_lab(mean_fusk))
