# Xenopus microbiome – Order-level stacked bar plot
# Plot 1: Antibiotic treatments, no Fusarium
# Plot 2: Fusarium treatments only

library(ggplot2)
library(tidyr)
library(dplyr)

# ── Load data ─────────────────────────────────────────────────────────────────
tax  <- read.table("~/Documents/GitHub/Gnotobiotic_Xenopus/order_taxonomy.txt",  header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
meta <- read.table("~/Documents/GitHub/Gnotobiotic_Xenopus/molly_xenopus_map.txt", header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)

# ── Clean up column names: extract just the order label ──────────────────────
extract_order <- function(x) {
  parts <- strsplit(x, ";")[[1]]
  o <- grep("^o__", parts, value = TRUE)
  if (length(o) == 0) return(x)
  sub("^o__", "", trimws(o))
}
colnames(tax) <- sapply(colnames(tax), extract_order)

# ── Identify top 10 orders by total read count across all samples ─────────────
order_totals <- colSums(tax)
top10 <- names(sort(order_totals, decreasing = TRUE))[1:10]

# ── Collapse non-top-10 into "Other" ─────────────────────────────────────────
tax_top <- tax[, top10]
tax_top$Other <- rowSums(tax[, !colnames(tax) %in% top10])
tax_top$SampleID <- rownames(tax_top)

# ── Convert to relative abundance (%) ────────────────────────────────────────
tax_top <- tax_top %>%
  mutate(total = rowSums(across(-SampleID))) %>%
  mutate(across(-c(SampleID, total), ~ . / total * 100)) %>%
  select(-total)

# ── Merge with metadata ───────────────────────────────────────────────────────
df <- left_join(tax_top, meta[, c("SampleID","Antibiotic","Fusarium","FusariumDose")],
                by = "SampleID") %>%
  filter(!is.na(Antibiotic))

# ── Reshape to long format ────────────────────────────────────────────────────
orders_all <- c(top10, "Other")
long <- df %>%
  pivot_longer(cols = all_of(orders_all),
               names_to = "Order", values_to = "RelAbund") %>%
  mutate(
    Order    = factor(Order, levels = rev(c(top10, "Other"))),
    SampleID = factor(SampleID, levels = sort(unique(SampleID)))
  )

# ── Color palette ─────────────────────────────────────────────────────────────
order_colors <- c(
  "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
  "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
  "#888888"
)
names(order_colors) <- c(top10, "Other")

# ── Shared theme ──────────────────────────────────────────────────────────────
bar_theme <- theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey92", color = NA),
    legend.position  = "right",
    legend.key.size  = unit(0.5, "cm"),
    panel.spacing    = unit(0.6, "lines"),
    plot.title       = element_text(face = "bold")
  )

# ── Plot 1: Antibiotic treatments, no Fusarium ───────────────────────────────
long_abx <- long %>%
  filter(Fusarium == "No") %>%
  mutate(
    Antibiotic = factor(Antibiotic,
                        levels = c("Control","Trimethoprim","Gentamycin",
                                   "Ampicillin","PSN","AllAB","Sulfamethazine"))
  )

p1 <- ggplot(long_abx, aes(x = SampleID, y = RelAbund, fill = Order)) +
  geom_bar(stat = "identity", width = 0.85) +
  facet_wrap(~ Antibiotic, scales = "free_x", ncol = 4) +
  scale_fill_manual(values = order_colors,
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Relative abundance of bacterial orders — antibiotic treatments (no Fusarium)",
    subtitle = "Top 10 orders; remaining taxa collapsed to 'Other'",
    x        = "Sample", y = "Relative abundance (%)", fill = "Order"
  ) +
  bar_theme

print(p1)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/xenopus_orders_antibiotics.pdf", p1, width = 16, height = 7, dpi = 150)


# ── Plot 2: Fusarium treatments only ─────────────────────────────────────────
# Grid layout: rows = antibiotic (No AB vs AllAB), columns = Fusarium dose
long_fusk <- long %>%
  filter(Fusarium == "Yes") %>%
  mutate(
    DoseLabel = factor(paste0("Fk ", as.integer(FusariumDose)),
                       levels = paste0("Fk ", c(100, 1000, 10000, 100000))),
    ABLabel   = ifelse(Antibiotic == "FuskNoAB", "No Antibiotic", Antibiotic)
  )

p2 <- ggplot(long_fusk, aes(x = SampleID, y = RelAbund, fill = Order)) +
  geom_bar(stat = "identity", width = 0.85) +
  facet_grid(ABLabel ~ DoseLabel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = order_colors,
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Relative abundance of bacterial orders — Fusarium treatments",
    subtitle = "Columns = Fusarium dose; rows = antibiotic treatment",
    x        = "Sample", y = "Relative abundance (%)", fill = "Order"
  ) +
  bar_theme +
  theme(strip.text.y = element_text(angle = 0))

print(p2)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/xenopus_orders_fusarium.pdf", p2, width = 16, height = 7, dpi = 150)


#Top 10 orders for reference ────────────────────────────────────────
print(data.frame(Order = top10, TotalReads = order_totals[top10]))
