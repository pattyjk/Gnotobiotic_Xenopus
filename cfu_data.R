########################################
## Xenopus laevis CFU data
## one plate count per treatment group per census date (no replication)
########################################
dir <- "~/Documents/GitHub/Gnotobiotic_Xenopus/"

library(ggplot2)
library(tidyr)
library(dplyr)
source(file.path(dir, "pub_theme.R"))   # theme_pub(), parse_group_header(), palettes, save_pub()

# Data (rows = dates, columns = groups; "Date " header has a trailing space)
raw <- read.delim(file.path(dir, "cfU-data.txt"), header = TRUE, check.names = FALSE,
                  strip.white = TRUE)
names(raw) <- trimws(names(raw))
raw$date <- as.Date(raw$Date, format = "%m/%d/%y")
raw$Date <- NULL
group_cols <- setdiff(names(raw), "date")
stopifnot(!anyNA(raw$date), length(group_cols) == 16)

grp <- parse_group_header(group_cols)
grp$group <- group_cols

# Reshape to long format
long <- pivot_longer(raw, cols = all_of(group_cols), names_to = "group", values_to = "CFU") %>%
  left_join(grp[, c("group", "Antibiotic", "Category", "Dose", "Label", "LabelExpr")], by = "group") %>%
  mutate(day = as.numeric(date - min(date)))

#---------------------------------------
# Statistical comparison: each group vs Control, paired by census date
# (Wilcoxon signed-rank, BH-adjusted; n = 5 dates per group)
#---------------------------------------
ctrl_vals <- long$CFU[long$group == "Ctrl"]
non_ctrl  <- setdiff(group_cols, "Ctrl")
p_vals <- sapply(non_ctrl, function(g) {
  grp_vals <- long$CFU[long$group == g]
  wilcox.test(grp_vals, ctrl_vals, paired = TRUE, exact = FALSE)$p.value
})
p_adj <- p.adjust(p_vals, method = "BH")
stat_result <- data.frame(
  group       = non_ctrl,
  label       = grp$Label[match(non_ctrl, grp$group)],
  p_raw       = round(p_vals, 4),
  p_adj_BH    = round(p_adj, 4),
  significant = ifelse(p_adj < 0.05, "*", "")
) %>% arrange(p_adj_BH)
print(stat_result)
# no sig difference from control, so no change in bacterial abundance

#---------------------------------------
# Supplementary figure: CFU over time (log10), faceted by category,
# Control repeated in every facet as a grey reference
#---------------------------------------
facet_levels <- c("Antibiotic", "Fusarium", "Antibiotic + Fusarium")
ctrl <- long %>% filter(Category == "Control") %>% select(day, CFU)
ctrl <- tidyr::crossing(ctrl, Facet = factor(facet_levels, levels = facet_levels))

trt <- long %>% filter(Category != "Control") %>%
  mutate(Facet = factor(as.character(Category), levels = facet_levels))

# colour per group: antibiotics -> abx_cols; Fusarium groups -> fk_dose_cols;
# AllAB + Fusarium distinguished by dashed line
col_key <- ifelse(trt$Category == "Antibiotic", as.character(trt$Antibiotic),
                  paste0("10^", round(log10(trt$Dose))))
trt$col_key <- col_key
key_vals <- c(abx_cols[setdiff(abx_levels, "Control")], fk_dose_cols)
key_labs <- c(setdiff(abx_levels, "Control"), paste0("Fk~", names(fk_dose_cols)))
names(key_labs) <- names(key_vals)
present <- intersect(names(key_vals), unique(trt$col_key))

facet_lab <- as_labeller(c("Antibiotic" = "Antibiotic",
                           "Fusarium" = "Fusarium",
                           "Antibiotic + Fusarium" = "AllAB + Fusarium"))

p_cfu <- ggplot(trt, aes(day, CFU)) +
  geom_line(data = ctrl, colour = "grey60", linewidth = 1.2, alpha = 0.6) +
  geom_point(data = ctrl, colour = "grey60", size = 1.4, alpha = 0.8) +
  geom_line(aes(colour = col_key, group = group,
                linetype = Category == "Antibiotic + Fusarium"), linewidth = 0.5) +
  geom_point(aes(colour = col_key), size = 1.1) +
  facet_wrap(~ Facet, nrow = 1, labeller = facet_lab) +
  scale_colour_manual(values = key_vals[present], breaks = present,
                      labels = parse(text = key_labs[present]), name = NULL) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
  scale_y_log10(limits = c(50, 1200), breaks = c(50, 100, 200, 500, 1000)) +
  scale_x_continuous(breaks = sort(unique(long$day))) +
  labs(x = "Days since 30 Oct 2024", y = "Bacterial load (CFU)",
       caption = paste0("Grey: Control. One plate count per group per date (no replication); ",
                        "dashed lines: AllAB pretreatment + Fusarium. ",
                        "Paired Wilcoxon vs Control, all BH-adjusted p >= ",
                        sprintf("%.2f", min(p_adj)), ".")) +
  theme_pub() +
  theme(legend.position = "right", legend.key.width = unit(4, "mm"),
        legend.text = element_text(size = 6.5),
        plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey25"),
        plot.caption.position = "plot",
        panel.spacing.x = unit(4, "mm")) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(linewidth = 0.8)))

print(p_cfu)
save_pub(p_cfu, file.path(dir, "CFU_faceted.pdf"), width_mm = 180, height_mm = 70)
