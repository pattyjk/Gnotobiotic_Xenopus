########################################
## Xenopus laevis survivorship analysis
## Kaplan-Meier curves reconstructed from weekly census counts
########################################
dir <- "~/Documents/GitHub/Gnotobiotic_Xenopus/"

library(ggplot2)
library(tidyr)
library(dplyr)
library(survival)
library(patchwork)
source(file.path(dir, "pub_theme.R"))   # theme_pub(), recode_abx(), parse_group_header(), palettes, save_pub()

#---------------------------------------
# Read census counts (rows = census dates, columns = treatment groups)
# NB: the "Date " header has a trailing space; check.names = FALSE keeps "Fk 10^3" etc.
#---------------------------------------
raw <- read.delim(file.path(dir, "xenous_survivorship.txt"), header = TRUE,
                  check.names = FALSE, strip.white = TRUE)
names(raw) <- trimws(names(raw))
raw$date <- as.Date(raw$Date, format = "%m/%d/%y")
raw$Date <- NULL
group_cols <- setdiff(names(raw), "date")
stopifnot(!anyNA(raw$date), length(group_cols) == 16)

# group metadata from the column headers (canonical names, category, Fusarium dose)
grp <- parse_group_header(group_cols)
grp$group <- group_cols
grp$Label[grp$header == "ALL AB/AF"] <- "AllAB"
grp$LabelExpr[grp$header == "ALL AB/AF"] <- "AllAB"
day0 <- min(raw$date)
days <- as.numeric(raw$date - day0)

#---------------------------------------
# DATA CHECK: census counts must be non-increasing.
# 'Fk 10^6' reads 14, 10, 11, 8, 8 (an increase on 11/12/24) - almost certainly a
# data-entry problem. We do NOT edit the raw value; for the KM reconstruction a
# cumulative minimum is used (i.e. the 11/12 count is treated as 10 alive).
# Please check the lab notebook and correct xenous_survivorship.txt if needed.
#---------------------------------------
for (g in group_cols) {
  if (any(diff(raw[[g]]) > 0)) {
    warning(sprintf("Census count increases over time in group '%s' (%s); using cumulative minimum for KM reconstruction - check raw data.",
                    g, paste(raw[[g]], collapse = ", ")), call. = FALSE)
  }
}

# Percent alive at each census (descriptive table, raw counts)
pct <- raw
for (col in group_cols) pct[[col]] <- raw[[col]] / raw[[col]][1] * 100

#---------------------------------------
# Expand census counts to one row per animal
# deaths are assigned to the census date at which the animal was first missing;
# animals alive at the last census are right-censored there
#---------------------------------------
expand_to_individual <- function(group_name, counts, time_days) {
  alive <- cummin(counts)                       # see DATA CHECK above
  died  <- -diff(alive)
  rows  <- data.frame(group = group_name, time = time_days[-1], status = 1, count = died)
  rows  <- rbind(rows, data.frame(group = group_name, time = max(time_days),
                                  status = 0, count = alive[length(alive)]))
  rows[rows$count > 0, ]
}

indiv <- bind_rows(lapply(group_cols, function(g) expand_to_individual(g, raw[[g]], days)))
indiv_exp <- indiv[rep(seq_len(nrow(indiv)), indiv$count), ]
indiv_exp$count <- NULL
indiv_exp <- left_join(indiv_exp, grp[, c("group", "Antibiotic", "Category", "Dose", "Label")],
                       by = "group")

#---------------------------------------
# KM fit and statistics
#---------------------------------------
km_fit <- survfit(Surv(time, status) ~ group, data = indiv_exp)
cat("\n-- Kaplan-Meier: animals, deaths, median survival (days) --\n")
print(summary(km_fit)$table[, c("records", "events", "median")])

# Log-rank test: all 16 groups
cat("\n-- Log-rank test: all groups --\n")
lr_all <- survdiff(Surv(time, status) ~ group, data = indiv_exp)
print(lr_all)
p_all <- lr_all$pvalue

# Log-rank test: pooled treatment category (Control / Antibiotic / Fusarium / Antibiotic + Fusarium)
cat("\n-- Log-rank test: treatment category (groups pooled) --\n")
lr_cat <- survdiff(Surv(time, status) ~ Category, data = indiv_exp)
print(lr_cat)
p_cat <- lr_cat$pvalue

# Log-rank test within each category, including Control as reference
p_within <- sapply(c("Antibiotic", "Fusarium", "Antibiotic + Fusarium"), function(cc) {
  sub_df <- indiv_exp[indiv_exp$Category %in% c("Control", cc), ]
  survdiff(Surv(time, status) ~ group, data = sub_df)$pvalue
})
cat("\n-- Log-rank test within category (+ Control) --\n")
print(round(p_within, 3))

# Pairwise log-rank: each treatment vs Control (Benjamini-Hochberg adjusted)
cat("\n-- Pairwise log-rank vs Ctrl (Benjamini-Hochberg adjusted) --\n")
non_ctrl <- setdiff(group_cols, "Ctrl")
p_vals <- sapply(non_ctrl, function(g) {
  sub_df <- indiv_exp[indiv_exp$group %in% c("Ctrl", g), ]
  tryCatch(survdiff(Surv(time, status) ~ group, data = sub_df)$pvalue,
           error = function(e) NA)
})
p_adj <- p.adjust(p_vals, method = "BH")
result <- data.frame(
  group = non_ctrl,
  label = grp$Label[match(non_ctrl, grp$group)],
  p_raw = round(p_vals, 4),
  p_adj = round(p_adj, 4),
  sig   = ifelse(p_adj < 0.05, "*", "")
) %>% arrange(p_adj)
print(result)
# no impact of treatment on survival (no pairwise comparison with Control has p_adj < 0.05)

#---------------------------------------
# KM step-curve data (add the day-0 point; one step per census)
#---------------------------------------
km_df <- data.frame(time = km_fit$time, surv = km_fit$surv,
                    group = rep(sub("^group=", "", names(km_fit$strata)), km_fit$strata))
km_df <- bind_rows(data.frame(time = 0, surv = 1, group = group_cols), km_df) %>%
  arrange(group, time) %>%
  left_join(grp[, c("group", "Antibiotic", "Category", "Dose", "LabelExpr")], by = "group") %>%
  mutate(pct = 100 * surv)

ctrl_df <- km_df %>% filter(group == "Ctrl") %>% select(time, pct)
census_breaks <- days

km_panel <- function(cat_name, colour_values, title_text, p_value, linetypes = NULL) {
  d <- km_df %>% filter(Category == cat_name)
  d$LabelExpr <- factor(d$LabelExpr, levels = names(colour_values))
  p <- ggplot(d, aes(time, pct, colour = LabelExpr)) +
    geom_step(data = ctrl_df, aes(time, pct), inherit.aes = FALSE,
              colour = "grey55", linewidth = 1.6, alpha = 0.45) +
    geom_step(linewidth = 0.55, direction = "hv") +
    geom_point(data = d %>% filter(time == max(time)), shape = 3, size = 1.1,
               show.legend = FALSE) +                   # + = censored at last census
    annotate("text", x = 0.5, y = 3, hjust = 0, vjust = 0, size = 6 / .pt,
             label = paste0("log-rank (+ Control) ", fmt_p(p_value))) +
    scale_colour_manual(values = colour_values,
                        labels = function(x) parse(text = x), name = NULL) +
    scale_x_continuous(breaks = census_breaks, limits = c(0, max(days) + 0.5),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Days since 30 Oct 2024", y = "Survival (%)", subtitle = title_text) +
    theme_pub() +
    theme(legend.position = "bottom", legend.key.width = unit(4, "mm"),
          legend.text = element_text(size = 6.5),
          legend.margin = margin(0, 0, 0, 0),
          plot.subtitle = element_text(face = "bold", colour = "black", size = 8.5)) +
    guides(colour = guide_legend(ncol = 2, byrow = TRUE,
                                 override.aes = list(linewidth = 0.9)))
  p
}

ab_lab <- grp$LabelExpr[grp$Category == "Antibiotic"]
ab_abx <- grp$Antibiotic[grp$Category == "Antibiotic"]
ab_cols <- setNames(abx_cols[ab_abx], ab_lab)
ab_cols <- ab_cols[order(match(ab_abx, abx_levels))]

fk_lab  <- grp$LabelExpr[grp$Category == "Fusarium"]
fk_cols <- setNames(fk_dose_cols[paste0("10^", log10(grp$Dose[grp$Category == "Fusarium"]))], fk_lab)

abfk_lab  <- grp$LabelExpr[grp$Category == "Antibiotic + Fusarium"]
abfk_cols <- setNames(fk_dose_cols[paste0("10^", log10(grp$Dose[grp$Category == "Antibiotic + Fusarium"]))], abfk_lab)

pA <- km_panel("Antibiotic", ab_cols, "Antibiotic", p_within["Antibiotic"])
pB <- km_panel("Fusarium", fk_cols, "Fusarium (spores/mL)", p_within["Fusarium"]) +
  theme(axis.title.y = element_blank())
pC <- km_panel("Antibiotic + Fusarium", abfk_cols, "AllAB + Fusarium", p_within["Antibiotic + Fusarium"]) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(linewidth = 0.9))) +
  theme(axis.title.y = element_blank())

p_km <- (pA | pB | pC) +
  plot_annotation(
    caption = paste0("Grey band: Control (n = ", raw$Ctrl[1], "). All 16 groups: log-rank ",
                     fmt_p(p_all), "; pooled categories: log-rank ", fmt_p(p_cat),
                     ". Deaths assigned to census date."),
    theme = theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey25"))) &
  theme(axis.title.x = element_text(size = 8))

print(p_km)
save_pub(p_km, file.path(dir, "xenopus_survival_km.pdf"), width_mm = 180, height_mm = 85)
