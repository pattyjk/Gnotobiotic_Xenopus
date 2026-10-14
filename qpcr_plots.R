############################################
## 16S qPCR: bacterial load on skin swabs
## Gnotobiotic Xenopus project
############################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- paths ----
dir <- "~/Documents/GitHub/Gnotobiotic_Xenopus/"
source(file.path(dir, "pub_theme.R"))
qpcr_file <- file.path(dir, "qpcr_data.txt")
map_file  <- file.path(dir, "molly_xenopus_map.txt")
out_dir   <- dir

# ---- data ----
qpcr <- read.delim(qpcr_file)
meta <- read.delim(map_file)

# NTC wells: Content == "NTC" (Cq NA = no amplification)
ntc <- filter(qpcr, Content == "NTC")
cat("NTC wells:", nrow(ntc), "| amplified:", sum(!is.na(ntc$Cq)), "\n")
qpcr <- filter(qpcr, Content == "Unkn", SampleID %in% meta$SampleID)

grp_levels <- c("Control", "Ampicillin", "AllAB", "PSN", "Sulfamethazine",
                "Gentamicin", "Trimethoprim", "Fk", "AllAB + Fk")
qpcr <- qpcr %>%
  select(-Antibiotic) %>%
  left_join(select(meta, SampleID, Antibiotic, Fusarium, FusariumDose), by = "SampleID") %>%
  mutate(Antibiotic = recode_abx(Antibiotic),
         Group = case_when(Fusarium == "Yes" & Antibiotic == "AllAB" ~ "AllAB + Fk",
                           Fusarium == "Yes" ~ "Fk",
                           TRUE ~ Antibiotic),
         Group = factor(Group, levels = grp_levels))

# ---- technical replicate agreement ----
reps <- qpcr %>%
  group_by(SampleID) %>% arrange(Plate, Well, .by_group = TRUE) %>%
  mutate(rep = paste0("rep", row_number())) %>% ungroup() %>%
  select(SampleID, rep, LogCopy, Plate) %>%
  pivot_wider(names_from = rep, values_from = c(LogCopy, Plate)) %>%
  mutate(layout = ifelse(Plate_rep1 == Plate_rep2, "Same plate", "Different plates"))

fit <- anova(lm(LogCopy ~ factor(SampleID), data = qpcr))
icc <- (fit$`Mean Sq`[1] - fit$`Mean Sq`[2]) / (fit$`Mean Sq`[1] + fit$`Mean Sq`[2])
mad <- mean(abs(reps$LogCopy_rep1 - reps$LogCopy_rep2))
r_by <- reps %>% group_by(layout) %>%
  summarise(n = n(), r = cor(LogCopy_rep1, LogCopy_rep2)) %>%
  mutate(lab = sprintf("%s (n = %d, r = %.2f)", layout, n, r))
lims <- range(c(reps$LogCopy_rep1, reps$LogCopy_rep2)) + c(-0.05, 0.05)

p_rep <- ggplot(reps, aes(LogCopy_rep1, LogCopy_rep2, colour = layout)) +
  geom_abline(linetype = 2, colour = "grey60", linewidth = 0.4) +
  geom_point(size = 1.4, alpha = 0.9) +
  coord_equal(xlim = lims, ylim = lims) +
  scale_colour_manual(values = c("Same plate" = "#332288", "Different plates" = "#CC6677"),
                      labels = setNames(r_by$lab, r_by$layout), name = NULL) +
  annotate("text", x = lims[2], y = lims[1], hjust = 1, vjust = 0, size = 2.6,
           label = sprintf("ICC = %.2f\nmean |\u0394| = %.2f log", icc, mad)) +
  labs(x = expression(Replicate~1~(log[10]~"16S copies")),
       y = expression(Replicate~2~(log[10]~"16S copies"))) +
  theme_pub() +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.text = element_text(size = 6.5))

# ---- per-sample load (mean of technical replicates) ----
load <- qpcr %>%
  group_by(SampleID, Group, Antibiotic, Fusarium, FusariumDose) %>%
  summarise(LogCopy = mean(LogCopy), .groups = "drop")
write.csv(load, file.path(out_dir, "qpcr_per_sample.csv"), row.names = FALSE)

# ---- statistics ----
ab_only <- filter(load, Fusarium == "No")
kw <- kruskal.test(LogCopy ~ Group, data = droplevels(ab_only))

contrast <- function(g1, g2) {
  x <- load$LogCopy[load$Group == g1]; y <- load$LogCopy[load$Group == g2]
  data.frame(contrast = paste(g1, "vs", g2), n1 = length(x), n2 = length(y),
             delta_log10 = round(median(x) - median(y), 3),
             fold = round(10^(median(x) - median(y)), 2),
             p = signif(wilcox.test(x, y, exact = FALSE)$p.value, 3))
}
stats <- bind_rows(
  lapply(setdiff(levels(droplevels(ab_only$Group)), "Control"), contrast, g2 = "Control"),
  contrast("Fk", "Control"),
  contrast("AllAB + Fk", "AllAB"),
  contrast("AllAB + Fk", "Fk"))
fk <- filter(load, Group == "Fk")
dose <- cor.test(log10(as.numeric(fk$FusariumDose)), fk$LogCopy, method = "spearman", exact = FALSE)
print(kw); print(stats)
cat(sprintf("Fk load vs log10 dose: rho = %.2f, p = %.2f\n", dose$estimate, dose$p.value))
write.csv(stats, file.path(out_dir, "qpcr_stats.csv"), row.names = FALSE)

# ---- load by treatment group ----
ctrl_med <- median(load$LogCopy[load$Group == "Control"])
grp_cols <- c(abx_cols[grp_levels[1:7]], "Fk" = "#E34A33", "AllAB + Fk" = unname(abx_cols["AllAB"]))
ytop <- max(load$LogCopy)
brackets <- data.frame(
  x1 = match(c("Control", "AllAB"), grp_levels),
  x2 = match(c("Fk", "AllAB + Fk"), grp_levels),
  y  = ytop + c(0.10, 0.32),
  lab = fmt_p(stats$p[match(c("Fk vs Control", "AllAB + Fk vs AllAB"), stats$contrast)]))

p_load <- ggplot(load, aes(Group, LogCopy)) +
  geom_hline(yintercept = ctrl_med, linetype = 2, colour = "grey60", linewidth = 0.4) +
  annotate("rect", xmin = 7.5, xmax = 9.5, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_boxplot(aes(colour = Group), outlier.shape = NA, width = 0.55, linewidth = 0.4,
               fill = NA, show.legend = FALSE) +
  geom_point(aes(fill = Group), shape = 21, colour = "black", stroke = 0.2, size = 1.6,
             position = position_jitter(width = 0.12, seed = 1), show.legend = FALSE) +
  geom_segment(data = brackets, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
  geom_segment(data = brackets, aes(x = x1, xend = x1, y = y, yend = y - 0.03), inherit.aes = FALSE, linewidth = 0.3) +
  geom_segment(data = brackets, aes(x = x2, xend = x2, y = y, yend = y - 0.03), inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = brackets, aes(x = (x1 + x2) / 2, y = y + 0.015, label = lab), vjust = 0, inherit.aes = FALSE, size = 2.5) +
  annotate("text", x = 8.5, y = min(load$LogCopy) - 0.1, label = "Fk-exposed",
           size = 2.5, fontface = "italic") +
  scale_colour_manual(values = grp_cols) + scale_fill_manual(values = grp_cols) +
  scale_x_discrete(labels = c(grp_levels[1:7], "Fk alone", "AllAB + Fk")) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  labs(x = NULL, y = expression(log[10]~"16S copies per reaction"),
       subtitle = paste0("Antibiotics vs Control (no Fk): Kruskal-Wallis ", fmt_p(kw$p.value))) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

# ---- load vs Fusarium dose ----
dose_df <- filter(load, Fusarium == "Yes") %>%
  mutate(Dose = as.numeric(FusariumDose),
         Pretreatment = factor(ifelse(Antibiotic == "AllAB", "AllAB", "None"), levels = c("None", "AllAB")))

p_dose <- ggplot(dose_df, aes(Dose, LogCopy, fill = Pretreatment, shape = Pretreatment)) +
  geom_hline(yintercept = ctrl_med, linetype = 2, colour = "grey60", linewidth = 0.4) +
  annotate("text", x = 1e2, y = ctrl_med, label = "Control median", vjust = -0.5, hjust = 0.1,
           size = 2.4, colour = "grey40") +
  geom_point(size = 1.8, stroke = 0.25, colour = "black",
             position = position_jitter(width = 0.06, height = 0, seed = 1)) +
  scale_x_log10(breaks = 10^(2:5), labels = scales::label_math(10^.x, format = log10)) +
  scale_shape_manual(values = c(None = 21, AllAB = 24)) +
  scale_fill_manual(values = c(None = "#E34A33", AllAB = unname(abx_cols["AllAB"]))) +
  labs(x = "Fusarium dose (spores/mL)", y = expression(log[10]~"16S copies per reaction"),
       fill = "Antibiotic\npretreatment", shape = "Antibiotic\npretreatment",
       caption = sprintf("Fk alone: Spearman \u03c1 = %.2f, %s", dose$estimate, fmt_p(dose$p.value))) +
  theme_pub() +
  theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey25"))

# ---- save ----
save_pub(p_rep,  file.path(out_dir, "qpcr_replicates"),     width_mm = 75,  height_mm = 75)
save_pub(p_load, file.path(out_dir, "qpcr_load_by_group"),   width_mm = 100, height_mm = 75)
save_pub(p_dose, file.path(out_dir, "qpcr_load_vs_dose"),    width_mm = 80,  height_mm = 65)

combined <- p_load / ((p_rep + theme(legend.position = "bottom", legend.direction = "vertical")) | p_dose) +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 10))
save_pub(combined, file.path(out_dir, "qpcr_combined"), width_mm = 140, height_mm = 160)
