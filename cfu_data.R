########################################
##Xenopus laevis CFU data
########################################
library()

data<-read.delim('~/Documents/GitHub/Gnotobiotic_Xenopus/cfU-data.txt', header=T)


library(ggplot2)
library(tidyr)
library(dplyr)

#Data
dates <- as.Date(c("2024-10-30", "2024-11-06", "2024-11-12",
                   "2024-11-19", "2024-11-30"))

raw <- data.frame(
  date      = dates,
  Gent      = c(660, 500, 486,  85, 245),
  Chlor     = c(816, 280, 217, 231, 287),
  Trim      = c(945, 850, 272, 224, 552),
  PSN       = c(867, 864, 200, 490, 392),
  ALL_AB_AF = c(825, 320, 324, 568, 210),
  Amp       = c(715, 640, 240, 688, 424),
  Strep     = c(936, 666, 639, 174, 246),
  Sulfa     = c(340, 360, 160, 217, 588),
  Ctrl      = c(756, 357, 133, 375, 410),
  Fk_1e2    = c(740, 609, 204, 408, 222),
  Fk_1e3    = c(804, 560, 430, 126, 147),
  Fk_1e4    = c(330, 261, 270, 696, 192),
  Fk_1e5    = c(517, 162, 756, 252, 306),
  Fk_1e6    = c(434, 740, 869, 232, 304),
  AB_Fk_1e3 = c(876, 297, 297, 472, 208),
  AB_Fk_1e6 = c(290, 603, 693, 440, 168)
)

label_map <- c(
  Gent      = "Gentamicin",
  Chlor     = "Chloramphenicol",
  Trim      = "Trimethoprim",
  PSN       = "PSN",
  ALL_AB_AF = "All AB/AF",
  Amp       = "Ampicillin",
  Strep     = "Streptomycin",
  Sulfa     = "Sulfamethoxazole",
  Ctrl      = "Control",
  Fk_1e2    = "Fk 10\u00b2",
  Fk_1e3    = "Fk 10\u00b3",
  Fk_1e4    = "Fk 10\u2074",
  Fk_1e5    = "Fk 10\u2075",
  Fk_1e6    = "Fk 10\u2076",
  AB_Fk_1e3 = "AB + Fk 10\u00b3",
  AB_Fk_1e6 = "AB + Fk 10\u2076"
)

group_cols <- names(raw)[-1]

#Reshape to long format
long <- pivot_longer(raw, cols = all_of(group_cols),
                     names_to = "group", values_to = "CFU") %>%
  mutate(
    label = label_map[group],
    category = case_when(
      group %in% c("Gent","Chlor","Trim","PSN","ALL_AB_AF","Amp","Strep","Sulfa") ~ "Antibiotic",
      group %in% c("Fk_1e2","Fk_1e3","Fk_1e4","Fk_1e5","Fk_1e6")               ~ "Antifungal (Fk)",
      group %in% c("AB_Fk_1e3","AB_Fk_1e6")                                     ~ "Combination",
      group == "Ctrl"                                                             ~ "Control"
    ),
    day = as.numeric(date - min(date))
  )

# ── Plot 1: Raw CFU over time — all groups ────────────────────────────────────
p1 <- ggplot(long, aes(x = date, y = CFU, color = label,
                       linetype = category, group = label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_x_date(date_labels = "%b %d", date_breaks = "1 week") +
  labs(
    title    = "CFU over time by treatment group",
    subtitle = "Raw colony forming unit counts, Oct 30 – Nov 30, 2024",
    x        = NULL, y = "CFU/mL", color = "Group", linetype = "Category"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right", legend.key.width = unit(1.5, "cm"))

print(p1)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/CFU_all_groups.pdf", p1, width = 12, height = 6, dpi = 150)

# ── Plot 2: Log10-transformed CFU (helps when counts span orders of magnitude) 
p2 <- p1 +
  scale_y_log10() +
  labs(
    title    = "CFU over time — log\u2081\u2080 scale",
    subtitle = "Log-transformed CFU counts"
  )

print(p2)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/CFU_log_scale.pdf", p2, width = 12, height = 6, dpi = 150)

# Faceted by category
p3 <- ggplot(long, aes(x = date, y = CFU, color = label, group = label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ category, ncol = 2, scales = "free_y") +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 weeks") +
  labs(
    title    = "CFU over time — faceted by treatment category",
    x        = NULL, y = "CFU", color = "Group"
  ) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"), legend.position = "right")

print(p3)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/CFU_faceted.pdf", p3, width = 12, height = 7, dpi = 150)

#CFU normalized to day 0 (fold change)
long_norm <- long %>%
  group_by(group) %>%
  mutate(CFU_norm = CFU / CFU[date == min(date)]) %>%
  ungroup()

p4 <- ggplot(long_norm, aes(x = date, y = CFU_norm, color = label,
                            linetype = category, group = label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_x_date(date_labels = "%b %d", date_breaks = "1 week") +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), "\u00d7")) +
  labs(
    title    = "CFU fold change relative to day 0",
    subtitle = "Values >1 indicate increase; <1 indicate decrease vs. baseline",
    x        = NULL, y = "Fold change (relative to Oct 30)",
    color = "Group", linetype = "Category"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right", legend.key.width = unit(1.5, "cm"))

print(p4)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/CFU_fold_change.pdf", p4, width = 12, height = 6, dpi = 150)


# Statistical comparison: Repeated-measures approach
ctrl_vals <- long$CFU[long$group == "Ctrl"]

non_ctrl <- unique(long$group[long$group != "Ctrl"])
p_vals <- sapply(non_ctrl, function(g) {
  grp_vals <- long$CFU[long$group == g]
  wilcox.test(grp_vals, ctrl_vals, paired = TRUE, exact = FALSE)$p.value
})
p_adj <- p.adjust(p_vals, method = "BH")

stat_result <- data.frame(
  group       = non_ctrl,
  label       = label_map[non_ctrl],
  p_raw       = round(p_vals, 4),
  p_adj_BH    = round(p_adj, 4),
  significant = ifelse(p_adj < 0.05, "*", "")
) %>% arrange(p_adj_BH)

print(stat_result)
#no sig difference from control, so no change in bacterial abundance

