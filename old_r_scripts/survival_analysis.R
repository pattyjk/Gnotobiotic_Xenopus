########################################
##Xenopus laevis survivorship analysis
########################################
#load packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(survival)

#read in data
#set better dates
dates <- as.Date(c("2024-10-30", "2024-11-06", "2024-11-12",
                   "2024-11-19", "2024-11-30"))
raw <- data.frame(
  date      = dates,
  Gent      = c(15, 10,  6,  5,  5),
  Chlor     = c(12,  7,  7,  7,  7),
  Trim      = c(15, 10,  8,  8,  8),
  PSN       = c(17, 12, 10, 10, 10),
  ALL_AB_AF = c(11, 10,  9,  8,  7),
  Amp       = c(11, 10, 10,  8,  8),
  Strep     = c(12,  9,  9,  6,  6),
  Sulfa     = c(10,  9,  8,  7,  7),
  Ctrl      = c(12,  7,  7,  5,  5),
  Fk_1e2    = c(10,  7,  6,  6,  6),
  Fk_1e3    = c(12, 10, 10,  7,  7),
  Fk_1e4    = c(10,  9,  9,  8,  8),
  Fk_1e5    = c(11,  9,  9,  9,  6),
  Fk_1e6    = c(14, 10, 11,  8,  8),
  AB_Fk_1e3 = c(12,  9,  9,  8,  8),
  AB_Fk_1e6 = c(10,  9,  9,  8,  7)
)

# Calculate percent survival relative to day 0
pct <- raw
group_cols <- names(raw)[-1]
for (col in group_cols) {
  pct[[col]] <- raw[[col]] / raw[[col]][1] * 100
}

#reshape to long format
long_raw <- pivot_longer(raw, cols = all_of(group_cols),
                         names_to = "group", values_to = "n")
long_pct <- pivot_longer(pct, cols = all_of(group_cols),
                         names_to = "group", values_to = "survival_pct")
long <- left_join(long_raw, long_pct, by = c("date", "group"))

# Assign group categories for faceting / coloring
long <- long %>%
  mutate(category = case_when(
    group %in% c("Gent","Chlor","Trim","PSN","ALL_AB_AF","Amp","Strep","Sulfa") ~ "Antibiotic",
    group %in% c("Fk_1e2","Fk_1e3","Fk_1e4","Fk_1e5","Fk_1e6")               ~ "Antifungal (Fk)",
    group %in% c("AB_Fk_1e3","AB_Fk_1e6")                                     ~ "Combination",
    group == "Ctrl"                                                             ~ "Control"
  ))

# ── Pretty labels for plotting ────────────────────────────────────────────────
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
  Fk_1e2    = "Fk 10²",
  Fk_1e3    = "Fk 10³",
  Fk_1e4    = "Fk 10⁴",
  Fk_1e5    = "Fk 10⁵",
  Fk_1e6    = "Fk 10⁶",
  AB_Fk_1e3 = "AB + Fk 10³",
  AB_Fk_1e6 = "AB + Fk 10⁶"
)
long$label <- label_map[long$group]

#Percent survival — all groups together
p1 <- ggplot(long, aes(x = date, y = survival_pct,
                       color = label, linetype = category, group = label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_x_date(date_labels = "%b %d", date_breaks = "1 week") +
  scale_y_continuous(limits = c(0, 110), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Xenopus laevis survivorship by treatment group",
    subtitle = "Percent survival relative to day 0 (Oct 30, 2024)",
    x        = NULL,
    y        = "Survival (%)",
    color    = "Group",
    linetype = "Category"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right",
        legend.key.width = unit(1.5, "cm"))

print(p1)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/xenopus_survival_all.pdf", p1, width = 11, height = 6, dpi = 150)

#plot, faceted by category
p2 <- ggplot(long, aes(x = date, y = survival_pct, color = label, group = label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ category, ncol = 2) +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 weeks") +
  scale_y_continuous(limits = c(0, 110), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Xenopus laevis survivorship — faceted by treatment category",
    subtitle = "Percent survival relative to day 0",
    x        = NULL,
    y        = "Survival (%)",
    color    = "Group"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold"))

print(p2)
ggsave("~/Documents/GitHub/Gnotobiotic_Xenopus/xenopus_survival_faceted.pdf", p2, width = 11, height = 7, dpi = 150)

####################################
##Kaplan-Meier style analysis
####################################
  expand_to_individual <- function(group_name, counts, obs_dates) {
    n_start <- counts[1]
    rows <- list()
    survivors <- seq_len(n_start)
    current_alive <- n_start
    
    for (i in seq_along(counts)) {
      if (i == 1) next
      died <- current_alive - counts[i]
      time_days <- as.numeric(obs_dates[i] - obs_dates[1])
      if (died > 0) {
        rows[[length(rows) + 1]] <- data.frame(
          group = group_name,
          time  = time_days,
          status = 1,  # event = death
          count = died
        )
      }
      current_alive <- counts[i]
    }
    # Remaining survivors are right-censored at last observation
    last_time <- as.numeric(obs_dates[length(obs_dates)] - obs_dates[1])
    rows[[length(rows) + 1]] <- data.frame(
      group = group_name,
      time  = last_time,
      status = 0,  # censored
      count = current_alive
    )
    bind_rows(rows)
  }
  
  indiv_list <- lapply(group_cols, function(g) {
    expand_to_individual(g, raw[[g]], dates)
  })
  indiv <- bind_rows(indiv_list)
  
  # Expand weighted rows to one row per animal
  indiv_exp <- indiv[rep(seq_len(nrow(indiv)), indiv$count), ]
  indiv_exp$count <- NULL
  
  # Add pretty labels and categories
  indiv_exp$label    <- label_map[indiv_exp$group]
  indiv_exp$category <- case_when(
    indiv_exp$group %in% c("Gent","Chlor","Trim","PSN","ALL_AB_AF","Amp","Strep","Sulfa") ~ "Antibiotic",
    indiv_exp$group %in% c("Fk_1e2","Fk_1e3","Fk_1e4","Fk_1e5","Fk_1e6")               ~ "Antifungal (Fk)",
    indiv_exp$group %in% c("AB_Fk_1e3","AB_Fk_1e6")                                     ~ "Combination",
    indiv_exp$group == "Ctrl"                                                             ~ "Control"
  )
  
  # Fit KM curves
  km_fit <- survfit(Surv(time, status) ~ group, data = indiv_exp)
  cat("\n── Kaplan-Meier median survival times (days) ──\n")
  print(summary(km_fit)$table[, c("records","events","median")])
  
  # Log-rank test: all groups vs control
  cat("\n── Log-rank test: all groups ──\n")
  lr_all <- survdiff(Surv(time, status) ~ group, data = indiv_exp)
  print(lr_all)
  
  # Pairwise log-rank: each treatment vs control
  if (requireNamespace("pairSurv", quietly = TRUE)) {
    # pairSurv::pairwise_survdiff if available
    library(pairSurv)
    cat("\n── Pairwise log-rank vs control ──\n")
    pw <- pairwise_survdiff(Surv(time, status) ~ group,
                            data = indiv_exp, p.adjust.method = "BH")
    print(pw)
  } else {
    # Manual pairwise vs Ctrl
    cat("\n── Pairwise log-rank vs Ctrl (Benjamini-Hochberg adjusted) ──\n")
    non_ctrl <- unique(indiv_exp$group[indiv_exp$group != "Ctrl"])
    p_vals <- sapply(non_ctrl, function(g) {
      sub_df <- indiv_exp[indiv_exp$group %in% c("Ctrl", g), ]
      tryCatch(
        survdiff(Surv(time, status) ~ group, data = sub_df)$pvalue,
        error = function(e) NA
      )
    })
    p_adj <- p.adjust(p_vals, method = "BH")
    result <- data.frame(
      group   = non_ctrl,
      label   = label_map[non_ctrl],
      p_raw   = round(p_vals, 4),
      p_adj   = round(p_adj, 4),
      sig     = ifelse(p_adj < 0.05, "*", "")
    ) %>% arrange(p_adj)
    print(result)
  }
  
#no impact of treatment on survival