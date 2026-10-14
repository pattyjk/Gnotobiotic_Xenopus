############################################
## Assemble the four main figures (Gnotobiotic Xenopus)
##   Figure 1  Bacterial load:       (a) culturable CFU over time, (b) 16S qPCR load
##   Figure 2  Community:            (a) Bray-Curtis PCoA, (b) richness by antibiotic, (c) richness vs Fk dose
##   Figure 3  Host outcomes:        (a) Kaplan-Meier survival, (b) SVL
##   Figure 4  Differential taxa:    (a) antibiotic LFC heatmap, (b) Fusarium forest plot
##   Figure 5  Order composition:    (a) group means by antibiotic, (b) Fusarium groups by dose
## Each analysis script is run in its own environment; its ggplot objects are reused here.
## Requires: everything the individual scripts need (vegan, survival, patchwork,
## phyloseq, ANCOMBC (CVXR < 1.9), tidyverse, scales, ggh4x).
############################################
library(ggplot2)
library(patchwork)

dir <- "~/Documents/GitHub/Gnotobiotic_Xenopus/"   # data + scripts + pub_theme.R
out_dir <- file.path(dir, "Main_figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(dir, "pub_theme.R"))

# Run a script with `dir` pointed at this folder; return its environment
run_script <- function(file) {
  lines <- readLines(file.path(dir, file))
  lines[grepl("^dir\\s*<-", lines)] <- sprintf('dir <- "%s"', dir)
  env <- new.env(parent = globalenv())
  message("Running ", file, " ...")
  eval(parse(text = lines), envir = env)
  env
}

cfu  <- run_script("cfu_data.R")
qpcr <- run_script("qpcr_plots.R")
comm <- run_script("r_analysis.R")
surv <- run_script("survival_analysis.R")
svl  <- run_script("xenopus_size.R")
da   <- run_script("differnetial_abundance.R")
bars <- run_script("stacked_bar.R")

# Panel letters sit in the outer margin so they never collide with titles/strips
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 11),
                   plot.tag.location = "plot", plot.tag.position = c(0, 1),
                   plot.tag.justification = c(0, 1))
W <- function(p) wrap_elements(full = p) +   # keep each panel's internal layout intact
  theme(plot.margin = margin(t = 14, r = 2, b = 2, l = 2))   # room for the panel letter

# ---- Figure 1: bacterial load ----
fig1 <- (W(cfu$p_cfu) / W(qpcr$p_load)) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "a") & tag_theme
save_pub(fig1, file.path(out_dir, "Figure1_bacterial_load"), width_mm = 180, height_mm = 155)

# ---- Figure 2: community composition and diversity ----
#   (a) PCoA  (b) richness by antibiotic  (c) richness vs Fusarium dose
fig2 <- (W(comm$p_all) / (W(comm$p_alpha_abx) | W(comm$p_alpha_fusk_pt))) +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "a") & tag_theme
save_pub(fig2, file.path(out_dir, "Figure2_community"), width_mm = 180, height_mm = 150)

# ---- Figure 3: host outcomes ----
fig3 <- (W(surv$p_km) / W(svl$p_svl)) +
  plot_layout(heights = c(85, 75)) +
  plot_annotation(tag_levels = "a") & tag_theme
save_pub(fig3, file.path(out_dir, "Figure3_host_outcomes"), width_mm = 180, height_mm = 160)

# ---- Figure 4: differential abundance ----
fig4 <- (W(da$p_ab) / W(da$p_fus)) +
  plot_layout(heights = c(da$h_ab, da$h_fus)) +
  plot_annotation(tag_levels = "a") & tag_theme
save_pub(fig4, file.path(out_dir, "Figure4_differential_taxa"),
         width_mm = 180, height_mm = da$h_ab + da$h_fus)

# ---- Figure 5: order-level composition (group means) ----
#   bar widths proportional to the number of groups in each panel; one shared legend
fig5 <- (bars$pA + bars$pB) +
  plot_layout(widths = c(nlevels(bars$mean_abx$GroupKey), nlevels(bars$mean_fusk$GroupKey)),
              guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "right",
        plot.tag = element_text(face = "bold", size = 11))
save_pub(fig5, file.path(out_dir, "Figure5_order_composition"), width_mm = 180, height_mm = 85)

message("Figures written to ", out_dir)
