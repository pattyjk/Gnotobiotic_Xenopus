############################################
## Shared publication style for Gnotobiotic Xenopus figures
## source("pub_theme.R") at the top of every plotting script
############################################
library(ggplot2)

# ---- treatment names: one canonical spelling everywhere ----
# Canonical antibiotic names (order = plotting order)
abx_levels <- c("Control", "Ampicillin", "AllAB", "PSN", "Sulfamethazine",
                "Gentamicin", "Trimethoprim", "Streptomycin", "Chloramphenicol")

# Map any spelling used in the raw files to the canonical name
recode_abx <- function(x) {
  key <- c(
    "Control" = "Control", "Ctrl" = "Control",
    "Ampicillin" = "Ampicillin", "Amp" = "Ampicillin",
    "AllAB" = "AllAB", "ALL AB/AF" = "AllAB", "All AB/AF" = "AllAB",
    "PSN" = "PSN",
    "Sulfamethazine" = "Sulfamethazine", "Sulfa" = "Sulfamethazine",
    "Sulfamethoxazole" = "Sulfamethazine",   # NOTE: confirm which sulfonamide was used
    "Gentamycin" = "Gentamicin", "Gentamicin" = "Gentamicin", "Gent" = "Gentamicin",
    "Trimethoprim" = "Trimethoprim", "Trim" = "Trimethoprim",
    "Streptomycin" = "Streptomycin", "Strep" = "Streptomycin",
    "Chloramphenicol" = "Chloramphenicol", "Chlor" = "Chloramphenicol",
    "FuskNoAB" = "None")
  out <- unname(key[trimws(as.character(x))])
  ifelse(is.na(out), as.character(x), out)
}

# Colour per antibiotic (Paul Tol "muted" + grey control); colour-blind safe
abx_cols <- c(
  "Control"         = "#7F7F7F",
  "Ampicillin"      = "#88CCEE",
  "AllAB"           = "#332288",
  "PSN"             = "#DDCC77",
  "Sulfamethazine"  = "#AA4499",
  "Gentamicin"      = "#117733",
  "Trimethoprim"    = "#882255",
  "Streptomycin"    = "#44AA99",
  "Chloramphenicol" = "#999933",
  "None"            = "#7F7F7F")

# ---- Fusarium (Fk) exposure: the fungus, NOT an antifungal ----
fk_doses <- c(1e2, 1e3, 1e4, 1e5, 1e6)

# Numeric dose (spores/mL) -> ordered factor "10^2" ... "10^6" (all levels kept,
# so 1e5 can never fall out as NA). Use with scale_*_discrete(labels = fk_labeller)
fk_dose_factor <- function(dose) {
  d <- suppressWarnings(as.numeric(as.character(dose)))
  lab <- ifelse(is.na(d) | d == 0, NA, paste0("10^", round(log10(d))))
  factor(lab, levels = paste0("10^", log10(fk_doses)))
}
# Parses "10^5" -> superscript on axes/legends
fk_labeller <- function(x) parse(text = ifelse(is.na(x), "NA", x))

# Sequential dose palette (light -> dark)
fk_dose_cols <- c("10^2" = "#FDD49E", "10^3" = "#FDBB84", "10^4" = "#FC8D59",
                  "10^5" = "#E34A33", "10^6" = "#B30000")

# Treatment category (for survival / CFU / qPCR panels)
category_levels <- c("Control", "Antibiotic", "Fusarium", "Antibiotic + Fusarium")
category_cols <- c("Control" = "#7F7F7F", "Antibiotic" = "#332288",
                   "Fusarium" = "#E34A33", "Antibiotic + Fusarium" = "#882255")

# Parse survival/CFU column headers ("Fk 10^3", "AB + Fk 10^6", "AB+ Fk 10^6", "Gent", ...)
parse_group_header <- function(h) {
  h <- trimws(h)
  fk  <- grepl("Fk", h)
  ab  <- grepl("^AB", h) | grepl("ALL AB", h, ignore.case = TRUE)
  dose <- suppressWarnings(as.numeric(sub(".*10\\^?([0-9]+).*", "1e\\1", ifelse(fk, h, NA))))
  data.frame(
    header     = h,
    Antibiotic = ifelse(fk & ab, "AllAB", ifelse(fk, "None", recode_abx(h))),
    Fusarium   = ifelse(fk, "Yes", "No"),
    Dose       = dose,
    Category   = factor(ifelse(fk & ab, "Antibiotic + Fusarium",
                         ifelse(fk, "Fusarium",
                         ifelse(recode_abx(h) == "Control", "Control", "Antibiotic"))),
                        levels = category_levels),
    # plain-text label; for superscripts use parse(text = sub("10\\^", "10^", ...)) or LabelExpr
    Label      = ifelse(fk & ab, paste0("AllAB + Fk 10^", round(log10(dose))),
                 ifelse(fk, paste0("Fk 10^", round(log10(dose))), recode_abx(h))),
    LabelExpr  = ifelse(fk & ab, paste0("AllAB~+~Fk~10^", round(log10(dose))),
                 ifelse(fk, paste0("Fk~10^", round(log10(dose))), recode_abx(h))),
    stringsAsFactors = FALSE)
}

# ---- theme ----
# Base 8 pt suits 85 mm (single) / 180 mm (double) column journal widths
theme_pub <- function(base_size = 8, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text         = element_text(colour = "black", size = base_size),
      axis.title        = element_text(colour = "black", size = base_size + 1),
      axis.line         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(1.5, "mm"),
      legend.title      = element_text(size = base_size + 0.5),
      legend.text       = element_text(size = base_size),
      legend.key.size   = unit(3.5, "mm"),
      legend.background = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = base_size + 0.5, margin = margin(2, 0, 2, 0)),
      plot.title        = element_text(face = "bold", size = base_size + 2, hjust = 0, margin = margin(0, 0, 3, 0)),
      plot.subtitle     = element_text(size = base_size, hjust = 0, colour = "grey25", margin = margin(0, 0, 4, 0)),
      plot.title.position = "plot",
      panel.grid        = element_blank(),
      plot.margin       = margin(4, 6, 4, 4))
}

# ---- save: vector PDF + 600-dpi PNG at journal size ----
save_pub <- function(p, file, width_mm = 85, height_mm = 70) {
  base <- sub("\\.(pdf|png)$", "", file)
  # cairo_pdf embeds fonts and handles Unicode (e.g. Greek letters); falls back to pdf()
  if (capabilities("cairo")) {
    ggsave(paste0(base, ".pdf"), plot = p, width = width_mm, height = height_mm,
           units = "mm", device = cairo_pdf)
  } else {
    ggsave(paste0(base, ".pdf"), plot = p, width = width_mm, height = height_mm,
           units = "mm", useDingbats = FALSE)
  }
  ggsave(paste0(base, ".png"), plot = p, width = width_mm, height = height_mm,
         units = "mm", dpi = 600, bg = "white")
  invisible(p)
}

# Format p-values for annotation
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.2g", p))
