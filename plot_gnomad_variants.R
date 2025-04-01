suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_gnomad_variants.R <variants.csv> <AF_cutoff> [comma-separated list of aa_pos: 1,17,32,113]")

}

csv_file <- args[1]
af_cutoff <- as.numeric(args[2])
highlight_input <- ifelse(length(args) >= 3, args[3], "")

# Parse highlight positions
highlight <- NULL
if (highlight_input != "") {
  highlight_positions <- as.numeric(unlist(str_split(highlight_input, ",")))
  highlight <- data.frame(aa_pos = highlight_positions)
}

# Load variant data
variants <- read_csv(csv_file, show_col_types = FALSE)

if (!"aa_pos" %in% names(variants)) {
  stop("Input CSV must contain 'aa_pos' column parsed from hgvsp")
}

# Define plotting attributes
variants <- variants %>%
  mutate(
    freq_status = ifelse(exome_af > af_cutoff, "AboveCutoff", "BelowCutoff"),
    plot_freq = ifelse(exome_af > af_cutoff, af_cutoff, exome_af),
    shape = ifelse(freq_status == "AboveCutoff", 17, 16),
    color = case_when(
      freq_status == "AboveCutoff" & exome_ac_hom > 0 ~ "#ff6500",     # triangle + hom → orange
      freq_status == "AboveCutoff" ~ "#1c663c",                        # triangle → dark green
      exome_ac_hom > 0 ~ "#ff6500",                                    # circle + hom → orange
      TRUE ~ "#2CA25F"                                                 # circle default green
    )
  )

prot_length <- max(variants$aa_pos, na.rm = TRUE)
y_max_af <- af_cutoff + af_cutoff * 0.1

# Build the Rainfall AF plot
p <- ggplot(variants, aes(x = aa_pos, y = plot_freq)) +
  geom_point(aes(shape = factor(shape), color = color), size = 2) +
  scale_shape_manual(values = c(`16` = 16, `17` = 17)) +
  scale_color_identity() +
  scale_y_reverse(limits = c(y_max_af, -0.0000005), expand = c(0, 0)) +
  xlim(-prot_length * 0.01, prot_length + prot_length * 0.01) +
  labs(
    x = "Amino Acid Position",
    y = "gnomAD freq",
    title = paste("gnomAD AF Rainfall Plot:", tools::file_path_sans_ext(basename(csv_file)))
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "none"
  )

# Add vertical highlight lines if residue numbers are provided
if (!is.null(highlight) && nrow(highlight) > 0) {
  p <- p + geom_vline(data = highlight, aes(xintercept = aa_pos), linetype = "dotted", color = "blue")
}

# View plot
#print (p)

# Save plot
output_file <- str_replace(csv_file, "\\.csv$", "_rainfall_plot.png")
ggsave(output_file, p, width = 8, height = 4, dpi = 300)
cat("Rainfall plot saved to:", output_file, "\n")
