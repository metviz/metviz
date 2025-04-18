# Load necessary libraries
library(ggplot2)
library(dplyr)
library(data.table)
library(plotly)

# Load ClinVar variant data
input_data <- fread("~/Path_clinvar.feb19.goldstars.pf", sep = "\t", header=TRUE, stringsAsFactors=FALSE)

# Define domain/region data
pfam_data <- list(
  Domain = data.frame(
    type = c("Domain", "Domain"),
    description = c("Inhibitor I9", "Peptidase S8"),
    start = c(77, 155),
    end = c(149, 461),
    color = c("#2dcf00", "#ff5353")
  ),
  Region = data.frame(
    type = c("Region"),
    description = c("C-terminal domain"),
    start = c(450),
    end = c(692),
    color = c("#5b5bff")
  )
)

# Combine domain and region data
domain_region_data <- rbind(pfam_data$Domain, pfam_data$Region)

# Function to process only pathogenic data
gene_data_pathogenic <- function(gene_symbol) {
  gene_data <- input_data %>%
    filter(GeneSymbol == gene_symbol, ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic")) %>%
    mutate(
      ClinicalSignificance = factor(ClinicalSignificance, levels = c("Likely pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic")),
      y_position = case_when(
        ClinicalSignificance == "Likely pathogenic" ~ 0.5,
        ClinicalSignificance == "Pathogenic/Likely pathogenic" ~ 1,
        ClinicalSignificance == "Pathogenic" ~ 1.5
      )
    ) %>%
    filter(!is.na(pos_aa)) %>%
    mutate(pos_aa = as.numeric(pos_aa))
  
  return(gene_data)
}

# Function to plot only pathogenic and likely pathogenic variants
create_gene_plot_pathogenic <- function(gene_symbol) {
  gene_data <- gene_data_pathogenic(gene_symbol)
  
  if (nrow(gene_data) == 0) {
    stop(paste("No pathogenic data found for gene:", gene_symbol))
  }
  
  max_pos <- max(c(gene_data$pos_aa, domain_region_data$end), na.rm = TRUE)
  
  p <- ggplot() +
    # Protein backbone
    geom_rect(aes(xmin = 1, xmax = max_pos, ymin = -0.1, ymax = 0), 
              fill = "grey80", alpha = 0.8, color = "black", size = 0.3) +
    
    # Domains and regions clearly below backbone
    geom_rect(data = domain_region_data,
              aes(xmin = start, xmax = end, ymin = -0.1, ymax = 0.01, fill = description),
              alpha = 1, color = "black", size = 0.3) +
    geom_text(data = domain_region_data,
              aes(x = (start + end)/2, y = -0.05, label = description),
              vjust = 0.5, size = 3, color = "black") +
    
    # Variants above zero line clearly separated
    geom_segment(data = gene_data,
                 aes(x = pos_aa, xend = pos_aa, y = 0, yend = y_position),
                 color = "gray70") +
    geom_point(data = gene_data,
               aes(x = pos_aa, y = y_position, color = ClinicalSignificance),
               size = 3) +
    
    scale_color_manual(values = c(
      "Likely pathogenic" = "red",
      "Pathogenic/Likely pathogenic" = "#CD5B45",
      "Pathogenic" = "darkred"
    )) +
  ## This is specific to PCSK9 (replace it for others)
    scale_fill_manual(values = c(
      "Inhibitor I9" = "#2dcf00",
      "Peptidase S8" = "#ff5353",
      "C-terminal domain" = "#5b5bff"
    )) +
    scale_x_continuous(breaks = seq(0, max_pos, by = 50), minor_breaks = NULL) +
    scale_y_continuous(breaks = c(0, 0.5, 1, 1.5), 
                       labels = c("Domain/Region", "(Likely pathogenic) 0.5",
                                  "(Pathogenic/Likely pathogenic) 1", "(Pathogenic) 1.5"),
                       limits = c(-0.5, 2)) +
    labs(title = paste("Pathogenicity Along", gene_symbol, "Gene"),
         x = "Amino Acid Position", y = "",
         color = "Clinical Significance", fill = "Domain/Region") +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 10, angle = 0, hjust = 1))
  
  p=ggplotly(p, tooltip = c("pos_aa", "ClinicalSignificance"))
  return(p)
}


# Example usage
plot <- create_gene_plot_pathogenic("PCSK9")
print(plot)
