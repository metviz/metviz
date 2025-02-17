# Load necessary libraries
suppressMessages(library(dplyr))
suppressMessages(library(shiny))
suppressMessages(library(Cairo)) 
suppressMessages(library(shinyjs))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(jsonlite))
suppressMessages(library(RCurl))
suppressMessages(library(XML))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(ggplot2))
suppressMessages(library(zoo))
suppressMessages(library(ggrepel))
suppressMessages(library(data.table))
suppressMessages(library(stringi))
suppressMessages(library(reshape2))
suppressMessages(library(bedr))
suppressMessages(library(plotly))
suppressMessages(library(httr))
suppressMessages(library(httr2))
theme_set(theme_cowplot(font_size=12)) 
suppressMessages(library(tidyr))

# --- Data Extraction & Feature Processing ---

extract_g2p <- function(g, u) {
  url <- paste0("https://g2p.broadinstitute.org/api/gene/", g, "/protein/", u, "/protein-features")
  print(url)
  
  # Read table from URL (assumes tab-delimited response)
  g2p_data <- check_file_complete <- function(url) {
    # Download the file temporarily to check line count
    temp_file <- tempfile()
    download.file(url, destfile = temp_file, method = "wget", extra = "-r -p --random-wait")
    #download.file(url, temp_file, quiet = TRUE)
    message(paste("Temporary file saved at:", temp_file))
    print(temp_file)
    # Count the number of lines in the file
    num_lines <- length(readLines(temp_file))
    print(num_lines)
    # Read the file using read.table
    g2p_data <- read.table(temp_file, sep = "	", header = TRUE, fileEncoding = "utf-8", stringsAsFactors = FALSE)
    
    # Check if the number of rows matches the number of lines minus the header
    if (nrow(g2p_data) == (num_lines - 1)) {
      message(paste("File is completely read.", nrow(g2p_data), "| Total lines:", num_lines -1))
    } else {
      warning(paste("File may be incomplete! Rows read:", nrow(g2p_data), "| Total lines:", num_lines))
    }
    
    # Clean up temporary file
    g2p_data$numlines=num_lines
    
    return(g2p_data)
  }
  
  g2p_data <- check_file_complete(url)
  g2p_data$protein_length=tail(g2p_data$residueId, n=1)
  
  # Extract features (using backticks for columns with special characters)
  features <- list(
    Alphafold_pLDDT = if ("AlphaFold.confidence..pLDDT." %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`AlphaFold.confidence..pLDDT.` != "", c("residueId", "AA", "AlphaFold.confidence..pLDDT.")] %>% 
      drop_na()
    } else data.frame(),
    
    Disulfide_bond = if ("Disulfide.bond..UniProt." %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`Disulfide.bond..UniProt.` != "", c("residueId", "AA", "Disulfide.bond..UniProt.")] %>% 
      drop_na()
    } else data.frame(),
    
    Compositional_bias = if ("Compositional.bias..UniProt." %in% colnames(g2p_data)) {
      g2p_data %>%
        filter(!is.na(`Compositional.bias..UniProt.`) & `Compositional.bias..UniProt.` != "") %>%
        select(residueId, AA, `Compositional.bias..UniProt.`)
    } else data.frame(),
    
    Domain = if ("Domain..UniProt." %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`Domain..UniProt.` != "", c("residueId", "AA", "Domain..UniProt.")] %>% 
      drop_na()
    } else data.frame(),
    
    Glycosylation = if ("Glycosylation..UniProt." %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`Glycosylation..UniProt.` != "", c("residueId", "AA", "Glycosylation..UniProt.")] %>% 
      drop_na()
    } else data.frame(),
    
    Modified_residue = if ("Modified.residue..UniProt." %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`Modified.residue..UniProt.` != "", c("residueId", "AA", "Modified.residue..UniProt.")] %>% 
      drop_na()
    } else data.frame(),
    
    Motif = if ("Motif..UniProt." %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`Motif..UniProt.` != "", c("residueId", "AA", "Motif..UniProt.")] %>% 
      drop_na()
    } else data.frame(),
    
    Region = if ("Region..UniProt." %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`Region..UniProt.` != "", c("residueId", "AA", "Region..UniProt.")] %>% 
      drop_na()
    } else data.frame(),
    
    Repeat = if ("Repeat..UniProt." %in% colnames(g2p_data)) {
      g2p_data %>%
        filter(!is.na(`Repeat..UniProt.`) & `Repeat..UniProt.` != "") %>%
        select(residueId, AA, `Repeat..UniProt.`)
    } else data.frame(),
    
    ASA = if ("Accessible.surface.area..Å..." %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`Accessible.surface.area..Å...` != "", c("residueId", "AA", "Accessible.surface.area..Å...")] %>% 
      drop_na()
    } else data.frame(),
    
    Hydropathy = if ("Hydropathy" %in% colnames(g2p_data)) {
      g2p_data[g2p_data$Hydropathy != "", c("residueId", "AA", "Hydropathy")] %>% 
      drop_na()
    } else data.frame(),
    
    Acetylation = if ("Acetylation" %in% colnames(g2p_data)) {
      g2p_data %>%
        filter(!is.na(`Acetylation`) & `Acetylation` != "") %>%
        select(residueId, AA, `Acetylation`)
    } else data.frame(),
    
    Disease_associated_PTM = if ("Disease.associated.PTMs" %in% colnames(g2p_data)) {
      g2p_data %>% 
        filter(!is.na(`Disease.associated.PTMs`) & `Disease.associated.PTMs` != "") %>%
        select (residueId, AA, "Disease.associated.PTMs")
    } else data.frame(),
    
    Methylation = if ("Methylation" %in% colnames(g2p_data)) {
      g2p_data %>% 
        filter(!is.na(`Methylation`) & `Methylation` != "") %>%
        select (residueId, AA, "Methylation")
    } else data.frame(),
    
    O_GalNAc = if ("O.GalNAc" %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`O.GalNAc` != "", c("residueId", "AA", "O.GalNAc")] %>% 
      drop_na()
    } else data.frame(),
    
    O_GlcNAc = if ("O.GlcNAc" %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`O.GlcNAc` != "", c("residueId", "AA", "O.GlcNAc")] %>% 
      drop_na()
    } else data.frame(),
    
    Phosphorylation = if ("Phosphorylation" %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`Phosphorylation` != "", c("residueId", "AA", "Phosphorylation")] %>% 
      drop_na()
    } else data.frame(),
    
    SNP_associated_PTM = if ("SNP.associated.PTMs" %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`SNP.associated.PTMs` != "", c("residueId", "AA", "SNP.associated.PTMs")] %>% 
      drop_na()
    } else data.frame(),
    
    Regulatory_sites = if ("Regulatory.sites" %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`Regulatory.sites` != "", c("residueId", "AA", "Regulatory.sites")] %>% 
      drop_na()
    } else data.frame(),
    
    Substrate_genes = if ("Substrate.genes" %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`Substrate.genes` != "", c("residueId", "AA", "Substrate.genes")] %>% 
      drop_na()
    } else data.frame(),
    
    SUMOylation = if ("SUMOylation" %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`SUMOylation` != "", c("residueId", "AA", "SUMOylation")] %>% 
      drop_na()
    } else data.frame(),
    
    Ubiquitination = if ("Ubiquitination" %in% colnames(g2p_data)) {
      g2p_data[g2p_data$`Ubiquitination` != "", c("residueId", "AA", "Ubiquitination")] %>% 
      drop_na()
    } else data.frame()
  )
  
  # Parse disulfide bond pairs if available
  if (nrow(features$Disulfide_bond) > 0) {
    disulfide_pairs <- features$Disulfide_bond %>%
      group_by(`Disulfide.bond..UniProt.`) %>%
      summarise(residueId1 = min(residueId),
                residueId2 = max(residueId)) %>%
      ungroup()
  } else {
    disulfide_pairs <- data.frame()
  }
  
  return(list(data = g2p_data, features = features, disulfide_pairs = disulfide_pairs))
}

# Function to get protein length
get_protein_length <- function(g2p_data) {
  max(g2p_data$residueId)
}
tail(g2p_data$residueId, n=1)
# --- Confidence Assignment ---
confidence_colors <- c(
  "Very low (0-50)" = "#FF6B6B",
  "Low (50-70)" = "#FFD93D",
  "Confident (70-90)" = "#4D96FF",
  "Very high (90-100)" = "#6BCB77"
)

classify_confidence <- function(pLDDT) {
  case_when(
    pLDDT >= 90 ~ "Very high (90-100)",
    pLDDT >= 70 ~ "Confident (70-90)",
    pLDDT >= 50 ~ "Low (50-70)",
    TRUE ~ "Very low (0-50)"
  )
}

# --- Input Data ---
gene_name <- "DDX3X" ## "PCSK9"   
uniprot_id <- "O00571" ## "Q8NBP7" 
gene_description <- "Proprotein convertase subtilisin/kexin type 9"

data_and_features <- extract_g2p(gene_name, uniprot_id)
#data_and_features <- extract_g2p("PCSK9", "Q8NBP7")

# --- Fetch Data ---
g2p_data <- data_and_features$data
features <- data_and_features$features
disulfide_pairs <- data_and_features$disulfide_pairs

# Add Confidence_Level factor to g2p_data
g2p_data$Confidence_Level <- factor(
  classify_confidence(g2p_data$AlphaFold.confidence..pLDDT.),
  levels = names(confidence_colors)
)

# Protein Length
protein_length <- get_protein_length(g2p_data)
print(protein_length)





# --- Plotting Functions ---

## 📌 Function to plot AlphaFold pLDDT
plot_pLDDT <- function(features, protein_length) {
  if (nrow(features$Alphafold_pLDDT) == 0) return(NULL)
  
  p <- ggplot(g2p_data, aes(x = residueId, y = AlphaFold.confidence..pLDDT., 
                            fill = Confidence_Level,
                            text = paste("Residue:", residueId,
                                         "<br>AA:", AA,
                                         "<br>pLDDT:", round(AlphaFold.confidence..pLDDT., 1),
                                         "<br>Confidence:", Confidence_Level))) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = confidence_colors) +
    labs(title = "Protein Features View",
         x = "Residue Position",
         y = "AlphaFold pLDDT") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
    #scale_x_continuous(expand = expansion(mult = 0.01))
  
  return(ggplotly(p, tooltip = "text") %>%
           layout(hoverlabel = list(bgcolor = "white"),
                  legend = list(orientation = "h", y = -0.2)))
}

## 📌 Function to plot Accessible Surface Area (ASA)
plot_ASA <- function(features, protein_length) {
  if (nrow(features$ASA) == 0) return(NULL)
  
  p <- ggplot(features$ASA, aes(x = residueId, y = Accessible.surface.area..Å..., 
                                text = paste("Residue:", residueId, 
                                             "<br>AA:", AA,
                                             "<br>ASA:", round(Accessible.surface.area..Å..., 1)))) +
    geom_bar(stat = "identity", fill = "#2ca02c", width = 1) +
    labs(title = "",
         x = "",
         y = "ASA") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot Disulfide Bonds
plot_disulfide_bonds <- function(disulfide_pairs, protein_length) {
  if (nrow(disulfide_pairs) == 0) return(NULL)
  
  disulfide_pairs$span <- abs(disulfide_pairs$residueId2 - disulfide_pairs$residueId1)
  disulfide_pairs <- disulfide_pairs[order(disulfide_pairs$span), ]
  disulfide_pairs$y_end <- seq(2, 2 + (nrow(disulfide_pairs) - 1) * 0.1, by = 0.1)
  
  p <- ggplot() +
    geom_segment(data = disulfide_pairs, aes(x = residueId1, y = 1.5, xend = residueId1, yend = y_end), color = "blue", linewidth = 1) +
    geom_segment(data = disulfide_pairs, aes(x = residueId2, y = 1.5, xend = residueId2, yend = y_end), color = "blue", linewidth = 1) +
    geom_segment(data = disulfide_pairs, aes(x = residueId1, y = y_end, xend = residueId2, yend = y_end), color = "blue", linewidth = 1) +
    labs(title = "",
         x = "",
         y = "Disulfide Bonds") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p))
}

## 📌 Function to plot Compositional Bias
plot_compositional_bias <- function(features, protein_length) {
  if (nrow(features$Compositional_bias) == 0) return(NULL)
  
  p <- ggplot(features$Compositional_bias, aes(x = residueId, y = 1,
                                               text = paste("Residue:", residueId,
                                                            "<br>AA:", AA,
                                                            "<br>Feature:", `Compositional.bias..UniProt.`))) +
    geom_point(color = "red", shape=1) +
    labs(title = "",
         x = "",
         y = "Compositional Bias") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot Domains
plot_domains <- function(features, protein_length) {
  if (nrow(features$Domain) == 0) return(NULL)
  
  p <- ggplot(features$Domain, aes(x = residueId, y = 1,
                                   text = paste("Residue:", residueId,
                                                "<br>AA:", AA,
                                                "<br>Domain:", `Domain..UniProt.`))) +
    geom_point(color = "blue", shape=15) +
    labs(title = "",
         x = "",
         y = "Domains") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}


## 📌 Function to plot Regions
plot_regions <- function(features, protein_length) {
  if (nrow(features$Region) == 0) return(NULL)
  
  p <- ggplot(features$Region, aes(x = residueId, y = 1,
                                   text = paste("Residue:", residueId,
                                                "<br>AA:", AA,
                                                "<br>Region:", `Region..UniProt.`))) +
    geom_point(color = "gray", shape=22) +
    labs(title = "",
         x = "",
         y = "Regions") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot Repeats
plot_repeats <- function(features, protein_length) {
  if (nrow(features$Repeat) == 0) return(NULL)
  
  p <- ggplot(features$Repeat, aes(x = residueId, y = 1,
                                   text = paste("Residue:", residueId,
                                                "<br>AA:", AA,
                                                "<br>Repeat:", `Repeat..UniProt.`))) +
    geom_point(color = "pink", shape=15) +
    labs(title = "",
         x = "",
         y = "Repeats") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot Hydropathy
plot_hydropathy <- function(features, protein_length) {
  if (nrow(features$Hydropathy) == 0 ) return(NULL)
  
  # Define colors for positive and negative hydropathy values
  features$Hydropathy$color <- ifelse(as.numeric(features$Hydropathy$Hydropathy) > 0, "red", "blue")
  
  p <- ggplot(features$Hydropathy, aes(x = residueId, y = as.numeric(Hydropathy),
                                       fill = color,
                                       text = paste("Residue:", residueId,
                                                    "<br>AA:", AA,
                                                    "<br>Hydropathy:", Hydropathy))) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_identity(guide = "none") +  # Use the predefined 'color' column for fill
    labs(title = "",
         x = "",
         y = "Hydropathy") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot Acetylation
plot_acetylation <- function(features, protein_length) {
  if (nrow(features$Acetylation) == 0) return(NULL)
  
  p <- ggplot(features$Acetylation, aes(x = residueId, y = 1,
                                        text = paste("Residue:", residueId,
                                                     "<br>AA:", AA,
                                                     "<br>Acetylation:", `Acetylation`))) +
    geom_point(color = "orange", shape=21,fill="pink") +
    labs(title = "",
         x = "",
         y = "Acetylation") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot Disease-associated PTMs
plot_disease_associated_PTM <- function(features, protein_length) {
  if (nrow(features$Disease_associated_PTM) == 0) return(NULL)
  
  p <- ggplot(features$Disease_associated_PTM, aes(x = residueId, y = 1,
                                                   text = paste("Residue:", residueId,
                                                                "<br>AA:", AA,
                                                                "<br>Disease-associated PTM:", `Disease.associated.PTMs`))) +
    geom_point(color = "brown", shape=25, fill="blue") +
    labs(title = "",
         x = "",
         y = "Disease-associated PTMs") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot Methylation
plot_methylation <- function(features, protein_length) {
  if (nrow(features$Methylation) == 0) return(NULL)
  
  p <- ggplot(features$Methylation, aes(x = residueId, y = 1,
                                        text = paste("Residue:", residueId,
                                                     "<br>AA:", AA,
                                                     "<br>Methylation:", `Methylation`))) +
    geom_point(color = "green", shape=13) +
    labs(title = "",
         x = "",
         y = "Methylation") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot O-GalNAc
plot_O_GalNAc <- function(features, protein_length) {
  if (nrow(features$O_GalNAc) == 0) return(NULL)
  
  p <- ggplot(features$O_GalNAc, aes(x = residueId, y = 1,
                                     text = paste("Residue:", residueId,
                                                  "<br>AA:", AA,
                                                  "<br>O-GalNAc:", `O.GalNAc`))) +
    geom_point(color = "darkcyan") +
    labs(title = "",
         x = "",
         y = "O-GalNAc", shape=14) +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot O-GlcNAc
plot_O_GlcNAc <- function(features, protein_length) {
  if (nrow(features$O_GlcNAc) == 0) return(NULL)
  
  p <- ggplot(features$O_GlcNAc, aes(x = residueId, y = 1,
                                     text = paste("Residue:", residueId,
                                                  "<br>AA:", AA,
                                                  "<br>O-GlcNAc:", `O.GlcNAc`))) +
    geom_point(color = "cyan", shape=14) +
    labs(title = "",
         x = "",
         y = "O-GlcNAc") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot Phosphorylation
# --- Phosphorylation Plot ---
plot_phosphorylation <- function(features, protein_length) {
  if (nrow(features$Phosphorylation) == 0) return(NULL)
  
  p <- ggplot(features$Phosphorylation, aes(x = residueId, y = 1,
                                            text = paste("Residue:", residueId,
                                                         "<br>AA:", AA,
                                                         "<br>Phosphorylation:", Phosphorylation))) +
    geom_point(color = "darkred", shape=19) +
    labs(title = "",
         x = "",
         y = "Phosphorylation") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}


## 📌 Function to plot SNP_associated_PTM
plot_SNP_associated_PTM <- function(features, protein_length) {
  if (nrow(features$SNP_associated_PTM) == 0) return(NULL)
  
  p <- ggplot(features$SNP_associated_PTM, aes(x = residueId, y = 1,
                                            text = paste("Residue:", residueId,
                                                         "<br>AA:", AA,
                                                         "<br>SNP_associated_PTM:", `SNP.associated.PTMs`))) +
    geom_point(color = "orange", shape=11) +
    labs(title = "",
         x = "",
         y = "SNP_associated_PTM") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot Regulatory_sites
plot_Regulatory_sites <- function(features, protein_length) {
  if (nrow(features$Regulatory_sites) == 0) return(NULL)
  
  p <- ggplot(features$Regulatory_sites, aes(x = residueId, y = 1,
                                               text = paste("Residue:", residueId,
                                                            "<br>AA:", AA,
                                                            "<br>Regulatory_sites:", `Regulatory.sites`))) +
    geom_point(color = "red", shape=10) +
    labs(title = "",
         x = "",
         y = "Regulatory_sites") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}


## 📌 Function to plot Substrate_genes
plot_substrate_genes <- function(features, protein_length) {
  if (nrow(features$Substrate_genes) == 0) return(NULL)
  
  p <- ggplot(features$Substrate_genes, aes(x = residueId, y = 1,
                                           text = paste("Residue:", residueId,
                                                        "<br>AA:", AA,
                                                        "<br>Substrate_genes:", Substrate.genes))) +
    geom_point(color = "cadetblue", shape=14) +
    labs(title = "",
         x = "",
         y = "Substrate_genes") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}


## 📌 Function to plot SUMOylation
plot_sumoylation <- function(features, protein_length) {
  if (nrow(features$SUMOylation) == 0) return(NULL)
  
  p <- ggplot(features$SUMOylation, aes(x = residueId, y = 1,
                                           text = paste("Residue:", residueId,
                                                        "<br>AA:", AA,
                                                        "<br>SUMOylation:", SUMOylation))) +
    geom_point(color = "darkgoldenrod", shape=15) +
    labs(title = "",
         x = "",
         y = "SUMOylation") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}
## 📌 Function to plot Ubiquitination
plot_ubiquitination <- function(features, protein_length) {
  if (nrow(features$Ubiquitination) == 0) return(NULL)
  
  p <- ggplot(features$Ubiquitination, aes(x = residueId, y = 1,
                                            text = paste("Residue:", residueId,
                                                         "<br>AA:", AA,
                                                         "<br>Ubiquitination:", Ubiquitination))) +
    geom_point(color = "darkblue", shape=18) +
    labs(title = "",
         x = "Residue Position",
         y = "Ubiquitination") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}



## 📌 Function to plot Glycosylation, Modified Residues, Motifs Combined Plot
plot_glycosylation_modified_motifs <- function(features, protein_length) {
  # Standardize column names for Glycosylation, Modified_residue, and Motif
  glycosylation_data <- if (nrow(features$Glycosylation) > 0) {
    features$Glycosylation %>%
      rename(feature_value = `Glycosylation..UniProt.`) %>%
      mutate(feature_type = "Glycosylation")
  } else data.frame(residueId = numeric(), AA = character(), feature_value = character(), feature_type = character())
  
  modified_residue_data <- if (nrow(features$Modified_residue) > 0) {
    features$Modified_residue %>%
      rename(feature_value = `Modified.residue..UniProt.`) %>%
      mutate(feature_type = "Modified Residue")
  } else data.frame(residueId = numeric(), AA = character(), feature_value = character(), feature_type = character())
  
  motif_data <- if (nrow(features$Motif) > 0 && !all(is.na(features$Motif))) {
    features$Motif %>%
      rename(feature_value = `Motif..UniProt.`) %>%
      mutate(feature_type = "Motif")
  } else data.frame(residueId = numeric(), AA = character(), feature_value = character(), feature_type = character())
  
  # Combine data
  combined_data <- bind_rows(glycosylation_data, modified_residue_data, motif_data)
  
  # Remove rows with NA values in residueId or feature_value
  combined_data <- combined_data %>%
    filter(!is.na(residueId) & !is.na(feature_value))
  
  # Return NULL if no data is available
  if (nrow(combined_data) == 0) return(NULL)
  
  # Create the plot
  p <- ggplot(combined_data, aes(x = residueId, y = 1,
                                 text = paste("Residue:", residueId,
                                              "<br>AA:", AA,
                                              "<br>Feature Type:", feature_type,
                                              "<br>Feature:", feature_value),
                                 color = feature_type,
                                 shape = feature_type)) +
    geom_point(size = 3) +
    labs(title = "",
         
         x = "",
         y = "Glycosylation") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # Reduce margins
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

# --- Create Plots Dynamically ---
plots <- list()
if (nrow(g2p_data) > 0) {
  plots$p1_interactive <- plot_pLDDT(features, protein_length)
}
if (nrow(features$ASA) > 0) {
  plots$p2_interactive <- plot_ASA(features, protein_length)
}
if (nrow(disulfide_pairs) > 0) {
  plots$p3_interactive <- plot_disulfide_bonds(disulfide_pairs, protein_length)
}
if (nrow(features$Compositional_bias) > 0) {
  plots$p4_interactive <- plot_compositional_bias(features, protein_length)
}
if (nrow(features$Domain) > 0) {
  plots$p5_interactive <- plot_domains(features, protein_length)
}
if (nrow(bind_rows(features$Glycosylation, features$Modified_residue, features$Motif)) > 0) {
  plots$p6_interactive <- plot_glycosylation_modified_motifs(features, protein_length)
}
if (nrow(features$Region) > 0) {
  plots$p7_interactive <- plot_regions(features, protein_length)
}
if (nrow(features$Repeat) > 0) {
  plots$p8_interactive <- plot_repeats(features, protein_length)
}
if (nrow(features$Hydropathy) > 0) {
  plots$p9_interactive <- plot_hydropathy(features, protein_length)
}
if (nrow(features$Acetylation) > 0) {
  plots$p10_interactive <- plot_acetylation(features, protein_length)
}
if (nrow(features$Disease_associated_PTM) > 0) {
  plots$p11_interactive <- plot_disease_associated_PTM(features, protein_length)
}
if (nrow(features$Methylation) > 0) {
  plots$p12_interactive <- plot_methylation(features, protein_length)
}
if (nrow(features$O_GalNAc) > 0) {
  plots$p13_interactive <- plot_O_GalNAc(features, protein_length)
}
if (nrow(features$O_GlcNAc) > 0) {
  plots$p14_interactive <- plot_O_GlcNAc(features, protein_length)
}
if (nrow(features$SNP_associated_PTM) > 0) {
  plots$p15_interactive <- plot_SNP_associated_PTM(features, protein_length)
}
if (nrow(features$Regulatory_sites) > 0) {
  plots$p16_interactive <- plot_Regulatory_sites(features, protein_length)
}
if (nrow(features$Substrate_genes) > 0) {
  plots$p17_interactive <- plot_substrate_genes(features, protein_length)
}
if (nrow(features$SUMOylation) > 0) {
  plots$p18_interactive <- plot_sumoylation(features, protein_length)
}
 if (nrow(features$Phosphorylation) > 0) {
   plots$p19_interactive <- plot_phosphorylation(features, protein_length)
}
 if (nrow(features$Ubiquitination) > 0) {
   plots$p20_interactive <- plot_ubiquitination(features, protein_length)
 }


# Remove NULL plots from the list
plots <- Filter(Negate(is.null), plots)
print(length(plots))
# Check if there are any plots left after filtering
if (length(plots) == 0) {
  stop("No data available to plot.")
}

# Normalize heights
total_height <- sum(c(1, rep(0.6, length(plots) - 1)))
normalized_heights <- c(1, rep(0.6, length(plots) - 1)) / total_height

# Arrange all plots in a vertical layout with custom heights
subplot(plots, nrows = length(plots), shareX = TRUE, heights = normalized_heights) %>%
  layout(showlegend = TRUE, title = list(
      text = "Protein Feature Visualization in R",
      x = 0.5, 
      xanchor = "center",
      y = 1,
      yanchor = "top",
      font = list(size = 20)
    ),
    xaxis = list(title = "Amino Acid Position")
    ,
    annotations = list(
      list(
        x = 0.5,
        y = -0.15,
        xref = "paper",
        yref = "paper",
        text = paste(gene_name, "-", uniprot_id, " (", protein_length, " aa)", sep = ""),
        showarrow = FALSE,
        xanchor = "center",
        yanchor = "top"
      )
    )
    
    )
