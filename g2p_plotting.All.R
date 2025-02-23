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
    # Count the number of lines in the file
    num_lines <- length(readLines(temp_file))
    # Read the file using read.table
    g2p_data <- read.table(temp_file, sep = "\t", quote = "", header = TRUE, fileEncoding = "utf-8", stringsAsFactors = FALSE)
    
    # Check if the number of rows matches the number of lines minus the header
    if (nrow(g2p_data) == (num_lines - 1)) {
      message(paste("File is completely read.", nrow(g2p_data), "| Total lines:", num_lines - 1))
    } else {
      warning(paste("File may be incomplete! Rows read:", nrow(g2p_data), "| Total lines:", num_lines))
    }
    
    # Clean up temporary file
    g2p_data$numlines <- num_lines
    
    return(g2p_data)
  }
  
  g2p_data <- check_file_complete(url)
  g2p_data$protein_length <- tail(g2p_data$residueId, n = 1)
  
  # Define a helper function to extract features
  extract_feature <- function(df, col_name) {
    if (col_name %in% colnames(df)) {
      df %>%
        filter(!is.na(!!sym(col_name)) &!!sym(col_name)!= "") %>%
        select(residueId, AA,!!sym(col_name))
    } else {
      data.frame()
    }
  }
  
  # Extract features using the helper function
  features <- list(
    Alphafold_pLDDT = extract_feature(g2p_data, "AlphaFold.confidence..pLDDT."),
    Disulfide_bond = extract_feature(g2p_data, "Disulfide.bond..UniProt."),
    Compositional_bias = extract_feature(g2p_data, "Compositional.bias..UniProt."),
    Domain = extract_feature(g2p_data, "Domain..UniProt."),
    Glycosylation = extract_feature(g2p_data, "Glycosylation..UniProt."),
    Modified_residue = extract_feature(g2p_data, "Modified.residue..UniProt."),
    Motif = extract_feature(g2p_data, "Motif..UniProt."),
    Region = extract_feature(g2p_data, "Region..UniProt."),
    Repeat = extract_feature(g2p_data, "Repeat..UniProt."),
    ASA = extract_feature(g2p_data, "Accessible.surface.area..Å..."),
    Hydropathy = extract_feature(g2p_data, "Hydropathy"),
    Acetylation = extract_feature(g2p_data, "Acetylation"),
    Disease_associated_PTM = extract_feature(g2p_data, "Disease.associated.PTMs"),
    Methylation = extract_feature(g2p_data, "Methylation"),
    O_GalNAc = extract_feature(g2p_data, "O.GalNAc"),
    O_GlcNAc = extract_feature(g2p_data, "O.GlcNAc"),
    Phosphorylation = extract_feature(g2p_data, "Phosphorylation"),
    SNP_associated_PTM = extract_feature(g2p_data, "SNP.associated.PTMs"),
    Regulatory_sites = extract_feature(g2p_data, "Regulatory.sites"),
    Substrate_genes = extract_feature(g2p_data, "Substrate.genes"),
    SUMOylation = extract_feature(g2p_data, "SUMOylation"),
    Ubiquitination = extract_feature(g2p_data, "Ubiquitination"),
    # Additional features from the provided file
    Amino_acid_properties = extract_feature(g2p_data, "Amino.acid.properties"),
    Secondary_structure_PDBe_SIFTS = extract_feature(g2p_data, "Secondary.structure..PDBe.SIFTS."),
    Secondary_structure_DSSP_3_state = extract_feature(g2p_data, "Secondary.structure..DSSP.3.state.."),
    Secondary_structure_DSSP_9_state = extract_feature(g2p_data, "Secondary.structure..DSSP.9.state.."),
    Phi_angle_degrees = extract_feature(g2p_data, "Phi.angle..degrees.."),
    Psi_angle_degrees = extract_feature(g2p_data, "Psi.angle..degrees.."),
    Active_site = extract_feature(g2p_data, "Active.site..UniProt."),
    Binding_site = extract_feature(g2p_data, "Binding.site..UniProt."),
    Chain = extract_feature(g2p_data, "Chain..UniProt."),
    Coiled_coil = extract_feature(g2p_data, "Coiled.coil..UniProt."),
    Cross_link = extract_feature(g2p_data, "Cross.link..UniProt."),
    DNA_binding = extract_feature(g2p_data, "DNA.binding..UniProt."),
    Initiator_methionine = extract_feature(g2p_data, "Initiator.methionine..UniProt."),
    Intramembrane = extract_feature(g2p_data, "Intramembrane..UniProt."),
    Lipidation = extract_feature(g2p_data, "Lipidation..UniProt."),
    Mutagenesis = extract_feature(g2p_data, "Mutagenesis..UniProt."),
    Non_adjacent_residues = extract_feature(g2p_data, "Non.adjacent.residues..UniProt."),
    Non_standard_residue = extract_feature(g2p_data, "Non.standard.residue..UniProt."),
    Non_terminal_residue = extract_feature(g2p_data, "Non.terminal.residue..UniProt."),
    Peptide = extract_feature(g2p_data, "Peptide..UniProt."),
    Propeptide = extract_feature(g2p_data, "Propeptide..UniProt."),
    Sequence_conflict = extract_feature(g2p_data, "Sequence.conflict..UniProt."),
    Sequence_uncertainty = extract_feature(g2p_data, "Sequence.uncertainty..UniProt."),
    Signal = extract_feature(g2p_data, "Signal..UniProt."),
    Site = extract_feature(g2p_data, "Site..UniProt."),
    Topological_domain = extract_feature(g2p_data, "Topological.domain..UniProt."),
    Transit_peptide = extract_feature(g2p_data, "Transit.peptide..UniProt."),
    Transmembrane = extract_feature(g2p_data, "Transmembrane..UniProt."),
    Zinc_finger = extract_feature(g2p_data, "Zinc.finger..UniProt."),
    Molar_mass_g_mol = extract_feature(g2p_data, "Molar.mass..g.mol."),
    Pocket_number_fpocket = extract_feature(g2p_data, "Pocket.number..fpocket.."),
    Druggability_score_fpocket = extract_feature(g2p_data, "Druggability.score..fpocket.."),
    Intra_chain_Hydrogen_bond_PDB = extract_feature(g2p_data, "Intra.chain.Hydrogen.bond..PDB."),
    Intra_chain_Hydrogen_bond_AlphaFold2 = extract_feature(g2p_data, "Intra.chain.Hydrogen.bond..AlphaFold2."),
    Intra_chain_Non_bonded_interaction_PDB = extract_feature(g2p_data, "Intra.chain.Non.bonded.interaction..PDB."),
    Intra_chain_Non_bonded_interaction_AlphaFold2 = extract_feature(g2p_data, "Intra.chain.Non.bonded.interaction..AlphaFold2."),
    Intra_chain_Disulfide_bond_PDB = extract_feature(g2p_data, "Intra.chain.Disulfide.bond..PDB."),
    Intra_chain_Disulfide_bond_AlphaFold2 = extract_feature(g2p_data, "Intra.chain.Disulfide.bond..AlphaFold2."),
    Intra_chain_Salt_bridge_PDB = extract_feature(g2p_data, "Intra.chain.Salt.bridge..PDB."),
    Intra_chain_Salt_bridge_AlphaFold2 = extract_feature(g2p_data, "Intra.chain.Salt.bridge..AlphaFold2."),
    Inter_chain_Hydrogen_bond_PDB = extract_feature(g2p_data, "Inter.chain.Hydrogen.bond..PDB."),
    Inter_chain_Non_bonded_interaction_PDB = extract_feature(g2p_data, "Inter.chain.Non.bonded.interaction..PDB."),
    Inter_chain_Disulfide_bond_DB = extract_feature(g2p_data, "Inter.chain.Disulfide.bond..PDB."),
    Inter_chain_Salt_bridge_PDB = extract_feature(g2p_data, "Inter.chain.Salt.bridge..PDB.")
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



# --- Input Data ---
gene_name <- "CACNA1A" ## "PCSK9" ## "DDX3X" ##    
uniprot_id <- "O00555" ## "Q8NBP7" ## "O00571" ## 
gene_description <- "Proprotein convertase subtilisin/kexin type 9"

data_and_features <- extract_g2p(gene_name, uniprot_id)
#data_and_features <- extract_g2p("PCSK9", "Q8NBP7")

# --- Fetch Data ---
g2p_data <- data_and_features$data
features <- data_and_features$features
disulfide_pairs <- data_and_features$disulfide_pairs


# --- Confidence Assignment ---
confidence_colors <- c(
  "Very low (0-50)" = "#FF6B6B",
  "Low (50-70)" = "#FFD93D",
  "Confident (70-90)" = "#4D96FF",
  "Very high (90-100)" = "#6BCB77"
)


# Add Confidence_Level factor to g2p_data
g2p_data$Confidence_Level <- factor(
  classify_confidence(g2p_data$AlphaFold.confidence..pLDDT.),
  levels = names(confidence_colors)
)

# Protein Length
protein_length <- get_protein_length(g2p_data)
print(protein_length)


## --- Plotting Functions ---(adding one at a time)




## 📌 Function to plot Amino_acid_properties
plot_amino_acid_properties <- function(features, protein_length) {
  if (nrow(features$Amino_acid_properties) == 0) return(NULL)
  
  # Convert Amino.acid.properties to factor to assign colors
  features$Amino_acid_properties$Amino.acid.properties <- as.factor(features$Amino_acid_properties$Amino.acid.properties)
  
  # Check the unique levels present in the dataset
  current_levels <- unique(features$Amino_acid_properties$Amino.acid.properties)
  print(current_levels)
  
  # Amino Acid Residue Properties
  amino_acid_residue_colors <- c(
    "Aliphatic" = "#FFD93D",
    "Aromatic" = "#FF6B6B",
    "Polar/Neutral" = "#6BCB77",
    "Positively-charged" = "#4D96FF",
    "Negatively-charged" = "#C70039",
    "Special, a very reactive sulfhdryl group" = "#696969",
    "Special, No backbone hydrogen" = "#696969",
    "default" = "#808080"
  )
  
  # Handle any missing or unexpected levels by assigning them a default color
  missing_levels <- setdiff(current_levels, names(amino_acid_residue_colors))
  if (length(missing_levels) > 0) {
    amino_acid_residue_colors <- c(amino_acid_residue_colors, setNames(rep("white", length(missing_levels)), missing_levels))
  }
  
  p <- ggplot(features$Amino_acid_properties, aes(x = residueId, y = 1, fill = `Amino.acid.properties`,
                                                  text = paste("Residue:", residueId,
                                                               "<br>AA:", AA,
                                                               "<br>Amino_Acid_Property:", `Amino.acid.properties`))) +
    geom_bar(stat = "identity", width = 1) +
    labs(title = "Amino Acid Properties", y = NULL) +  # Add title back
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length)) +
    scale_fill_manual(values = amino_acid_residue_colors)
  
  return(ggplotly(p, tooltip = "text"))
}

## 📌 Function to plot Secondary_structure_DSSP_9_state
plot_secondary_structure_DSSP_9_state <- function(features, protein_length) {
  if (nrow(features$Secondary_structure_DSSP_9_state) == 0) return(NULL)
  
  # Check the unique levels present in the dataset
  current_levels <- unique(features$Secondary_structure_DSSP_9_state$Secondary.structure..DSSP.9.state..)
  print(current_levels)  # Debugging: View the unique levels
  
  # Convert the column to a factor to ensure correct mapping
  features$Secondary_structure_DSSP_9_state$Secondary.structure..DSSP.9.state.. <- 
    as.factor(features$Secondary_structure_DSSP_9_state$Secondary.structure..DSSP.9.state..)
  
  # Define colors for each DSSP 9-state type
  ss_colors <- c(
    "H (α-helix)"          = "magenta",
    "G (3₁₀-helix)"        = "magenta",
    "I (π-helix)"          = "magenta",
    "E (parallel sheets)"  = "yellow",
    "B (beta bridge)"      = "yellow",
    "T (turn)"             = "paleturquoise",  
    "S (bend)"             = "paleturquoise",
    "C (loop/coil)"        = "white",
    "P (polyproline helix)"= "palegreen"
  )
  
  # Handle any missing or unexpected levels by assigning them a default color
  missing_levels <- setdiff(current_levels, names(ss_colors))
  if (length(missing_levels) > 0) {
    ss_colors <- c(ss_colors, setNames(rep("white", length(missing_levels)), missing_levels))
  }
  
  # Create the plot
  p <- ggplot(features$Secondary_structure_DSSP_9_state, 
              aes(x = residueId, 
                  y = 1, 
                  fill = Secondary.structure..DSSP.9.state..,
                  text = paste("Residue:", residueId,
                               "<br>AA:", AA,
                               "<br>Secondary Structure:", Secondary.structure..DSSP.9.state..))) +
    geom_bar(stat = "identity", width = 1) +
    #labs(title = "Secondary Structure (DSSP 9-state)", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length)) +
    scale_fill_manual(values = ss_colors)  # Apply the custom color mapping
  
  return(ggplotly(p, tooltip = "text"))
}



## 📌 Function to plot Secondary_structure_PDBe_SIFTS
plot_secondary_structure_PDBe_SIFTS <- function(features, protein_length) {
  if (nrow(features$Secondary_structure_PDBe_SIFTS) == 0) return(NULL)
  
  # Convert only the relevant column to factor
  features$Secondary_structure_PDBe_SIFTS$Secondary.structure..PDBe.SIFTS. <- 
    as.factor(features$Secondary_structure_PDBe_SIFTS$Secondary.structure..PDBe.SIFTS.)
  
  # Ensure all possible levels are present, even if not in the data
  levels(features$Secondary_structure_PDBe_SIFTS$Secondary.structure..PDBe.SIFTS.) <- 
    c("Helix", "Beta strand", "Turn", "Other")

    # Define custom colors for secondary structures
  ## "C (loop/coil)" "H (helix)"     "B (strand)" 
  ss_colors <- c(
    "Helix" = "magenta",
    "Beta strand" = "yellow",
    "Turn" = "paleturquoise"  # Pale blue
   )
  # Dynamically handle any other levels present in the data
  # If any unknown levels are present, we will color them as "Other"
  levels_present <- levels(features$Secondary_structure_PDBe_SIFTS$Secondary.structure..PDBe.SIFTS.)
  missing_levels <- setdiff(levels_present, names(ss_colors))
  
  # Add missing levels to the color mapping as "Other"
  if (length(missing_levels) > 0) {
    ss_colors <- c(ss_colors, setNames(rep("white", length(missing_levels)), missing_levels))
  }
  
  p <- ggplot(features$Secondary_structure_PDBe_SIFTS, 
              aes(x = residueId, 
                  y = 1, 
                  fill = Secondary.structure..PDBe.SIFTS.,
                  text = paste("Residue:", residueId,
                               "<br>AA:", AA,
                               "<br>Secondary Structure PDBe SIFTS:", Secondary.structure..PDBe.SIFTS.))) +
    geom_bar(stat = "identity", width = 1) +
   # labs(title = "Secondary Structure (PDBe SIFTS)", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length)) +
    scale_fill_manual(values = ss_colors) 
  
  return(ggplotly(p, tooltip = "text"))
}


## 📌 Function to plot Secondary_structure_DSSP_3_state
plot_secondary_structure_DSSP_3_state <- function(features, protein_length) {
  if (nrow(features$Secondary_structure_DSSP_3_state) == 0) return(NULL)
  
  # Check the unique levels present in the dataset
  current_levels <- unique(features$Secondary_structure_DSSP_3_state$Secondary.structure..DSSP.3.state..)
  print(current_levels)  # Debugging: View the unique levels
  
  # Convert the column to a factor to ensure correct mapping
  features$Secondary_structure_DSSP_3_state$Secondary.structure..DSSP.3.state.. <- 
    as.factor(features$Secondary_structure_DSSP_3_state$Secondary.structure..DSSP.3.state..)
  
  # Define colors for each DSSP 3-state type
  ss_colors1 <- c(
    "H (helix)"   = "magenta",
    "B (strand)"  = "yellow",
    "C (loop/coil)" = "paleturquoise"
  )
  
  # Handle any missing or unexpected levels by assigning them a default color
  missing_levels <- setdiff(current_levels, names(ss_colors1))
  if (length(missing_levels) > 0) {
    ss_colors1 <- c(ss_colors1, setNames(rep("white", length(missing_levels)), missing_levels))
  }
  
  # Create the plot
  p <- ggplot(features$Secondary_structure_DSSP_3_state, 
              aes(x = residueId, 
                  y = 1, 
                  fill = Secondary.structure..DSSP.3.state..,
                  text = paste("Residue:", residueId,
                               "<br>AA:", AA,
                               "<br>Secondary Structure DSSP3:", Secondary.structure..DSSP.3.state..))) +
    geom_bar(stat = "identity", width = 1) +
    #labs(title = "Secondary Structure (DSSP 3-state)", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length)) +
    scale_fill_manual(values = ss_colors1)  # Apply the custom color mapping
  
  return(ggplotly(p, tooltip = "text"))
}

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
    # labs(title = "Protein Features View",
    #      x = "Residue Position",
    #      y = "AlphaFold pLDDT") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"), # Reduce margins
      legend.position="none"
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
    #scale_x_continuous(expand = expansion(mult = 0.01))
  
  return(ggplotly(p, tooltip = "text") %>%
           layout(hoverlabel = list(bgcolor = "white")))
}

## 📌 Function to plot Accessible Surface Area (ASA)
plot_ASA <- function(features, protein_length) {
  if (nrow(features$ASA) == 0) return(NULL)
  
  p <- ggplot(features$ASA, aes(x = residueId, y = Accessible.surface.area..Å..., 
                                text = paste("Residue:", residueId, 
                                             "<br>AA:", AA,
                                             "<br>ASA:", round(Accessible.surface.area..Å..., 1)))) +
    geom_bar(stat = "identity", fill = "green", width = 1, color = "darkgreen") +  # Add darkgreen outline
    # labs(title = "",
    #      x = "",
    #      y = "ASA") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
                                                            "<br>Compositional_bias:", `Compositional.bias..UniProt.`))) +
    geom_point(color = "red", shape=18) +
    # labs(title = "",
    #      x = "",
    #      y = "Compositional Bias") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
    # labs(title = "",
    #      x = "",
    #      y = "Domains") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
    # labs(title = "",
    #      x = "",
    #      y = "Regions") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
    # labs(title = "",
    #      x = "",
    #      y = "Repeats") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
    # labs(title = "",
    #      x = "",
    #      y = "Hydropathy") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
    # labs(title = "",
    #      x = "",
    #      y = "Acetylation") +  # Remove Y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y = element_blank(),  # Remove Y-axis title
      axis.text.y = element_blank(),   # Remove Y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove Y-axis ticks
      panel.background = element_blank(),  # Blank background
      panel.grid.major = element_blank(),  # No grid lines
      panel.grid.minor = element_blank(),  # No grid lines
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position="none"
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  # Reduce margins
      legend.position="none"
    ) +
    scale_x_continuous(expand = expansion(mult = 0.01), limits = c(1, protein_length))
  
  return(ggplotly(p, tooltip = "text"))
}

# --- Create Plots Dynamically ---
plots <- list()
# Execute the function if data exists
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
if (nrow(features$Amino_acid_properties) > 0) {
  plots$p20_interactive <- plot_amino_acid_properties(features, protein_length)
}
if (nrow(features$Secondary_structure_PDBe_SIFTS) > 0) {
  plots$p_secondary_structure_pdbe_sifts <- plot_secondary_structure_PDBe_SIFTS(features, protein_length)
}
if (nrow(features$Secondary_structure_DSSP_3_state) > 0) {
  plots$p_secondary_structure_dssp_3_state <- plot_secondary_structure_DSSP_3_state(features, protein_length)
}

if (nrow(features$Secondary_structure_DSSP_9_state) > 0) {
  plots$p_secondary_structure_dssp_9_state <- plot_secondary_structure_DSSP_9_state(features, protein_length)
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

#Arrange all plots in a vertical layout with custom heights
subplot(plots, nrows = length(plots), shareX = TRUE, heights = normalized_heights) %>%
  layout(title = list(
      text = "Protein Feature Visualization in R",
      x = 0.5,
      xanchor = "center",
      y = 1,
      yanchor = "top",
      font = list(size = 20)
    ),
    xaxis = list(title =  paste("Amino Acid Position\n", gene_name, "-", uniprot_id, " (", protein_length, " aa)", sep = "") ),
    theme(legend.position = "none") 
    )


# features$Secondary_structure_DSSP_3_state
# levels_present <- unique(features$Secondary_structure_PDBe_SIFTS$Secondary.structure..PDBe.SIFTS.)
