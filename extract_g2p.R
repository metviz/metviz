# Function to extract G2P data from broadinstitute Genomics 2 Proteins (G2P) Portal
extract_g2p <- function(g, u) {
  url <- paste0("https://g2p.broadinstitute.org/api/gene/", g, "/protein/", u, "/protein-features")
  print(url)
  
  g2p_data <- read.table(url, sep = '\t', header = TRUE, fileEncoding = "utf-8")
  return(g2p_data)
}

# Fetch data using HGNC Gene ID and corresponding Uniprot ID
g2p_data <- extract_g2p("GATA4", "P43694")
