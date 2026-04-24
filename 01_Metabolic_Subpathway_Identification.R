# ==============================================================================
# 01_Metabolic_Subpathway_Identification.R
# Description: KEGG metabolic pathway download, network construction, and 
#              subpathway identification using LTE algorithm
# ==============================================================================

library(igraph)
library(KEGGgraph)
library(KEGGREST)
library(clusterProfiler)
library(org.Hs.eg.db)

# --- Step 1: Define 84 KEGG metabolic pathways ---
meta <- c("00010","00020","00030","00040","00051","00052","00053","00500","00520","00620","00630",
          "00640","00650","00562","00190","00910","00920","00061","00062","00071","00100","00120",
          "00140","00561","00564","00565","00600","00590","00591","00592","01040","00230","00240",
          "00250","00260","00270","00280","00290","00310","00220","00330","00340","00350","00360",
          "00380","00400","00410","00430","00440","00450","00470","00480","00510","00513","00512",
          "00515","00514","00532","00534","00533","00531","00563","00601","00603","00604","00511",
          "00730","00740","00750","00760","00770","00780","00785","00790","00670","00830","00860",
          "00130","00900","00232","00524","00980","00982","00983")

# --- Step 2: Download KEGG KGML files ---
kgml_dir <- "./data/raw/KEGG_KGML"
dir.create(kgml_dir, showWarnings = FALSE, recursive = TRUE)

for (pathway_id in meta) {
  kgml_file <- file.path(kgml_dir, paste0("hsa", pathway_id, ".xml"))
  if (!file.exists(kgml_file)) {
    tryCatch({
      keggGet(paste0("hsa", pathway_id), "kgml", file = kgml_file)
    }, error = function(e) {
      message(paste("Failed to download:", pathway_id))
    })
  }
}

# --- Step 3: Parse KGML and extract subpathways using LTE algorithm ---
setwd(kgml_dir)
idx <- paste0("hsa", meta, ".xml")
subpathway <- list()

for (i in seq_along(idx)) {
  pathway_list <- list()
  pathway <- parseKGML2Graph(idx[i], expandGenes = TRUE, genesOnly = TRUE)
  Nodes <- nodes(pathway)
  Edges <- edges(pathway)
  Edges <- Edges[sapply(Edges, length) > 0]
  
  if (length(Edges) == 0) {
    pathway_list <- list(Nodes)
    names(pathway_list) <- paste(idx[i], 1, sep = "-")
    subpathway <- append(subpathway, pathway_list)
    next
  }
  
  res <- lapply(seq_along(Edges), function(t) {
    name <- names(Edges)[t]
    len <- length(Edges[[t]])
    do.call(rbind, lapply(seq_len(len), function(n) c(name, Edges[[t]][n])))
  })
  result <- data.frame(do.call(rbind, res), stringsAsFactors = FALSE)
  
  a <- paste(result$X1, result$X2, sep = "-")
  b <- paste(result$X2, result$X1, sep = "-")
  if (length(intersect(a, b)) > 0) {
    result <- result[-which(a %in% intersect(a, b)), ]
  }
  
  net <- graph_from_data_frame(result, directed = FALSE)
  fc <- fastgreedy.community(net)
  
  if (modularity(fc) < 0.3) {
    pathway_list <- list(Nodes)
    names(pathway_list) <- paste(idx[i], 1, sep = "-")
    subpathway <- append(subpathway, pathway_list)
    next
  }
  
  for (j in seq_along(fc)) {
    pathway_list <- append(pathway_list, fc[j])
  }
  
  for (k in sort(seq_along(pathway_list), decreasing = TRUE)) {
    if (length(pathway_list[[k]]) < 5) pathway_list[k] <- NULL
  }
  
  if (length(pathway_list) == 0) {
    pathway_list <- list(Nodes)
    names(pathway_list) <- paste(idx[i], 1, sep = "-")
    subpathway <- append(subpathway, pathway_list)
    next
  }
  
  for (l in seq_along(pathway_list)) {
    names(pathway_list)[l] <- paste(idx[i], l, sep = "-")
  }
  subpathway <- append(subpathway, pathway_list)
}

# --- Step 4: Convert Entrez IDs to Gene Symbols ---
subpathway_d <- data.frame()
subpathway_gene <- data.frame()

for (i in seq_along(subpathway)) {
  name <- names(subpathway)[i]
  geneID <- unlist(strsplit(subpathway[[i]], ":"))[2 * seq_along(subpathway[[i]])]
  genesymbol <- bitr(geneID, fromType = "ENTREZID", toType = c("SYMBOL"), OrgDb = "org.Hs.eg.db")
  temp_data <- cbind(data.frame(name = name), genesymbol)
  subpathway_d <- rbind(subpathway_d, temp_data)
  temp_data2 <- data.frame(name = name, gene = paste(genesymbol$SYMBOL, collapse = ", "))
  subpathway_gene <- rbind(subpathway_gene, temp_data2)
}

# --- Step 5: Save results ---
dir.create("./data/processed", showWarnings = FALSE, recursive = TRUE)
write.csv(subpathway_d, "./data/processed/subpathway_gene_symbol.csv", row.names = FALSE)
write.csv(subpathway_gene, "./data/processed/subpathway_gene.csv", row.names = FALSE)
saveRDS(subpathway, "./data/processed/subpathway_list.rds")

cat("Subpathway identification complete. Results saved.\n")
