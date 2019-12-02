#' This function is used to plot a graph that shows all the
#' residue modification positons in the selected protein sequences.
#' Input two interested proteins, the plot will show the modified positon
#' for each selected protein and comparing to see whether these two proteins
#' have the same modified position
#'
#' @param  protein1 A protein Uniprot
#' @param  protein2 A protein Uniprot
#'
#' @return Returns a plot to show all the protein residue positions
#' @examples
#' PlotResModification("Q9UHB7", "Q9UKV5")
#'
#' @exportClass
#' @name PlotResModif
#' @import UniprotR
#' @import stringi
#' @import scales
#' @import pheatmap
#'


citation("UniprotR");
citation("stringi");
citation("scales");
citation("pheatmap");

if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if(!requireNamespace("UniprotR", quietly = TRUE)) {
  install.packages("UniprotR")
}

if(!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}

if(!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}

if(!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(UniprotR);
library(stringr);
library(scales);
library(pheatmap);

PlotResModification <- function(protein1 = "Q9UHB7", protein2 = "Q9UKV5") {
  pro1_r <- as.character(UniprotR::GetPTM_Processing(protein1)$Modified.residue);
  pro1_res <- stringr::str_extract_all(pro1_r, stringr::regex("\\d+\\s[A-Za-z]."));
  pro1_res <- unlist(pro1_res);

  pro2_r <- as.character(UniprotR::GetPTM_Processing(protein2)$Modified.residue);
  pro2_res <- stringr::str_extract_all(pro2_r, stringr::regex("\\d+\\s[A-Za-z]."));
  pro2_res <- unlist(pro2_res);

  pro1_res_mat <- matrix(1, nrow = 1,
                         ncol = length(pro1_res),
                         dimnames = list(as.character(protein1), c(pro1_res)));
  pro2_res_mat <- matrix(1, nrow = 1,
                         ncol = length(pro2_res),
                         dimnames = list(as.character(protein2), c(pro2_res)));

  # https://stackoverflow.com/questions/5738773/r-how-to-merge-two-matrix-according-to-their-column-and-row-names
  # authors: Chase and Joris Meys
  # According to above code to merge two matrix together
  # and convert list to matrix

  tab <- merge(pro1_res_mat, pro2_res_mat, by = "row.names", all = TRUE);
  tab[is.na(tab)] <- 0;

  df <- as.matrix(tab[-1]);
  rownames(df) <- tab[,1];

  col<- colorRampPalette(c("white", "purple"))(256);
  f <- pheatmap::pheatmap(df, scale = "none", col = col,
                cluster_cols = FALSE, cluster_rows = FALSE,
                legend = FALSE);
  return(df);
}







