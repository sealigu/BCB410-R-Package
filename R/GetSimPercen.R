#' This function is used to calculate the similarities between
#' two protein sequences with a percentage output.
#'
#' @param  protein1 A protein Uniprot
#' @param  protein2 A protein Uniprot
#'
#' @return Return a percentage number tells the
#' similarities between selected protein sequences
#' @examples
#' GetSimilarPercentage("Q9UHB7", "Q9UKV5")
#'
#' @name GetSimPercen
#' @rdname GetSimPercen
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

GetSimilarPercentage <- function(protein1, protein2){
  pro1_seq <- as.character(UniprotR::GetSequences(protein1)$Sequence);
  pro2_seq <- as.character(UniprotR::GetSequences(protein2)$Sequence);

  if(identical(pro1_seq, "NA") | identical(pro2_seq, "NA")){
    return("One of the protein sequence is N/A.");
  }

  #By using the Levenshtein distance algorithm
  #for measuring the difference between sequences

  pro1len <- nchar(pro1_seq);
  pro2len <- nchar(pro2_seq);

  #check the sequence length is not 0
  if(identical(pro1len, 0)){
    return(pro2len);
  }
  if(identical(pro2len,0)){
    return(pro1len);
  }

  #create a new empty matrix
  tab <- matrix(0, nrow = pro1len+1, ncol = pro2len+1);
  #initial first row of the matrix
  for(i in 1:(pro1len+1)){
    tab[i, 1] <- i;
  }
  #initial first col of the matrix
  for(j in 1:(pro2len+1)){
    tab[1, j] <- j;
  }

  for(i in 2:pro1len+1){
    ch1 <- substr(pro1_seq, i-1, i-1);
    for(j in 2:pro2len+1){
      ch2 <- substr(pro2_seq, j-1, j-1);
      if(identical(ch1, ch2)){
        hit = 0;
      }else{
        hit = 1;
      }

      tab[i,j] <- min(tab[i-1, j]+1, tab[i, j-1]+1, tab[i-1, j-1]+hit);
    }
  }

  tab_val <- tab[pro1len+1, pro2len+1];
  #using percent() function from library(scales) to convert decimal to percentage
  return(percent(1-tab_val/max(pro1len, pro2len)));
}
