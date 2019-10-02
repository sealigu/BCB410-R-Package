# This .R file is used to analyze the indivudual ubiquitin protein
#
# The first function is used to plot a graph that shows all the
# interacted proteins with the two input proteins.
#
# The second function is used to calculate the difference between
# two protein sequences with a percentage output.
#
#
#
#
#

citation("UniprotR");
citation("stringi");
citation("scales");
library(UniprotR);
library(stringr);
library(scales);
library(pheatmap);

Q9UHB7_seq <- as.character(GetSequences("Q9UHB7")$Sequence);
Q9UHB7_seq

Q9UHB7inter <- as.character(GetProteinInteractions("Q9UHB7")$Interacts.with);
Q9UHB7inter
Q9UHB7sp <- strsplit(Q9UHB7inter, ";"); #create list
Q9UHB7sp <- unlist(Q9UHB7sp);
Q9UHB7sp

Q9UKV5_seq <- as.character(GetSequences("Q9UKV5")$Sequence);
Q9UKV5_seq

Q9UKV5inter <- as.character(GetProteinInteractions("Q9UKV5")$Interacts.with);
Q9UKV5inter
Q9UKV5sp <- strsplit(Q9UKV5inter, ";");
Q9UKV5sp <- unlist(Q9UKV5sp);
Q9UKV5sp;

Q9UHB7_PTM <- GetPTM_Processsing("Q9UHB7");
as.character(Q9UHB7_PTM$Modified.residue);
str_extract_all(Q9UHB7_PTM$Modified.residue, regex("\\d+\\s[A-Za-z]."));


#################################################################
# Example PlotResModification("Q9UHB", "Q9UKV5")
PlotResModification <- function(protein1, protein2) {
  pro1_r <- as.character(GetPTM_Processsing(protein1)$Modified.residue);
  pro1_res <- str_extract_all(pro1_r, regex("\\d+\\s[A-Za-z]."));
  pro1_res <- unlist(pro1_res);

  pro2_r <- as.character(GetPTM_Processsing(protein2)$Modified.residue);
  pro2_res <- str_extract_all(pro2_r, regex("\\d+\\s[A-Za-z]."));
  pro2_res <- unlist(pro2_res);

  pro1_res_mat <- matrix(1, nrow = 1,
                         ncol = length(pro1_res),
                         dimnames = list(as.character(protein1), c(pro1_res)));
  pro2_res_mat <- matrix(1, nrow = 1,
                         ncol = length(pro2_res),
                         dimnames = list(as.character(protein2), c(pro2_res)));

  #https://stackoverflow.com/questions/5738773/r-how-to-merge-two-matrix-according-to-their-column-and-row-names
  #According to above code to merge two matrix together
  tab <- merge(pro1_res_mat, pro2_res_mat, by = "row.names", all = TRUE);
  tab[is.na(tab)] <- 0;

  df <- as.matrix(tab[-1]);
  rownames(df) <- tab[,1];

  col<- colorRampPalette(c("white", "purple"))(256);
  f <- pheatmap(df, scale = "none", col = col);
  return(f);
}
l <- PlotResModification("Q9UHB7", "Q9UKV5");
l

#################################################################
# Example PlotProteinInteractions("Q9UHB7", "Q9UKV5")
PlotProteinInteractions <- function(protein1, protein2) {
  pro1_obj <- as.character(GetProteinInteractions(protein1)$Interacts.with);
  pro2_obj <- as.character(GetProteinInteractions(protein2)$Interacts.with);

  pro1_list <- strsplit(pro1_obj, ";"); #list
  pro2_list <- strsplit(pro2_obj, ";");
  pro1_sp <- unlist(pro1_list); #character
  pro2_sp <- unlist(pro2_list);

  #change "Itself" into its uniprot id
  for(i in 1:length(pro1_sp)){
    if(identical(pro1_sp[i], "Itself")){
      pro1_sp[i] <- as.character(protein1);
    }
  }

  for(j in 1:length(pro2_sp)){
    if(identical(pro2_sp[j], "Itself")){
      pro2_sp[j] <- as.character(protein2);
    }
  }

  pro1_mat <- matrix(1, nrow =  1,
                     ncol = length(pro1_sp),
                     dimnames = list(as.character(protein1), c(pro1_sp)));
  pro2_mat <- matrix(1, nrow = 1,
                     ncol = length(pro2_sp),
                     dimnames = list(as.character(protein2), c(pro2_sp)));

  #https://stackoverflow.com/questions/5738773/r-how-to-merge-two-matrix-according-to-their-column-and-row-names
  #According to above code to merge two matrix together

  df <- merge(pro1_mat, pro2_mat, by = "row.names", all = TRUE);
  df[is.na(df)] <- 0; #replacing NA's with 0
  tab <- as.matrix(df[-1]);
  rownames(tab) <- df[,1];

  col<- colorRampPalette(c("white", "blue"))(256);
  final_tab <- pheatmap(tab, scale = "none", col = col);
  return(final_tab);
}

q <- PlotProteinInteractions("Q9UHB7", "Q9UKV5");
q

##############################################################
# Example GetSimilarPercentage("Q9UHB7", "Q9UKV5");
GetSimilarPercentage <- function(protein1, protein2){
  pro1_seq <- as.character(GetSequences(protein1)$Sequence);
  pro2_seq <- as.character(GetSequences(protein2)$Sequence);

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
p <- GetSimilarPercentage("Q86Y01", "Q9NYG5");
p
