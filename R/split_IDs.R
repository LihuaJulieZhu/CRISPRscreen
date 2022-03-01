#' @title Split the gene ID column into multiple columns using delimiter
#' specified by split
#' @param all_gRNAs Results from merge_gRNAs
#' @param split Used to extract gene symbol, gene_id, and gRNA from c
#' olumn_to_split of gRNA_count_files. Set split to ":" for extracting gene symbol,
#' gene ID, and gRNA from : separated string, e.g.,
#' TMEM192:NM_001100389.1:TGGGTCGTCTTCAATACTCT.
#' @param column_to_split The name of the column that contains the 
#' annotation. Default to "gene".
#' @param new_column_names Default to c("Symbol", "ID", "gRNA") if 
#' the annotation in column_to_split is in the form of Symbol:ID:gRNA such as
#' MEM192:NM_001100389.1:TGGGTCGTCTTCAATACTCT, where split is ":"
#'
#' @author Lihua Julie Zhu
#'
split_IDs <- function(all_gRNAs, split = ":",
                     column_to_split = "gene",
         new_column_names = c("Symbol", "ID", "gRNA"))
{
  geneIDs <- do.call(rbind, lapply(as.character(all_gRNAs[, column_to_split]),
                              function(id) {unlist(strsplit(id, split))}))
  all_gRNAs <- cbind(geneIDs, all_gRNAs)
  if (ncol(geneIDs) == length(new_column_names))
  {
    colnames(all_gRNAs)[1:ncol(geneIDs)] <- new_column_names
  }
  else
  {
    cat("Length of new_column_names is ", length(new_column_names),
        ". However, the number of extracted IDs is ", ncol(geneIDs), "!")
    stop("Length of new_column_names is not the same as the number of extracted IDs!")
  }
  n.gRNAs <- as.data.frame(table(all_gRNAs[,1]))
  colnames(n.gRNAs) <- c(colnames(all_gRNAs)[1], "n.gRNAs")
  all_gRNAs <- merge(n.gRNAs, all_gRNAs)
  all_gRNAs
}
