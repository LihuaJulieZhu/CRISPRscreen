#' @title Merge multiple gRNA count files from different sample groups
#' @param gRNA_count_files A vector of file names. Each file contains a list of
#' gRNAs with count information.
#' @param sample_names Sample names corresponding to each file in
#' gRNA_count_files with the same order
#' @param sep delimiter in the gRNA_count_files. Default to space
#' @param header TRUE or FALSE indicating whether the files specified in
#' gRNA_count_files contain header or not. Default to FALSE
#' 
#' @author Lihua Julie Zhu
#' @importFrom utils read.table
#'
#'
merge_gRNA_files <- function(gRNA_count_files,
                             sample_names , 
                             sep =" ", 
                             header = FALSE) {
  if (missing(gRNA_count_files)) {
    stop("gRNA_count_files is a required parameter!")
  }
  if (missing(sample_names)) {
    stop("sample_names is a required parameter!")
  }
  gRNAs <- read.table(gRNA_count_files[1], sep = sep, header = header,
                      stringsAsFactors = FALSE)
  colnames(gRNAs) <- c("gene", sample_names[1])
  for (i in 2:length(gRNA_count_files)){
    x <- read.table(gRNA_count_files[i], sep = sep, header = header,
                    stringsAsFactors = FALSE)
    colnames(x) <- c("gene", sample_names[i])
    gRNAs <- merge(gRNAs, x, all = TRUE)
  }

  for (i in 2:(length(gRNA_count_files) + 1)) {
      gRNAs[is.na(gRNAs[,i]), i] <- 0
  }
  gRNAs
}
