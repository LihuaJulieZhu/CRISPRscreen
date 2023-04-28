#' @title A work flow function for merging multiple groups, statistical analysis,
#'and filtering results based on p-value/adjusted p-values and log odds ratio.
#'
#' @param gRNA_count_files A vector of file names. Each file contains a list of
#' gRNAs with count information.
#' @param sample_names Sample names corresponding to each file in
#' gRNA_count_files with the same order
#' @param sep delimiter in the gRNA_count_files. Default to space
#' @param header TRUE or FALSE indicating whether the files specified in
#' gRNA_count_files contain header or not. Default to FALSE
#' @param comparisons_group1 The denominator for calculating Log odds ratio
#' @param comparisons_group2 The numerator for calculating log odds ratio
#' @param multiAdjMethod A vector of character strings containing the names of 
#' the multiple testing procedures for which adjusted p-values are to be 
#' computed. This vector should include any of the following: 
#' "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH", 
#' and "TSBH". Please type ?multtest::mt.rawp2adjp for details
#' @param min_total_count gRNAs with total count summed over all groups less 
#' than min_total_count are excluded from further statistical analysis.
#' @param maxP maximum p-value cutoff. Please note that if multiAdjMethod is set,
#' this is the maximum adjusted p-value
#' @param gene_col For calculating the number of significant gRNAs per gene.
#' Default to 2, i.e., the second column contains the gene symbol.
#' @param min_odds_ratio gRNAs with odds_ratio >= min_odds_ratio will be kept in
#' the filtered results, default to 1
#' @param max_odds_ratio gRNAs with odds_ratio <= max_odds_ratio will be kept in
#' the filtered results if provided by users
#' @param ngRNAs_per_gene The designed number of gRNAs per gene in the 
#' screening library
#' @param split Used to extract gene symbol, gene_id, and gRNA from c
#' olumn_to_split of gRNA_count_files. Set split to ":" for extracting gene symbol,
#' gene ID, and gRNA from : separated string, e.g.,
#' TMEM192:NM_001100389.1:TGGGTCGTCTTCAATACTCT.
#' @param column_to_split The name of the column that contains the 
#' annotation. Default to "gene".
#' @param new_column_names Default to c("Symbol", "ID", "gRNA") if 
#' the annotation in column_to_split is in the form of Symbol:ID:gRNA such as
#' MEM192:NM_001100389.1:TGGGTCGTCTTCAATACTCT, where split is ":"
#' @param drug_target a data frame with drug target information. The one
#' included in this package is from drug bank 2018
#' (https://go.drugbank.com/targets). Users can input their #' own drug
#' target data frame. The latest release can be found at
#' https://go.drugbank.com/releases/latest (access approval required)
#' @param output_file The output file name
#' @param ... Set additional parameters for the funciton is_druggable 
#' @return A list of two objects with the first one being a data frame
#' containing all gRNAs >=  min_total_count, and the second one being a list
#' containing four data frames: genes with at least one gRNA, two, three,
#' and four gRNAs being significant respectively. Significance is determined
#' using the filtering criteria set by users. The filtering criteria
#' includes multiAdjMethod, min_total_count, maxP, min_odds_ratio,
#' and max_odds_ratio. Default settings will be used if users do not modify
#' these parameters.
#' @author Lihua Julie Zhu
#' @import dplyr
#' @export
#'
#' @examples
#' gRNA_count_files<- c(system.file("extdata", c("A.txt", "B.txt"),
#'             package = "CRISPRscreen"))
#' drug_target <- readRDS(system.file("extdata", "drug_targets_2018.RDS",
#'             package = "CRISPRscreen"))
#' results <- analyze_CRISPR_screen_data(gRNA_count_files,
#'             sample_names = c("A", "B"),
#'             comparisons_group1 = "B",
#'             comparisons_group2 = "A",
#'             multiAdjMethod = "BH",
#'             min_total_count = 1, drug_target = drug_target,
#'             output_file = "allgRNAs_stats_results.xlsx")
#'
analyze_CRISPR_screen_data <- function(gRNA_count_files,
                                       sample_names = c("A", "B", "C"),
                                       sep =" ",
                                       header = FALSE,
                                       comparisons_group1 = c("B", "C"),
                                       comparisons_group2 = c("A", "A"),
                                       multiAdjMethod,
                                       min_total_count = 6,
                                       maxP = 0.05,
                                       gene_col = 2,
                                       min_odds_ratio = 1,
                                       max_odds_ratio = 1000000,
                                       ngRNAs_per_gene = 4,
                                       split = ":",
                                       column_to_split = "gene",
                                       new_column_names = c("Symbol", "ID", "gRNA"),
                                       drug_target,
                                       output_file = "allgRNAs_stats_results.xlsx",
                                       ...
)
{
  if(missing(comparisons_group1) || missing(comparisons_group2) ||
     length(comparisons_group1) != length(comparisons_group2 )
  )
    stop("Parameters comparisons_group1 and comparisons_group2 are required and their lengths need to be equal!\n")
  if (length(setdiff(c(comparisons_group1, comparisons_group2), sample_names )) > 0)
    stop("The names of comparisons_group1 and comparisons_group2 need to be consistent with sample_names")
  if (length(output_file) != length(comparisons_group1))
    stop("The number of output_file needs to be the same as the number of comparisons specified in comparisions_group1!\n")
  
  if(missing(drug_target)) 
    drug_target <- readRDS(system.file("extdata", "drug_targets_2018.RDS",
                               package = "CRISPRscreen"))
  all_gRNAs <- merge_gRNA_files(gRNA_count_files,
                                sample_names = sample_names,
                                sep = sep,
                                header = header)
  
  annotation <- split_IDs(all_gRNAs,
                          split = split,
                          column_to_split = column_to_split,
                          new_column_names = new_column_names)
  
  stats_results <- lapply(1:length(comparisons_group1), function(i) {
    col_count1 <- which(colnames(all_gRNAs) == comparisons_group1[i])
    col_count2 <- which(colnames(all_gRNAs) == comparisons_group2[i])
    # if (length(comparisons_group1) == 2)
    # unwanted_col <- sym(comparisons_group1[-i])
    DE_gRNAs(all_gRNAs,
             col_count1 = col_count1,
             col_count2 = col_count2,
             min_total_count = min_total_count,
             multiAdjMethod = multiAdjMethod) # %>%
    # select(!(!!unwanted_col))
  })
  
  stats_results <- Reduce(function(x,y) merge(x,y, all = TRUE), stats_results)
  
  rm(all_gRNAs)
  
  n.col <- length(new_column_names) + 2
  stats_results <- merge(annotation[, 1:n.col], stats_results, all.x = TRUE)
  stats_results <- is_druggable(stats_results, drug_target, ...)
  write_xlsx(stats_results, path = output_file)
  
  rm(annotation)
  
  filtered_results <- lapply(1:length(comparisons_group1), function(i) {
    this_comparison <- paste(comparisons_group1[i], "vs",
                             comparisons_group2[i], sep = ".")
    if (exists(deparse(substitute(min_odds_ratio))))
      filter_gRNAs(stats_results,
                   comparison = this_comparison,
                   multiAdjMethod = multiAdjMethod,
                   maxP = maxP, gene_col = gene_col,
                   min_odds_ratio = min_odds_ratio,
                   ngRNAs_per_gene = ngRNAs_per_gene,
                   output_file = output_file[i])
    else if (exists(deparse(substitute(max_odds_ratio))))
      filter_gRNAs(stats_results,
                   comparison = this_comparison,
                   multiAdjMethod = multiAdjMethod,
                   maxP = maxP, gene_col = gene_col,
                   max_odds_ratio = max_odds_ratio,
                   ngRNAs_per_gene = ngRNAs_per_gene,
                   output_file = output_file[i])
  })
  
  list(stats_results = stats_results, filtered_results = filtered_results)
}
