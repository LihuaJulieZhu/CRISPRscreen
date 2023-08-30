#' Generate Gene-Level Statistical Summary for Each Gene
#'
#' @param x the first list item (named as stats_results) in the
#' analysis results from analyze_CRISPR_screening_data
#' @param comparisons_group1 The numerator for calculating Log odds ratio
#' @param comparisons_group2 The denominator for calculating log odds ratio
#' @param combine_odds_ratio_method The method to summarize gRNA-level odds ratio
#' to gene-level odds ratio. Supported methods are "median", "mean", "max", and "min"
#' @param combine_p_value_method The method to summarize gRNA-level p.value
#' to gene-level p-avalue. Supported methods are Supported methods are
#' "min", "median", "mean", and "max"
#' @param odds_ratio_threshold odds ratio threshold for selecting the gRNAs.
#' It is used to compute n.gRNA.ge.odds_ratio_threshold for each gene
#' which equals the number of gRNAs with odds.ratio >= odds_ratio_threshold,
#' default to 2
#' @param min_gRNAs_ge_odds_ratio_threshold For computing
#' diff.expressed.adjp0.05 and
#' diff.expressed.adjp0.01. Default 3, meaning that
#' at least 3 gRNAs with odds ratio >= odds_ratio_threshold and adjusted
#' p.value < 0.05 and < 0.01 respectively.
#' @param multi_adj_method A vector of character strings containing the names of
#' the multiple testing procedures for which adjusted p-values are to be
#' computed. This vector should include any of the following:
#' "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH",
#' and "TSBH". Please type ?multtest::mt.rawp2adjp for details
#'
#' @return A data frame with gene-level summary including Symbol,
#'  p.value.gene, odds.ratio.gene, n.gRNAs.ge.odds_ratio_threshold,
#'  diff.expressed.adjp0.05,
#'  diff.expressed.adjp0.01, and log2.odds.ratio. It can be used
#'  to generate a plot with n.gRNAs.ge.odds_ratio_threshold as Y-axis
#'  and log.odds.ratio as the x-axis
#'
#' @import dplyr
#' @importFrom stats p.adjust
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
#'             multi_adj_method = "BH",
#'             min_total_count = 1, drug_target = drug_target,
#'             output_file = "allgRNAs_stats_results.xlsx")
#'
#' geneResults <- get_gene_level_summary(results$stats_results,
#'      comparisons_group1 = "B",
#'      comparisons_group2 = "A",
#'      combine_odds_ratio_method = "median",
#'      combine_p_value_method= "min")
#'  
#' 
get_gene_level_summary <- function(x,
                                comparisons_group1,
                                comparisons_group2,
                                combine_odds_ratio_method =
                                  c("median", "mean", "max", "min"),
                                combine_p_value_method =
                                  c("min", "median", "mean", "max"),
                                odds_ratio_threshold = 2,
                                min_gRNAs_ge_odds_ratio_threshold = 3,
                                multi_adj_method = "BH")
{
  combine_odds_ratio_method <- match.arg(combine_odds_ratio_method)
  combine_p_value_method <- match.arg(combine_p_value_method)
  odds.ratio.column.name <- paste0(paste(comparisons_group1,
                                    comparisons_group2,
                             sep =".vs."), ".odds.ratio")
  p.value.column.name <- paste0(paste(comparisons_group1,
                                         comparisons_group2,
                                         sep =".vs."), ".p.value")
  odds.ratio.column <-
    which(colnames(x) == odds.ratio.column.name)

  comp1.count.col <- which(colnames(x) == paste0(comparisons_group1,
                                           ".count"))
  comp1.total.col <- which(colnames(x) == paste0(comparisons_group1,
                                                 ".total"))
  comp2.count.col <- which(colnames(x) == paste0(comparisons_group2,
                                                 ".count"))
  comp2.total.col <- which(colnames(x) == paste0(comparisons_group2,
                                                 ".total"))

  x[x[,odds.ratio.column] == Inf, odds.ratio.column] <-
    (x[x[,odds.ratio.column]== Inf, comp1.count.col] + 1)/
    (x[x[,odds.ratio.column] == Inf, comp1.total.col]+ 1)/
    ((x[x[,odds.ratio.column] == Inf, comp2.count.col] + 1)/
     (x[x[,odds.ratio.column] == Inf,comp2.total.col] + 1))

  x[x[,odds.ratio.column] == 0, odds.ratio.column] <-
    (x[x[,odds.ratio.column]== 0, comp1.count.col] + 1)/
    (x[x[,odds.ratio.column] == 0, comp1.total.col]+ 1)/
    ((x[x[,odds.ratio.column] == 0, comp2.count.col] + 1)/
       (x[x[,odds.ratio.column] == 0,comp2.total.col] + 1))

  x <- x %>%  group_by(Symbol) %>%
    mutate(p.value.gene  = match.fun(
      combine_p_value_method)(!!as.name(p.value.column.name)),
         odds.ratio.gene = match.fun(
           combine_odds_ratio_method)(!!as.name(odds.ratio.column.name)),
         n.gRNAs.ge.odds.ratio.threshold =
           sum(as.numeric(!!as.name(odds.ratio.column.name) >=
                            odds_ratio_threshold)))
  x2 <- x %>%
    select(Symbol,
           p.value.gene,
           odds.ratio.gene,
          n.gRNAs.ge.odds.ratio.threshold) %>%
    unique %>%
    mutate(adjusted.p.value.gene =
           p.adjust(p.value.gene, method = multi_adj_method),
         diff.expressed.adjp0.05 =
           adjusted.p.value.gene < 0.05 &&
           n.gRNAs.ge.odds.ratio.threshold >=
           min_gRNAs_ge_odds_ratio_threshold,
         diff.expressed.adjp0.01 =
           adjusted.p.value.gene < 0.01 &&
           n.gRNAs.ge.odds.ratio.threshold >=
           min_gRNAs_ge_odds_ratio_threshold,
         log2.odds.ratio = log2(odds.ratio.gene))


    writexl::write_xlsx(x2,
                    path =
                      paste0(paste(comparisons_group1,
                                          comparisons_group2,
                                          sep =".vs."),
                             "_geneLevel_stats_summary.xls"))
    x2
  }
