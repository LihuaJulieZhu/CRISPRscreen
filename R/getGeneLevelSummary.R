#' Generate Gene-Level Statistical Summary for Each Gene
#'
#' @param x the first list item (named as stats_results) in the
#' analysis results from analyze_CRISPR_screening_data
#' @param comparisons_group1 The denominator for calculating Log odds ratio
#' @param comparisons_group2 The numerator for calculating log odds ratio
#' @param combine_odds_ratio_method The method to summarize gRNA-level odds ratio
#' to gene-level odds ratio. Supported methods are "median", "mean", "max", and "min"
#' @param combine_p_value_method The method to summarize gRNA-level p.value
#' to gene-level p-avalue. Supported methods are Supported methods are
#' "min", "median", "mean", and "max"
#' @param odds.ratio.threshold odds.ratio threshold for selecting the gRNAs.
#' It is used to compute n.gRNA.ge.odds.ratio.threshold for each gene
#' which equals the number of gRNAs with odds.ratio >= odds.ratio.threshold,
#' default to 2
#' @min.gRNAs.ge.odds.ratio.threshold For computing
#' diff.expressed.adjp0.05 and
#' diff.expressed.adjp0.01. Default 3, meaning that
#' at least 3 gRNAs with odds ratio >= odds.ratio.threshold and adjusted
#' p.value < 0.05 and < 0.01 respectively.
#' @param multiAdjMethod A vector of character strings containing the names of
#' the multiple testing procedures for which adjusted p-values are to be
#' computed. This vector should include any of the following:
#' "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH",
#' and "TSBH". Please type ?multtest::mt.rawp2adjp for details
#'
#' @return A data frame with gene-level summary including Symbol,
#'  p.value.gene, odds.ratio.gene, n.gRNAs.ge.odds.ratio.threshold,
#'  diff.expressed.adjp0.05,
#'  diff.expressed.adjp0.01, and log2.odds.ratio. It can be used
#'  to generate a plot with n.gRNAs.ge.odds.ratio.threshold as Y-axis
#'  and log.odds.ratio as the x-axis
#'
#' @import dplyr
#' @export
#'
#' @examples
getGeneLevelSummary <- function(x,
                                comparisons_group1,
                                comparisons_group2,
                                combine_odds_ratio_method =
                                  c("median", "mean", "max", "min"),
                                combine_p_value_method =
                                  c("min", "median", "mean", "max"),
                                odds.ratio.threshold = 2,
                                min.gRNAs.ge.odds.ratio.threshold = 3,
                                multiAdjMethod = "BH")
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
                            odds.ratio.threshold)))
  x2 <- x %>%
    select(Symbol,
           p.value.gene,
           odds.ratio.gene,
          n.gRNAs.ge.odds.ratio.threshold) %>%
    unique %>%
    mutate(adjusted.p.value.gene =
           p.adjust(p.value.gene, method = multiAdjMethod),
         diff.expressed.adjp0.05 =
           adjusted.p.value.gene < 0.05 &&
           n.gRNAs.ge.odds.ratio.threshold >=
           min.gRNAs.ge.odds.ratio.threshold,
         diff.expressed.adjp0.01 =
           adjusted.p.value.gene < 0.01 &&
           n.gRNAs.ge.odds.ratio.threshold >=
           min.gRNAs.ge.odds.ratio.threshold,
         log2.odds.ratio = log2(odds.ratio.gene))


    writexl::write_xlsx(x2,
                    path =
                      paste0(paste(comparisons_group1,
                                          comparisons_group1,
                                          sep =".vs."),
                             "_geneLevel_stats_summary.xls"))
    x2
  }
#   ###### new style of graph
# genes_to_label <- c("WNK1", "MYC")
# library(ggrepel)
# library(ggplot2)
#
# ggplot(data = x2, aes(x = log2.odds.ratio,
#                       y = n.gRNAs.ge2.odds.ratio)) +
#   geom_point() +
#   theme_bw() +
#   geom_point(data = x2[x2$Symbol %in% genes_to_label,],
#              aes(x=log2.odds.ratio,
#                  y=n.gRNAs.ge2.odds.ratio),
#              color = "orange",
#              size = 6) +
#   geom_text_repel(aes(x=log2.odds.ratio,
#                       y=n.gRNAs.ge2.odds.ratio, label = Symbol),
#                   color = "orange",
#                   max.overlaps =  80,
#                   data = x2[x2$Symbol %in% genes_to_label,],
#                   size = 6) +
#   theme(text=element_text(size=24, family="Arial"),
#         axis.text=element_text(size=24, family="Arial", face="bold"),
#         axis.title=element_text(size=24,family="Arial",face="bold"),
#         legend.text = element_text(size=24, family = "Arial"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(colour = "black", size = 1.5)) +
#   theme(legend.title = element_blank()) +  guides(colour = guide_legend(override.aes = list(size=12)))
