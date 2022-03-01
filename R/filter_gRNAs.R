#' @title filter gRNAs
#' @param gRNAs Results from DE_gRNAs
#' @param comparison The comparison name. For example, use "B vs A" for comparing
#' sample B and A
#' @param multiAdjMethod A vector of character strings containing the names of 
#' the multiple testing procedures for which adjusted p-values are to be 
#' computed. This vector should include any of the following: 
#' "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH", 
#' and "TSBH". Please type ?multtest::mt.rawp2adjp for details
#' @param maxP maximum p-value cutoff. Please note that if multiAdjMethod is set,
#' this is the maximum adjusted p-value
#' @param min_odds_ratio gRNAs with odds_ratio >= min_odds_ratio will be kept in
#' the filtered results, default to 1
#' @param max_odds_ratio gRNAs with odds_ratio <= max_odds_ratio will be kept in
#' the filtered results
#' @param ngRNAs_per_gene The designed number of gRNAs per gene in the 
#' screening library
#' @param gene_col For calculating the number of significant gRNAs per gene.
#' Default to 2, i.e., the second column contains the gene symbol.
#' @param output_file The output file name
#'  @return A list of four data frames: genes with at least one gRNA, two, three,
#' and four gRNAs being significant respectively. Significance is determined
#' using the filtering criteria set by users. The filtering criteria
#' includes multiAdjMethod, maxP, min_odds_ratio,
#' and max_odds_ratio. Default settings will be used if users do not modify
#' these parameters.
#' @author Lihua Julie Zhu
#' 
#' @importFrom writexl write_xlsx
#'
filter_gRNAs <- function(gRNAs,
                         comparison,
                         multiAdjMethod = "BH",
                         maxP = 0.05,
                         min_odds_ratio = 1,
                         max_odds_ratio,
                         ngRNAs_per_gene = 4,
                         gene_col = 2,
                         output_file = "enrichedInC.xlsx"
                         )
{
  if (missing(gRNAs)) {
    stop("gRNAs is a required parameter!")
  }
  if (missing(max_odds_ratio) && missing(min_odds_ratio))
  {
    stop("Please specify max_odds_ratio or min_odds_ratio to obtain enriched or depleted gRNAs!\n")
  }
  if (!missing(max_odds_ratio) && !missing(min_odds_ratio))
  {
    stop("Please specify either max_odds_ratio or min_odds_ratio to obtain enriched or depleted gRNAs!\n")
  }
  if (multiAdjMethod != "none") {
    if(!missing(comparison) && comparison != "")
        p_col <- sym(paste(comparison, multiAdjMethod, "adjusted.p.value", sep ="."))
    else
        p_col <- sym(paste(multiAdjMethod, "adjusted.p.value", sep ="."))  # called from volcano_plot
  }
  else if (!missing(comparison) && comparison != "")
        p_col <- sym(paste(comparison, "p.value", sep ="."))
  else
        p_col <- sym("p.value")
  if (!missing(comparison) && comparison != "")
      odds_col <- sym(paste(comparison, "odds.ratio", sep = "."))
  else
      odds_col <- sym("log_odds_ratio")  # called from volcano_plot
  if (!missing(min_odds_ratio))
    gRNAs %>%
      filter(!!odds_col > min_odds_ratio) %>%
      filter(!!p_col < maxP) -> x1
  if (!missing(max_odds_ratio))
    gRNAs %>%
    filter(!!odds_col < max_odds_ratio) %>%
    filter(!!p_col < maxP) -> x1

  y <- as.data.frame(table(x1[,gene_col]))

  temp <- lapply(1:ngRNAs_per_gene, function(i) {
    x1[x1[,gene_col] %in% y[y[,2] >= i ,1],]
  })

  sheet.names <- unlist(lapply(1:ngRNAs_per_gene, function(i) {
     paste("min", i, "gRNAs", sep ="")
  }))
  names(temp) <- sheet.names
  write_xlsx(temp, path = output_file)
  temp
}
