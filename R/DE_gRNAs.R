#' @title Differential gRNA analysis
#' @param gRNAs Results from merge_gRNA_files
#' @param col_count1 the column used as the numerator for calculating odds ratio
#' @param col_count2 the column used as the denominator for calculating odds ratio
#' @param min_total_count gRNAs with total observed counts summed over all groups
#' less than min_total_count are excluded before statistical
#' analysis
#' @param multiAdjMethod A vector of character strings containing the names of 
#' the multiple testing procedures for which adjusted p-values are to be 
#' computed. This vector should include any of the following: 
#' "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH", 
#' and "TSBH". Please type ?multtest::mt.rawp2adjp for details. Default to "BH"
#' #'
#' @author Lihua Julie Zhu
#' @importFrom multtest mt.rawp2adjp
#' @importFrom stringr str_c
#' @importFrom methods is
#' @importFrom stats fisher.test
#'

DE_gRNAs <- function(gRNAs,
                     col_count1 = 3,
                     col_count2 = 2,
                     min_total_count = 6,
                     multiAdjMethod = "BH")
{
  if (missing(gRNAs)) {
    stop("gRNAs is required!")
  }
  gRNAs <- unique(gRNAs)
  total1 = sum(as.numeric(gRNAs[, col_count1]))
  total2 = sum(as.numeric(gRNAs[, col_count2]))
  group1 = colnames(gRNAs)[col_count1]
  group2 = colnames(gRNAs)[col_count2]
  gRNAs = gRNAs[rowSums(gRNAs[, c(col_count1, col_count2)]) >= min_total_count, ]
  gRNAs = cbind(gRNAs, c(total1), c(total2))
  colnames(gRNAs)[dim(gRNAs)[2]] = str_c(c(group2, "total"),
                                                 collapse  = ".")
  colnames(gRNAs)[dim(gRNAs)[2] - 1] = str_c(c(group1,
                                                     "total"), collapse = ".")
  unique.pair = unique(cbind(as.numeric(gRNAs[, col_count1]),
                             as.numeric(gRNAs[, col_count2])))
  colnames(unique.pair) = c(colnames(gRNAs)[col_count1],
                            colnames(gRNAs)[col_count2])
  pvalue = do.call(rbind, lapply(1:dim(unique.pair)[1], function(i) {
                temp = matrix(c(unique.pair[i, 1], unique.pair[i, 2],
                total1 - unique.pair[i, 1], total2 - unique.pair[i, 2]), nrow = 2,
                dimnames = list(c(group1, group2), c("yes", "no")))
                r1 = fisher.test(temp)
    c(unique.pair[i, ], r1$p.value, r1$est)
  }))
  colnames(pvalue) = c(str_c(c(group1, "count"), collapse = "."),
                       str_c(c(group2, "count"), collapse = "."),
                       "p.value",
                       str_c(c(group1, "vs", group2,"odds.ratio"), collapse = "."))
  colnames(gRNAs)[col_count1] = str_c(c(group1, "count"),
                                          collapse = ".")
  colnames(gRNAs)[col_count2] = str_c(c(group2, "count"),
                                          collapse = ".")
  gRNAs = merge(gRNAs, pvalue)
  if (multiAdjMethod != "none") {
    procs = c(multiAdjMethod)
    ptemp = gRNAs[, dim(gRNAs)[2] - 1]
    res <- mt.rawp2adjp(ptemp, procs)
    adjp = unique(res$adjp)
    colnames(adjp)[1] =  "p.value"
    colnames(adjp)[2] = str_c(c(group1, "vs", group2, multiAdjMethod, 
                                "adjusted.p.value"),
                              collapse  = ".")
    temp = merge(gRNAs, adjp, all.x = TRUE)
    colnames(temp)[1] = str_c(c(group1, "vs", group2, "p.value"), collapse =".")
    temp[is.na(temp[, dim(temp)[2]] ),dim(temp)[2]] = 1
    temp
  }
  else {
    colnames(gRNAs)[dim(gRNAs)[2] - 1] = str_c(c(group1, "vs", group2, "p.value"), 
                                               collapse =".")
    gRNAs
  }
}
