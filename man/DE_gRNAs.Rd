% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DE_gRNAs.R
\name{DE_gRNAs}
\alias{DE_gRNAs}
\title{Differential gRNA analysis}
\usage{
DE_gRNAs(
  gRNAs,
  col_count1 = 3,
  col_count2 = 2,
  min_total_count = 6,
  multi_adj_method = "BH"
)
}
\arguments{
\item{gRNAs}{Results from merge_gRNA_files}

\item{col_count1}{the column used as the numerator for calculating odds ratio}

\item{col_count2}{the column used as the denominator for calculating odds ratio}

\item{min_total_count}{gRNAs with total observed counts summed over all groups
less than min_total_count are excluded before statistical
analysis}

\item{multi_adj_method}{A vector of character strings containing the names of
the multiple testing procedures for which adjusted p-values are to be
computed. This vector should include any of the following:
"Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH",
and "TSBH". Please type ?multtest::mt.rawp2adjp for details. Default to "BH"
#'}
}
\description{
Differential gRNA analysis
}
\author{
Lihua Julie Zhu
}
