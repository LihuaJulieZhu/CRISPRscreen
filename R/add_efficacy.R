#' @title Add gRNA efficacy to each gRNA
#'
#' @param gRNAs DNAStringSet object containing a set of gRNAs. Please note the
#' sequences must contain PAM appended after gRNAs, e.g.,
#' ATCGAAATTCGAGCCAATCCCGG where ATCGAAATTCGAGCCAATCC is the gRNA and CGG is
#' the PAM
#' @param rule_set Specify a rule set scoring system for calculating gRNA
#' efficacy. Please note that Root_RuleSet2_2016 requires the following python
#' packages with specified verion and python 2.7.  1. scikit-learn 0.16.1 2.
#' pickle 3. pandas 4. numpy 5. scipy
##' @param BSgenomeName BSgenome object. Please refer to available.genomes in
#' BSgenome package. For example,
#' \itemize{
#' \item{BSgenome.Hsapiens.UCSC.hg19} - for hg19,
#' \item{BSgenome.Mmusculus.UCSC.mm10} - for mm10
#' \item{BSgenome.Celegans.UCSC.ce6} - for ce6
#' \item{BSgenome.Rnorvegicus.UCSC.rn5} - for rn5
#' \item{BSgenome.Drerio.UCSC.danRer7} - for Zv9
#' \item{BSgenome.Dmelanogaster.UCSC.dm3} - for dm3
#' }
#' @param output_dir the directory where the off target analysis and reports
#' will be written to
#' @param overwrite overwrite the existing files in the output directory or
#' @param ... Please see additional parameters in offTargetAnalysis of 
#' CRISPRseek pacakage
#' not, default FALSE
#' @author Lihua Julie Zhu
#' @seealso CRISPRseek
#' @importFrom CRISPRseek offTargetAnalysis
#' @return gRNA efficacy
#' @author Lihua Julie Zhu
#' @export
#' @examples
#'  library(CRISPRseek)
#' 	library("BSgenome.Hsapiens.UCSC.hg19")
#' 	output_dir <- getwd()
#'
add_efficacy <- function(gRNAs,
                         rule_set = "CRISPRscan",
                         BSgenomeName,
                         output_dir,
                         overwrite = FALSE, ...)
{
  if (rule_set == "CRISPRscan")
    featureWeightMatrixFile <- system.file("extdata", "Morenos-Mateo.csv",
                                         package = "CRISPRseek")
  else
      featureWeightMatrixFile = system.file("extdata", "DoenchNBT2014.csv",
                                        package = "CRISPRseek")
  res <-  offTargetAnalysis(gRNAs,
                            gRNAoutputName = names(gRNAs),
                            rule.set = rule_set,
                            findgRNAs = FALSE,
                            exportAllgRNAs = "no",
                            annotatePaired = FALSE,
                            annotateExon = FALSE,
                            useEfficacyFromInputSeq = FALSE,
                            BSgenomeName = BSgenomeName,
                            max.mismatch = 0,
                            featureWeightMatrixFile = featureWeightMatrixFile,
                            outputDir = output_dir,
                            overwrite = overwrite, ...)
  res
}
