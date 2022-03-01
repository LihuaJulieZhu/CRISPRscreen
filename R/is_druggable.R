#' @title Annotate genes with drug information
#' @param gRNAs  a data frame with minimum one column named as the first
#' element in match_col_name such as the output from
#' analyze_CRISPR_screen_data
#' @param drug_target a data frame with drug target information. The one
#' included in this package is from drug bank 2018
#' (https://go.drugbank.com/targets). Users can input their #' own drug
#' target data frame. The latest release can be found at
#' https://go.drugbank.com/releases/latest (access approval required)
#' @param match_col_name Only Gene and Symbol supported for the first element of
#' match_col_name. The first element must be a valid column name in gRNAs data
#' frame and the second element must be a valid column name in drug_target
#' data frame
#' @return an annotated data frame with the last column is_druggable added and the
#' rest of the columns are from gRNAs
#' @author Lihua Julie Zhu
#' @import dplyr
#' @export
#'
#' @examples
#'
#' x <- readRDS(system.file("extdata", "example_data.RDS",
#'             package = "CRISPRscreen"))
#' drug_target <- readRDS(system.file("extdata", "drug_targets_2018.RDS",
#'             package = "CRISPRscreen"))
#' x <- is_druggable(x[1:10,], drug_target,
#'          match_col_name = c("Symbol", "Gene Name"))
#'

is_druggable <- function(gRNAs, drug_target,
                         match_col_name = c("Symbol", "Gene Name"))
{
  drug_target %>%
    mutate(`Gene Name` = as.character(toupper(`Gene Name`))) -> drug_target
  if (match_col_name[1] == "Symbol") {
      gRNAs %>% mutate(Symbol = as.character(toupper(Symbol))) %>%
         mutate(is.druggable = Symbol %in% drug_target$`Gene Name`) -> gRNAs
  }
  else if (match_col_name[1] == "Gene") {
      gRNAs %>% mutate(Symbol = as.character(toupper(Gene))) %>%
         mutate(is.druggable = Symbol %in% drug_target$`Gene Name`) -> gRNAs
  }
  else {
      stop("Only Gene and Symbol supported for match_col_name 1!")
  }
  gRNAs
}
