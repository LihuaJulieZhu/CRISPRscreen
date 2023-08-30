#' @title Create Customized Volcano Plot
#'
#' @param x a data frame such as the output from analyze_CRISPR_screen_data
#' with minimum of three columns: Gene Identifier, gene symbol,
#' odds_ratio, and adj_pvalue/FDR
#' @param adj_pvalue_col the column where adj_pvalues or FDRs are stored
#' @param odds_ratio_col the column where odds_ratios are stored
#' @param adj_pvalue_cutoff default to 0.05 which is used for drawing horizontal
#'  line and specifying enriched/depleted gRNAs/genes
#' @param odds_ratio_cutoff default to 2 which used for drawing vertical lines
#' and specifying enriched/depleted gRNAs/genes
#' @param genes_to_label the genes that will be labeled in the figure
#' @param label_color the color used for highlighting the gene label specified by
#' genes_to_label
#' @param genes_to_highlight the list of genes that will be highlighted in
#'  the figure
#' @param highlight_color the color used for highlighting the genes specified
#' by genes_to_highlight
#' @param gene_col the column containing the gene symbol
#' @param inverse_odds_ratio indicates whether to reverse the odds_ratio
#' @param alpha for specifying the percentage of data to plot
#' @param max.overlaps for setting max.overlaps in the function geom_text_repel
#' @param colors for specifying colors for depleted, enriched, and unchanged
#' gRNAs/genes respectively
#' @param legend_color_labels for labeling the color legend corresponding to
#' the colors specified in colors
#' @param legend_size_labels for labeling the size legend corresponding to the
#' size specified in point_size
#' @param plot_unit Indicates whether plot the data at gRNA or gene level
#' gRNA level means plotting one point per gRNA while gene level means plotting
#' one point per gene with FDR and log_odds_ratio summarized across gRNAs using
#' functions summarize_FDR_fun and summarize_log_odds_ratio_fun
#' @param point_size only applicable when plot_unit is set to gene. Default
#' setting of c(0.2, 1, 2, 3, 4) to plot point size 0.2, 1, 2, 3, and 4 for
#'  genes with 0, 1, 2, 3, and 4 gRNAs significantly increased/decreased.
#'  Please note that any gene with gRNAs altered significantly in both
#'  directions are removed
#' @param summarize_FDR_fun only applicable when plot_unit is set to gene
#' default to min indicating the minimum FDR of enriched gRNAs for each gene
#'  is used for plotting the enriched genes and the minimum FDR of depleted
#'  gRNAs for each gene is used for plotting depleted ones.
#' @param summarize_log_odds_ratio_fun only applicable when plot_unit is set
#' to gene, default to max indicating that the maximum odds_ratio of each gene
#' will be used for plotting enriched gene and the minimum odds_ratio will be
#' used for plotting depleted gene
#' @param highlight_genes_4gRNAs only applicable if plot_unit is gene.
#' It is used for indicating whether to use different color for the genes
#' with 4 gRNAs signficantly altered in the same direction.
#' @param label_size size of the genes to highlight in the plot
#' @author Lihua Julie Zhu
#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#' @importFrom stats median
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#'
#' x <- readRDS(system.file("extdata", "example_data.RDS",
#'             package = "CRISPRscreen"))
#' genes_to_label = c("GPAA1")
#'
#' # Plot one point for each gene by setting plot_unit to "gene"
#' volcano_plot(x, adj_pvalue_col=12, odds_ratio_col=11,
#'        odds_ratio_cutoff = 1.5,
#'        genes_to_label = genes_to_label,
#'        genes_to_highlight = genes_to_label,  plot_unit = "gene",
#'        gene_col = 2, highlight_color = "red", label_color ="blue",
#'       inverse_odds_ratio = FALSE)
#'
#' #  Plot one point for each gRNA by setting plot_unit to "gRNA"
#' volcano_plot(x, adj_pvalue_col=12, odds_ratio_col=11,
#'        genes_to_label = genes_to_label,
#'        genes_to_highlight = genes_to_label,  plot_unit = "gRNA",
#'        gene_col = 2, highlight_color = "red", label_color ="blue",
#'       inverse_odds_ratio = FALSE)
#'
volcano_plot <- function(x,
                         adj_pvalue_col,
                         odds_ratio_col,
                         adj_pvalue_cutoff = 0.05,
                         odds_ratio_cutoff = 1.5,
                         genes_to_label,
                         label_color = "purple",
                         genes_to_highlight,
                         highlight_color = "purple",
                         gene_col = 1,
                         inverse_odds_ratio,
                         alpha = 0.5,
                         max.overlaps = 80,
                         colors = c("purple", "blue", "gray", "orange"),
                         legend_color_labels = c("4 gRNAs Altered", "Depleted",
                                                 "Enriched", "No Change" ),
                         legend_size_labels = c("0 gRNA Altered",
                                                "1 gRNA Altered",
                                                paste(2:4, "gRNAs Altered")),
                         plot_unit = c("gene", "gRNA"),
                         point_size = c(1, 2, 4, 6, 9),
                         summarize_FDR_fun = "min",
                         summarize_log_odds_ratio_fun = "max",
                         highlight_genes_4gRNAs = TRUE,
                         label_size = 4)
{
  if (!highlight_genes_4gRNAs)
  {
      colors <- colors[1:3]
      legend_color_labels <- legend_color_labels[1:3]
  }
  plot_unit <- match.arg(plot_unit)
  colnames(x)[adj_pvalue_col] <- "BH.adjusted.p.value"
  colnames(x)[odds_ratio_col] <- "log_odds_ratio"
  colnames(x)[gene_col] <- "Symbol"

  plot_unit <- match.arg(plot_unit)

  if (plot_unit == "gene")
  {
          up.gRNAs <- filter_gRNAs(x, min_odds_ratio = odds_ratio_cutoff,
                        multi_adj_method = "BH", gene_col = gene_col,
                        maxP = adj_pvalue_cutoff, output_file = "temp.xlsx")
          down.gRNAs <- filter_gRNAs(x,
                         multi_adj_method = "BH", gene_col = gene_col,
                         max_odds_ratio = 1/odds_ratio_cutoff,
                         maxP = adj_pvalue_cutoff, output_file = "temp.xlsx")
          unlink("temp.xlsx")

          x.up <- x %>% filter(log_odds_ratio >= odds_ratio_cutoff &
                               BH.adjusted.p.value < adj_pvalue_cutoff)
          x.down <- x %>% filter(log_odds_ratio <= 1/odds_ratio_cutoff &
                                 BH.adjusted.p.value < adj_pvalue_cutoff)
          x.nc <- x %>% filter(BH.adjusted.p.value >= adj_pvalue_cutoff |
                               (log_odds_ratio < odds_ratio_cutoff &
                                 log_odds_ratio > 1/odds_ratio_cutoff))

          genes.sig <- union(x.up$Symbol, x.down$Symbol)

          x.nc <- x.nc %>% filter(!Symbol %in% genes.sig)

          genes.incons <- intersect(x.up$Symbol, x.down$Symbol)
          genes.incons <- setdiff(genes.incons,
                                  c(genes_to_label,
                                    genes_to_highlight))
          x.up <- x.up %>% filter(!Symbol %in% genes.incons)
          x.down <- x.down %>% filter(!Symbol %in% genes.incons)


          if (summarize_FDR_fun == "min")
          {
             x.up %>%
                 group_by(Symbol) %>%
                 summarise(minP = min(BH.adjusted.p.value, na.rm=TRUE)) -> minP
             x.up <- inner_join(x.up, minP)

             x.down %>%
                 group_by(Symbol) %>%
                 summarise(minP = min(BH.adjusted.p.value, na.rm=TRUE)) -> minP
             x.down <- inner_join(x.down, minP)

             x.nc %>%
                 group_by(Symbol) %>%
                 summarise(minP = min(BH.adjusted.p.value, na.rm=TRUE)) -> minP
             x.nc <- inner_join(x.nc, minP)
          }

          if (summarize_log_odds_ratio_fun == "max")
          {
             x.up %>% group_by(Symbol) %>%
                 summarise(maxOR = max(log_odds_ratio, na.rm=TRUE)) -> maxOR
             x.up <- inner_join(x.up, maxOR)

             x.down %>% group_by(Symbol) %>%
                 summarise(maxOR = min(log_odds_ratio, na.rm=TRUE)) -> maxOR
             x.down <- inner_join(x.down, maxOR)

             x.nc %>% group_by(Symbol) %>%
                 summarise(maxOR = median(log_odds_ratio, na.rm=TRUE)) -> maxOR
             x.nc <- inner_join(x.nc, maxOR)
          }

          x.up %>% select(Symbol, maxOR, minP) %>% unique %>%
               mutate(n.sig.gRNAs = point_size[1]) -> up
          x.down %>% select(Symbol, maxOR, minP) %>% unique %>%
               mutate(n.sig.gRNAs = point_size[1]) -> down

          x.nc %>% select(Symbol, maxOR, minP) %>% unique %>%
               mutate(n.sig.gRNAs = point_size[1]) -> nc
          #return(list(x.up, x.down, x.nc, up, down, nc, up.gRNAs, down.gRNAs))
          up[up[,1] %in% up.gRNAs[[1]]$Symbol,]$n.sig.gRNAs <- point_size[2]
          up[up[,1] %in% up.gRNAs[[2]]$Symbol,]$n.sig.gRNAs <- point_size[3]
	  up[up[,1] %in% up.gRNAs[[3]]$Symbol,]$n.sig.gRNAs <- point_size[4]
          up[up[,1] %in% up.gRNAs[[4]]$Symbol,]$n.sig.gRNAs <- point_size[5]

          down[down[,1] %in% down.gRNAs[[1]]$Symbol,]$n.sig.gRNAs <- point_size[2]
          down[down[,1] %in% down.gRNAs[[2]]$Symbol,]$n.sig.gRNAs <- point_size[3]
          down[down[,1] %in% down.gRNAs[[3]]$Symbol,]$n.sig.gRNAs <- point_size[4]
          down[down[,1] %in% down.gRNAs[[4]]$Symbol,]$n.sig.gRNAs <- point_size[5]

          x <- rbind(up, down, nc)
          odds_ratio_col = 2
          adj_pvalue_col = 3
          colnames(x) <-
            c("Symbol", "log_odds_ratio", "BH.adjusted.p.value", "n.sig.gRNAs")
          #return(list(x,x.up, x.down, up.gRNAs, down.gRNAs, genes.incons))
  }
  else
  {
          x <- x %>% mutate(n.sig.gRNAs = 1)
          #legend_size_labels <- legend_size_labels[1]
          #point_size <- point_size[1]
  }

  x <- x %>% filter(!log_odds_ratio %in% c(0, Inf))

  #x[x[,odds_ratio_col] == 0, odds_ratio_col] <-
           #rnorm(length(x[x[,odds_ratio_col] == 0, odds_ratio_col]), 0.1, 0.01)
          #min(x[x[, odds_ratio_col] != 0, odds_ratio_col]) * 20, 0.01)
  #x[x[,odds_ratio_col] == Inf, odds_ratio_col] <-
           #rnorm(length(x[x[,odds_ratio_col] == Inf, odds_ratio_col]),
  #        max(x[x[, odds_ratio_col] != Inf, odds_ratio_col]) /20, 0.1)

  if (inverse_odds_ratio)
  {
        x[,odds_ratio_col] <- log2(1/x[,odds_ratio_col])
   }
  else
  {
    x[,odds_ratio_col] <- log2(x[,odds_ratio_col])
  }

  x$diffexpressed <- "No change"
  x$diffexpressed[x[, odds_ratio_col] >= log2(odds_ratio_cutoff) &
                    x[, adj_pvalue_col] < adj_pvalue_cutoff] <- "Enriched"
  x$diffexpressed[x[, odds_ratio_col] <= -log2(odds_ratio_cutoff) &
                    x[, adj_pvalue_col] < adj_pvalue_cutoff] <- "Depleted"

  if (highlight_genes_4gRNAs)
  {
     x$diffexpressed[x$n.sig.gRNAs == point_size[5]] <- "4 gRNAs Altered"
  }
  colnames(x)[adj_pvalue_col] <- "p_val_adj"
  if (plot_unit == "gRNA")
      p <- ggplot(data= x,
                      aes(x=log_odds_ratio,
                          y=-log10(p_val_adj), col= diffexpressed))
  else
      p <- ggplot(data= x,
                aes(x=log_odds_ratio, y=-log10(p_val_adj), col= diffexpressed,
                    size = n.sig.gRNAs))
  p <- p +
         geom_point(alpha = alpha) +
         theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank()) +
         geom_text_repel(aes(x=log_odds_ratio,
                    y=-log10(p_val_adj), label = Symbol),
                    color = label_color,
                    max.overlaps =  max.overlaps,
                    data = x[x[, odds_ratio_col] > log2(odds_ratio_cutoff) &
                             x[, adj_pvalue_col] < adj_pvalue_cutoff &
                             x$Symbol %in% genes_to_label &
                             !x$Symbol %in% genes_to_highlight,], size = 2) +
         geom_text_repel(aes(x=log_odds_ratio, y=-log10(p_val_adj),
                         label = Symbol), color = highlight_color,
                         max.overlaps =  max.overlaps,
                         data = x[x[, odds_ratio_col] > log2(odds_ratio_cutoff) &
                                  x[, adj_pvalue_col] < adj_pvalue_cutoff &
                                  x$Symbol %in% genes_to_highlight,],
                         size = label_size) +
         scale_colour_manual(values = colors, labels = legend_color_labels) +
         geom_vline(xintercept=c(-log2(odds_ratio_cutoff),
                                 log2(odds_ratio_cutoff)), col=  "grey") +
         geom_hline(yintercept=-log10(adj_pvalue_cutoff), col=  "grey") +
         theme(legend.title = element_blank())

      if (plot_unit == "gRNA")
        p
      else if(length(legend_size_labels) > 1)
        p + scale_size(name   = "# of significant gRNAs",
               breaks = point_size,
               labels = legend_size_labels)
     else
       p + scale_size(name   = "# of significant gRNAs",
                      breaks = point_size) + 
                      theme(legend.position ="none")
}
