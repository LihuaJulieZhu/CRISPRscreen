#' @title Create the distribution of log odds ratio and summarized FDR for given genes
#' with or without pathway annotation
#'
#' @param x a data frame such as the output from analyze_CRISPR_screen_data
#' with minimum of three columns: Gene Identifier, gene symbol,
#' odds_ratio, and adj_pvalue/FDR
#' @param genes_to_highlight the list of genes that will be highlighted in
#'  the figure. When plot_facet is set to by_gene or none, genes_to_highlight
#'  needs to be specified as a character vector.
#' @param category_anno If plot_facet is set to by_category, then
#' category_anno needs to be specified as a data frame with category specified in
#' anno_category_col and gene symbol specified in anno_gene_col.
#' @param anno_category_col If plot_facet is set to by_category, then
#' category_anno needs to be specified with category specified in
#' anno_category_col and gene symbol specified in anno_gene_col.
#' @param anno_gene_col If plot_facet is set to by_category, then
#' category_anno needs to be specified with category specified in
#' anno_category_col and gene symbol specified in anno_gene_col.
#' @param highlight_color_up the color used for highlighting the
#' genes specified by genes_to_highlight and up
#' @param highlight_color_down the color used for highlighting the
#' genes specified by genes_to_highlight and down
#' @param gene_col the column containing the gene symbol
#' @param alpha for specifying the percentage of data to plot in histogram,
#' default to 0.5 (50 percent)
#' @param odds_ratio_col the column where odds_ratios are stored
#' @param FDR_col the column where adj_pvalues or FDRs are stored
#' @param bins the number of bins used for histogram generation using
#' geom_histogram, default to 100
#' @param inverse_odds_ratio indicates whether to reverse the odds_ratio
#' @param plot_facet either by by_gene, by_category, or none for facet plotting
#' If plot_facet is set to by_category, then category_anno needs to be
#' specified with category specified in anno_category_col and gene symbol
#' specified in anno_gene_col.
#' @param gene_label_size the size of gene label
#' @param dashline_color the color for drawing the dashed line for log odds
#' ratio of 0
#' @param FDR_label_size the size of FDR label
#' @param cutoff specify where to draw to vertical dash line to indicate the
#' cutoff. Default to log odds ratio of 0.
#' @param include_FDR Indicates whether to include FDR on the right side of
#' the plot for each gene. Default to TRUE.
#' @param Indicates whether the genes_to_highlight is ordered.
#' Default to TRUE, i.e., genes_to_highlight will be plotted from top to bottom.
#' @author Lihua Julie Zhu
#' @import patchwork
#' @import ggplot2
#' @import dplyr
#' @importFrom stats density
#' @return a list of ggplot object of size 2
#' @export
#'
#' @examples
#'
#' # add category on the right side of the bottom figure by setting
#' # plot_facet to by_category
#' # load the example category data frame
#'
#' category_anno <- readRDS(system.file("extdata", "example_category_ann.RDS",
#'             package = "CRISPRscreen"))
#'
#' x <- readRDS(system.file("extdata", "example_data.RDS",
#'             package = "CRISPRscreen"))
#'
#' p <- frequency_plot(x = x, category_anno = category_anno[1:6,],
#'                    genes_to_highlight = category_anno[1:6,2],
#'                    plot_facet = "by_category")
#' p[[1]]/p[[2]]
#'
#' # Set plot_facet to by_gene
#' genes_to_highlight <- category_anno[1:6, 2]
#' p <- frequency_plot(x = x, genes_to_highlight = genes_to_highlight,
#'                    plot_facet = "by_gene")
#'  p[[1]]/p[[2]]

frequency_plot <-function(x, genes_to_highlight,
			 category_anno,
       anno_category_col = 1,
       anno_gene_col = 2,
       highlight_color_up = "red",
       highlight_color_down = "blue",
       gene_col = 2,
       alpha = 0.5,
       odds_ratio_col = 11,
       FDR_col = 12,
       bins = 100,
       inverse_odds_ratio = FALSE,
       plot_facet = c("by_gene", "by_category", "none"),
       gene_label_size = 12,
       dashline_color = "purple",
       FDR_label_size = 4,
       cutoff = 0,
       include_FDR = TRUE,
       gene_list_ordered = TRUE
)
{
  if (missing(x) || class(x) != "data.frame")
     stop(paste("x is required as a data.frame object with gene symbol, odds_ratio or fold_change",
       "and FDR specified in gene_col, odds_ratio_col, and FDR_col respectively!"))
  plot_facet <- match.arg(plot_facet)

  if (missing(category_anno) &&  plot_facet == "by_category") {
     stop(paste("When plot_facet is set to by_category, category_anno needs to be specified",
        "with category specified in anno_category_col and gene symbol specified in anno_gene_col!"))
  }
  if (missing(genes_to_highlight) && (plot_facet == "by_gene" || plot_facet == "none")) {
     stop(paste("When plot_facet is set to by_gene or none, genes_to_highlight needs to be",
        "specified as a character vector!"))
  }
  colnames(x)[odds_ratio_col] <- "log_odds_ratio"
  colnames(x)[gene_col] <- "Gene"
  colnames(x)[FDR_col] <- "FDR"

  if (plot_facet == "by_gene") {
     category_anno <- cbind(Category = "", Gene = genes_to_highlight)
  }
   #x[x[,odds_ratio_col] == 0, odds_ratio_col] <-
   #	rnorm(length(x[x[,odds_ratio_col] == 0, odds_ratio_col]), 0.1, 0.01)

   x <- x %>% filter(log_odds_ratio != 0) %>%
       filter(log_odds_ratio != Inf)

  #x[x[,odds_ratio_col] == Inf, odds_ratio_col] <-
  #      rnorm(length(x[x[,odds_ratio_col] == Inf, odds_ratio_col]),
  #        max(x[x[, odds_ratio_col] != Inf, odds_ratio_col]) /20, 0.1)

  if (inverse_odds_ratio) {
        x[,odds_ratio_col] <- log2(1/x[,odds_ratio_col])
   }
  else {
    x[,odds_ratio_col] <- log2(x[,odds_ratio_col])
  }

   myhjust <- abs(min(x$log_odds_ratio) - cutoff)/(max(x$log_odds_ratio) -
                                                     min(x$log_odds_ratio))
   max.x.pos <- max(x$log_odds_ratio)
   min.x.pos <- min(x$log_odds_ratio)

   p1 <- ggplot(data= x,
         aes(x = log_odds_ratio)) +
         #geom_histogram(alpha = alpha, bins = bins) +
     scale_x_continuous(limits = c(min.x.pos, max.x.pos)) +
          geom_density(color = "grey") +
         theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.title.x = element_text(angle = 0, hjust = myhjust))
   temp <- density(x$log_odds_ratio, adjust = 1)
   temp <- data.frame(x = temp$x, y = temp$y)

if (plot_facet == "by_category" || plot_facet == "by_gene") {
      anno <- category_anno[, c(anno_category_col, anno_gene_col)]
      colnames(anno) <- c("Category", "Gene")
      anno <- merge(anno, x, by = "Gene", all.x = FALSE)

      anno %>%  filter(!is.na(FDR)) %>%
         group_by(Gene) %>%
         summarise(minP = min(FDR, na.rm=TRUE)) -> minP
      anno <- merge(anno, minP)

      if (gene_list_ordered)
      {
          anno <- anno %>% filter(Gene %in% genes_to_highlight)
          anno <- merge(anno,
                        cbind(Gene = genes_to_highlight,
                             gene_orders = 1:length(genes_to_highlight))) %>%
            arrange(gene_orders) %>%
            mutate(Gene = factor(Gene, levels = rev(genes_to_highlight)))
      }
      anno <- anno %>% mutate(highlight_color = ifelse(log_odds_ratio >= 0,
                                               highlight_color_up,
                                               highlight_color_down))
      anno$highlight_color <- as.factor(anno$highlight_color)
      p2 <- anno %>% ggplot( aes(x = log_odds_ratio, y = Gene)) +
          scale_x_continuous(limits = c(min.x.pos, max.x.pos)) +
          geom_point(aes(x = log_odds_ratio, y = Gene,  color = highlight_color)) +
          geom_vline(xintercept  = cutoff, color = dashline_color, linetype= 2) +
          #geom_line(aes(x = log_odds_ratio, y = Gene, color = highlight_color), linetype = 1 ) +
          geom_tile(aes(fill = log_odds_ratio,  color = highlight_color) ) +
          scale_colour_manual(values=c("blue"="blue", "red" = "red", "green"="green",
                                       "purple" = "purple",
                                     "yellow"="yellow","orange"="orange"),
                            breaks=c("blue","red", "green","purple","yellow","orange")) +
          facet_grid(Category ~ ., scales = "free", space = "free") +
          theme(legend.position = "none", strip.text.y = element_text(angle = 0, hjust = 0) ,
             plot.title=element_text(hjust = 1),
             axis.title.x=element_text(angle = 0, hjust = myhjust),
             strip.placement = "outside",
             plot.caption = element_text(hjust = myhjust)) +
             labs(x = "<-      ->",
                caption = "Depleted   Enriched", size = gene_label_size,
                family="Arial", face="bold")
      if (include_FDR)
           p2 <- p2 +  ggtitle("FDR     ") +
             geom_text(aes(x = max.x.pos, y = Gene,
                           label = signif(minP, digits=3)),
                    size = FDR_label_size, vjust = "inward",
                    hjust = "inward")

     if( plot_facet == "by_gene")
            p2 <- p2 + theme(strip.background = element_blank(),
                           strip.text = element_blank())
       #p1/p2 +  plot_layout(heights=  c(1, 10))

       list(p1, p2)
   }
   else {
      p2 <- ggplot(data= temp,
        aes(x = x, y = 1)) +
        geom_tile(aes(fill = y)) +
        scale_fill_gradient(low = "white", high = "black") +
        geom_vline(xintercept  = 0, color = dashline_color, linetype= 2) +
        theme_grey() +
        labs(x = "") +
        scale_y_discrete(expand = c(0,0)) +
        coord_fixed(ratio=1) +
        theme(legend.position = "none", plot.title = element_text(hjust =1),
           axis.ticks = element_blank(), axis.text.x = element_blank(),
           axis.title = element_text(size = gene_label_size),
           axis.title.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
           axis.title.x = element_text(angle = 0, hjust = myhjust))

      p <- list(length = length(genes_to_highlight))
      for (i in 1:length(genes_to_highlight)) {
         highlight_color <- ifelse(x[x[, gene_col] ==
                                      genes_to_highlight[i],
                                    odds_ratio_col] >= 1, highlight_color_up,
                                   highlight_color_down)
         p[[i]] <- p2 + geom_vline(xintercept = x[x[, gene_col] ==
                                                    genes_to_highlight[i],
               odds_ratio_col],
               color = highlight_color) +
               labs(y = genes_to_highlight[i]) +
               geom_text(x = max(temp[,1]), y = 1, size = FDR_label_size,
                   label = signif(min(x[x[, gene_col] ==
                                          genes_to_highlight[i], FDR_col]),
                                  digits=3),
                   vjust = "inward", hjust = "inward") +
           scale_colour_manual(values=c("blue"="blue", "red" = "red", "green"="green",
                                        "purple" = "purple",
                                        "yellow"="yellow","orange"="orange"),
                               breaks=c("blue","red", "green","purple","yellow","orange"))

      }
      p[[1]] <- p[[1]] + ggtitle("FDR       ")
      p[[length(genes_to_highlight)]] <-
          p[[length(genes_to_highlight)]] +
          theme(plot.caption = element_text(hjust = 0.5 ,
                                            size= gene_label_size)) +
          labs(x = "<-      ->", caption = "Depleted   Enriched")

      fig <- paste("p[[1]] ", paste(unlist(lapply(1:length(genes_to_highlight),
          function(i) { paste("/ p[[2]][", i, "]", sep ="")})), collapse= " "))
      list(p1, p, fig)
   }
}
