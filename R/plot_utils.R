#### PCA related ####

#' Plots PC1 by PC2 using ggplot
#'
#' Creates a PCA plot using ggplot2
#'
#' @param transformed_data transformed (e.g., log, vst, rlog) expression data
#' @param sample_map df of sample annotations
#' @param leg_row_num how many rows the leg should be
#' @param gene_num number of genes (ranked by variance) to use
#' @param return_data whether to return the pca data
#'
#' @return ggplot object
#' @export
#'
#' @examples
pca_plotter <- function(transformed_data, sample_map,leg_row_num=3, gene_num=Inf, return_data=FALSE){

  # Filter out any samples not listed in sample_map
  cur_subset_mat <- data.frame(transformed_data)
  cur_subset_mat <- transformed_data[,colnames(transformed_data) %in% rownames(sample_map)]
  #return(cur_subset_mat)
  # Taken from DESeq2 plotPCA function
  # Calculates the row wise variance
  rv <- matrixStats::rowVars(as.matrix(cur_subset_mat))


  # select the gene_num genes by variance
  # the seq_len thing looks weird, but it was in DESeq2 function
  # so leaving it.
  select <- order(rv, decreasing=TRUE)[seq_len(min(gene_num, length(rv)))]
  #return(select)
  # perform a PCA on the data in assay(x) for the selected genes
  fnmt_pcomp <- prcomp(t((cur_subset_mat)[select,]))

  var_exp <- (fnmt_pcomp$sdev^2)/sum(fnmt_pcomp$sdev^2)

  # Puts first 3 pcs into df
  plot_data <- data.frame(pc1=fnmt_pcomp$x[,1],
                          pc2=fnmt_pcomp$x[,2],
                          pc3=fnmt_pcomp$x[,3])
  # sorting because paranoia

  plot_data <- plot_data[sort(rownames(plot_data)),]
  # Getting PCs
  plot_data <- data.frame(pc1=fnmt_pcomp$x[,1],
                          pc2=fnmt_pcomp$x[,2],
                          pc3=fnmt_pcomp$x[,3])
  # merges sample metadata into df by rowname
  plot_data <- merge(plot_data, sample_map, by=0)
  # gets rid of extraneous column
  rownames(plot_data) <- plot_data$Row.names
  plot_data$Row.names <- NULL
  # changes metadata to column name to group
  colnames(plot_data)[4] <- 'group'
  plot_data$group <- as.factor(plot_data$group)
  #eturn(plot_data)
  # This just makes the labels for axises
  axlab.1 <- paste("PC1 (", signif(var_exp[1]*100, digits=4),"%)", sep="")
  axlab.2 <- paste("PC2 (", signif(var_exp[2]*100, digits=4), "%)", sep="")
  axlab.3 <- paste("PC3 (", signif(var_exp[3]*100, digits=4), "%)", sep="")

  if(return_data){
    return(plot_data)
  }
  # And here comes the plot!
  plt1 <- ggplot(data=plot_data,
                 aes(x=pc1,
                     y=pc2,
                     fill=group,
                     colour=group,
                     shape=group,
                     label=row.names(plot_data))) +
    # The geom point aes specificies colouring by group
    # and changes point shape by group as well
    geom_point(size = rel(1.95), aes(shape=factor(group), colour=factor(group))) +
    # geom_point(size = rel(1.5)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    # gets a pretty colour set
    # stat_ellipse(alpha=.15, geom = "polygon") +
    # scale_colour_brewer(palette="Set1") +
    # scale_fill_brewer(palette="Set1") +
    labs(x=axlab.1, y=axlab.2) + theme_bw()

  # This all just setting the themes the way I like it

  plt2 <- plt1 + theme(plot.margin = unit(c(1,1,1,1), "cm"),panel.background = element_blank(), axis.title.y=element_text(size=rel(1.75), face="bold", margin=margin(0,7.5,0,0)), axis.title.x=element_text(size=rel(1.75), face="bold",margin=margin(7.5,0,0,0)),axis.text.y=element_text(size=rel(1.5),colour="black"),axis.text.x=element_text(size=rel(1.5), colour="black"), legend.title=element_blank(),legend.key = element_blank(),legend.text=element_text(size=rel(1.25)),legend.position = 'bottom',panel.border=element_rect(fill=NA,colour="black", size=rel(1.9))) + guides(colour= guide_legend(override.aes = list(size=rel(3.75))))
  #title=element_text(size=22,


  # this just splits the legend into two rows
  # when there is more than 3 groups because of
  # ugly formatting
  # if(length(levels(factor(plot_data$group))) > 3){
  #   plt2 <- plt2 + guides(col=guide_legend(nrow = 2))
  # }
  group_num <- length(levels(factor(plot_data$group)))
  if (group_num > 6){
    plt2 <- plt2 + scale_shape_manual(values = seq(1,group_num))
  }
  plt2 <- plt2 + guides(col=guide_legend(nrow = leg_row_num))
  plt2 <- plt2 + scale_x_continuous(breaks = pretty(plot_data$pc1, n=7)) + scale_y_continuous(breaks = pretty(plot_data$pc2, n=7))


  return(plt2)
}




#### ggplot2 theme and helper functions ####

#' Increases ggplot2 theme size
#'
#' Functions which increases the relative size of
#' fonts used in ggplot2. It also makes some other
#' adjustments like increasing margin size, ticks,
#' etc.
#' @param rel_size
#'
#' @import ggplot2
#'
#' @return
#' @keywords internal
#'
#' @examples
theme_large_font <- function(rel_size=1.5){
  out_theme <- theme(plot.title = element_text(size=rel(rel_size),
                                               margin = margin(t = 0, r = 0, b = .5, l = 0, unit = 'line')),
                     # title = element_text(vjust = .5, size=rel(rel_size)),
                     title = element_text( size=rel(rel_size)),
                     axis.title.x = element_text(margin = margin(t = 1.5, r = 0, b = 0, l = 0, unit = 'line')),
                     axis.title.y = element_text(margin = margin(t = 0, r = .5, b = 0, l = 0, unit = 'line')),
                     axis.text = element_text(size=rel(rel_size)),
                     axis.text.x = element_text( vjust = -.75),
                     axis.text.y = element_text( hjust = 1),
                     strip.text = element_text(size=rel(rel_size)),
                     ## Legend stuff
                     legend.position = 'bottom',
                     legend.text=element_text(size=rel(rel_size)),
                     legend.background = element_blank(),
                     legend.title = element_blank(),
                     legend.key = element_blank(),
                     ## Axis stuff
                     axis.line = element_line(colour='black', size = rel(rel_size )),
                     axis.ticks = element_line(colour='black', size = rel(rel_size )),
                     axis.ticks.length = unit(.2, "cm"),
                     ## Panel stuff
                     #panel.grid.major = element_line(colour='gray', size = rel(rel_size - .25)),
                     panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()
  )
  return(out_theme)
}

#' Returns expanded coordinates of a ggplot
#'
#' This returns coordinates which will expand
#' the plotting area when used with
#'
#' @param plt ggplot object
#' @param prop proportion to exand the x & y limits
#'
#' @return list with x & y limits as vector
#' @export
#'
#' @examples
#' \dontrun{
#' xy_lims <- increase_xy_lims(plt, .15)
#' plt + coord_cartesian(xlim = xy_lims$x_lims, ylim = xy_lims$y_lims)
#' }
#'
increase_xy_lims <- function(plt, prop=.25){

  lims <- list(c(min(plt$data$x), max(plt$data$x)), c(min(plt$data$y), max(plt$data$y))
  ) %>%
    lapply(., function(x){
      temp_lims <-  x
      temp_lims[1] <- temp_lims[1] - abs(mean(x)*prop)
      temp_lims[2] <- temp_lims[2] + abs(mean(x)*prop)
      temp_lims

    }) %>% setNames(., c('x_lims', 'y_lims'))
  return(lims)
}

