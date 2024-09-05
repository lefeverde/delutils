#### Table plots and functions ####

#' Makes a flextable that can be (sometimes, atleast) used as-is
#'
#' This function creates a flextable using customizations that I generally like.
#' It returns a flextable object which can be further customized.
#'
#' @param dat_tbl data.frame or tibble to be turned into a flextable
#'
#' @return flextable object
#' @export
#'
#' @import flextable
#' @importFrom officer fp_border
#'
#' @examples
make_flextable <- function(dat_tbl){
  std_border <-
    fp_border(color="gray", width = 1)

  ft_tbl <- dat_tbl %>%
    flextable(.) %>%
    align(., align = 'center', part='all') %>%
    font(fontname='Calibri', part='all') %>%
    fontsize(size=12, part = 'all') %>%
    color(., i=1, color='black', part = 'header') %>%
    bold(., part='all') %>%
    border_inner_h(., border = std_border) %>%
    border_inner_v(., border = std_border) %>%
    border_outer(., border = std_border) %>%
    theme_booktabs %>%
    autofit(., -.1, -.1)
  return(ft_tbl)
}


#' Colors a vector of DNA/peptide sequences via html
#'
#' @param sequences vector of DNA or peptide sequences
#'
#' @importFrom crayon make_style
#'
#' @return
#' @export
#'
#' @examples
#' #' library(knitr)
#' library(kableExtra)
#'
#' # Example 1: DNA sequences
#' dna_sequences <- c("ACTG", "GCTA", "TACG")
#' colored_dna_sequences <- color_dna_sequence(dna_sequences)
#'
#' # Creating a table with kable for DNA sequences
#' dna_data <- data.frame(
#'   DNA_Sequence = dna_sequences,
#'   Colored_DNA_Sequence = colored_dna_sequences
#' )
#' dna_data %>%
#'   kable('html', escape = FALSE) %>%
#'   kable_styling('striped', full_width = FALSE)
#'
#'
#' # Example 2: Peptide sequences
#' peptide_sequences <- c("MKTAY", "LLGTP", "IVAGL")
#' colored_peptide_sequences <- color_dna_sequence(peptide_sequences)
#'
#' # Creating a table with kable for peptide sequences
#' peptide_data <- data.frame(
#'   Peptide_Sequence = peptide_sequences,
#'   Colored_Peptide_Sequence = colored_peptide_sequences
#' )
#' peptide_data %>%
#'   kable('html', escape = FALSE) %>%
#'   kable_styling('striped', full_width = FALSE)
#'
color_dna_sequence <- function(sequences){


  # Define custom styles
  whiter <- make_style(rgb(1, 1, 1))  # Pure white color
  dark_grey_bg <- make_style(rgb(0.5, 0.5, 0.5), bg = TRUE)  # Dark grey background

  # All the IUPAC ambiguity letters minus N
  dark_grey_bg_letters <- c("M", "R", "W", "S", "Y", "K", "V", "H", "D", "B")

  # Custom color scheme
  # Define custom HTML styles instead of crayon
  color_scheme <- c(
    A = '<span style="color:black; background-color:rgb(255, 128, 128); padding:0;">A</span>',
    C = '<span style="color:black; background-color:rgb(128, 255, 128); padding:0;">C</span>',
    G = '<span style="color:black; background-color:rgb(128, 255, 255); padding:0;">G</span>',
    T = '<span style="color:black; background-color:rgb(255, 204, 128); padding:0;">T</span>',
    U = '<span style="color:black; background-color:rgb(255, 255, 128); padding:0;">U</span>',
    M = '<span style="color:white; background-color:rgb(128, 128, 128); padding:0;">M</span>',
    R = '<span style="color:white; background-color:rgb(128, 128, 128); padding:0;">R</span>',
    W = '<span style="color:white; background-color:rgb(128, 128, 128); padding:0;">W</span>',
    S = '<span style="color:white; background-color:rgb(128, 128, 128); padding:0;">S</span>',
    Y = '<span style="color:white; background-color:rgb(128, 128, 128); padding:0;">Y</span>',
    K = '<span style="color:white; background-color:rgb(128, 128, 128); padding:0;">K</span>',
    V = '<span style="color:white; background-color:rgb(128, 128, 128); padding:0;">V</span>',
    H = '<span style="color:white; background-color:rgb(128, 128, 128); padding:0;">H</span>',
    D = '<span style="color:white; background-color:rgb(128, 128, 128); padding:0;">D</span>',
    B = '<span style="color:white; background-color:rgb(128, 128, 128); padding:0;">B</span>',
    N = '<span style="color:white; background-color:grey; padding:0;">N</span>'
  )

  # Function to color DNA sequences
  color_dna <- function(dna_string) {
    for (nuc in names(color_scheme)) {
      dna_string <- gsub(nuc, color_scheme[[nuc]], dna_string, fixed = TRUE)
    }
    return(dna_string)
  }
  return(color_dna(sequences))

}





#### Venn diagrams ####
#' Venn maker
#'
#' @param set_list named vector list
#' @param tit plot title
#'
#' @return
#' @keywords internal
#'
#' @import eulerr
#'
#' @examples
venn_maker <- function(set_list, tit=''){
  if(is.null(names(set_list))){
    stop(call. = TRUE, 'set_list needs to be named')
  }
  euler(set_list) %>%
    plot(.,
         fills=list(fill=c('red', 'blue'),
                    alpha=.5),
         col=c('red', 'blue'),
         quantities = list(cex = 1.125),
         fontsize = 14,
         text_args = list(font = 20),
         legend=list(cex=1.5, alpha=1))
}





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
    return(list(plot_data, axlab.1, axlab.2, axlab.3))
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

