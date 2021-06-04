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
