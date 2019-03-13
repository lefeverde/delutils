#' rbinds a named list of data.frames
#'
#'
#' rbinds together data.frames in which the list name
#' is turned into a column and bound to the data.frame
#' becoming the first column of the data.frame
#'
#' @param df_list named list of data.frames
#' @param sum_col_name colname of new column bound to data.frame
#'
#' @return a single data.frame row bound data.frame
#' @export
#'
#' @examples
#' library(dplyr)
#' library(lefutils)
#' l <- list(data.frame(A=seq(1,5)), data.frame(A=seq(6, 7))) %>%
#'  setNames(., c('Cond1', 'Cond2'))
#'  rbind_named_df_list(l)
#'
#' # Doesn't work because list contains matrix
#' \dontrun{
#' l <- list(data.frame(A=seq(1,5)),
#'  data.frame(A=seq(6, 7)),
#'  matrix(data=NA, ncol=1, nrow=5)) %>%
#'  setNames(., c('Cond1', 'Cond2', 'Cond3'))
#'  rbind_named_df_list(l)
#' }
rbind_named_df_list <- function(df_list, sum_col_name='col_from_list'){
  # Makes sure list has names, everything is
  # either a data.frame or NULL
  nms <- names(df_list)
  if(is.null(nms)){
    stop('df_list needs to be a named list')
  }
  item_class_checker <- sapply(df_list, simplify = TRUE, function(x){
    class(x) %in% c('data.frame', 'NULL')
  }) %>%
    all(.)
  if(! item_class_checker){
    stop('df_list contains items which are not data.frame')
  }

  out_df <- list()
  empty_dfs <- NULL
  for(x in nms){
    temp_pth <- data.frame(df_list[[x]])
    # Checks if cur df is null or empty
    if(nrow(temp_pth) > 0 | is.null(temp_pth)){
      temp_pth <- cbind(x, temp_pth)
      row.names(temp_pth) <- NULL
      names(temp_pth)[1] <- sum_col_name
      out_df[[x]] <- temp_pth
    } else {
      empty_dfs <- c(empty_dfs, x)
      next
    }
  }
  if(!is.null(empty_dfs)){
    warning_str <- paste('Some items contained empty or otherwise disagreeable data.frames and were not included in the bound data.frame. These include:\n',
                         paste0(head(empty_dfs), sep='\n'))
    warning(warning_str)
  }
  # TODO wrap this trycatch when data.frames have
  # different columns, s.t., an error w/ offending
  # items is spit out
  out_df <- do.call(rbind, out_df)
  row.names(out_df) <- NULL
  return(out_df)
}

