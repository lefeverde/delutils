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
#' library(delutils)
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
  item_class_checker <-
    vapply(df_list, function(x){
      any(class(x) %in% c('data.frame', 'NULL'))
    }, logical(1)) %>%
    all(.)
  if(! item_class_checker){
    stop('df_list contains items which are not data.frame or tibble')
  }



  out_df <- list()
  empty_dfs <- NULL
  for(x in nms){
    temp_pth <- df_list[[x]]
    valid_conditions <- c(
      nrow(temp_pth) >= 1, # more than 1 rows
      !is.null(temp_pth) # df is not NULL
    ) %>% all
    # Check valid conditions
    if(valid_conditions){
      temp_pth <- cbind(x,
                        temp_pth,
                        stringsAsFactors = FALSE)
      row.names(temp_pth) <- NULL
      names(temp_pth)[1] <- sum_col_name
      out_df[[x]] <- temp_pth
    } else {
      empty_dfs <- c(empty_dfs, x)
      next
    }
  }
  if(!is.null(empty_dfs)){
    warning_str <-
      paste0(head(empty_dfs), collapse ='\n') %>%
      paste0('Some items contained empty or otherwise disagreeable data.frames and were not included in the bound data.frame.\nThese include:\n', . )
    warning(warning_str, sep='')
  }

  n_cols <- lapply(out_df, function(x){
    colnames(x) %>%
      length
  }) %>% do.call(c, .)

  if(length(unique(n_cols)) > 1){
    suspect_dfs <-  names(out_df)[which(n_cols != mode(n_cols))]
    e_message <-
      paste0(suspect_dfs, collapse='\n') %>%

      paste0('data.frames have different number of columns. check:\n', . )
    stop(e_message, sep='')
  }

  # TODO wrap this trycatch when data.frames have
  # different columns, s.t., an error w/ offending
  # items is spit out
  out_df <- do.call(rbind, c(out_df, stringsAsFactors=FALSE))
  row.names(out_df) <- NULL
  return(out_df)
}

