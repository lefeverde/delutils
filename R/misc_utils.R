#' Removes un-informative cols from a data.frame
#'
#' @param dat a data.frame
#' @param n_levels thresholded number of items to keep column
#'
#' @return a data.frame without un-informative columns
#' @export
#'
#' @examples
#' library(lefutils)
#' dat <- data.frame(A=rep(1,4), B=seq(1,4), C=c('f1', 'f2', 'f1', 'f1'))
#' cleaned_dat <- clean_by_levels(dat)
clean_by_levels <- function(dat, n_levels=1){
  if(! 'data.frame' %in% class(dat)){
    emsg <- paste0('dat needs to be a data.frame not a ', class(dat))
    stop(emsg)
  }

  stopifnot
  if( !is.numeric(n_levels)){
    if(!is.integer(n_levels)){
      emsg <- paste0('n_levels needs to be an integer >= 1')
      stop(emsg)
    }
  }
  if( (n_levels %% 1 != 0) |  (n_levels < 1)){
    emsg <- paste0('n_levels needs to be an integer >= 1')
    stop(emsg)
  }

  nitems_per_col <- sapply(dat, simplify = TRUE, function(x){
    length(unique(x))
  })
  cols_to_keep <- which(nitems_per_col > n_levels)
  cleaned_df <- dat[,cols_to_keep, drop=FALSE]
  return(cleaned_df)
}

#' Finds rows in a data.frame which contain NA values
#'
#' @param df data.frame or tibble
#'
#' @return indexes of rows with NA values
#' @export
#'
#' @examples
#'
find_na_rows <- function(df){
  if(!("data.frame" %in% class(df))){
    stop("df needs to be either a data.frame or tibble")
  }

  weird_rows <-
    apply(df, 2, function(x){
      which(is.na(x))
    }) %>% unlist

  weird_rows <-
    unique(weird_rows)
  return(weird_rows)
}

