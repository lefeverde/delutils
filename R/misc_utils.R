##### General data wrangling #####


#' Removes un-informative cols from a data.frame
#'
#' @param dat a data.frame
#' @param n_levels thresholded number of items to keep column
#'
#' @return a data.frame without un-informative columns
#' @export
#'
#' @examples
#' library(delutils)
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


#' Convert Column Names to Lowercase, Replace Whitespace, and Remove Leading/Trailing Whitespace
#'
#' This function converts the column names of a data frame to lowercase,
#' replaces whitespace (including tabs) with underscores, and removes any leading or trailing whitespace.
#'
#' @param dat A data frame.
#' @return A data frame with cleaned column names.
#' @examples
#' # Example data frame
#' my_df <- data.frame("First Name\t" = c("John", "Jane", "Alice"),
#'                     "Last Name  " = c("Doe", "Smith", "Johnson"),
#'                     "Age " = c(25, 30, 35))
#'
#' # Clean column names: convert to lowercase, replace whitespace, and remove leading/trailing whitespace
#' cleaned_df <- clean_colnames(my_df)
#' colnames(cleaned_df)
#' @importFrom stringr str_replace_all str_trim
#' @export
clean_colnames <- function(dat) {
  if(!("data.frame" %in% class(dat))){
    stop("df needs to be either a data.frame or tibble")
  }
  fixed_names <- colnames(dat)


  fixed_names <- colnames(dat) %>%
    str_trim(side = "both") %>%
    gsub("\\s+", "_", .) %>%
    tolower()

  colnames(dat) <-fixed_names
  return(dat)
}

#' Finds rows in a data.frame which contain NA values
#'
#' @param df data.frame or tibble
#'
#' @return indexes of rows with NA values
#' @export
#'
#' @examples
#' library(delutils)
#' dat <- data.frame(A=c(NA,2:4), B=c(1:3, NA), C=1:4)
#' find_na_rows(dat)
#'
find_na_rows <- function(df){
  if(!("data.frame" %in% class(df))){
    stop("df needs to be either a data.frame or tibble")
  }

  weird_rows <-
    lapply(df, function(x){
      which(is.na(x))
    }) %>% unlist(use.names = FALSE)

  weird_rows <-
    unique(weird_rows)
  return(weird_rows)
}

#' Finds rows in a data.frame which contain Inf values
#'
#' @param df data.frame or tibble
#'
#' @return A vector of unique row indices containing Inf values
#' @export
#'
#' @examples
#' # Example usage:
#' dat <- data.frame(A = c(Inf, 2:4), B = c(1:3, Inf), C = 1:4)
#' find_inf_rows(dat)
find_inf_rows <- function(df){
  if(!("data.frame" %in% class(df))){
    stop("df needs to be either a data.frame or tibble")
  }

  inf_rows <-
    lapply(df, function(col){
      which(is.infinite(col))
    }) %>% unlist(use.names = FALSE)

  inf_rows <- unique(inf_rows)
  return(inf_rows)
}


#' Uniquefy by variance
#'
#' Uniquefies an matrix by the row variance. This
#' ranks variance for each as specifified in the
#' map data.frame. The row.names of the expression
#' matrix and map need to have at least some in common.
#' Function first identifies which row.names are shared
#' and then subsets.
#' Idea for using max var to de-duplicate probes came
#' form here: https://www.biostars.org/p/51756/#51875
#'
#'
#' @param mat expression matrix
#' @param map gene map
#' @param dup_col col with potential duplicates
#'
#' @return filtered expression matrix
#' @export
#'
#' @examples
uniquefy_by_variance <- function(mat, map, dup_col=NULL){

  if(!any(class(mat) %in% c('data.frame','matrix'))){
    stop('expr_mat needs to be either a data.frame or matrix of gene expression value.')
  }

  if(is.null(dup_col)){
    dup_col <- names(map)[1]
    warning(paste0('Using map ', dup_col, ' as dup_col'))

  }

  if(!any(row.names(mat) %in% row.names(map))){
    stop('None of the row.names matched between mat and map.')
  }
  # Sorts and keeps only rows shared between map and mat
  common_ids <- intersect(row.names(mat), row.names(map))
  mat <- mat[common_ids,]
  map <- map[common_ids,dup_col]

  # Gets row vars and filters
  rv_mat <- data.frame(row.names = common_ids,
                       group=map,
                       vals=matrixStats::rowVars(mat))
  rv_mat <- uniquefy_base(rv_mat, 'group', 'vals')
  filt_mat <- mat[rv_mat,]
  return(filt_mat)

}


#' Return a unique data.frame by max value
#'
#' @param dat data.frame with 1st column as string and second as numeric
#' @param group_col string denoting group column (col with duplicates)
#' @param val_col string denoting value column
#'
#' @return unique data.frame
#' @export
#'
#' @examples
#' set.seed(42)
#' dat <- data.frame(group=c('a', 'a', 'a', 'b', 'c', 'c'), vals=runif(6))
#' dat$vals[6] <- -4
uniquefy_by_abs_max <- function(dat, group_col=NULL, val_col=NULL){
  #TODO make this workhorse functions which
  # and split the variance and other uniquefiers
  # into compatible functions
  if(is.null(group_col)){
    group_col <- names(dat)[1]
  }
  if(is.null(val_col)){
    val_col <- names(dat)[2]
  }
  temp_dat <- data.frame(dat)
  temp_dat[val_col] <- abs(temp_dat[val_col])
  rows_to_keep <- uniquefy_base(temp_dat, group_col, val_col)
  dat <- dat[rows_to_keep,]
  return(dat)
}





#' Base function which does the uniquefying
#'
#' This is meant to be an internal function
#' which is the commonly used across the
#' different uniquefiers. It works by
#' taking uniquefying based on largest value
#' and returning a logical vector
#'
#' @param dat data frame
#' @param group_col group_col string denoting group column (col with duplicates)
#' @param val_col string denoting value column
#'
#' @return logical vector
#' @keywords internal
#'
#' @examples
uniquefy_base <- function(dat, group_col, val_col){

  ordered_dat <- dat[order(dat[group_col], -dat[val_col]),]
  ordered_dat <- ordered_dat[!duplicated(ordered_dat[group_col]),]
  logical_idxs <- row.names(dat) %in% row.names(ordered_dat)

  return(logical_idxs)
}

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

#' de_rbind_df
#'
#' @param df data.frame
#' @param col_to_split columnt to split data.frame by
#'
#' @return
#' @export
#'
#' @keywords internal
#'
#' @examples
de_rbind_df <- function(df, col_to_split=1){

  e_message <- 'col_to_split needs to be one of a "character", "numeric", or "integer" which matches a column in the data.frame.'
  if(!class(col_to_split) %in% c("character", "numeric", "integer")){
    stop(e_message)
  }
  df_nms <- colnames(df)
  df_idxs <- seq(1, ncol(df))
  valid_col <- col_to_split %in% df_nms || col_to_split %in% df_idxs
  if(!valid_col){
    stop(e_message)
  }

  splits_vec <- df[, col_to_split]
  df[,col_to_split] <- NULL
  split_list <- split(df, splits_vec)
  return(split_list)



}

