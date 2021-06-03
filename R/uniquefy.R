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


#### TODO Functions ####

