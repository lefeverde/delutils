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


