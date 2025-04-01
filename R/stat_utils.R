#### Statisitc related wranglers ####

#' Creates a much nicer formatted design matrix than the default model.matrix
#'
#' I made this function because the way the model.matrix function
#' creates a desing matrix, can cause issues with differential
#' gene expression tools, such as edgeR, DESeq2,
#' and limma. If the intercept column is included,
#' it first removes the parens and parens from the
#' intercept column and then removing the variable
#' name from the colnames and replaces them with
#' factor level vs reference level syntax. If the
#' design is specified without an intercept, variable
#' concatentation is simply removed.
#'
#' @param ... formula with outcome as first variable
#' @param data data.frame
#' @param show_warnings logical whether to show wanrings
#'
#' @return model.matrix
#' @export
#'
#' @importFrom limma strsplit2
#' @importFrom stats model.matrix
#'
#'
#' @examples
better_model_matrix <- function(..., data, show_warnings=TRUE){
  m <- model.matrix(..., data=data)

  # Getting variables (columns) from formula
  mod_str <- deparse(...)

  var_str <-
    strsplit2(mod_str, '\\+') %>%
    drop(.) %>%
    trimws(.) %>%
    gsub('~', '', .) #%>%


  # cleaning model matrix
  if(var_str[1] == "0"){
    # checks if intercept is 1st term
    outcome_var <- var_str[2]
  } else {
    outcome_var <- var_str[1]
  }
  if(show_warnings){
    w <- paste0('\n', outcome_var, ' used as outcome variable.\nIf this is incorrect, reorder the formula so that the outcome variable comes first.'   )
    warning(w)
  }
  df_str <- deparse(substitute(data))
  colnames(m) <-
    gsub(outcome_var, '', colnames(m)) %>%
    # gsub(df_str, '', .) %>%
    gsub('\\(', '', .) %>%
    gsub('\\)', '', .)

  model_has_intercept <-
    colnames(m) %>%
    grepl('Intercept', ., ignore.case = TRUE) %>%
    any(.)

  if(model_has_intercept){
    colnames_to_prettify <- levels(data[,outcome_var])
    ref_level <- colnames_to_prettify[1]
    colnames_to_prettify <- colnames_to_prettify[-1]
    stopifnot(all(colnames_to_prettify %in% colnames(m)))
    pretty_colnames <- sapply(colnames_to_prettify, function(x){
      paste0(x, '_vs_', ref_level)
    })
    col_idxs <- match(colnames_to_prettify, colnames(m))
    colnames(m)[col_idxs] <- pretty_colnames

  }
  return(m)
}







#' Selects the top n genes by variance
#'
#' @param expr_mat expression mat
#' @param n_genes number of genes to extract
#'
#' @return gene matrix
#' @export
#'
#' @importFrom matrixStats rowVars
#'
#' @examples
select_genes_by_variance <- function(expr_mat, n_genes=5000){
  rv_idx <- rowVars(expr_mat) %>%
    order(., decreasing = TRUE)
  rv_idx <- rv_idx[1:n_genes]
  filt_expr_mat <- expr_mat[rv_idx,]
  return(filt_expr_mat)
}


#' Gene expression filter by quantile
#'
#' This function first removes rows in which
#' the max value equals the min value and
#' then filters the resulting genes by row
#' medians. Rows with a median greater than
#' the resulting value of the quantile are
#' retained.
#'
#' @param expr_mat an expression matrix of class matrix
#' @param quant_val quantile threshold
#'
#' @return filtered by quantile expresion matrix
#' @export
#'
#' @importFrom stats quantile
#' @importFrom Biobase rowMax rowMedians
#'
#' @examples
gene_quantile_filter <- function(expr_mat, quant_val=.05){
  rm <- rowMax(expr_mat)
  qthresh <- quantile(rm, quant_val)
  filt_mat <- expr_mat[rm > qthresh,]
  filt_mat <- filt_mat[rowMedians(filt_mat) > qthresh,]
  return(filt_mat)
}









#' Creates matrix of overlap ratios
#'
#' This function creates a matrix of the ratio
#' of overlapping genesets. It takes a list of genes
#' or alternatively 2 gene lists
#'
#' @param l1 gene list
#' @param l2 optional 2nd gene list
#'
#' @return top diagonal matrix of overlap ratios
#' @export
#'
#' @examples
#' \dontrun{
#' # returns a self (l1 vs l1) overlap matrix
#' genes <- paste0('gene', seq(1, 100))
#' set.seed(42)
#' n_genes_to_sample <- rnbinom(10, 10, .5)
#' l1 <- lapply(seq_along(n_genes_to_sample), function(i){
#'  n <- n_genes_to_sample[i]
#'  sample(genes, n)
#' }) %>%
#'  setNames(., paste0('L1_', LETTERS[1:10]))
#' omat <- get_overlap_matix(l1)
#' }
#'
get_overlap_matix <- function(l1, l2=NULL){
  if(is.null(l2)){
    l2 <- l1
  }

  l1_nms <- names(l1)
  l2_nms <- names(l2)


  w <- matrix(NA, nrow=length(l1_nms), ncol=length(l2_nms))
  rownames(w) <- l1_nms
  colnames(w) <- l2_nms

  for (i in 1:nrow(w)) {
    for (j in 1:ncol(w)) {
      w[i,j] = overlap_ratio(l1[[i]], l2[[j]])
    }
  }
  return(w)
}


#' Extract genes sets from a pathway results data.frame
#'
#' pathway results data.frame with a column
#' of gene sets seperated by a delimiter
#'
#' @param path_df data.frame with gene sets
#' @param delim delimeter seperating genes (default: '/')
#'
#' @return list with names as pathway ids
#' @export
#'
#' @examples
get_geneSets <- function(path_df, delim="/"){
  nms <- path_df$ID
  gene_sets <-  strsplit(path_df$geneID, delim, fixed=TRUE)
  names(gene_sets) <- nms
  return(gene_sets)
}


#' get results from limma fit object
#'
#' Takes a limma fit object created using either \link[limma]{eBayes} or
#' \link[limma]{treat} and returns a long data.frame of the results.
#' If the `fit` object contains `lods`, then results are obtained using
#' \link[limma]{topTable}. Else, if the fit object contains a `treat.lfc`
#' topTreat is used. If no coefficients are passed, it tries to guess
#' which are the relevant ones checking which are made up of 1's and 0's.
#'
#' @param limma_fit fit object produced by limma
#' @param coefs coefficients to return
#' @param print_summary logical, whether print summary of the results
#'
#' @return data.frame of results in long(ish) format
#' @export
#'
#' @import tibble dplyr
#' @importFrom limma topTable
#'
#'
#' @examples
get_limma_results <- function(limma_fit, coefs=NULL, print_summary=TRUE){
  #TODO add tests
  #TODO think about handling levels with whitespace
  #TODO maybe automagically figure type of gene annots
  d <- limma_fit$design
  if(is.null(coefs)){
    d <- limma_fit$design
    one_hot_cols <- sapply(seq(1,ncol(d)), function(i){
      all(d[,i] %in% c(1,0))
    })
    coefs <- colnames(d)[one_hot_cols]
    coefs <- coefs[!grepl('intercept', coefs, ignore.case = T)]
  }

  res_list <- lapply(coefs, function(x){
    if('lods' %in% names(limma_fit)){
      cur_res <-
        topTable(limma_fit, coef = x, number = Inf, confint=TRUE)
    }
    if('treat.lfc' %in% names(limma_fit)){
      cur_res <-
        topTreat(limma_fit, coef = x, number = Inf, confint=TRUE)
    }
    cur_res <- cur_res %>%
      tibble::rownames_to_column(., 'ensembl_gene_id') %>%
      cbind(coefficient=x, .)
  }) %>% do.call(rbind, .)

  if(print_summary){
    expression_summaries(res_list) %>%
      print
  }
  return(res_list)
}





#' Returns DESeq2 results in long format
#'
#' @param dds_object DDS object after running DESeq
#' @param use_shrinkage logical whether to perform lfcshrinkage
#' @param gene_annots Gene annotations which if not given are assumed
#'
#' @return data.frame of results in long(ish) format
#' @export
#'
#' @importFrom SummarizedExperiment mcols
#' @importFrom DESeq2 resultsNames results
#'
#' @examples
get_deseq2_results <- function(dds_object, gene_annots=NULL, use_shrinkage=FALSE){
  if(is.null(gene_annots)){
    gene_annots <-
      mcols(dds_object) %>%
      data.frame(.) %>%
      dplyr::select_if(., is.character)

    warning_msg <- paste0('gene annots not given.\nUsing: ',
                          paste(names(gene_annots), collapse = ', '),
                          ' for gene annotations')
    warning(warning_msg)

  }

  # Get results by coefs which
  # have `_vs_` in them
  results_list <- list()
  res_names <- DESeq2::resultsNames(dds_object)
  res_names <- res_names[grepl('_vs_', res_names)]

  for(i in  res_names){
    # Checks whether to get results via shrinkage
    if(use_shrinkage){
      temp_res <- lfcShrink(dds_object,
                            coef = i,
                            type = 'normal',
                            parallel=TRUE)
    } else if(!use_shrinkage){
      temp_res <- results(dds_object, name = i)
    }
    temp_res <-
      temp_res %>%
      data.frame(.) %>%
      rownames_to_column(., var = "ensembl_gene_id") %>%
      data.frame(coefficient=i, .)

    results_list[[i]] <- temp_res
  }
  names(results_list) <- NULL
  res_df <- do.call(rbind,results_list)
  rn <- res_df$ensembl_gene_id
  res_df$row <- NULL
  sorted_annots <- gene_annots[match(rn, row.names(gene_annots)),,drop=FALSE]
  out_df <- data.frame(ensembl_gene_id=rn, sorted_annots, res_df)
  row.names(out_df) <- NULL
  return(out_df)
}



#' Returns summaries of gene expression
#'
#'
#' This just summarizes the number of significant up & down
#' regulated genes per group
#'
#' @param long_res_table res df in long format
#' @param log2_thresh threshold to genes by lfc
#' @param logFC_col str of logFC_col
#' @param padj_col str of padj_col
#'
#' @return table summary
#' @export
#'
#' @examples
expression_summaries <- function(long_res_table, log2_thresh=0, logFC_col = NULL, padj_col=NULL){
  #TODO fix hardcoding of coefficient
  if(is.null(logFC_col)){
    potential_cols <- c('log2FoldChange','logFC')
    logFC_col <- na.omit(match(potential_cols, colnames(long_res_table)) )
    colnames(long_res_table)[logFC_col] <- 'log2FoldChange'
  }
  if(is.null(padj_col)){
    potential_cols <- c('padj','adj.P.Val')
    logFC_col <- na.omit(match(potential_cols, colnames(long_res_table)) )
    colnames(long_res_table)[logFC_col] <- 'padj'
  }
  group_sums <-
    long_res_table %>%
    dplyr::filter(padj < .05) %>%
    dplyr::filter(.,abs(log2FoldChange) >= log2_thresh) %>%
    group_by(., coefficient) %>%
    summarise(.,
              up=sum(log2FoldChange > 0),
              down=sum(log2FoldChange < 0))
  return(group_sums)
}



#### Statistical functions ####



#' Function performs lm over accross a number of genes
#'
#' @param expr_mat matrix of expression values
#' @param p_data matrix with sample data
#' @param model_formula a string of the form: y ~ x + o
#'
#' @return lm list
#' @export
#'
#' @examples
basic_linear_regression <- function(expr_mat, p_data, model_formula){
  if(class(model_formula) != 'character'){
    stop('model_formula needs to be a string of the form y ~ + x1 + x2...')
  }
  expr_mat <- data.frame(expr_mat[row.names(p_data),], check.names = F)
  l <- lapply(expr_mat, function(x){
    tmp <- summary(lm(as.formula(model_formula), data=p_data))
  })
  out_df <- lmlist_to_df(l)
  # nms <- names(l)
  # x_list <- lapply(l, function(x){
  #   x <- x$coefficients[2,]
  # })
  # out_df <- do.call(rbind, x_list)
  # row.names(out_df) <- nms
  return(out_df)
}

#' Performs logistic regression with numerous Y vals
#'
#' @param expr_mat matrix of expression values
#' @param p_data matrix with sample data
#' @param model_formula a string of the form: y ~ x + o
#'
#' @return lm list
#' @export
#'
#' @examples
basic_logistic_regression <- function(expr_mat, p_data, model_formula){
  ### Args:
  ### expr_mat: matrix of expression values
  ### p_data: matrix with sample data
  ### model_formula: needs to be y ~ x + o
  ### y is outcome var, x is gene, o is other covariates

  if(class(model_formula) != 'character'){
    stop('model_formula needs to be a string of the form y ~ + x1 + x2...')
  }
  expr_mat <- data.frame(expr_mat, check.names = F)
  l <- lapply(expr_mat, function(x){
    tmp <- summary(glm(as.formula(model_formula), family = 'binomial', data=p_data))
    tmp <- coefficients(tmp)
    if(nrow(tmp) == 1){
      temp_row <- matrix(rep(NA, 4), nrow = 1)
      colnames(temp_row) <- colnames(tmp)
      x <- temp_row
    } else {
      x <- tmp[2,]
    }
  })
  res_logreg <- data.frame(do.call(rbind, l), check.names = F)
  colnames(res_logreg) <- c("Estimate","Std_Error","z_or_t_value", "p_value" )
  row.names(res_logreg) <- names(l)
  return(res_logreg)
}




#' Compare the sample contributions according to their annotation level across the components.
#'
#' Wilcoxon or Kruskal-Wallis tests are performed depending on the number of levels in the considered annotation.
#' @title Comparison of distributions of sample groups
#' @param A A matrix of dimensions 'samples x components' containing the sample contributions
#' @param annot A matrix of dimensions 'samples x variables' containing the sample annotations
#' @param colAnnot The name of the column of \code{annot} to be considered
#' @return A vector of p-values
#' @author Anne Biton
#' @seealso \code{wilcox.test}, \code{kruskal.test}
#' @keywords external
#' @import foreach
#' @import doParallel
wilcoxOrKruskalOnA <- function (A, colAnnot, annot) {

  comp <- NULL
  A <- A[rownames(annot),]

  res_tests <- foreach::foreach(comp=as.list(A),.combine = c, .errorhandling = "stop") %dopar% {
    annotComp <- data.frame(comp=comp)
    annotComp[[colAnnot]] <- as.factor(annot[[colAnnot]])

    if (length(levels(annotComp[[colAnnot]])) == 2)
      res.test <- wilcox.test(comp~eval(as.name(colAnnot)), data = annotComp, na.action = "na.omit")
    else
      res.test <- kruskal.test(comp~eval(as.name(colAnnot)), data = annotComp, na.action = "na.omit")
    return(res.test$p.value)
  }
  return(unlist(res_tests))

}

#### Distance functions ####

#' Calculates the overlap coefficient of 2  sets
#'
#'
#' @param x set vector
#' @param y set vector
#'
#' @return overlap ratio/coefficient
#' @export
#'
#'
overlap_ratio <- function(x, y, check_sets=TRUE) {
  if(check_sets){
    stopifnot("x contains duplicated values"=!any(unclass(table(x)) > 1))
    stopifnot("y contains duplicated values"=!any(unclass(table(y)) > 1))
  }
  numer <- length(intersect(x, y))
  denom <- min(length(x), length(y))

  return(numer/denom)
}

#' Calculates the jaccard distance between 2 named sets
#'
#' @param x set vector
#' @param y set vector
#'
#' @return jaccard distance
#' @export
#'
#' @examples
jaccard <- function(x, y, check_sets=TRUE){
  if(check_sets){
    stopifnot("x contains duplicated values"=!any(unclass(table(x)) > 1))
    stopifnot("y contains duplicated values"=!any(unclass(table(y)) > 1))
  }

  numer <- length(intersect(x, y))
  denom <- length(union(x, y))

  return(numer/denom)
}

#' Calculates the Dice coefficient
#'
#' @param x
#' @param y
#' @param check_sets
#'
#' @return
#' @export
#'
#' @examples
dice <- function(x, y, check_sets=TRUE){
  if(check_sets){
    stopifnot("x contains duplicated values"=!any(unclass(table(x)) > 1))
    stopifnot("y contains duplicated values"=!any(unclass(table(y)) > 1))
  }
  numer <- 2*length(intersect(x, y))
  denom <- length(x) + length(y)

  return(numer/denom)
}

#' Order DNA sequences by Levenshtein edit distance
#'
#' This function takes a vector of DNA sequences, computes pairwise distances
#' between them, and orders the sequences based on their similarity using
#' hierarchical clustering.
#'
#' @param dna_strings A character vector containing DNA sequences (A, C, G, T).
#'
#' @return A character vector of DNA sequences ordered by similarity.
#'
#' @examples
#' dna_strings <- c("ACTG", "ACCG", "ATGC", "AGTC")
#' ordered_dna <- order_by_similarity(dna_strings)
#' print(ordered_dna)
#'
#' @importFrom Biostrings DNAStringSet stringDist
#' @export
order_by_similarity <- function(dna_strings) {
  dna_set <- DNAStringSet(dna_strings)
  dist_matrix <- stringDist(dna_set)
  hc <- hclust(as.dist(dist_matrix))
  ordered_indices <- hc$order

  # Return the DNA strings in the new order
  ordered_dna_strings <- dna_strings[ordered_indices]

  return(ordered_dna_strings)
}


#### TODO Functions ####

#' Split names by delim character
#'
#' @param df a data.frame
#' @param delim character to split names by
#' @param split_item index of split string to keep
#' @param remove_space should whitespace be replaced with "_" default: TRUE
#'
#' @return vector of names
#' @keywords internal
#'
#' @importFrom limma strsplit2
#'
#' @examples
.split_fix_names <- function(df, delim=':', split_item=1, remove_space=TRUE){
  nms <- strsplit2(names(df), delim)[,split_item]
  if(remove_space){
    nms <- gsub(' ', '_',nms)
  }
  return(nms)
}



#' Title
#'
#' @param n_genes
#' @param m_samples
#' @param seed
#' @param groups
#'
#' @return
#' @keywords internal
#'
#' @import DESeq2
#'
#' @examples
.create_test_dds <- function(n_genes=500, m_samples=60, seed=42, groups=NULL){
  set.seed(seed)
  dds <- makeExampleDESeqDataSet(n=n_genes, m = m_samples)
  mcols(dds) <- NULL
  mcols(dds) <- data.frame(gene_id=names(dds),
                           row.names=names(dds))
  if(is.null(groups)){
    groups <-
      c(rep('A', m_samples/3),
        rep('B', m_samples/3),
        rep('C', m_samples/3)) %>%
      factor(.)
  }
  colData(dds) <- DataFrame(condition=groups,
                            row.names = colnames(dds))

  return(dds)

}

#' Turn summary object from lm into a data.frame
#'
#' @param lmlist lm summary object
#' @param coef the coefficient to extract
#'
#' @return data.frame with coefficient
#' @keywords internal
#'
#' @examples
.lmlist_to_df <- function(lmlist, coef=2){

  rnms <- try(limma::strsplit2(names(lmlist), ' ')[,2], silent = T) # gets rid of extranous name
  if(class(rnms) == 'try-error'){
    rnms <- names(lmlist)
  }
  cnms <- c("Estimate","Std_Error","z_or_t_value", "p_value" )

  coef_list <- lapply(lmlist, function(x){
    x <- x$coefficients[coef,]
  })
  out_df <- data.frame(do.call(rbind, coef_list))
  row.names(out_df) <- rnms
  colnames(out_df) <- cnms
  return(out_df)

}

#' Gets lfc from a pathway res df
#'
#' This returns the lfc from a genelist as a
#' named vector of genes by lfc
#'
#' @param gene_list list of gene vectors
#' @param res_df long res df table
#'
#' @return
#' @keywords internal
#'
#' @examples
.get_lfc <- function(gene_list, res_df){
  ol <- lapply(gene_list, function(x){
    lfc <- res_df$logFC[match(x, res_df$external_gene_name)] %>%
      setNames(., x)
  })
  names(ol) <- names(gene_list)
  return(ol)
}


#' Split a named vector of DEGs by fold change
#'
#' @param lfc_list DEG vector of FCs
#'
#' @return
#' @keywords internal
#'
#' @examples
.split_by_reg <- function(lfc_list){
  split_degs <- list(
    sort(lfc_list[lfc_list > 0, drop=F], decreasing = TRUE),
    sort(lfc_list[lfc_list < 0, drop=F]))
  names(split_degs) <- c('upreg', 'downreg')
  return(split_degs)
}








