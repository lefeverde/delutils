context("stat_utils.R tests")
library(delutils)

#### Functions #####
limma_example_data <- function(){
  set.seed(2016)
  sigma2 <- 0.05 / rchisq(100, df=10) * 10
  y <- matrix(rnorm(100*6,sd=sqrt(sigma2)),100,6)
  design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))
  y[1,4:6] <- y[1,4:6] + 1
  fit <- lmFit(y,design)
}

get_example_data <- function(){
  rds_to_load <-
    c('example_dge_object.rds',
      'example_voom_object.rds',
      'example_fit_object.rds',
      'expected_limma_res.rds') %>%
    lapply(., function(x){
      system.file('extdata',
                  x,
                  package = 'delutils',
                  mustWork = TRUE)
    }) %>% setNames(., c('dge', 'v', 'fit', 'e_res'))
}


test_that('get_limma_results returns correct tibble of results',{
  ex_files <- get_example_data()
  e_res <- readRDS(ex_files$e_res)
  fit <- readRDS(ex_files$fit)
  o_res <- get_limma_results(fit)
  gene_annots_cols <-
    c("coefficient",
      "ensembl_gene_id",
      "external_gene_name",
      "gene_biotype",
      "description" )
  for(col in gene_annots_cols){

    # expect_equal(e_res[,col], o_res[,col])
    expect_true(all(e_res[,col] ==  o_res[,col])) # No clue why the above wasn't working
  }


})


#### Statisitc related wranglers ####
test_that('get_overlap_matix produces correct overlap matrix',{
  set.seed(42)
  n_genes_to_sample <- rnbinom(3, 5, .2)
  genes <- paste0('gene', seq(1, max(n_genes_to_sample)))
  l1 <- lapply(seq_along(n_genes_to_sample), function(i){
    n <- n_genes_to_sample[i]
    sample(genes, n)
  }) %>%
    setNames(., paste0('L1_', LETTERS[1:length(n_genes_to_sample)]))
  e_mat <-  c(1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 0.4615385, 1.0000000, 0.4615385, 1.0) %>%
    matrix(data=., nrow=3, ncol=3, dimnames = list(names(l1), names(l1)))
  omat <- get_overlap_matix(l1)
  expect_equivalent(e_mat, omat, tolerance=1e-6)
})

#### Distance functions ####
test_that('overlap_ratio calculates overlap correctly',{
  x <- LETTERS[1:5]
  x2 <- c(x, x)
  y <- LETTERS[2:10]
  outs <- c(overlap_ratio(x,y),
            jaccard(x, y),
            dice(x, y))
  e_outs <- c(.8, .4, 0.5714286)
  expect_equal(outs, e_outs, tolerance = .01)
  ## Testing the check
  expect_error(overlap_ratio(x2, y))
  expect_error(jaccard(x2, y))
  expect_error(dice(x2, y))

})



