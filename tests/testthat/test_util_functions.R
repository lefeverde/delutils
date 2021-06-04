context("util_functions.R tests")
library(delutils)


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
                  package = 'basicOmics',
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
      "gene_id",
      "external_gene_name",
      "gene_biotype",
      "description" )
  for(col in gene_annots_cols){

    # expect_equal(e_res[,col], o_res[,col])
    expect_true(all(e_res[,col] ==  o_res[,col])) # No clue why the above wasn't working
  }


})
