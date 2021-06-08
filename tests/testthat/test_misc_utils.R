
library(delutils)
context("misc_utils")

create_example_dataframe <- function(){
  set.seed(42)
  mat <- data.frame(group=c('a', 'a', 'a', 'b', 'c', 'c'), vals=runif(6) * 100)
  mat$vals[6] <- -400
  mat$vals <- round(mat$vals)
  row.names(mat) <- paste0('r_', seq(1,6))
  return(mat)
}


#### clean_by_levels ####

test_that('clean_by_levels correctly checks input for correct type', {
  expect_error(clean_by_levels(dat=matrix()), 'dat needs to be a data.frame not a matrix')
  expect_error(clean_by_levels(dat=list()), 'dat needs to be a data.frame not a list')
  expect_error(clean_by_levels(dat=data.frame(), n_levels='foo'), 'n_levels needs to be an integer >= 1')
  expect_error(clean_by_levels(dat=data.frame(), n_levels=-1), 'n_levels needs to be an integer >= 1')
})

test_that('clean_by_levels cleans a data.frame', {
  in_dat <- data.frame(A=rep(1,4), B=seq(1,4), C=c('f1', 'f2', 'f1', 'f1'))
  ex_dat <- data.frame(B=seq(1,4), C=c('f1', 'f2', 'f1', 'f1'))
  expect_equal(clean_by_levels(dat=in_dat), ex_dat)

})


#### rbind_named_df_list ####
test_that('rbind_named_df_list returns a bound data.frame',{
  el <- data.frame(
    col_from_list=paste0('Cond',
                         c(rep(2, 5), rep(3,2))),A=seq(1,7))

  l <- list(data.frame(A=seq(1,5)),
            data.frame(A=seq(6, 7))) %>%
    setNames(., c('Cond2', 'Cond3'))
  expect_equal(rbind_named_df_list(l), el)
  # Making sure it works with tibbles
  l2 <- lapply(l, tibble::as_tibble)
  expect_equal(rbind_named_df_list(l2), el)
})

test_that('rbind_named_df_list checks for non data.frames', {
  l <- list(matrix(data=NA, ncol=1, nrow=5),
            data.frame(A=seq(1,5)),
            data.frame(A=seq(6, 7))) %>%
    setNames(., c('Cond1', 'Cond2', 'Cond3'))
  expect_error(rbind_named_df_list(l),"df_list contains items which are not data.frame or tibble" )
})



test_that('rbind_named_df_list can (gracefully) handle empty data.frames or NULL objects',{
  el <- data.frame(col_from_list=paste0('Cond', c(rep(2, 5), rep(3,2))),A=seq(1,7))
  l_empty <- list(data.frame(),
                  data.frame(A=seq(1,5)),
                  data.frame(A=seq(6, 7))) %>%
    setNames(., c('Cond1', 'Cond2', 'Cond3'))
  l_null <- list(NULL,
                 data.frame(A=seq(1,5)),
                 data.frame(A=seq(6, 7))) %>%
    setNames(., c('Cond1', 'Cond2', 'Cond3'))
  expect_equal(expect_warning(rbind_named_df_list(l_empty), 'Cond1'), el)
  expect_equal(expect_warning(rbind_named_df_list(l_null), 'Cond1'), el)
})

test_that('rbind_named_df_list complains loudly about data.frames with mismatched columns',{
  l <- list(data.frame(A=seq(1,5), B=seq(8,12)),
            data.frame(A=seq(6, 7))) %>%
    setNames(., c('Cond2', 'Cond3'))
  expect_error(rbind_named_df_list(l))
})

#### uniquefy ####

test_that('uniquefy (base) function returns a logical vector of uniquefied rows by largest value data.frame',{
  mat <-  create_example_dataframe()
  e_vec <- c(FALSE,TRUE,FALSE,TRUE,TRUE,FALSE)
  out_vec <- uniquefy_base(mat, 'group', 'vals')
  expect_equal(out_vec, e_vec)

})

test_that('uniquefy by abs max returns a uniquefied (by abs(max)) data.frame',{
  mat <-  create_example_dataframe()
  e_mat <- data.frame(group=c('a', 'b', 'c'),
                      vals=c(94,83,-400),
                      row.names = c('r_2','r_4','r_6'))
  out_mat <- uniquefy_by_abs_max(mat)

  expect_equal(row.names(out_mat), row.names(e_mat))
  expect_equal(out_mat$group, e_mat$group)
  expect_equal(out_mat$vals, e_mat$vals)
})

test_that('uniquefy by variance returns a uniquefied (by variance) data.frame of expression data',{


  map <- data.frame(row.names=paste0('gene', seq(1,6)),
                    duplicated_ids=c('a','a','b','b','b','c'))
  gene_vars <- c(2, 1, .5, 1, 4, 1)
  set.seed(42)
  mat <- lapply(gene_vars, function(x){
    rnorm(10, 0, x)
  }) %>% do.call(rbind, .)
  row.names(mat) <- row.names(map)
  e_vec <-  c(TRUE,FALSE,FALSE,FALSE,TRUE,TRUE)
  e_mat <- mat[e_vec,]
  out_mat <- uniquefy_by_variance(mat, map, 'duplicated_ids')
  expect_equal(out_mat, e_mat)
})
