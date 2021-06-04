context('uniquefy scripts')
library(delutils)

create_example_dataframe <- function(){
  set.seed(42)
  mat <- data.frame(group=c('a', 'a', 'a', 'b', 'c', 'c'), vals=runif(6) * 100)
  mat$vals[6] <- -400
  mat$vals <- round(mat$vals)
  row.names(mat) <- paste0('r_', seq(1,6))
  return(mat)
}

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
