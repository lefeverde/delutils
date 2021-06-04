context("misc_utils")
library(delutils)

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
