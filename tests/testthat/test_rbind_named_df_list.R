context("rbind_named_df_list")
library(delutils)



options(stringsAsFactors = FALSE)

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
