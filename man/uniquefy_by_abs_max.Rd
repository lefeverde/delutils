% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc_utils.R
\name{uniquefy_by_abs_max}
\alias{uniquefy_by_abs_max}
\title{Return a unique data.frame by max value}
\usage{
uniquefy_by_abs_max(dat, group_col = NULL, val_col = NULL)
}
\arguments{
\item{dat}{data.frame with 1st column as string and second as numeric}

\item{group_col}{string denoting group column (col with duplicates)}

\item{val_col}{string denoting value column}
}
\value{
unique data.frame
}
\description{
Return a unique data.frame by max value
}
\examples{
set.seed(42)
dat <- data.frame(group=c('a', 'a', 'a', 'b', 'c', 'c'), vals=runif(6))
dat$vals[6] <- -4
}
