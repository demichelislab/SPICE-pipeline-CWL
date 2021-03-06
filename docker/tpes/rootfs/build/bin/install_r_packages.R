#!/usr/bin/env Rscript

packages   <- c(
  'docopt',
  'TPES'
)
repository <- 'https://cloud.r-project.org/'
install.packages(packages, repos = repository)
rr         <- as.integer(!all(sapply(packages, function (x) require(x, character.only = TRUE))))
quit(save = 'no', status = rr)
