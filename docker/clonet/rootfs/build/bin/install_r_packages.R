#!/usr/bin/env Rscript

# To be removed when ipgrah is updated to a version > 1.1.2
igraph_deps  <- c('irlba', 'magrittr', 'Matrix', 'pkgconfig')
packages     <- c('docopt', 'dgof', 'pso', 'reshape2', 'sets', igraph_deps)
repository   <- 'https://cloud.r-project.org/'
download_dir <- file.path(getwd(), 'tmp')
dir.create(download_dir)
install.packages(packages, repos = repository)
download.packages('igraph', download_dir, repos = repository)
rr         <- as.integer(!all(sapply(packages, function (x) require(x, character.only = TRUE))))
quit(save = 'no', status = rr)
