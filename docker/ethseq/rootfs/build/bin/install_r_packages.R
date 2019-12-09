#!/usr/bin/env Rscript

packages      <- c('docopt', 'MASS', 'geometry', 'data.table', 'plot3D')
bioc_packages <- c('SNPRelate', 'gdsfmt')
repository    <- 'https://cloud.r-project.org/'
install.packages(c(packages, 'BiocManager'), repos = repository)
BiocManager::install(bioc_packages)
install.packages('EthSEQ')
to_check_packages <-  c(packages, bioc_packages, 'EthSEQ')
rr <- as.integer(!all(sapply(to_check_packages, function (x) require(x, character.only = TRUE))))
quit(save = 'no', status = rr)
