#!/usr/bin/env Rscript

bioc_packages <- c('DNAcopy', 'cghFLasso')
repository    <- 'https://cloud.r-project.org/'
install_path  <- '/installed/lib/R/library'
dir.create(install_path, recursive = TRUE)
source('http://bioconductor.org/biocLite.R')
.libPaths(install_path)
biocLite(bioc_packages)
slmseg_archive    <- list.files('tmp', pattern = '^SLMSeg_[0-9\\.]+\\.tar\\.gz$',
                                full.names = TRUE, recursive = TRUE)
n_slmseg_archives <- length(slmseg_archive)
if (n_slmseg_archives != 1) {
  if (n_slmseg_archives == 0) {
    cat('No SLMSeg archive found\n')
  } else {
    cat('Found more than one SLMSeg archive:\n', paste0(slmseg_archive, collapse = '\n'))
  }
  stop('no', 1)
}
install.packages(slmseg_archive, type = 'source')
to_check_packages <-  c(bioc_packages, 'SLMSeg')
rr <- as.integer(!all(sapply(to_check_packages, function (x) require(x, character.only = TRUE))))
quit(save = 'no', status = rr)
