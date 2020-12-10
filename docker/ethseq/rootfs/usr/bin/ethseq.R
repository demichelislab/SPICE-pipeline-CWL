#!/usr/bin/env Rscript

`%+%` <- paste0
`%_%` <- paste
nl    <- function (ll, n = 1) ll %+% Reduce(paste0, rep('\n', n))

die <- function (msg, status = 1) {
  cat(msg)
  quit(save = 'no', status = status)
}

silent_lib_loader <- function (...) {
  libs <- c(...)
  for (l in libs) {
    suppressMessages(library(l, quietly = TRUE, character.only = TRUE,
                             warn.conflicts = FALSE))
  }
}

parse_arg <- function (arguments, arg_name, fun) {
  argument <- getElement(arguments, arg_name)
  tryCatch(fun(argument), error = function (e) {
    die(sprintf('Error while trying to validate "%s" argument:\n\t%s\n',
                arg_name, e$message))
  })
}

parse_arg_pred <- function (arguments, arg_name, pred, error, converter = identity) {
  parse_arg(arguments, arg_name, function (x) {
              res <- converter(x)
              if (pred(res)) {
                res
              } else {
                msg <- error
                if (is.function(error)) {
                  msg <- error(res)
                }
                stop(msg)
              }
  })
}

check_arg_in_list <- function (...) {
  function (aa) !anyNA(pmatch(aa, c(...), duplicates.ok = TRUE))
}

check_if_positive <- function (x) !is.na(x) && x > 0
conv_to_numeric   <- function (x) suppressWarnings(as.numeric(x))
conv_to_integer   <- function (x) suppressWarnings(as.integer(x))

make_main <- function (docstring, libs = c(), main_function = function () {}) {
  main_template <- function (...) {
    do.call(silent_lib_loader, as.list(c('docopt', libs)))
    main_env <- environment(main_function)
    main_env$cli_args <- docopt::docopt(docstring)
    main_function(...)
  }

  function (...) {
    tryCatch(main_template(...), error = function (e) {
      die(e$message)
    })
  }
}

docstring <- '
Usage:
  ethseq [options] MODEL_FILE VCF_INPUT OUTPUT_FOLDER

-h --help          Shows help informations
MODEL_FILE         The file containing the model to use for ethnicity annotation.
VCF_INPUT          The file containing the genotypes of the samples to analyze.
OUTPUT_FOLDER      The folder where the output will be created.

options:
  --num_threads N  Specifies the name for the sample to be used in the
                     output [default: 1]
  --enable_plot    If specified a graphical report is produced to show the
                     position of the sample in the PCA space of the genotypes.
  --log_file PATH  The path to the output file
'

run_ethseq <- function (vcf_dat, output_path, model_path, enable_plot = FALSE, threads = 1) {
  res <- ethseq.Analysis(target.vcf = vcf_dat,
                         out.dir = output_path,
                         model.gds = model_path,
                         space = '3D',
                         composite.model.call.rate = 0.99,
                         run.genotype = FALSE,
                         cores = threads)
  files_to_delete <- paste0(c('Aggregated', 'Target'), '.gds')
  if (!enable_plot) {
    files_to_delete <- c(files_to_delete, 'Report.pdf')
  }
  for (ff in file.path(output_path, files_to_delete)) {
    if (file.exists(ff)) {
      file.remove(ff)
    }
  }
  res
}

main_function <- function () {
  if (!is.null(cli_args$log_file)) {
    sink(cli_args$log_file)
  }
  greater_than_0    <- 'The value must be greater than 0'
  num_threads       <- parse_arg_pred(cli_args, 'num_threads', check_if_positive,
                                      greater_than_0, converter = conv_to_integer)

  res <- run_ethseq(
    cli_args$VCF_INPUT,
    cli_args$OUTPUT_FOLDER,
    cli_args$MODEL_FILE,
    enable_plot = cli_args$enable_plot,
    num_threads
  )
  quit('no', !res)
}

(make_main(docstring, c('EthSEQ'), main_function))()
