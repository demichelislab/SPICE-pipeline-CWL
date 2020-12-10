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

check_if_positive    <- function (x) !is.na(x) && x > 0
check_if_zero_to_one <- function (x) !is.na(x) && x >= 0 && x <= 1
conv_to_numeric      <- function (x) suppressWarnings(as.numeric(x))
conv_to_integer      <- function (x) suppressWarnings(as.integer(x))

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
  tpes [options] SAMPLE_NAME SEGFILE SNV_PILEUP_FILE PLOIDY OUTPUT_FOLDER

-h --help        Shows help informations
SAMPLE_NAME      The name of the sample whose purity have to be inferred.
SEGFILE          The file containing copy number segments for the sample.
SNV_PILEUP_FILE  The file containing the coverage of the SNVs.
PLOIDY           The ploidy of the sample.
OUTPUT_FOLDER    The folder where the output will be created.

options:
  --report_filename NAME                 The name of the file where to write the pdf report.
                                            If specified the report is generated otherwise only
                                            the text output is produced.
  --output_filename NAME                 The name of the file where to write the result of the
                                            analysis. [default: tpes.tsv]
  --reference_mapping_bias RMB           The allelic fraction corresponding to the
                                           reference mapping bias (see TPES
                                           manual). [default: 0.47].
  --max_allelic_fraction MAX_AF          The maximum allelic fraction to consider when
                                           filtering SNVs (see TPES manual).
                                           [default: 0.55]
  --min_coverage MIN_COV                 The minimum coverage to consider when
                                           filtering SNVs (see TPES manual).
                                           [default: 10]
  --min_alternative_reads MIN_ALT_READS  The minimum number of alternative reads
                                           to consider when filtering SNVs (see
                                           TPES manual). [default: 5]
  --min_snvs MIN_SNVs                    The minimum number of SNVs that a
                                           sample needs to have to be able to
                                           infer purity (see TPES manual).
                                           [default: 10]
  --log_file PATH                        The path to the output file
'

run_tpes <- function (
    sample_name,
    segfile,
    snv_read_counts,
    ploidy,
    reference_mapping_bias = 0.47,
    max_allelic_fraction   = 0.55,
    min_coverage           = 10,
    min_alternative_reads  = 5,
    min_snvs               = 10,
    run_report             = FALSE
  ) {
  tpes_function <- if (run_report) TPES_report else TPES_purity
  ploidy_table  <- data.frame(
    sample = sample_name,
    ploidy = ploidy
  )
  tpes_function(
    sample_name,
    segfile,
    snv_read_counts,
    ploidy_table,
    RMB         = reference_mapping_bias,
    maxAF       = max_allelic_fraction,
    minCov      = min_coverage,
    minAltReads = min_alternative_reads,
    minSNVs     = min_snvs
  )
}

main_function <- function () {
  if (!is.null(cli_args$log_file)) {
    sink(cli_args$log_file)
  }
  greater_than_0         <- 'The value must be greater than 0'
  zero_to_one            <- 'The value must be between 0 and 1 included'
  ploidy                 <- parse_arg_pred(cli_args, 'PLOIDY', function (x) is.na(x) || check_if_positive(x),
                              greater_than_0, converter = conv_to_numeric)
  reference_mapping_bias <- parse_arg_pred(cli_args, 'reference_mapping_bias', check_if_zero_to_one,
                              zero_to_one, converter = conv_to_numeric)
  max_allelic_fraction   <- parse_arg_pred(cli_args, 'max_allelic_fraction', check_if_zero_to_one,
                              zero_to_one, converter = conv_to_numeric)
  min_coverage           <- parse_arg_pred(cli_args, 'min_coverage', check_if_positive,
                              greater_than_0, converter = conv_to_numeric)
  min_alternative_reads  <- parse_arg_pred(cli_args, 'min_alternative_reads', check_if_positive,
                              greater_than_0, converter = conv_to_integer)
  min_snvs               <- parse_arg_pred(cli_args, 'min_snvs', check_if_positive,
                              greater_than_0, converter = conv_to_integer)
  output_folder_info     <- file.info(cli_args$OUTPUT_FOLDER)

  if (is.na(output_folder_info$isdir)) {
    if (!dir.create(cli_args$OUTPUT_FOLDER, recursive = TRUE)) {
      die(sprintf('Error while trying to create "%s" folder', cli_args$OUTPUT_FOLDER))
    }
  } else {
    if (!output_folder_info$isdir) {
      die(sprintf('The path specified as output folder (%s) already exists and is a file',
                  cli_args$OUTPUT_FOLDER))
    }
  }

  segfile    <- read.delim(cli_args$SEGFILE)
  snv_pileup <- read.delim(cli_args$SNV_PILEUP_FILE)

  setwd(cli_args$OUTPUT_FOLDER)

  res <- run_tpes(
    cli_args$SAMPLE_NAME,
    segfile,
    snv_pileup,
    ploidy,
    reference_mapping_bias = reference_mapping_bias,
    max_allelic_fraction   = max_allelic_fraction,
    min_coverage           = min_coverage,
    min_alternative_reads  = min_alternative_reads,
    min_snvs               = min_snvs
  )

  if (!is.null(cli_args$report_filename)) {
    pdf(cli_args$report_filename)
    run_tpes(
      cli_args$SAMPLE_NAME,
      segfile,
      snv_pileup,
      ploidy,
      reference_mapping_bias = reference_mapping_bias,
      max_allelic_fraction   = max_allelic_fraction,
      min_coverage           = min_coverage,
      min_alternative_reads  = min_alternative_reads,
      min_snvs               = min_snvs,
      run_report             = TRUE
    )
    dev.off()
  }

  write.table(
    res,
    cli_args$output_filename,
    quote     = FALSE,
    sep       = '\t',
    na        = '',
    row.names = FALSE
  )

}

(make_main(docstring, c('TPES'), main_function))()
