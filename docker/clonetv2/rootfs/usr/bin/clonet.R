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
  clonet [options] SEGFILE NORMAL_SNPS_PILEUP TUMOR_SNPS_PILEUP OUTPUT_FOLDER

-h --help           Shows help informations
SEGFILE             The file containing copy number segments for the sample.
NORMAL_SNPS_PILEUP  The file containing the pileup on the SNPs of the normal
                      sample
TUMOR_SNPS_PILEUP   The file containing the pileup on the SNPs of the normal
                      sample
OUTPUT_FOLDER       The folder where the output will be created.

options:
  --num_threads THREADS                          The number of threads to use [default: 1]
  --plot_filename NAME                          The name of the file where to write the pdf plot
                                                   If specified the plot is generated otherwise only
                                                   the text output is produced.
  --snv_table PATH                               The name of the file that contains the SNVs.
                                                   If this parameter is passed the program computes
                                                   the clonality of the SNVs otherwise it is not
                                                   computed.
  --beta_table_filename FILENAME                 The name of the output file containing the beta
                                                   table. [default: beta_table.tsv]
  --ploidy_table_filename FILENAME               The name of the output file containing the ploidy
                                                   table. [default: ploidy_table.tsv]
  --admixture_table_filename FILENAME            The name of the output file containing the
                                                   admixture table. [default: admixture_table.tsv]
  --scna_clonality_table_filename FILENAME       The name of the output file containing the
                                                   CN clonality table. [default: scna_clonality_table.tsv]
  --allele_specific_cna_table_filename FILENAME  The name of the output file containing the allele specific
                                                  copy number table. [default: allele_specific_scna_table.tsv]
  --snv_clonality_table_filename FILENAME        The name of the output file containing the SNV clonality
                                                   table. [default: snv_clonality_table.tsv]
  --log_file PATH                                The path to the output file.
'

read_vep_table <- function (filename) {
  vep_lines       <- readLines(filename)
  is_heading_line <- grepl('^#', vep_lines)
  dl              <- vep_lines[!is_heading_line]
  lr              <- c(sub('^#', '', tail(vep_lines[is_heading_line], 1)), dl)
  list(vep_data      = fread(text = lr, na.strings = c('-', '.')),
       raw_lines     = vep_lines,
       heading_lines = which(is_heading_line))
}

get_snv_containing_lines <- function (vep_data) {
  vep_keys          <- with(vep_data, paste(Location, Allele, sep = ':'))
  vep_only_snv_keys <- vep_keys[vep_data$VARIANT_CLASS == 'SNV']
  match(vep_keys, vep_only_snv_keys)
}

main_function <- function () {
  if (!is.null(cli_args$log_file)) {
    sink(cli_args$log_file)
  }
  greater_than_0     <- 'The value must be greater than 0'
  zero_to_one        <- 'The value must be between 0 and 1 included'
  num_threads        <- parse_arg_pred(cli_args, 'num_threads', check_if_positive,
                              greater_than_0, converter = conv_to_integer)
  output_folder_info <- file.info(cli_args$OUTPUT_FOLDER)

  if (is.na(output_folder_info$isdir)) {
    if (!dir.create(cli_args$OUTPUT_FOLDER, recursive = TRUE)) {
      die(sprintf('Error while trying to create "%s" folder', cli_args$OUTPUT_FOLDER))
    }
  } else {
    if (!output_folder_info$isdir) {
      die(sprintf('The path specified as output folder (%s) already exists but is a file',
                  cli_args$OUTPUT_FOLDER))
    }
  }


  segfile            <- fread(cli_args$SEGFILE)
  normal_snps_pileup <- fread(cli_args$NORMAL_SNPS_PILEUP)
  tumor_snps_pileup  <- fread(cli_args$TUMOR_SNPS_PILEUP)
  snv_data           <- if (!is.null(cli_args$snv_table)) { read_vep_table(cli_args$snv_table) }
  snv_table          <- snv_data$vep_data

  setwd(cli_args$OUTPUT_FOLDER)

  beta_table                <- compute_beta_table(
    segfile,
    tumor_snps_pileup,
    normal_snps_pileup,
    n_cores = num_threads
  )
  ploidy_table              <- compute_ploidy(beta_table)
  admixture_table           <- compute_dna_admixture(
    beta_table   = beta_table,
    ploidy_table = ploidy_table,
    n_cores      = num_threads
  )
  scna_clonality_table      <- compute_scna_clonality_table(
    beta_table      = beta_table,
    ploidy_table    = ploidy_table,
    admixture_table = admixture_table,
    n_cores         = num_threads
  )
  allele_specific_cna_table <- compute_allele_specific_scna_table(
    beta_table      = beta_table,
    ploidy_table    = ploidy_table,
    admixture_table = admixture_table,
    n_cores         = num_threads
  )

  if (!is.null(snv_table)) {
    output_file <- cli_args$snv_clonality_table_filename
    snv_clonality_table <- compute_snv_clonality(
      sample_id       = unique(beta_table$sample),
      snv_read_count  = snv_table,
      beta_table      = beta_table,
      ploidy_table    = ploidy_table,
      admixture_table = admixture_table,
      n_cores         = num_threads
    )
    snv_clonality_cols <- setdiff(colnames(snv_clonality_table), colnames(snv_table))
    writeLines(
      with(
        snv_data,
        c(raw_lines[head(heading_lines, -1)],
          paste0(c(raw_lines[tail(heading_lines, 1)], snv_clonality_cols), collapse = '\t'))
      ),
      output_file
    )
    output_data <- data.table(
      snv_table,
      snv_clonality_table[
        get_snv_containing_lines(snv_data$vep_data),
        snv_clonality_cols
      ]
    )
    fwrite(
      output_data,
      output_file,
      append = TRUE,
      sep = '\t',
      quote = FALSE,
      na = '-'
    )
  }

  if (!is.null(cli_args$plot_filename)) {
    pdf(cli_args$plot_filename)
    print(
      check_ploidy_and_admixture(
        beta_table      = beta_table,
        ploidy_table    = ploidy_table,
        admixture_table = admixture_table
      )
    )
    dev.off()
  }

  table_names <- c(
    'beta_table',
    'ploidy_table',
    'admixture_table',
    'scna_clonality_table',
    'allele_specific_cna_table'
  )

  for (tn in table_names) {
    if (exists(tn)) {
      fwrite(get(tn), cli_args[[tn %+% '_filename']], sep = '\t')
    }
  }

}

(make_main(docstring, c('data.table', 'CLONETv2'), main_function))()
