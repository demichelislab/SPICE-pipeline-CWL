#!/usr/bin/env Rscript

`%+%` <- paste0
`%_%` <- paste
nl    <- function (ll, n = 1) ll %+% Reduce(paste0, rep('\n', n))

add_zero_to_single_digit_chromosome <- function (x) {
  sub('^(chr)?(\\d)$', '\\10\\2', x)
}

get_main_chrs <- function(chr_names, strict = FALSE){
  postfix <- if (strict) '' else '(_.*)?'
  grepl(paste0('^(chr)?([01]?[1-9]|10|2[0-4]|[XY]|MT?|Un?)', postfix, '$'), chr_names)
}

chr_order <- function (chromosomes) {
  is_main_chr <- get_main_chrs(chromosomes)
  order(
    !is_main_chr,
    local({
      x <- sub('^chr([^_]+).*$',  '\\1',  chromosomes[is_main_chr])
      x <- sub('^(\\d)$',        '0\\1',                         x)
      x <- sub('^M$',               'Z',                         x)
      c(sub('^Un.*$',         'ZZ',   x), chromosomes[!is_main_chr])
    }),
    c(sub('^[^_]+_?', '', chromosomes[is_main_chr]),
      rep('', sum(!is_main_chr)))
  )
}

order_rank <- function (x, order_fun) {
  x_unq <- unique(x)
  match(x, x_unq[order_fun(x_unq)])
}

chr_rank  <- function (x) order_rank(x, chr_order)

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

check_if_between    <- function (x, min_value = -Inf, max_value = Inf) !is.na(x) && x >= min_value && x <= max_value
check_if_positive    <- function (x) !is.na(x) && x > 0
check_if_zero_to_one <- function (x) check_if_between(x, 0, 1)
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
  pipeline_utils pileup_to_genotype SNPS_PILEUP OUTPUT_FOLDER [--sample_name SAMPLE_NAME --snps_list_file SNPS_LIST --min_depth_of_coverage DEPTH --heterozygous_perc PERC --output_filename FILENAME]
  pipeline_utils pileup_to_coverage --sample_name SAMPLE_NAME SNVS_PILEUP OUTPUT_FOLDER [--output_filename FILENAME]
  pipeline_utils vcf_snvs_to_bed VCF_FILE OUTPUT_FOLDER [--output_filename FILENAME]
  pipeline_utils mutect_only_snvs VCF_FILE OUTPUT_FOLDER [--output_filename FILENAME]
  pipeline_utils merge_vep_snv_coverage VEP_ANNOTATED_VARIANTS SNV_COVERAGE_NORMAL_FILE SNV_COVERAGE_TUMOR_FILE OUTPUT_FOLDER [--output_filename FILENAME]

-h --help                      Shows help informations
SNPS_PILEUP                    The file containing the pileup of SNPs
SNVS_PILEUP                    The file containing the pileup of SNVs
VCF_FILE                       The file containing SNVs and indels
VEP_ANNOTATED_VARIANTS         The variants annotated by VEP in tabular format.
SNV_COVERAGE_NORMAL_FILE       The file containing the coverage of SNVs positions
                                 in the normal sample
SNV_COVERAGE_TUMOR_FILE        The file containing the coverage of SNVs positions
                                 in the tumor sample
OUTPUT_FOLDER                  The folder where the output will be created.
--sample_name SAMPLE_NAME      The name of the analyzed sample that is written
                                 to the output file.
--snps_list_file FILE          The SNPs whose genotype needs to be written to the output
                                 file.
--min_depth_of_coverage DEPTH  The minimum depth of coverage of the SNPs to be considered
                                 reliable (SNPs with less coverage than this will have
                                 unknown genotype in the output) [default: 10]
--heterozygous_perc PERC       The threshold that will be used to define the genotype.
                                 A SNP with allelic fraction (AF) less than heterozygous_perc
                                 will have "0/0" genoytpe. A SNP with AF between heterozygous_perc
                                 and 1 - heterozygous_perc will have "1/0" genoytpe. Finally
                                 a SNP with AF greater or equal to 1 - heterozygous_perc
                                 will have "1/1" genotype. [default: 0.2]
--output_filename FILE         The name of the file where the output will be produced.
                                 The default value depends on the tool that is called.
'

create_folder_if_not_exists <- function (folder) {
  output_folder_info <- file.info(folder)
  if (is.na(output_folder_info$isdir)) {
    if (!dir.create(folder, recursive = TRUE)) {
      die(sprintf('Error while trying to create "%s" folder', folder))
    }
  } else {
    if (!output_folder_info$isdir) {
      die(sprintf('The path specified as output folder (%s) already exists but is a file',
                  folder))
    }
  }
}

read_vcf_data <- function (filename) {
  vcf_lines       <- readLines(filename)
  is_heading_line <- grepl('^#', vcf_lines)
  dl              <- vcf_lines[!is_heading_line]
  lr              <- c(sub('^#', '', tail(vcf_lines[is_heading_line], 1)), dl)
  list(vcf_data      = fread(text = lr, na.strings = c('-', '.')),
       raw_lines     = vcf_lines,
       heading_lines = which(is_heading_line))
}

vcf_is_snv <- function (variants) {
  nucleotides_dna <- c('a', 'c', 'g', 't')
  with(
    variants,
    tolower(REF) %in% nucleotides_dna & tolower(ALT) %in% nucleotides_dna
  )
}

vep_is_snv <- function (variants) {
  variants$VARIANT_CLASS == 'SNV'
}

snps_file_match_fields  <- c('chr', 'pos', 'rsid', 'ref', 'alt')

pileup_to_genotype <- function (
  snps_pileup_path,
  output_folder,
  output_filename       = 'genotype.vcf',
  one_alternative       = TRUE,
  sample_name           = NULL,
  snps_list_file        = NULL,
  min_depth_of_coverage = 10,
  heterozygous_perc     = 0.2
) {
  sample_name        <- if (is.null(sample_name)) gsub('^([^.]*).*$', '\\1', basename(snps_pileup_path)) else sample_name
  output_file        <- paste0(output_folder, '/', output_filename)
  header             <- '##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' %+% sample_name
  genotype_intervals <- c(0, heterozygous_perc,  1 - heterozygous_perc, 1)
  genotype_map       <- c('0/0', '0/1', '1/1')
  unknown_genotype   <- './.'
  empty_fields       <- c('qual', 'filter', 'info', 'format')
  has_list_file      <- !is.null(snps_list_file)
  snps_list          <- NULL
  create_folder_if_not_exists(output_folder)
  snps_pileup        <- fread(
    snps_pileup_path,
    select = c(snps_file_match_fields, 'af', 'cov')
  )[
    cov >= min_depth_of_coverage
  ][,
    genotype := fcoalesce(
      genotype_map[findInterval(af, genotype_intervals, rightmost.closed = TRUE)],
      unknown_genotype
    )
  ][
    , c(snps_file_match_fields, 'genotype'), with = FALSE
  ]
  if (one_alternative) {
    snps_pileup$alt <- substr(snps_pileup$alt, 1, 1)
  }
  setkeyv(snps_pileup, snps_file_match_fields)
  output            <- snps_pileup
  if (has_list_file) {
    snps_list <- local({
      lll       <- readLines(snps_list_file)
      fread(
        text = grep('^#', lll, value = TRUE, invert = TRUE),
        col.names = c(snps_file_match_fields, 'a', 'b', 'c')
      )[
        , snps_file_match_fields, with = FALSE
      ]
    })
    setkeyv(snps_list, snps_file_match_fields)
    matching_rows <- merge(
      snps_list,
      snps_pileup,
      by              = snps_file_match_fields,
      all.x           = TRUE,
      allow.cartesian = TRUE
    )[
      , genotype := fcoalesce(genotype, unknown_genotype)
    ]
    output <- matching_rows
  }
  output <- unique(output, by = snps_file_match_fields)[
    order(chr_rank(chr), pos)
  ][
    , (empty_fields) := .(NA, NA, NA, 'GT')
  ]
  setcolorder(output, c(snps_file_match_fields, empty_fields, 'genotype'))
  writeLines(header, output_file)
  fwrite(output, output_file, sep = '\t', append = TRUE, na = '.', quote = FALSE)
}

vcf_snvs_to_bed <- function (
  vcf_file,
  output_folder,
  output_filename = 'vcf_regions.bed'
) {
  create_folder_if_not_exists(output_folder)
  output_file <- paste0(output_folder, '/', output_filename)
  vcf_content <- read_vcf_data(vcf_file)
  variants    <- with(vcf_content, vcf_data[vcf_is_snv(vcf_content$vcf_data)])
  output      <- copy(variants[, c('CHROM', 'POS')])[
    , start := POS - 1L
  ][
    , c('CHROM', 'start', 'POS')
  ][
    order(chr_rank(CHROM), POS)
  ]
  fwrite(output, output_file, sep = '\t', col.names = FALSE)
}

nucleotide_to_index <- function(num_rows, nucleotide, columns) {
  cbind(seq_len(num_rows), match(nucleotide, columns))
}

read_snv_coverage_data <- function (coverage_file) {
  nucleotide_cols <- c('A', 'C', 'G', 'T')
  ccc             <- fread(coverage_file)
  ccc_nuc_cov_df  <- as.data.frame(ccc[, ..nucleotide_cols])
  nr              <- nrow(ccc)
  cn              <- colnames(ccc)
  ref_cov_idx     <- nucleotide_to_index(nr, ccc$ref, nucleotide_cols)
  alt_cov_idx     <- nucleotide_to_index(nr, ccc$alt, nucleotide_cols)
  data.table(
    ccc[, -..nucleotide_cols],
    reference_cov   = ccc_nuc_cov_df[ref_cov_idx],
    alternative_cov = ccc_nuc_cov_df[alt_cov_idx]
  )
}

pileup_to_coverage <- function (
  snvs_pileup_path,
  output_folder,
  sample_name,
  output_filename = 'snv_coverage.tsv'
) {
  create_folder_if_not_exists(output_folder)
  output_file <- paste0(output_folder, '/', output_filename)
  output      <- with(
    read_snv_coverage_data(snvs_pileup_path),
    data.table(
      sample    = sample_name,
      chr       = chr,
      start     = pos,
      end       = pos,
      ref.count = reference_cov,
      alt.count = alternative_cov
    )
  )
  fwrite(
    output,
    output_file,
    sep   = '\t',
    quote = FALSE
  )
}

mutect_only_snvs <- function (
  vcf_file,
  output_folder,
  output_filename = 'vcf_only_snvs.vcf'
) {
  create_folder_if_not_exists(output_folder)
  output_file <- paste0(output_folder, '/', output_filename)
  vcf_content <- read_vcf_data(vcf_file)
  snv_lines   <- local({
    lll <- with(vcf_content, which(vcf_is_snv(vcf_data)))
    lll[vcf_content$vcf_data[lll][, order(chr_rank(CHROM), POS)]]
  })
  output <- with(vcf_content,
    raw_lines[
      c(heading_lines, snv_lines + length(heading_lines))
    ]
  )
  writeLines(output, output_file)
}

read_snv_coverage_data_to_vep <- function (coverage_file) {
  with(
    read_snv_coverage_data(coverage_file),
    data.table(
      Location = paste0(chr, ':', pos),
      Allele   = alt,
      rc_ref   = reference_cov,
      rc_alt   = alternative_cov,
      af       = af,
      cov      = cov
    )
  )
}

merge_vep_snv_coverage <- function (
  vep_annotated_variants,
  snv_coverage_normal_file,
  snv_coverage_tumor_file,
  output_folder,
  output_filename = 'vep_and_coverage.txt'
) {
  create_folder_if_not_exists(output_folder)
  output_file       <- paste0(output_folder, '/', output_filename)
  merge_fields      <- c('Location', 'Allele')
  snv_coverage_data <- local({
    ddd_n <- read_snv_coverage_data_to_vep(snv_coverage_normal_file)
    ddd_t <- read_snv_coverage_data_to_vep(snv_coverage_tumor_file)
    merge(
      ddd_n,
      ddd_t,
      by       = merge_fields,
      suffixes = c('_normal', '_tumor'),
      sort = FALSE
    )
  })
  vcf_content       <- read_vcf_data(vep_annotated_variants)
  output <- merge(
    vcf_content$vcf_data,
    snv_coverage_data,
    by    = merge_fields,
    all.x = TRUE,
    sort = FALSE
  )
  setcolorder(
    output,
    c(
      colnames(vcf_content$vcf_data),
      setdiff(colnames(snv_coverage_data), merge_fields)
    )
  )
  with(
    vcf_content,
    writeLines(
      c(
        head(raw_lines[heading_lines], -1),
        '#' %+% paste0(colnames(output), collapse = '\t')
      ),
     output_file
    )
  )
  fwrite(
    output,
    output_file,
    append = TRUE,
    sep   = '\t',
    na    = '-',
    quote = FALSE
  )
}

main_function <- function () {
  if (!is.null(cli_args$log_file)) {
    sink(cli_args$log_file)
  }
  greater_than_0        <- 'The value must be greater than 0'
  zero_to_one           <- 'The value must be between 0 and 1 included'
  min_depth_of_coverage <- parse_arg_pred(cli_args, 'min_depth_of_coverage', check_if_positive,
                             greater_than_0, converter = conv_to_integer)
  heterozygous_perc     <- parse_arg_pred(cli_args, 'heterozygous_perc', function (x) check_if_between(x, 0, 0.5),
                             'The value must be between 0 and 0.5 included', converter = conv_to_numeric)
  output_filename       <- cli_args$output_filename

  if (cli_args$pileup_to_genotype) {
    if (is.null(output_filename)) {
      output_filename <- 'genotype.vcf'
    }
    pileup_to_genotype(
      cli_args$SNPS_PILEUP,
      cli_args$OUTPUT_FOLDER,
      output_filename       = output_filename,
      sample_name           = cli_args$sample_name,
      snps_list_file        = cli_args$snps_list_file,
      min_depth_of_coverage = min_depth_of_coverage,
      heterozygous_perc     = heterozygous_perc
    )
  } else if (cli_args$pileup_to_coverage) {
    if (is.null(output_filename)) {
      output_filename <- 'snv_coverage.tsv'
    }
    pileup_to_coverage(
      cli_args$SNVS_PILEUP,
      cli_args$OUTPUT_FOLDER,
      sample_name     = cli_args$sample_name,
      output_filename = output_filename
    )
  } else if (cli_args$vcf_snvs_to_bed) {
    if (is.null(output_filename)) {
      output_filename <- 'vcf_regions.bed'
    }
    vcf_snvs_to_bed(
      cli_args$VCF_FILE,
      cli_args$OUTPUT_FOLDER,
      output_filename = output_filename
    )
  } else if (cli_args$mutect_only_snvs) {
    if (is.null(output_filename)) {
      output_filename <- 'vcf_only_snvs.vcf'
    }
    mutect_only_snvs(
      cli_args$VCF_FILE,
      cli_args$OUTPUT_FOLDER,
      output_filename = output_filename
    )
  } else if (cli_args$merge_vep_snv_coverage) {
    if (is.null(output_filename)) {
      output_filename <- 'variants_with_coverage.vcf'
    }
    merge_vep_snv_coverage(
      cli_args$VEP_ANNOTATED_VARIANTS,
      cli_args$SNV_COVERAGE_NORMAL_FILE,
      cli_args$SNV_COVERAGE_TUMOR_FILE,
      cli_args$OUTPUT_FOLDER,
      output_filename = output_filename
    )
  }
}

(make_main(docstring, c('data.table'), main_function))()
