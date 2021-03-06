cwlVersion: v1.1
class: CommandLineTool
baseCommand: [gatk4, CollectHsMetrics]

label: Computes the capture statistics for hybrid capture sequencing data

doc: |-
  Picard is a set of command line tools for manipulating high-throughput
  sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.
  CollectHsMetrics is used to capture several metrics useful to verify the
  quality of target-capture sequencing experiments.

requirements:
  SchemaDefRequirement:
    types:
      - $import: "../types/hsmetrics_output_map.yaml"
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  LoadListingRequirement:
    loadListing: shallow_listing
  InitialWorkDirRequirement:
    listing: |
      ${
        return [
          {
            class: "Directory",
            basename: inputs.output_directory_name,
            listing: [],
            writable: true
          }
        ]
      }

hints:
  DockerRequirement:
    dockerPull: demichelislab/gatk4:latest

inputs:
  bam_file:
    doc: The BAM file that contains the reads of the sample to analyze.
    type: File
    secondaryFiles: [^.bai, .bai]
    inputBinding:
      prefix: "--INPUT"
  reference_genome_fasta_file:
    doc: The file containing the reference genome in FASTA format.
    type: File
    secondaryFiles: [^.fai, .fai]
    inputBinding:
      prefix: "--REFERENCE_SEQUENCE"
  kit_bait_interval_file:
    doc: |-
      The interval list file containing the regions that are captured by the
      kit.
    type: File
    inputBinding:
      prefix: "--BAIT_INTERVALS"
  kit_target_interval_file:
    doc: |-
      The interval list file containing the regions that are targeted by the
      kit.
    type: File
    inputBinding:
      prefix: "--TARGET_INTERVALS"
  output_filename:
    doc: The name of the output file.
    type: string
    default: "hsmetrics.txt"
    inputBinding:
      prefix: "--OUTPUT"
      valueFrom: $(inputs.output_directory_name)/$(self)
  lenient_validation:
    doc: |-
      If true, the validation strategy used by Picard HsMetrics will be the
      lenient strategy.
    type: boolean
    default: true
    inputBinding:
      prefix: "--VALIDATION_STRINGENCY"
      valueFrom: "LENIENT"
  per_target_coverage_output_name:
    doc: |-
      The name of the file that will contain the coverage per target region.
      The table will be generated only when this option is provided.
    type: string?
    inputBinding:
      prefix: "--PER_TARGET_COVERAGE"
      valueFrom: $(inputs.output_directory_name)/$(self)
  output_directory_name:
    doc: |-
      The name of the folder where the picard HsMetrics output will be stored.
    type: string
    default: "hsmetrics"
  log_to_file:
    doc: |-
      If true, the output generated by the tool will be redirected to a file.
      Otherwise the output will be printed on the output.
    type: boolean
    default: true
  redirect_stdout_to_stderr:
    doc: |-
      If true, it includes the stderr output along with the stdout in the log
      file.
    type: boolean
    default: true
  log_filename:
    doc: The name of the output file that will contain the output.
    type: string
    default: "picard.log"

arguments:
  - valueFrom: "$(inputs.log_to_file ? '2> ' + inputs.log_filename + (inputs.redirect_stdout_to_stderr ? ' 1>&2' : '') : '')"
    shellQuote: false
    position: 99999

outputs:
  output:
    doc: The output directory produced by Picard HsMetrics.
    type: Directory?
    outputBinding:
      glob: $(inputs.output_directory_name)
  output_map:
    doc: |-
      A data structure that will allow easy access to the various outputs
      generated by picard HsMetrics.
    type: out_map:hsmetrics_output_map?
    outputBinding:
      glob: $(inputs.output_directory_name)
      outputEval: |
        ${
          var filt_fun = function (filt) { return function (ff) { return ff.basename.endsWith(filt) }; },
              files    = self[0].listing;
          var result   = {
            table: files.filter(filt_fun(inputs.output_filename))[0]
          };
          if (inputs.per_target_coverage_output_name != null) {
            result.per_target_coverage = files.filter(filt_fun(inputs.per_target_coverage_output_name))[0];
          }
          return result;
        }
  log_file:
    doc: |-
      The log file, if enabled, that captures the output produced by the tool.
    type: File?
    outputBinding:
      glob: $(inputs.log_filename)

$namespaces:
  out_map: "../types/hsmetrics_output_map.yaml#"
