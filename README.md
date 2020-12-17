# CWL SPICE-pipeline

## Overview

The SPICE pipeline is designed to generate a comprehensive view of the genomic landscape of matched tumor and normal samples by leveraging allele-specific information from high quality next generation sequencing data. Specifically, using targeted DNA data (e.g. WES), the pipeline first applies tumor purity and ploidy correction, generates a quantitative measure of aneuploidy (asP), and then calls a series of genomic aberrations, including allele-specific copy number aberrations, SNVs with copy number corrected allelic fractions and indels. To allow for reliable analyses, the SPICE pipeline includes several QC tools.
The pipeline is written using the Common Workflow Language (CWL) {REF CWL}, a standard specification for the description of computational workflows that enables easily portable and scalable pipelines. Using one of the many available CWL implementations, it is possible to run SPICE on a variety of architectures (from single machines to clusters or cloud services) to easily scale up as needed. In order to enable ease of use and reproducible analyses, the tools that are used in the pipeline are ready on Docker Hub as containers.
To run the pipeline, it is sufficient to create a single configuration file per tumor/normal pair, where the user provides the required options (e.g. BAM files, reference genome).


## Main steps

This section includes a brief description for each of the main analysis tools
included in the pipeline.

### CLONET

  Analyzes genomic data from next-generation sequencing experiments.
  CLONETv2 offers a set of functions to compute allele specific copy number
  and clonality from segmented data leveraging heterozygous SNPs position
  pileups. The package also calculates the clonality of single nucleotide
  variants (SNVs) given read counts at mutated positions.

### CNVkit

  CNVkit is a Python library and command-line software toolkit to infer and
  visualize copy number from high-throughput DNA sequencing data. It is designed
  for use with hybrid capture, including both whole-exome and custom target
  panels, and short-read sequencing platforms such as Illumina and Ion Torrent.

### ETHSeq

  EthSEQ provides an automated pipeline, implemented as R package, to annotate
  the ethnicity of individuals from WES data inspecting differential SNPs
  genotype profiles while exploiting variants covered by the specific assay. 

### MuTect2

  MuTect2 Call somatic short mutations via local assembly of haplotypes. Short
  mutations include single nucleotide variant (SNVs) and insertion and deletion
  (indel) alterations.

### PaCBAM

  PaCBAM is a C command line tool for the complete characterization of genomic
  regions and single nucleotide positions from next-generation sequencing data.
  PaCBAM implements a fast and scalable multi-core computational engine,
  generates exhaustive output files for downstream analysis, introduces an
  innovative on-the-fly read duplicates filtering strategy and provides
  comprehensive visual reports. 

### Picard HSMetrics

  Picard is a set of command line tools for manipulating high-throughput
  sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.
  CollectHsMetrics is used to capture several metrics useful to verify the
  quality of target-capture sequencing experiments.

### SPIA

  SPIA allows for the verification of two or more DNA samples deriving from the
  same or different individuals.

### TPES

  A bioinformatics tool for the estimation of the tumor purity from sequencing
  data. It uses the set of putative clonal somatic single nucleotide variants
  within copy number neutral segments to call tumor cellularity.

### VEP

  VEP determines the effect of variants (SNPs, insertions, deletions, CNVs or
  structural variants) on genes, transcripts, and protein sequence, as well as
  regulatory regions.

## Configuration

Following there is a complete reference of the configuration options of the
pipeline. These options needed in order to be able to run the pipeline.

### Mandatory parameters

#### bam_file_normal

The BAM file of the normal sample. The folder should contain BAM index (`.bai`
file) along to the `.bam` file.

**Example**:

```yaml
bam_file_normal:
  class: File
  path: path/to/normal.bam
```

#### bam_file_tumor

The BAM file of the tumor sample. The folder should contain BAM index (`.bai`
file) along to the `.bam` file.


**Example**:

```yaml
bam_file_tumor:
  class: File
  path: path/to/tumor.bam
```

#### reference_genome_fasta_file

The file containing the reference genome used to align the BAM files in FASTA
format. The folder should contain the `.fai` index and the dictionary (`.dict`)
file.

Example:

```yaml
reference_genome_fasta_file:
  class: File
  path: path/to/reference_genome.fasta
```

#### kit_target_bed_file

The BED file containing the regions targeted by the capture kit.

**Example**:

```yaml
kit_target_bed_file:
  class: File
  path: path/to/target_regions.bed
```

#### kit_bait_bed_file

The BED file containing the regions captured by the capture kit.

**Example**:

```yaml
kit_bait_bed_file:
  class: File
  path: path/to/bait_regions.bed
```

#### kit_target_interval_file

Contains the regions that are targeted by the kit in the GATK interval list
format.

**Example**:

```yaml
kit_bait_interval_file:
  class: File
  path: path/to/target_regions.interval_list
```

#### kit_bait_interval_file

Contains the regions that are captured by the kit in the GATK interval list
format.

**Example**:

```yaml
kit_bait_interval_file:
  class: File
  path: path/to/bait_regions.interval_list
```

#### snps_in_kit_vcf_file

The list of SNPs that are contained in the regions targeted by the sequencing
kit. The VCF must contain only the SNPs and SNPs with more than one alternative
allele must be either removed or changed.

**Example**:

```yaml
snps_in_kit_vcf_file:
  class: File
  path: path/to/snps_in_kit.vcf
```

#### ethseq_snps_vcf_file

The file with the SNPs that are included in the EthSEQ model.

**Example**:

```yaml
ethseq_snps_vcf_file:
  class: File
  path: path/to/snps_in_ethseq_model.vcf
```

#### ethseq_snps_gds_file

The GDS file containing the model of the SNPs used for ethnicity inference.

**Example**:

```yaml
ethseq_snps_gds_file:
  class: File
  path: path/to/ethseq_model.gds
```

#### spia_snps_vcf_file

The file with the SNPs that are used by SPIA to compute the genotype distance.

**Example**:

```yaml
spia_snps_vcf_file:
  class: File
  path: path/to/spia_snps.vcf
```

#### sample_sex

The sex of the patient from which the sample was collected. Either "m" or "f".

**Example**:

```yaml
sample_sex: m
```

#### vep_reference_genome_version

The name of the reference genome that VEP have to to use for annotation.

**Example**:

```yaml
vep_reference_genome_version: GRCh38
```

#### vep_data_directory

The folder where VEP can find the annotation database.

**Example**:

```yaml
vep_data_directory:
  class: Directory
  path: path/to/vep_data_directory
```

### Optional parameters

#### cnvkit_accessible_regions_bed

An optional file containing the regions of the genome that are accessible
(meaning that can be sequenced).

**Example**:

```yaml
accessible_regions_bed:
  class: File
  path: path/to/accessible/regions.bed
```

#### threads

The number of parallel threads that will be used for computation in tools that
support parallel computation.

**Example**:

```yaml
threads: 5
```

#### create_reports

If true enables the creation of graphical reports. Only the tools that support
this type of output will include such outputs files. By default the option is
set to `false`

**Example**:

```yaml
create_reports: true
```

#### log_to_file

If true, the output generated by each tool will be redirected to a file.
Otherwise the output will be printed on the output. By default the options is
set to `true`.

**Example**:

```yaml
log_to_file: false
```

### The complete example configuration

```yaml
bam_file_normal:
  class: File
  path: path/to/normal.bam
bam_file_tumor:
  class: File
  path: path/to/tumor.bam
reference_genome_fasta_file:
  class: File
  path: path/to/reference_genome.fasta
kit_target_bed_file:
  class: File
  path: path/to/target_regions.bed
kit_bait_bed_file:
  class: File
  path: path/to/bait_regions.bed
kit_bait_interval_file:
  class: File
  path: path/to/target_regions.interval_list
kit_bait_interval_file:
  class: File
  path: path/to/bait_regions.interval_list
snps_in_kit_vcf_file:
  class: File
  path: path/to/snps_in_kit.vcf
ethseq_snps_vcf_file:
  class: File
  path: path/to/snps_in_ethseq_model.vcf
ethseq_snps_gds_file:
  class: File
  path: path/to/ethseq_model.gds
spia_snps_vcf_file:
  class: File
  path: path/to/spia_snps.vcf
sample_sex: m
vep_reference_genome_version: GRCh38
vep_data_directory:
  class: Directory
  path: path/to/vep_data_directory
accessible_regions_bed:
  class: File
  path: path/to/accessible/regions.bed
threads: 5
create_reports: true
log_to_file: false
```

## Running the pipeline

The pipeline is expected to be run on a linux operating system (no tests were
run on different operating systems) using one of the available implementations
that can be found [here](https://www.commonwl.org/#Implementations). Multiple
versions of the CWL standard are available. The SPICE CWL pipeline is based on
CWL version `1.1`. Before running the pipeline make sure that the selected
implementation supports version `1.1` of the specification.

In order to run the pipeline you need the CWL files where the pipeline is
described (these are available in this repository in the cwl folder). The
configuration needs to be saved to a `.yaml` file. In the examples below we show
how to run the pipeline using `cwltool` (the official CWL implementation).

### Standard run

If run in this way the pipeline launches all tools within docker containers.

```bash
cwltool path/to/workflows/pipeline.cwl path/to/parameters.yaml
```

### Running without containers

In order to run without using containers add the `--no-container` option to the
command line as shown below. In order to run in this way all the tools used by
the pipeline must be available as commands (the executables need to be in one of
the folders included in the `$PATH` environment variable).

```bash
cwltool --no-container path/to/workflows/pipeline.cwl path/to/parameters.yaml
```

### Running with other containers

In order to run using a different container runtime just use the specific
option. For example to run using singularity just add the `--singularity`
option like shown below. The `cwltool` runner supports other runtimes but only
docker and singularity have been tested.

```bash
cwltool --singularity path/to/workflows/pipeline.cwl path/to/parameters.yaml
```

## Outputs

The pipeline will create a folder with the output of each tool. Below is
represented the tree of folders created as output of the pipeline. Each
subfolder of the `data/` folder contains the output of the corresponding tool.
If option `log_to_file` is set to `true` the pipeline will create a folder named
`logs` with a log file for each step that is part of the pipeline.

```ini
data/
|-- clonet           # Purity, ploidy, corrected log2, allele specific cn and
|                      # clonality
|-- cnvkit           # Copy number segments
|-- ethseq           # Ethnicity information
|-- hsmetrics_normal # QC metrics normal sample
|-- hsmetrics_tumor  # QC metrics normal sample
|-- mutect2          # SNVs and indel calls
|-- snps_pileups     # Pileup of the SNPs that are within kit regions.
|-- snvs_coverage    # Coverage data in tumor and normal samples of SNV
|                      # positions
|-- spia             # Genotype distance between normal and tumor sample.
|-- tpes             # Purity based on SNVs.
`-- vep              # Annotation for SNVs and indels.
logs/                # Contains log for each step if option log_to_file is set
                       # to true
```
## Fundings

This project is funded by the ERC (ERC-CoG-2014-648670), to F. Demichelis.
