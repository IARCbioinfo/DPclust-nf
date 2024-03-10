# DPclust-nf
## Nextflow pipelines to run DPclust on sequencing data

## Description
Nextflow pipeline to run tumor subclone detection software DPclust, perform postprocessing and run mutationTimeR to date amplifications

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
- [DPclust](https://github.com/Wedge-lab/dpclust). R package to cluster mutations to detect tumor subclones
- [dpclust3p](https://github.com/Wedge-lab/dpclust3p). R package with useful functions to format DPclust inputs 
- [MutationTimeR](https://github.com/gerstung-lab/MutationTimeR). R package to time amplification events


## Input
  | Type      | Description     |
  |-----------|---------------|
  | CNV_file    | Path to tab-separated file with copy number variants for all samples, in PURPLE format (see e.g. our purple-nf pipeline) |
  | CNV_summary_file    | Path to tab-separated file with copy number summary for all samples (in PURPLE output format) |
  | input_file | Path to Input file (tab-separated values) with 4 columns: sampleID, vcf, normal, and tumor | 
  
## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --param1    |            xx | ...... |
| --param2    |            xx | ...... |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --vcf_folder   | . | Path to folder with variants for all samples in VCF format |
| --bam_folder    | . | Path to folder with BAM/CRAM files for tumor and normal samples |
| --output_folder    | . | Path to output folder |
| --cpu    | 2 | Number of cpus used |
| --mem    | 8 | Memory used in Gb |
| --ext    | cram | Extension of alignment files, cram or bam |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |

## Usage
  ```
  nextflow run iarcbioinfo/DPclust-nf [-with-docker] [OPTIONS]
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | DPclust_inputs    | Folder with inputs to DPclust formatted using dpclust3p: a folder for each sample with a DPmasterfile summarizing the files necessary to run DPclust, formatted CNVs and allele frequency tables |
  | results/DPclust    | Folder with outputs from DPclust: a folder for each sample with the location (CCF and number of alterations) of each subclone |
  | results/MutationTimeR    | Folder with results from mutationtimeR: timing of CNVs and SNVs |


## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Nicolas Alcala*    |            alcalan@iarc.who.int | Developer to contact for support |
