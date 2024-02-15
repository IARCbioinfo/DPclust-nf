#! /usr/bin/env nextflow

// Copyright (C) 2024 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  DPclust-nf V1.0: nextflow pipeline to run DPclust     "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/DPclust-nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--CNV_file                      FILE                      Path to tab-separated file with copy number variants for all samples"
    log.info "--CNV_summary_file              FILE                      Path to tab-separated file with copy number summary for all samples (in PURPLE output format)"
    log.info "--input_file                    FILE                      Path to Input file (tab-separated values) with 4 columns: sampleID, vcf, normal, and tumor"
    log.info ""
    log.info "Optional arguments:"
    log.info "--vcf_folder                    FOLDER                    Path to folder with variants for all samples in VCF format (default: .)"
    log.info "--bam_folder                    FOLDER                    Path to folder with BAM/CRAM files for tumor and normal samples (default: .)"
    log.info '--output_folder                 FOLDER                    Output folder (default: .).'
    log.info '--cpu                           INTEGER                   Number of cpu used by bwa mem and sambamba (default: 2).'
    log.info '--mem                           INTEGER                   Size of memory used for mapping (in GB) (default: 8).'
    log.info ""
    log.info "Flags:"
    log.info "--help                                                    Display help"
    log.info ""
    exit 0
} else {
/* Software information */
    log.info "CNV_file         = ${params.CNV_file}"
    log.info "CNV_summary_file = ${params.CNV_summary_file}"
    log.info "input_file       = ${params.input_file}"

    log.info "vcf_folder       = ${params.vcf_folder}"
    log.info "bam_folder       = ${params.bam_folder}"
    log.info "output_folder    = ${params.output_folder}"
    log.info "cpu              = ${params.cpu}"
    log.info "mem              = ${params.mem}"
    log.info "help:              ${params.help}"
}

params.CNV_file = null
params.CNV_summary_file = null
params.vcf_folder = './'
params.input_file = null
params.bam_input_folder = null
params.bam_folder = "./"
params.cpu = 2
params.mem = 8
params.fai = "${projectDir}/data/hs38DH.fa.fai"
params.ign = "${projectDir}/data/ign_file"

params.ext = "cram"
ext_ind    = ".crai"
if(params.ext=="bam"){ ext_ind=".bai"}

process preprocDPclust{
    cpus params.cpu
    memory params.mem+'GB'
    tag {sample}

    input:
    path file_C
    path file_s
    path fai
    path ign
    tuple val(sample), file(bam), file(bai), file(vcf)

    output:
        tuple val(sample), path(vcf) , path("${sample}*.txt")
    publishDir "${params.output_folder}/DPclust_inputs/${sample}", mode: "copy", pattern: '{*.txt}'
    
   shell:
    '''
    ln -s !{projectDir}/bin/preproc_dpclust_master_file.R .
    Rscript !{projectDir}/bin/preproc_for_DPclust.R -c !{file_C} -s !{file_s} -n !{sample} -v !{vcf} -t !{bam} -f !{fai} -g !{ign}
    '''
}

process runDPclust{
    cpus params.cpu
    memory params.mem+'GB'
    tag {sample}

    input:
    tuple val(sample), path(vcf), path(DPinput)

    output:
    tuple val(sample), path(vcf), path("${sample}"), path(DPinput)
    publishDir "${params.output_folder}/results/DPclust", mode: "copy", pattern: "${sample}"

    shell:
    '''
    ln -s !{projectDir}/bin/dpclust_pipeline.R .
    Rscript !{projectDir}/bin/run_DPclust.R
    mv DPclust !{sample}
    '''
}

process DPclust_postproc{
    cpus params.cpu
    memory params.mem+'GB'
    tag {sample}

    input:
    tuple val(sample), path(vcf), path(DPoutput), path(DPinput)

    output:
    tuple val(sample), path(vcf), path("${sample}_cluster_locations.tsv")
    tuple path("${sample}*.Rdata"), path("${sample}*.pdf")
    publishDir "${params.output_folder}/results/DPclust/${sample}", mode: "copy", pattern: '*{.pdf,.Rdata,cluster_locations.tsv}'

    shell:
    '''
    Rscript !{projectDir}/bin/postproc_DPclust.R
    mv Cluster_CCF.pdf !{sample}_Cluster_CCF.pdf
    mv dataset_fracs_list.Rdata !{sample}_dataset_fracs_list.Rdata
    mv cluster_locations.tsv !{sample}_cluster_locations.tsv
    '''
}

process mutationtimeR{
    cpus params.cpu
    memory params.mem+'GB'
    tag {sample}

    input:
    path CNV
    path CNVsummary
    tuple val(sample), path(vcf), path(cluster_locations)

    output:
    path "${sample}*tim*"
    publishDir "${params.output_folder}/results/MutationTimeR", mode: "copy"

    shell:
    '''
    Rscript !{projectDir}/bin/run_mutationtimeR.R -n !{sample} -v !{vcf} -c !{CNV} -s !{CNVsummary} -k !{cluster_locations}
    touch done
    '''
}

workflow{
  // create channel with information about individual samples
  bams = Channel.fromPath("${params.input_file}")
                .splitCsv(header: true, sep: '\t', strip: true)
                        .map{ row -> [ row.sampleID , file(params.bam_folder+row.tumor), 
                        file(params.bam_folder+row.tumor+ext_ind), file(params.vcf_folder + row.vcf)] }

  preprocDPclust(params.CNV_file,params.CNV_summary_file,params.fai,params.ign, bams) | runDPclust | DPclust_postproc

  mutationtimeR(params.CNV_file,params.CNV_summary_file,DPclust_postproc.out[0] )

}