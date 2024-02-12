#! /usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO

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
    log.info "nextflow run iarcbioinfo/template-nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--<OPTION>                      <TYPE>                      <DESCRIPTION>"
    log.info ""
    log.info "Optional arguments:"
    log.info "--<OPTION>                      <TYPE>                      <DESCRIPTION>"
    log.info ""
    log.info "Flags:"
    log.info "--<FLAG>                                                    <DESCRIPTION>"
    log.info ""
    exit 0
} else {
/* Software information */
log.info "input_file = ${params.input_file}"
log.info "output_folder    = ${params.output_folder}"
log.info "help:                               ${params.help}"
}

params.CNV_file = null
params.CNV_summary_file = null
params.vcf_folder = '.'
params.input_file = null
params.bam_input_folder = null
params.bam_folder = ""
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
        tuple val(sample), path("${sample}*.txt")
    publishDir "${params.output_folder}/DPclust_inputs", mode: "copy"
    
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
    tuple val(sample), path(DPinput)

    output:
    tuple val(sample), path("${sample}")
    publishDir "${params.output_folder}/results/DPclust", mode: "copy"

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
    tuple val(sample), path(DPoutput)

    output:
    tuple val(sample), path("${sample}*.pdf") , path("${sample}*.Rdata")
    publishDir "${params.output_folder}/results/DPclust", mode: "copy"

    shell:
    '''
    Rscript !{projectDir}/bin/postproc_DPclust.R
    mv Cluster_CCF.pdf !{sample}_Cluster_CCF.pdf
    mv dataset_fracs_list.Rdata !{sample}_dataset_fracs_list.Rdata
    '''
}

process mutationtimeR{
    cpus params.cpu
    memory params.mem+'GB'
    tag {sample}

    input:
    tuple path(DPinput)

    output:
    path "results*"
    publishDir "${params.output_folder}/results/MutationTimeR", mode: "copy"

    shell:
    '''
    Rscript !{projectDir}/bin/run_mutationtimeR.R -i !{DPinput}
    '''
}

workflow {
  // create channel with information about individual samples
  bams = Channel.fromPath("${params.input_file}")
                .splitCsv(header: true, sep: '\t', strip: true)
                        .map{ row -> [ row.sampleID , file(params.bam_folder+row.tumor), 
                        file(params.bam_folder+row.tumor+ext_ind), file(params.vcf_folder + row.vcf)] }

  preprocDPclust(params.CNV_file,params.CNV_summary_file,params.fai,params.ign, bams) | runDPclust | DPclust_postproc
}
