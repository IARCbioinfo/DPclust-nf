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
log.info "help:                               ${params.help}"
}

params.CNV_file = null
params.CNV_summary_file = null
params.vcf_folder = '.'
params.bam_input_file = null
params.fai = ${baseDir}/data/hs38DH.fa.fai
params.ign = ${baseDir}/data/ign_file

params.ext = "cram"
ext_ind    = ".crai"
if(params.ext=="bam"){ ext_ind=".bai"}

// create channel with information about individual samples
bams = Channel.fromPath(params.bam_input_file).splitCsv(header: true, sep: '\t', strip: true)
                        .map{ row -> [ row.sample , file(row.tumor), file(row.tumor+ext_ind), file(row.vcf) ] }

// create channels with files

process preproc_DPclust {
    memory params.mem+'GB'
    cpus params.cpu
    tag { sample }

    input:
    path file_C
    path file_s
    path fai
    path ign
    string sample, file bam, file vcf from bams
  
    output:
    file '${sample}*.txt'

    shell:
    '''
    Rscript ${baseDir}/bin/preproc_DPclust.R -C ${file_C} -s ${file_s} -l ${sample} -v ${vcf} -t ${bam}
    '''
}

workflow main {
  preproc_DPclust(params.CNV_file,params.CNV_summary_file,params.fai,params.ign, bams)
}
