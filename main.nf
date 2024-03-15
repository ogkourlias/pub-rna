#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { prefetch; fastq_dump } from './modules/sra_query.nf'
include { convertBAMToFASTQ } from './modules/align.nf'

workflow {
    ids = Channel.fromPath(params.input_txt).splitText().distinct().flatten()
    prefetch(ids)
    fastq_dump(prefetch.output)
    convertBAMToFASTQ(fastq_dump.output)
    }