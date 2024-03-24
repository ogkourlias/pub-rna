#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { prefetch; fastq_dump } from './modules/sra_query.nf'
include { convertBAMToFASTQ; fastqcQualityControl; alignWithSTAR; sortBAM; markDuplicates; QCwithRNASeqMetrics; QCwithMultipleMetrics; identifyAlternativeSplicingSitesrMATS; identifyAlternativeSplicingSitesLeafCutter; convertBAMToCRAM } from './modules/align.nf'

workflow {
    ids = Channel.fromPath(params.input_txt).splitText().distinct().flatten()
    prefetch(ids)
    fastq_dump(prefetch.output)
    convertBAMToFASTQ(fastq_dump.output)
    fastqcQualityControl(convertBAMToFASTQ.output)
    alignWithSTAR(convertBAMToFASTQ.output)
    sortBAM(alignWithSTAR.output.bam_file)
    markDuplicates(sortBAM.output)
    QCwithRNASeqMetrics(markDuplicates.out.bamFile)
    QCwithMultipleMetrics(markDuplicates.out.bamFile)
    identifyAlternativeSplicingSitesrMATS(markDuplicates.out.bamFile)
    identifyAlternativeSplicingSitesLeafCutter(markDuplicates.out.bamFile)
    convertBAMToCRAM(markDuplicates.out.bamFile)
    }