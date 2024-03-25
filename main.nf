#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { prefetch; fastq_dump } from './modules/sra_query.nf'
include { convertBAMToFASTQ; fastqcQualityControl; alignWithSTAR; sortBAM; markDuplicates; QCwithRNASeqMetrics; QCwithMultipleMetrics; identifyAlternativeSplicingSitesrMATS; identifyAlternativeSplicingSitesLeafCutter; convertBAMToCRAM } from './modules/align.nf'
include { indexBam; encodeConvert; splitNCigarReads; AddOrReplaceReadGroups; baseRecalibrator; applyBQSR; haplotypeCaller; genomicsDBImport; jointGenotype } from './modules/genotype.nf'

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
    indexBam(markDuplicates.out.bamFile)
    encodeConvert(markDuplicates.out.bamFile, indexBam.output)
    splitNCigarReads(encodeConvert.output.bam_file, encodeConvert.output.index_file) 
    AddOrReplaceReadGroups(splitNCigarReads.output.bam_file, splitNCigarReads.output.index_file)
    baseRecalibrator(AddOrReplaceReadGroups.output.bam_file, AddOrReplaceReadGroups.output.index_file, params.knownSitesIndex)
    applyBQSR(AddOrReplaceReadGroups.output.bam_file, AddOrReplaceReadGroups.output.index_file, baseRecalibrator.output.table_file)
    haplotypeCaller(applyBQSR.output.bam_file, applyBQSR.output.index_file)
    gvcfs = haplotypeCaller.output.collect()
    chrs = Channel.from( 1..22 )
    genomicsDBImport(gvcfs, chrs)
    jointGenotype(genomicsDBImport.output)
    }