#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { prefetch; fastq_dump; fasterq_dump } from './modules/sra_query.nf'
include { convertBAMToFASTQ; fastqcQualityControl; alignWithSTAR; sortBAM; markDuplicates; QCwithRNASeqMetrics; QCwithMultipleMetrics; identifyAlternativeSplicingSitesrMATS; identifyAlternativeSplicingSitesLeafCutter; convertBAMToCRAM } from './modules/align.nf'
include { indexBam; encodeConvert; splitNCigarReads; AddOrReplaceReadGroups; baseRecalibrator; applyBQSR; haplotypeCaller; genomicsDBImport; jointGenotype } from './modules/genotype.nf'
include { fastqRm; bamRm; encodedRm; recaldRm; gdbRm; bqsrRm } from './modules/rm.nf'

workflow {
    ids = Channel.fromPath(params.input_txt).splitText().distinct().flatten().map { it.trim() }
    // Downloader
    fasterq_dump(ids)
    // Alignment
    fastqcQualityControl(fasterq_dump.out.fastq_dir)
    alignWithSTAR(fasterq_dump.out.fastq_dir)
    // Cleanup checkpoint. Fastq no longer needed.
    fastqRm(fasterq_dump.out.work_dir,
    fastqcQualityControl.out.work_dir,
    alignWithSTAR.out.work_dir)
    // Cleanup end.
    sortBAM(alignWithSTAR.output.bam_file)
    markDuplicates(sortBAM.output.bam_file)
    QCwithRNASeqMetrics(markDuplicates.out.bamFile)
    QCwithMultipleMetrics(markDuplicates.out.bamFile)
    identifyAlternativeSplicingSitesrMATS(markDuplicates.out.bamFile)
    identifyAlternativeSplicingSitesLeafCutter(markDuplicates.out.bamFile)

    // Genotyping
    encodeConvert(markDuplicates.out.bamFile)
    splitNCigarReads(encodeConvert.output.bam_file, encodeConvert.output.index_file) 
    // Cleanup checkpoint. BAM no longer needed.
    bamRm(alignWithSTAR.out.work_dir,
    sortBAM.out.work_dir,
    markDuplicates.out.work_dir,
    QCwithRNASeqMetrics.out.work_dir,
    QCwithMultipleMetrics.out.work_dir,
    identifyAlternativeSplicingSitesrMATS.out.work_dir,
    identifyAlternativeSplicingSitesLeafCutter.out.work_dir,
    encodeConvert.out.work_dir,
    splitNCigarReads.out.work_dir)
    // Cleanup end
    AddOrReplaceReadGroups(splitNCigarReads.output.bam_file)
    baseRecalibrator(AddOrReplaceReadGroups.output.bam_file, params.knownSitesIndex)
    applyBQSR(AddOrReplaceReadGroups.output.bam_file, baseRecalibrator.output.table_file)
    // Cleanup checkpoint. Encoded BAM no longer needed.
    encodedRm(encodeConvert.out.work_dir,
    splitNCigarReads.out.work_dir,
    AddOrReplaceReadGroups.out.work_dir,
    applyBQSR.out.work_dir)
    // Cleanup end
    haplotypeCaller(applyBQSR.output.bam_file)
    bqsrRm(applyBQSR.output.work_dir, haplotypeCaller.output)
    gvcfs = haplotypeCaller.output.collect()
    chrs = Channel.from( 1..22 )
    genomicsDBImport(gvcfs, chrs)
    jointGenotype(genomicsDBImport.output.gdb_path)
    }

    workflow GenotypeGVCFs{
    gvcfs = Channel.fromPath("asdf/").collect()
    chrs = Channel.from( 1..2 )
    genomicsDBImport(gvcfs, chrs)
    jointGenotype(genomicsDBImport.output.gdb_path)
    }

    workflow expanded{
    ids = Channel.fromPath(params.input_txt).splitText().distinct().flatten().map { it.trim() }
    // Downloader
    fasterq_dump(ids)
    // Alignment
    fastqcQualityControl(fasterq_dump.out.fastq_dir)
    alignWithSTAR(fasterq_dump.out.fastq_dir)
    // Cleanup checkpoint. Fastq no longer needed.
    fastqRm(fasterq_dump.out.work_dir,
    fastqcQualityControl.out.work_dir,
    alignWithSTAR.out.work_dir)
    // Cleanup end.
    sortBAM(alignWithSTAR.output.bam_file)
    markDuplicates(sortBAM.output.bam_file)
    QCwithRNASeqMetrics(markDuplicates.out.bamFile)
    QCwithMultipleMetrics(markDuplicates.out.bamFile)
    identifyAlternativeSplicingSitesrMATS(markDuplicates.out.bamFile)
    identifyAlternativeSplicingSitesLeafCutter(markDuplicates.out.bamFile)

    // Genotyping
    encodeConvert(markDuplicates.out.bamFile)
    splitNCigarReads(encodeConvert.output.bam_file, encodeConvert.output.index_file) 
    // Cleanup checkpoint. BAM no longer needed.
    bamRm(alignWithSTAR.out.work_dir,
    sortBAM.out.work_dir,
    markDuplicates.out.work_dir,
    QCwithRNASeqMetrics.out.work_dir,
    QCwithMultipleMetrics.out.work_dir,
    identifyAlternativeSplicingSitesrMATS.out.work_dir,
    identifyAlternativeSplicingSitesLeafCutter.out.work_dir,
    encodeConvert.out.work_dir,
    splitNCigarReads.out.work_dir)
    // Cleanup end
    AddOrReplaceReadGroups(splitNCigarReads.output.bam_file)
    baseRecalibrator(AddOrReplaceReadGroups.output.bam_file, params.knownSitesIndex)
    applyBQSR(AddOrReplaceReadGroups.output.bam_file, baseRecalibrator.output.table_file)
    // Cleanup checkpoint. Encoded BAM no longer needed.
    encodedRm(encodeConvert.out.work_dir,
    splitNCigarReads.out.work_dir,
    AddOrReplaceReadGroups.out.work_dir,
    applyBQSR.out.work_dir)
    // Cleanup end
    haplotypeCaller(applyBQSR.output.bam_file)
    bqsrRm(applyBQSR.output.work_dir, haplotypeCaller.output)
    gvcfs = haplotypeCaller.output.collect()
    chrs = Channel.from( 1..22 )
    genomicsDBImport(gvcfs, chrs)
    jointGenotype(genomicsDBImport.output.gdb_path)
    
    // Genotype QC
    fixGTAnnot(jointGenotype.out)
    splitMultiAllelicVariants(fixGTAnnot.output)
    filterVariants(splitMultiAllelicVariants.output)
    getMetrics(splitMultiAllelicVariants.output, filterVariants.output.filteredVcf)
    concatCHRFiles(filterVariants.output.filteredVcf.collect())
    convertToPlinkFormat(concatCHRFiles.output.sortedVcf)
    calculateMissingness(convertToPlinkFormat.output, splitMultiAllelicVariants.output)
    createHetFile(convertToPlinkFormat.output)
    findMissingSamples(calculateMissingness.output.missing)
    filterMissingSamples(findMissingSamples.output.filteredSamples, convertToPlinkFormat.output)
    findHetSamples(createHetFile.output)
    filterHetSamples(findHetSamples.output.failedHetSamples, filterMissingSamples.output)
    // createMetricsFile(filterHetSamples.output)
    filterLowAltFreq(filterHetSamples.output)
    filterRelated(filterLowAltFreq.output)
    beagle(filterRelated.output.bed, filterRelated.output.bim, filterRelated.output.fam)
    popProject(beagle.output.bed, beagle.output.bim, beagle.output.fam)
    // targetPCA(filterRelated.out.bed, filterRelated.out.bim, filterRelated.out.fam)
    finalPCA(beagle.out.bed, beagle.out.bim, beagle.out.fam)  
    }
