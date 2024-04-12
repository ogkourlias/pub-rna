workflow GETFASTQ {
    take:
      ids
    main:
      fasterq_dump(ids)
    emit:
      fasterq_dump.out.fastq_dir
      fasterq_dump.out.work_dir
}

workflow ALIGN {
  take:
    fastq_files
    fasterq_dump_work
  main:
    fastqcQualityControl(fastq_files)
    alignWithSTAR(fastq_files)
    sortBAM(alignWithSTAR.output.bam_file)
    markDuplicates(sortBAM.output.bam_file)
    QCwithRNASeqMetrics(markDuplicates.out.bamFile)
    QCwithMultipleMetrics(markDuplicates.out.bamFile)
    identifyAlternativeSplicingSitesrMATS(markDuplicates.out.bamFile)
    identifyAlternativeSplicingSitesLeafCutter(markDuplicates.out.bamFile)
  emit:
    markDuplicates.out.bamFile
    alignWithSTAR_work.out.work_dir
    markDuplicates.out.work_dir
    QCwithRNASeqMetrics.out.work_dir
    QCwithMultipleMetrics.out.work_dir
    identifyAlternativeSplicingSitesrMATS.out.work_dir
    identifyAlternativeSplicingSitesLeafCutter.out.work_dir
}

workflow TOGVCF {
  take:
    bam_files
    alignWithSTAR_work
    markDuplicates_work
    QCwithRNASeqMetrics_work
    QCwithMultipleMetrics_work
    identifyAlternativeSplicingSitesrMATS_work
    identifyAlternativeSplicingSitesLeafCutter_work
  main:
    indexBam(bam_files)
    encodeConvert(bam_files, indexBam.output.index_file)

    bamRm(alignWithSTAR_work,
    markDuplicates_work,
    QCwithRNASeqMetrics_work,
    QCwithMultipleMetrics_work,
    identifyAlternativeSplicingSitesrMATS_work,
    identifyAlternativeSplicingSitesLeafCutter_work,
    encodeConvert.out.work_dir)

    splitNCigarReads(encodeConvert.output.bam_file, encodeConvert.output.index_file) 
    AddOrReplaceReadGroups(splitNCigarReads.output.bam_file, splitNCigarReads.output.index_file)
    baseRecalibrator(AddOrReplaceReadGroups.output.bam_file, AddOrReplaceReadGroups.output.index_file, params.knownSitesIndex)
    applyBQSR(AddOrReplaceReadGroups.output.bam_file, AddOrReplaceReadGroups.output.index_file, baseRecalibrator.output.table_file)
    
    encodedRm(encodeConvert.out.work_dir,
    splitNCigarReads.out.work_dir,
    AddOrReplaceReadGroups.out.work_dir)

    haplotypeCaller(applyBQSR.output.bam_file, applyBQSR.output.index_file)
    gvcfs = haplotypeCaller.output.gvcf_file.collect()
  emit:
    gvcfs
}


workflow TOVCF {
  take:
    gvcfs
  main:
    chrs = Channel.from( 1..22 )
    genomicsDBImport(gvcfs, chrs)
    jointGenotype(genomicsDBImport.output.gdb_path)
}