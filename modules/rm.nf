#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastqRm {
  errorStrategy 'finish'
  time '2h'
  memory '4 GB'
  cpus 1
  
  input:
  path fasterq_dump_work
  path fastqcQualityControl_work
  path alignWithSTAR_work


  
  script:
  """
  sleep 5
  rm -rf ${fasterq_dump_work}/
  rm -rf ${fastqcQualityControl_work}/
  """
}

process bamRm {
  errorStrategy 'finish'
  time '2h'
  memory '4 GB'
  cpus 1
  
  input:
  path alignWithSTAR_work
  path sortBAM_work
  path markDuplicates_work
  path QCwithRNASeqMetrics_work
  path QCwithMultipleMetrics_work
  path identifyAlternativeSplicingSitesrMATS_work
  path identifyAlternativeSplicingSitesLeafCutter_work
  path encodeConvert_work
  path splitNCigarReads_work



  
  script:
  """
  sleep 5
  rm -rf ${alignWithSTAR_work}/
  rm -rf ${sortBAM_work}/
  rm -rf ${markDuplicates_work}/
  rm -rf ${QCwithRNASeqMetrics_work}/
  rm -rf ${QCwithMultipleMetrics_work}/
  rm -rf ${identifyAlternativeSplicingSitesrMATS_work}/
  rm -rf ${identifyAlternativeSplicingSitesLeafCutter_work}/
  rm -rf ${encodeConvert_work}/
  """
}

process encodedRm {
  errorStrategy 'finish'
  time '2h'
  memory '4 GB'
  cpus 1
  
  input:
  path encodeConvert_work
  path splitNCigarReads_work
  path AddOrReplaceReadGroups_work
  path applyBQSR_work


  
  script:
  """
  sleep 5
  rm -rf ${encodeConvert_work}/
  rm -rf ${splitNCigarReads_work}/
  rm -rf ${AddOrReplaceReadGroups_work}/
  """
}

process recaldRm {
  errorStrategy 'finish'
  time '2h'
  memory '4 GB'
  cpus 1
  
  input:
  path baseRecalibrator_work
  path applyBQSR_work
  path haplotypeCaller_work


  
  script:
  """
  sleep 5
  rm -rf ${baseRecalibrator_work}/
  rm -rf ${applyBQSR_work}/
  """
}

process gdbRm {
  errorStrategy 'finish'
  time '2h'
  memory '4 GB'
  cpus 1
  
  input:
  path genomicsDBImport_work
  path jointGenotype_work


  
  script:
  """
  sleep 5
  rm -rf ${haplotypeCaller_work}/
  rm -rf ${genomicsDBImport_work}/
  rm -rf ${jointGenotype_work}/
  """
}


process bqsrRm {
  errorStrategy 'finish'
  time '2h'
  memory '2 GB'
  cpus 1
  
  input:
  path bqsr_work
  path haplotypeCaller_out


  
  script:
  """
  sleep 5
  rm -rf ${bqsr_work}/bqsr/*.bam
  """
}

