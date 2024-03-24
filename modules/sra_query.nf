nextflow.enable.dsl=2

process fastq_dump {
  scratch true
  
  errorStrategy 'finish'
  time '6h'
  memory '8 GB'
  cpus 1
  
  input:
  path sra_path

  output:
  path "${sra_path.SimpleName}"
  
  script:
  """
  mkdir ${sra_path.SimpleName}
  fastq-dump --gzip --split-3 ${sra_path.SimpleName} -O ${sra_path.SimpleName} 
  """
}

process prefetch {
  scratch true
  
  errorStrategy 'finish'
  time '6h'
  memory '8 GB'
  cpus 1

  input:
  val id

  output:
  path "SRR*/*"

  script:
  """
  prefetch ${id}
  """ 
}


//workflow {
//    ids = Channel.fromPath(params.input_txt).splitText().distinct().flatten()
//    prefetch(ids)
//    fastq_dump(prefetch.output)
//    }