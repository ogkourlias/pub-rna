nextflow.enable.dsl=2

process fastq_dump {
  maxRetries 6
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '8 GB'
  cpus 1
  
  input:
  path sra_path

  output:
  path "${sra_path.SimpleName}", emit: fastq_dir
  val task.workDir, emit: work_dir

  script:
  """
  mkdir ${sra_path.SimpleName}
  fastq-dump --gzip --split-3 ${sra_path.SimpleName} -O ${sra_path.SimpleName} 
  """
}

process prefetch {
  maxRetries 6
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  memory '8 GB'
  cpus 1

  input:
  val id

  output:
  path "${id.trim()}*/*", emit: sra_path
  val task.workDir, emit: work_dir
  
  script:
  """
  prefetch ${id}
  """ 
}

process fasterq_dump {
  maxForks 20
  maxRetries 6
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  memory '16 GB'
  time '6h'
  cpus 4

  input:
  val id

  output:
  path "${id}", emit: fastq_dir
  val task.workDir, emit: work_dir
  
  script:
  """
  mkdir ${id}
  fasterq-dump ${id} -e 4 --include-technical -O ${id}
  gzip ${id}/*
  """ 
}

//workflow {
//    ids = Channel.fromPath(params.input_txt).splitText().distinct().flatten()
//    prefetch(ids)
//    fastq_dump(prefetch.output)
//    }