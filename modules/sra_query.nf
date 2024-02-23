nextflow.enable.dsl=2

process fastq_dump {
  debug true
  containerOptions '--bind /groups/'
  publishDir "${params.out_dir}", mode: 'move'
  errorStrategy 'finish'
  time '6h'
  memory '8 GB'
  cpus 1
  
  input:
  path sra_path

  output:
  path "*.fastq.gz"
  
  script:
  """
  fastq-dump --gzip --split-3 ${sra_path.SimpleName}
  rm -rf sra_path ${sra_path}
  """
}

process prefetch {
  debug true
  containerOptions '--bind /groups/'
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


workflow {
    ids = Channel.fromPath(params.input_txt).splitText().distinct().flatten()
    prefetch(ids)
    fastq_dump(prefetch.output)
    }