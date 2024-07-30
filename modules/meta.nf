nextflow.enable.dsl=2

process createBatches {
//  maxRetries 2
//  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '2h'
  memory '8 GB'
  cpus 1

  input:
  path in_dir
  
  output:
  path 'batches/*.txt'
  
  script:
  """
  createBatches.py \
  ${params.name} \
	${params.exp} \
	${params.gte} \
	${params.in_dir}/chrCHR.vcf.gz \
	${params.genelist} \
	${params.annotation} \
	jobtemplate.sh \
	${params.batch_size} \
	./
  """
}

process batchMap {
//  maxRetries 2
//  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '4 GB'
  cpus 1

  input:
  path batch
  
  output:
  path '*.txt'
  
  script:
  """
  NUMBER=\$(echo "${batch.SimpleName}" | tr -dc '0-9')
  java -Xmx2g -jar ${baseDir}/bin/MbQTL-1.4.5-SNAPSHOT-jar-with-dependencies.jar \
        -m mbqtl \
        -a ${params.annotation} \
        -e ${params.exp} \
        -g ${params.gte} \
        -v ${params.in_dir}/${batch.SimpleName}.vcf.gz \
        --chr \$NUMBER \
        --perm ${params.perms} \
        --minobservations ${params.minob} \
        -gl ${batch} \
        --replacemissinggenotypes \
        -o ${batch.BaseName} \
        --mingenotypecount ${params.mingt} \
        --fisherzmeta
  """
}

process combineOut {
  storeDir "${params.out_dir}/meta"
//  maxRetries 2
//  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '4 GB'
  cpus 1

  input:
  path batch
  
  output:
  path "${params.name}-merged-TopEffects.txt"
  
  script:
  """
  merge.py ./ ${params.name}-merged-TopEffects.txt
  """
}

process qval {
  storeDir "${params.out_dir}/meta"
//  maxRetries 2
//  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '16 GB'
  cpus 1

  input:
  path batch
  
  output:
  path "${params.name}-qval-merged-TopEffects.txt"
  
  script:
  """
  getQVal.R ${params.name}-merged-TopEffects.txt ${params.name}-qval-merged-TopEffects.txt
  """
}

process sig {
  storeDir "${params.out_dir}/meta"
//  maxRetries 2
//  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '16 GB'
  cpus 1

  input:
  path effects
  
  output:
  path "${params.name}-sig-qval-merged-TopEffects.txt"
  
  script:
  """
  countSig.py ${effects} ${params.name}-sig-qval-merged-TopEffects.txt
  """
}


workflow {
  vcfs = Channel.fromPath("${params.in_dir}/*.vcf.gz").collect()
  createBatches(vcfs)
  batchMap(createBatches.out.flatten())
  combineOut(batchMap.out.collect())
  qval(combineOut.out)
  sig(qval.out)
}