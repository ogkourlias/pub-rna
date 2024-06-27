nextflow.enable.dsl=2

process cramToBam {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '4 GB'
  cpus 1
  
  output:
  path "${cram_file.SimpleName}.bam"
  val task.workDir, emit: work_dir

  input:
  path cram_file

  script:
  """
  samtools view -b  -T ${params.referenceGenome} -o ${cram_file.SimpleName}.bam ${cram_file}
  """
}

process indexBam {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '6 GB'
  cpus 1

  
  output:
  path "indexBam/${bam_file.SimpleName}.bai", emit: index_file
  val task.workDir, emit: work_dir

  input:
  path bam_file

  script:
  """
  mkdir indexBam
  gatk --java-options "-Xmx4g" BuildBamIndex \
  I=${bam_file} \
  O=indexBam/${bam_file.SimpleName}.bai
  """
}

process encodeConvert {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '6 GB'
  cpus 1
  
  input:
  path sample_bam

  output:
  path "fixed/${sample_bam.SimpleName}.bam", emit: bam_file
  path "fixed/${sample_bam.SimpleName}.bai", emit: index_file
  val task.workDir, emit: work_dir

  script:
  """
  gatk --java-options "-Xmx4g" BuildBamIndex \
  I=${sample_bam} \
  O=${sample_bam.SimpleName}.bai
  mkdir fastqc unzipped fixed
  fastqc ${sample_bam} --outdir fastqc
  unzip fastqc/*.zip -d unzipped
  num=\$(sed '6!d' unzipped/${sample_bam.SimpleName}_fastqc/fastqc_data.txt | cut -d " " -f 4)
  if [ echo "\$num < 1.8" | bc ]  
  then
      gatk --java-options "-Xmx4g" FixMisencodedBaseQualityReads \
      -I ${sample_bam} \
      -O fixed/${sample_bam.SimpleName}.bam
  else
      mv ${sample_bam} fixed/${sample_bam.SimpleName}.bam
  fi
  mv ${sample_bam.SimpleName}.bai fixed/${sample_bam.SimpleName}.bai
  """
}

process splitNCigarReads {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '16h'
  memory '16 GB'
  cpus 1

  input:
  path sample_bam
  path index_file

  output:
  path "split/${sample_bam.SimpleName}.bam", emit: bam_file
  path "split/${index_file.SimpleName}.bai", emit: index_file
  val task.workDir, emit: work_dir

  script:
  """
  mkdir split
  gatk --java-options "-Xmx12g" SplitNCigarReads \
  -R ${params.referenceGenome}\
  -I ${sample_bam} \
  -O split/${sample_bam.SimpleName}.bam \
  """
}

process AddOrReplaceReadGroups {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '6 GB'
  cpus 1

  input:
  path split_bam
  
  output:
  path "readgroup/${split_bam.SimpleName}.bam", emit: bam_file
  val task.workDir, emit: work_dir

  script:
  """
  mkdir readgroup
  gatk --java-options "-Xmx4g" AddOrReplaceReadGroups \
  I=${split_bam} \
  O=readgroup/${split_bam.SimpleName}.bam \
  RGID=4 \
  RGLB=lib1 \
  RGPL=ILLUMINA \
  RGPU=unit1 \
  RGSM=${split_bam.SimpleName}
  """
}


process baseRecalibrator {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path split_bam
  path vcf_index
  
  output:
  path "recal/${split_bam.SimpleName}.table", emit: table_file
  val task.workDir, emit: work_dir

  script:
  """
  mkdir recal
  gatk --java-options "-Xmx16g" BaseRecalibrator \
  -I ${split_bam} \
  -R ${params.referenceGenome} \
  --known-sites ${params.knownSites} \
  -O recal/${split_bam.SimpleName}.table
  """
}

process applyBQSR {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '6 GB'
  cpus 1

  input:
  path sample_bam
  path table_file

  output:
  path "bqsr/${sample_bam.SimpleName}.bqsr.bam", emit: bam_file
  val task.workDir, emit: work_dir

  script:
  """
  mkdir bqsr
  gatk --java-options "-Xmx4g" ApplyBQSR \
  -I ${sample_bam} \
  -R ${params.referenceGenome} --bqsr-recal-file ${table_file} \
  -O bqsr/${sample_bam.SimpleName}.bqsr.bam 
  """
}

process haplotypeCaller {
  storeDir  "${params.out_dir}/${sample_bam.SimpleName}/gvcf"
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '24h'
  memory '48 GB'
  cpus 4

  input:
  path sample_bam
  
  output:
  file "${sample_bam.SimpleName}.gvcf.gz"

  script:
  """
  gatk --java-options "-Xmx42g" HaplotypeCaller \
  -R ${params.referenceGenome}\
  -I ${sample_bam} \
  -O ${sample_bam.SimpleName}.gvcf.gz \
  -ERC GVCF \
  -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 \
  -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22
  """
}

process indexGvcf {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path gvcf_file
  
  output:
  path "${gvcf_file.SimpleName}.gvcf.gz.tbi"
  val task.workDir, emit: work_dir

  script:
  """
  gatk --java-options "-Xmx6g" IndexFeatureFile \
  -I ${gvcf_file} \
  -O ${gvcf_file.SimpleName}.gvcf.gz.tbi
  """
}


process indexJointGvcf {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path gvcf_file
  
  output:
  path "${gvcf_file.SimpleName}.gvcf.gz.tbi", emit: tbi_file
  path "${gvcf_file}", emit: gvcf_file
  val task.workDir, emit: work_dir

  script:
  """
  gatk --java-options "-Xmx6g" IndexFeatureFile \
  -I ${gvcf_file} \
  -O ${gvcf_file.SimpleName}.gvcf.gz.tbi
  """
}

process combineGvcf {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '16h'
  memory '12 GB'
  cpus 1

  input:
  path gvcfs
  path gvcf_idx
  val i
  
  output:
  path "chr${i}.gvcf.gz"
  path "chr${i}.gvcf.gz.tbi"
  val task.workDir, emit: work_dir

  script:
  """
  for GVCF in ${gvcfs}
  do
    bcftools view \$GVCF --regions chr${i} -o \$GVCF.chr${i}.gvcf.gz -Oz
    tabix \$GVCF.chr${i}.gvcf.gz
    echo \$GVCF.chr${i}.gvcf.gz >> gvcf.list
  done
  gatk --java-options "-Xmx10g" CombineGVCFs \
  -R ${params.referenceGenome} \
  -V gvcf.list \
  -O chr${i}.gvcf.gz
  """
}

process jointGenotype {
  storeDir "${params.out_dir}/genotypes"
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '72h'
  memory '18 GB'
  cpus 1

  input:
  path gvcf
  
  output:
  path "${gvcf.SimpleName}.vcf.gz"

  script:
  """
  gatk --java-options "-Xmx16g" GenotypeGVCFs \
  -R ${params.referenceGenome}\
  -V gendb://${gvcf} \
  -O ${gvcf.SimpleName}.vcf.gz
  """
}

process chrSplit {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path gvcf_dir
  val i
  
  output:
  path "chr${i}.gvcf.gz"
  val task.workDir, emit: work_dir

  script:
  """
  bcftools view ${gvcf_file} --regions chr${i} -o chr${i}.gvcf.gz -Oz
  """
}

process combineChrGvcf {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  publishDir "${chr_dir}", mode: 'move'
  time '16h'
  memory '12 GB'
  cpus 1

  input:
  path chr_dir
  
  output:
  path "combined/combined.vcf.gz"
  val task.workDir, emit: work_dir

  script:
  """
  mkdir combined
  for VCF in ${chr_dir}/*.vcf.gz
  do
    echo \$VCF >> vcf.list
  done
  gatk --java-options "-Xmx10g" CombineGVCFs \
  -R ${params.referenceGenome} \
  -V vcf.list \
  -O combined/combined_${chr_dir}.vcf.gz
  """
}

process genomicsDBImport_copy {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '36h'
  memory '48 GB'
  cpus 1

  input:
  path gvcfs
  val i

  output:
  path "chr${i}.gdb", emit: gdb_path
  val task.workDir, emit: work_dir
  
  script:
  """
  for GVCF in ${gvcfs}
  do
    tabix \$GVCF
    bcftools view \$GVCF --regions chr${i} -o \$GVCF.chr${i}.gvcf.gz -Oz
    echo "\${GVCF%%.*}\t\$GVCF" >> cohort.sample_map
  done
  gatk --java-options "-Xmx20g" GenomicsDBImport \
  --genomicsdb-workspace-path chr${i}.gdb \
  --batch-size 200 \
  --sample-name-map cohort.sample_map \
  -L chr${i} \
  --reader-threads 4
  """
}

process genomicsDBImport {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '72h'
  memory '32 GB'
  cpus 4

  input:
  path gvcfs
  val i

  output:
  path "chr${i}.gdb", emit: gdb_path
  
  script:
  """
  for GVCF in ${gvcfs}
  do
    tabix \$GVCF
    echo "\${GVCF%%.*}\t\$GVCF" >> cohort.sample_map
  done
  gatk --java-options "-Xmx28g" GenomicsDBImport \
  --genomicsdb-workspace-path chr${i}.gdb \
  --batch-size 200 \
  --sample-name-map cohort.sample_map \
  -L chr${i} \
  --reader-threads 4
  """
}


//workflow {
//    cram_files = Channel.fromPath("${params.cram_dir}/*/*.cram")
//    cramToBam(cram_files)
//    indexBam(cramToBam.output)
//    encodeConvert(cramToBam.output, indexBam.output)
//   splitNCigarReads(encodeConvert.output.bam_file, encodeConvert.output.index_file) 
//    AddOrReplaceReadGroups(splitNCigarReads.output.bam_file, splitNCigarReads.output.index_file)
//    baseRecalibrator(AddOrReplaceReadGroups.output.bam_file, AddOrReplaceReadGroups.output.index_file, params.knownSitesIndex)
//    applyBQSR(AddOrReplaceReadGroups.output.bam_file, AddOrReplaceReadGroups.output.index_file, baseRecalibrator.output.table_file)
//    haplotypeCaller(applyBQSR.output.bam_file, applyBQSR.output.index_file)
//    gvcfs = haplotypeCaller.output.collect()
//    chrs = Channel.from( 1..22 )
//    genomicsDBImport(gvcfs, chrs)
//    jointGenotype(genomicsDBImport.output)
//}
