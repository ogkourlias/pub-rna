nextflow.enable.dsl=2

process convertBAMToFASTQ {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '8 GB'
  cpus 1
  
  input:
  path sample_dir
  
  output:
  path "fastq_output/${sample_dir}"

  
  shell:
  '''
  PAIRFLAG=0
  FASTQFLAG=0
  for file in !{sample_dir}/*; do
    if [[ "$file" == *"_2"* ]]; then
      PAIRFLAG=1
    fi
    if [[ "$file" == *".fastq.gz" ]]; then
      FASTQFLAG=1
    fi
  done
  mkdir -p fastq_output/!{sample_dir}
  if [[ $FASTQFLAG -eq 1 ]]; then
    cp !{sample_dir}/*.gz fastq_output/!{sample_dir}
  elif [[ $PAIRFLAG -eq 1 && $FASTQFLAG -eq 0 ]]; then
    samtools sort -n !{sample_dir}/!{sample_dir}_1.bam -o sorted_!{sample_dir}_1.bam
    samtools sort -n !{sample_dir}/!{sample_dir}_2.bam -o sorted_!{sample_dir}_2.bam
    samtools fastq -1 fastq_output/!{sample_dir.baseName}_1.fastq.gz -2 fastq_output/!{sample_dir.baseName}_2.fastq.gz !{sample_dir}.bam
  elif [[ $PAIRFLAG -eq 0 && $FASTQFLAG -eq 0 ]]; then
    samtools sort -n !{sample_dir}/!{sample_dir}.fastq.gz -o sorted_!{sample_dir}.fastq.gz
    samtools fastq -0 fastq_output/!{sample_dir.baseName}.fastq.gz !{sample_dir}/!{sample_dir}.fastq.gz
  fi
  '''
}


process fastqcQualityControl {
  publishDir "${params.out_dir}/${sample_dir}", mode: 'move'
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path sample_dir

  output:
  file "fastqc/*.zip"
  val task.workDir, emit: work_dir

  shell:
  '''
  mkdir fastqc
  for FILE in !{sample_dir}/*.fastq.gz; do
    fastqc $FILE --outdir fastqc --5
  done
  '''
}

process alignWithSTAR {
  
  publishDir "${params.out_dir}/${sample_dir}/star", mode: 'move', pattern: "*.{gz}"
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '48 GB'
  cpus 4

  input:
  path sample_dir

  output:
  path "${sample_dir}.bam", emit: bam_file
  path "*.gz"
  val task.workDir, emit: work_dir

  shell:
  '''
  # Get path of the first file in the input directory
  firstFile=$(ls -1 "!{sample_dir}" | sort | head -n 1)

  # Determine allowed number of mismatches based on read length
  readLength=$(samtools view !{sample_dir}/${firstFile} |head -n1 |awk '{print $10}'|tr -d "\\n" |wc -m)

  if [ $readLength -ge 90 ]; then
    numMism=4
  elif [ $readLength -ge 60 ]; then
    numMism=3
  else
    numMism=2
  fi

  # Set different --readFilesIn values for paired and single end
  if [[ -f !{sample_dir}/!{sample_dir}_2.fastq.gz ]]; then
     let numMism=$numMism*2
     readFilesInArgument="--readFilesIn !{sample_dir}/!{sample_dir}_1.fastq.gz !{sample_dir}/!{sample_dir}_2.fastq.gz"
  else
     readFilesInArgument="--readFilesIn !{sample_dir}/!{sample_dir}.fastq.gz"
  fi

  # Run the STAR command
  STAR --runThreadN 4 \
  --outFileNamePrefix !{sample_dir} \
  --outSAMtype BAM Unsorted \
  --genomeDir !{params.refDir} \
  --genomeLoad NoSharedMemory \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNmax ${numMism} \
  --twopassMode Basic \
  --quantMode GeneCounts \
  --outSAMunmapped Within \
  --readFilesCommand zcat \
  ${readFilesInArgument}

  mv !{sample_dir}Aligned.out.bam !{sample_dir}.bam

  # Gzip all output files
  gzip *.tab
  gzip *.out
'''
}

process sortBAM {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path sample_bam

  output:
  path "${sample_bam.SimpleName}.sorted.bam", emit: bam_file
  val task.workDir, emit: work_dir

  script:
  """
  samtools sort ${sample_bam} -o ${sample_bam.SimpleName}.sorted.bam
  """
}

process markDuplicates {
  publishDir "${params.out_dir}/${bam_file.SimpleName}/mark_duplicates", mode: 'move', pattern: "*.{gz}"
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '12 GB'
  cpus 1

  input:
  path bam_file

  output:
  path "${bam_file.SimpleName}.duplicates.bam", emit: bamFile
  path "${bam_file.SimpleName}_duplicates.txt.gz"
  val task.workDir, emit: work_dir

  script:
  """
  java -Xmx10g -jar /usr/bin/picard.jar MarkDuplicates \
      I=${bam_file} \
      O=${bam_file.SimpleName}.duplicates.bam \
      M=${bam_file.SimpleName}_duplicates.txt

  gzip ${bam_file.SimpleName}_duplicates.txt
  """
}

process QCwithRNASeqMetrics {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '12 GB'
  cpus 1

  publishDir "${params.out_dir}/${bam_file.SimpleName}", mode: 'move'

  input:
  path bam_file

  output:
  path "${bam_file.SimpleName}_rnaseqmetrics.gz"
  path "${bam_file.SimpleName}.chart.pdf.gz"
  val task.workDir, emit: work_dir
  
  script:
  """
  java -Xmx10g -jar /usr/bin/picard.jar CollectRnaSeqMetrics \
  I=${bam_file} \
  O=${bam_file.SimpleName}_rnaseqmetrics \
  CHART_OUTPUT=${bam_file.SimpleName}.chart.pdf \
  REF_FLAT=${params.refFlat} \
  STRAND=NONE \
  RIBOSOMAL_INTERVALS=${params.ribosomalIntervalList}

  gzip ${bam_file.SimpleName}.chart.pdf
  gzip ${bam_file.SimpleName}_rnaseqmetrics
  """
}

process QCwithMultipleMetrics {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '12 GB'
  cpus 1

  publishDir "${params.out_dir}/${bam_file.SimpleName}", mode: 'move'

  input:
  path bam_file

  output:
  path "${bam_file.SimpleName}/multiple_metrics*"
  val task.workDir, emit: work_dir
  
  script:
  """
  mkdir ${bam_file.SimpleName}
  java -Xmx10g -jar /usr/bin/picard.jar CollectMultipleMetrics I=${bam_file} \
  O=${bam_file.SimpleName}/multiple_metrics \
  R=${params.referenceGenome} \
  PROGRAM=CollectAlignmentSummaryMetrics \
  PROGRAM=QualityScoreDistribution \
  PROGRAM=MeanQualityByCycle \
  PROGRAM=CollectInsertSizeMetrics 

  gzip ${bam_file.SimpleName}/*
  """
}

process identifyAlternativeSplicingSitesrMATS {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "${params.out_dir}/${bam_file.SimpleName}", mode: 'move'
  
  input:
  path bam_file
  
  output:
  path "${bam_file.SimpleName}/*.txt.gz"
  val task.workDir, emit: work_dir
  
  shell:
  '''
  # 1. Check if the BAM file is derived from single or paired end reads
  numFASTQFiles=$(samtools view -H !{bam_file} | \
  grep "@PG" | \
  tr ' ' '\n' | \
  grep -oE '(.fq.gz|.fastq.gz|.fq|.fastq)($)' | \
  wc -l)

  if [ "${numFASTQFiles}" -eq 1 ]; 
  then
    end="single"
  elif [ "${numFASTQFiles}" -gt 1 ];
  then
    end="paired"
  else
    exit 1
  fi

  # 2. Check the read length
  readLength=$(samtools view !{bam_file} |head -n1 |awk '{print $10}'|tr -d "\\n" |wc -m)

  # 3. Create config file
  echo !{bam_file} > config.txt

  # 4. Run rMATS command
  python /opt/rmats_turbo_v4_1_2/rmats.py --b1 config.txt \
  --gtf !{params.gtfAnnotationFile} \
  --readLength ${readLength} \
  --od !{bam_file.SimpleName} \
  --tmp rmats_tmp \
  --task both \
  -t ${end}  \
  --statoff

  # 5. Gzip all output files
  gzip !{bam_file.SimpleName}/*.txt
  '''
}

process identifyAlternativeSplicingSitesLeafCutter {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "${params.out_dir}/${bam_file.SimpleName}", mode: 'move'

  input:
  path bam_file
  
  output:
  path "${bam_file.SimpleName}/*.junc.gz"
  val task.workDir, emit: work_dir
  
  shell:
  '''
  # 1. Index BAM file
  samtools index !{bam_file}

  # 2. Run regtools command
  mkdir !{bam_file.SimpleName}
  regtools junctions extract -s XS -a 8 -m 50 -M 500000 !{bam_file} -o !{bam_file.SimpleName}/!{bam_file.SimpleName}.junc 

  # 3. Gzip the resulting junctions file
  gzip !{bam_file.SimpleName}/!{bam_file.SimpleName}.junc
  '''
}

process convertBAMToCRAM {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '6h'
  memory '10 GB'
  cpus 1

  publishDir "${params.out_dir}/${bam_file.SimpleName}", mode: 'move'

  input:
  path bam_file
  
  output:
  path "${bam_file.SimpleName}"
  val task.workDir, emit: work_dir
  
  script:
  """
  mkdir cram
  samtools view -T ${params.referenceGenome} -C -o cram/${bam_file.SimpleName}.cram ${bam_file}
  """
}

def extractSampleName(String path) {
    // Split the path by ';' to handle multiple paths
    def paths = path.split(';')

    // Initialize a set to store unique sample names
    Set<String> sampleNameSet = new HashSet<>()
    
    // Extract and process the sample names for each path
    paths.each { pathPart ->
        // Extract the sample name and remove _1 or _2 if present
        def sampleNamePart = pathPart.tokenize("/").last().replaceAll(/\.(fastq|fq|bam|gz)/, '').replaceAll(/_(1|2)$/, '')
        
        // Add the cleaned sample name to the set
        sampleNameSet << sampleNamePart
    }
    
    // Join the unique sample names with '_'
    def sampleName = sampleNameSet.join("_")
    
    return sampleName
}

def checkIfSampleIsProcessed(String folderName, String sampleName) {
    sampleName = extractSampleName(sampleName)
    
    // An array containing all the folders that are expected for a succesful pipeline run for a sample
    def expectedFolders = [
        folderName + '/fastqc/' + sampleName,
        folderName + '/star/' + sampleName,
        folderName + '/multiple_metrics/' + sampleName,
        folderName + '/rna_seq_metrics/' + sampleName,
        folderName + '/rmats/' + sampleName,
        folderName + '/mark_duplicates/' + sampleName,
        folderName + '/cram/' + sampleName,
    ];

    // An array containing the expected number of files in order of the folders above
    def expectedNumberOfFiles = [1, 5, 8, 2, 36, 1, 1];

    // Loop through expected folders and number of expected files
    for (int i = 0; i < expectedNumberOfFiles.size; i++) {
        
        // Return false if the expected folder does not exist
        if (!new File(expectedFolders[i]).exists()) {
           return false;
        }

        // Return flase if the number of items in folder does not match expected number of items
        if (new File(expectedFolders[i]).list().length < expectedNumberOfFiles[i]){
           return false;
        }
      }
      return true;
}

process removeWorkDirs {
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '1h'
  memory '1 GB'
  cpus 1

  input:
  val bamToFastqWorkDir
  val fastQCWorkDir
  val alignWithStarWorkDir
  val sortBamWorkDir
  val markDuplicatesWorkDir
  val QCwithRNASeqMetricsWorkDir
  val QCwithMultipleMetricsWorkDir
  val identifyAlternativeSplicingSitesrMATSWorkDir
  val identifyAlternativeSplicingSitesLeafCutterWorkDir
  val convertBAMToCRAMWorkDir
  
  script:
  """
  sleep 5
  rm -r ${bamToFastqWorkDir} || echo 'Failed to remove work directory'
  rm -r ${fastQCWorkDir} || echo 'Failed to remove work directory'
  rm -r ${alignWithStarWorkDir} || echo 'Failed to remove work directory'
  rm -r ${sortBamWorkDir} || echo 'Failed to remove work directory'
  rm -r ${markDuplicatesWorkDir} || echo 'Failed to remove work directory'
  rm -r ${QCwithRNASeqMetricsWorkDir} || echo 'Failed to remove work directory'
  rm -r ${QCwithMultipleMetricsWorkDir} || echo 'Failed to remove work directory'
  rm -r ${identifyAlternativeSplicingSitesrMATSWorkDir} || echo 'Failed to remove work directory'
  rm -r ${identifyAlternativeSplicingSitesLeafCutterWorkDir} || echo 'Failed to remove work directory'
  rm -r ${convertBAMToCRAMWorkDir} || echo 'Failed to remove work directory'
  """
}

workflow {
    // Load list with sample paths from the input text file
    String sample_dirs = new File(params.sampleFile).text
    String[] sample_dirsArray = sample_dirs.split('\n')

    // Create a list of sample names extracted from the sample paths
    List<String> sampleNamesList = sample_dirsArray.collect { extractSampleName(it) }
    String[] sampleNamesArray = sampleNamesList as String[]
    
    // Create channels of the sample paths and sample names
    sample_dirsChannel = Channel.of(sample_dirsArray)
    sampleNamesChannel = Channel.of(sampleNamesArray)

    // Remove samples from the channels that are already in the output folder
    filteredPathsChannel = sample_dirsChannel.filter { !checkIfSampleIsProcessed(params.out_dir, it) }
    filteredNamesChannel = sampleNamesChannel.filter { !checkIfSampleIsProcessed(params.out_dir, it) }

    // Run pipeline
    convertBAMToFASTQ(filteredPathsChannel, filteredNamesChannel)
    fastqcQualityControl(convertBAMToFASTQ.out.fastqPath, convertBAMToFASTQ.out.sampleName)
    alignWithSTAR(convertBAMToFASTQ.out.fastqPath, convertBAMToFASTQ.out.sampleName)
    sortBAM(alignWithSTAR.out.bamFile, alignWithSTAR.out.sampleName)
    markDuplicates(sortBAM.out.bamFile, sortBAM.out.sampleName)
    QCwithRNASeqMetrics(markDuplicates.out.bamFile, markDuplicates.out.sampleName)
    QCwithMultipleMetrics(markDuplicates.out.bamFile, markDuplicates.out.sampleName)
    identifyAlternativeSplicingSitesrMATS(markDuplicates.out.bamFile, markDuplicates.out.sampleName)
    identifyAlternativeSplicingSitesLeafCutter(markDuplicates.out.bamFile, markDuplicates.out.sampleName)
    convertBAMToCRAM(markDuplicates.out.bamFile, markDuplicates.out.sampleName)

    // Remove all work directories for sample
    removeWorkDirs(
      convertBAMToFASTQ.out.workDir,
      fastqcQualityControl.out.workDir,
      alignWithSTAR.out.workDir,
      sortBAM.out.workDir,
      markDuplicates.out.workDir,
      QCwithRNASeqMetrics.out.workDir,
      QCwithMultipleMetrics.out.workDir,
      identifyAlternativeSplicingSitesrMATS.out.workDir,
      identifyAlternativeSplicingSitesLeafCutter.out.workDir,
      convertBAMToCRAM.out.workDir,
    )
}