"""
File:         main.nf
Created:      2024/04/11
Last Changed: 2024/04/25
Author:       Peter Riesebos & Orfeas Gkourlias
"""

nextflow.enable.dsl=2

process concatCHRFiles {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '4 GB'
    cpus 1

    input:
    path vcfsPath

    output:
    path "${params.project}.sorted.concat.vcf.gz", emit: sortedVcf

    script:
    """
    mkdir -p  tmp_sort
    for i in {21..22}; do echo "chr\$i.filtered.vcf.gz" >> vcflist.txt; done
    bcftools concat -f vcflist.txt -Oz -o "${params.project}.sorted.concat.vcf.gz"
    """
}

process fixGTAnnot {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcf

    output:
    path "${vcf.SimpleName}.gtfixed.sorted.concat.vcf.gz", emit: unfilteredVcf

    script:
    """
    # 1. fix GT field annotation
    dotrevive.py -i ${vcf} -o ${vcf.SimpleName}.gtfixed.sorted.concat.vcf.gz
    """
}

process splitMultiAllelicVariants {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcf

    output:
    path "${vcf.SimpleName}.split.gtfixed.sorted.concat.vcf.gz", emit: splitVcf

    script:
    """
    bcftools norm -m -any ${vcf} -Oz -o ${vcf.SimpleName}.split.gtfixed.sorted.concat.vcf.gz
    """
}


process filterVariants {
    publishDir "${params.outDir}/variant_filter", mode: "move", pattern: "*.{log.gz}"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcf

    output:
    path "${vcf.SimpleName}.filtered.vcf.gz", emit: filteredVcf
    path "${vcf.SimpleName}.filtered.log.gz", emit: filteredLog

    script:
    """
    mkdir -p  ${params.outDir}/${vcf.SimpleName}
    mkdir -p  ${params.outDir}/${vcf.SimpleName}/variant_filter
    mkdir -p  ${params.outDir}/${vcf.SimpleName}/get_metrics
    mkdir -p ${params.outDir}/${vcf.SimpleName}/missingness
    mkdir -p  ${params.outDir}/${vcf.SimpleName}/het/
    mkdir -p ${params.outDir}/${vcf.SimpleName}/related
    mkdir -p  ${params.outDir}/${vcf.SimpleName}/pop_pca
    mkdir -p  ${params.outDir}/${vcf.SimpleName}/final

    # 1. Define command arguments
    commandArguments="--input ${vcf} \
    --output no_multi_allelic \
    --call_rate 0.5 \
    --filtered_depth 5 \
    --genotype_quality 10 \
    --minor_allele_frequency 0.01 \
    --no_snv_vqsr_check \
    --no_indel_vqsr_check \
    --remove_non_pass_snv \
    --remove_non_pass_indel \
    --replace_poor_quality_genotypes \
    --output ${vcf.SimpleName}"

    # 2. Run the custom VCF filter script
    custom_vcf_filter.py \${commandArguments}
    
    mv ${vcf.SimpleName}-filtered.vcf.gz ${vcf.SimpleName}.filtered.vcf.gz
    mv ${vcf.SimpleName}-filtered.log.gz ${vcf.SimpleName}.filtered.log.gz
    """
}

process getMetrics {
    publishDir "${params.outDir}//metrics", mode: "move", pattern: "*.{gz,log}"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path unfilteredVcf
    path filteredVcf

    output:
    path "${filteredVcf.SimpleName}.varcount.txt.gz", emit: variant_count

    script:
    """
    mkdir -p  ${params.outDir}/${filteredVcf.SimpleName}/get_metrics
    zcat ${unfilteredVcf} | grep -v '^#' | wc -l | xargs -I {} echo -e "variant count:\t{}" > ${filteredVcf.SimpleName}.varcount.txt
    zcat ${filteredVcf} | grep -v '^#' | wc -l | xargs -I {} echo -e "variant count filtered:\t{}" >> ${filteredVcf.SimpleName}.varcount.txt
    gzip ${filteredVcf.SimpleName}.varcount.txt
    """
}

process convertToPlinkFormat {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcf

    output:
    path "*.bed", emit: bed
    path "*.bim", emit: bim
    path "*.fam", emit: fam

    script:
    """
    plink2 --vcf ${vcf} --make-bed --out ${vcf.SimpleName}
    """
}

process calculateMissingness {
    storeDir "${params.outDir}/missingness"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam
    path vcf

    output:
    path "${vcf.SimpleName}.preFilter.smiss"
    path "${vcf.SimpleName}.postFilter.smiss", emit: missing

    script:
    """
    mkdir -p ${params.outDir}/${vcf.SimpleName}/missingness
    plink2 --vcf ${vcf} --missing --out ${vcf.SimpleName}.preFilter
    plink2 --bfile ${bed.SimpleName} --missing --out ${vcf.SimpleName}.postFilter
    """
}

// if less than 50 samples use --bad-freqs. Needs an alternative solution?
process createHetFile {
    storeDir "${params.outDir}/het"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam

    output:
    path "*.het"

    script:
    """
    mkdir -p  ${params.outDir}/${bed.SimpleName}/het
    plink2 --bfile ${bed.SimpleName} --het cols=hom,het,nobs,f --bad-freqs --out ${bed.SimpleName}
    """
}

// Create a txt file with samples where missingness => 50%
process findMissingSamples {
    storeDir "${params.outDir}/missingness"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path smiss

    output:
    path "${smiss.SimpleName}-filtered-samples.txt", emit: filteredSamples

    script:
    """
    mkdir -p  ${params.outDir}/missingness
    filterMissingness.py ${smiss} ${smiss.SimpleName}-filtered-samples.txt --threshold ${params.missing}
    sleep 15s
    """
}

// Filter these samples. Output -> filtered plink files
process filterMissingSamples {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path filteredSamples
    path bed
    path bim
    path fam

    output:
    path "${bed.SimpleName}.keep.bed", emit: bed
    path "${bed.SimpleName}.keep.bim", emit: bim
    path "${bed.SimpleName}.keep.fam", emit: fam

    script:
    """
    plink2 --bfile ${bed.SimpleName} --keep ${filteredSamples} --make-bed --out ${bed.SimpleName}.keep
    """
}

process findHetSamples {
    storeDir "${params.outDir}/het/"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path het

    output:
    path "*.png"
    path "*Failed.txt"
    path "*FailedSamplesOnly.txt", emit: failedHetSamples

    script:
    """
    mkdir -p  ${params.outDir}/${het.SimpleName}/het/
    heterozygosityCheck.py ${het} ./
    """
}

process filterHetSamples {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path failedHetSamples
    path bed
    path bim
    path fam

    output:
    path "${bed.SimpleName}.bed", emit: bed
    path "${bed.SimpleName}.bim", emit: bim
    path "${bed.SimpleName}.fam", emit: fam

    script:
    """
    plink2 --bfile ${bed.SimpleName}.keep \
       --remove ${failedHetSamples} \
       --make-bed \
       --out ${bed.SimpleName}
    """
}

process createMetricsFile {
    storeDir "${params.outDir}/metrics"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam

    output:
    path "metrics_matrix.tsv", emit: metrics_matrix

    script:
    """
    CombineQCFiles.py metrics_matrix.tsv
    """
}


process filterLowAltFreq {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 4

    input:
    path bed
    path bim
    path fam

    output:
    path "${bed.SimpleName}.over10.bed"
    path "${bed.SimpleName}.over10.bim"
    path "${bed.SimpleName}.over10.fam"
    
    script:
    """
    plink2 --sample-counts 'cols=homref,het,homalt,missing' --bfile ${bed.SimpleName} --out ${bed.SimpleName}
    nonalt.py -i ${bed.SimpleName}.scount -o ${bed.SimpleName}-keep.txt
    plink2 --keep ${bed.SimpleName}-keep.txt --bfile ${bed.SimpleName} --make-bed --out ${bed.SimpleName}.over10
    """
}

process filterRelated {
    publishDir "${params.outDir}/final", mode: 'copy', pattern: "*.{txt,pdf,png,log}"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 4

    input:
    path bed
    path bim
    path fam

    output:
    path "${bed.SimpleName}.RelatednessCheck.bed", emit: bed
    path "${bed.SimpleName}.RelatednessCheck.bim", emit: bim
    path "${bed.SimpleName}.RelatednessCheck.fam", emit: fam
    path "${bed.SimpleName}.related.kin0"
    path "${bed.SimpleName}.RelatednessPassedSamples.txt"

    script:
    """
    mkdir -p ${params.outDir}/${bed.SimpleName}/related
    # 1. Do relatedness check
    plink2 --bfile ${bed.SimpleName}.over10 \
        --make-king-table \
        --king-table-filter ${params.kingTableFilter} \
        --out ${bed.SimpleName}.related \
        --threads 4 \

    # 2. Create sample list of non-related samples
    find_related_samples.R --kin_file ${bed.SimpleName}.related.kin0 --target_bed ${bed} --out ${bed.SimpleName}

    # 3. Remove samples that are not on the list created above
    plink2 --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --keep ${bed.SimpleName}.RelatednessPassedSamples.txt \
        --make-bed \
        --out ${bed.SimpleName}.RelatednessCheck
    """
}


process popProject {
    storeDir "${params.outDir}/pop_pca"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '24 GB'
    cpus 4

    input:
    path bed
    path bim
    path fam
    
    output:
    path '*.png'
    path '*.pdf'
    path '*.txt'
    path bed, emit: bed
    path bim, emit: bim
    path fam, emit: fam
    path 'PopAssignResults.txt', emit: pops
    
    script:
    """
    mkdir -p  ${params.outDir}/${bed.SimpleName}/pop_pca
    project_samples_to_superpop.R --ref_bed ${params.refPath} \
        --target_bed ${bed} --ref_pop ${params.refPop}
    """
}

process popSplit {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '24 GB'
    cpus 4

    input:
    path bed
    path bim
    path fam
    path assign
    
    output:
    path bed, emit: bed
    path bim, emit: bim
    path fam, emit: fam
    path 'PopAssignResults.txt', emit: pops
    path 'pops/*', emit: pop_dirs
    
    script:
    """
    pop.py -i ${assign} -o assigned.txt
    """
}

process targetPCA {
  storeDir "${params.outDir}/target/" 
  

  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path bed
  path bim
  path fam
  
  output:
  path "${bed.SimpleName}.toImputation.bed", emit: bed
  path "${bed.SimpleName}.toImputation.bim", emit: bim
  path "${bed.SimpleName}.toImputation.fam", emit: fam
  path "*.txt"
  path "*.pdf"
  path "*.png"
  path "*.log"

  
  script:
  """
  mkdir -p  ${params.outDir}/target
  # 1. Do PCA and find outliers
  target_pca.R --target_bed ${bed} --outlier_threshold ${params.populationOutlierThreshold} --out ${bed.SimpleName}

  # 2. Remove outlier samples
  plink2 --bed ${bed} \
    --bim ${bim} \
    --fam ${fam} \
    --output-chr 26 \
    --keep ${bed.SimpleName}.SamplesToInclude.txt \
    --make-bed \
    --threads 4 \
    --out ${bed.SimpleName}.toImputation
  """
}


process finalSNPandGenotypeQC {
  publishDir "${params.outDir}/final_snp_qc/", mode: 'move', pattern: "*.{txt,pdf,png,bed,bim,fam,gz,log}"
  

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path bed
  path bim
  path fam
  
  output:
  path "${bed.SimpleName}.bed"
  path "${bed.SimpleName}.bim"
  path "${bed.SimpleName}.fam"
  
  script:
  """
  mkdir -p  ${params.outDir}/final_snp_qc/
  plink2 --bed ${bed} \
    --bim ${bim} \
    --fam ${fam} \
    --maf ${params.maf} \
    --geno ${params.geno} \
    --mind ${params.mind} \
    --hwe ${params.hwe} \
    --autosome \
    --make-bed \
    --out chrAll_QC \
    --output-chr 26 \
    --not-chr 0 25-26 \
    --set-all-var-ids @:#[b38]\\\$r,\\\$a \
    --new-id-max-allele-len 10 truncate \
    --threads 4
  """
}

process PCA {
  publishDir "${params.outDir}/before_target/", mode: 'move', pattern: "*.{txt,pdf,png}"
  

  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path bed
  path bim
  path fam
  
  output:
  path '*.png'
  path '*.pdf'
  path '*.txt'
  
  script:
  """
  final_pca.R --target_bed ${bed} 
  """
}

process finalPCA {
  publishDir "${params.outDir}/final/", mode: 'move', pattern: "*.{txt,pdf,png}"
  

  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path bed
  path bim
  path fam
  
  output:
  path '*.png'
  path '*.pdf'
  path '*.txt'
  
  script:
  """
  mkdir -p  ${params.outDir}/${bed.SimpleName}/final/
  final_pca.R --target_bed ${bed} 
  """
}


process pop_filter {
errorStrategy 'ignore'
publishDir "${params.outDir}/${bed.SimpleName}/", mode: 'copy'
time '6h'
memory '8 GB'
cpus 4

input:
path bed
path bim
path fam
path pop_dir

output:
path "target/${pop_dir}/", emit: pop_dir
// path "${pop_dir.SimpleName}.toImputation.bed", emit: bed
// path "${pop_dir.SimpleName}.toImputation.bim", emit: bim
// path "${pop_dir.SimpleName}.toImputation.fam", emit: fam

script:
"""
mkdir -p target/${pop_dir}
plink2 --bfile ${bed.SimpleName} --keep ${pop_dir}/${pop_dir.SimpleName}.txt --make-bed --out ${pop_dir.SimpleName}/${bed.SimpleName}_${pop_dir.SimpleName}
if [ \$(cat ${pop_dir}/${pop_dir.SimpleName}.txt | wc -l) -gt 29 ]
then
    plink2 --bfile ${bed.SimpleName} --keep ${pop_dir}/${pop_dir.SimpleName}.txt --make-bed --out ${pop_dir.SimpleName}/${bed.SimpleName}_${pop_dir.SimpleName}

    # 1. Do PCA and find outliers
    target_pca.R --target_bed ${pop_dir.SimpleName}/${bed.SimpleName}_${pop_dir.SimpleName}.bed --outlier_threshold ${params.populationOutlierThreshold} --out target/${pop_dir}/${pop_dir.SimpleName}

    # 2. Remove outlier samples
    plink2 --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --output-chr 26 \
        --keep target/${pop_dir}/${pop_dir.SimpleName}.SamplesToInclude.txt \
        --make-bed \
        --threads 4 \
        --out target/${pop_dir}/${pop_dir.SimpleName}.toImputation
    final_pca.R --target_bed target/${pop_dir}/${pop_dir.SimpleName}.toImputation.bed
    mv ./*.pdf ./*.png target/${pop_dir}

else
    mv ${pop_dir.SimpleName}/${bed.SimpleName}_${pop_dir.SimpleName}.bed target/${pop_dir}/${pop_dir.SimpleName}.toImputation.bed
    mv ${pop_dir.SimpleName}/${bed.SimpleName}_${pop_dir.SimpleName}.bim target/${pop_dir}/${pop_dir.SimpleName}.toImputation.bim
    mv ${pop_dir.SimpleName}/${bed.SimpleName}_${pop_dir.SimpleName}.fam target/${pop_dir}/${pop_dir.SimpleName}.toImputation.fam
fi
"""
}

process pop_merge {
  publishDir "${params.outDir}/${proj_bed.SimpleName}/final", mode: "move", pattern: "*.{txt,pdf,png,bed,bim,fam,gz,log}"
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path pop_dirs
  path proj_bed
  
  output:
  path "merged.vcf.gz"
  path "merged.bed"
  path "merged.bim"
  path "merged.fam"
  path "*.pdf"
  path "*.txt"
  path "*.png"
  
  script:
  """
  for pop in ${pop_dirs}
  do
    plink2 --bfile \$pop/\$pop.toImputation --export vcf bgz --out \$pop
    tabix \$pop.vcf.gz
  done
  bcftools merge ./*.vcf.gz -Oz -o merged.vcf.gz
  plink2 --vcf merged.vcf.gz --make-bed --out merged
  final_pca.R --target_bed merged.bed
  """
}

process beagle {
    storeDir "${params.outDir}/${vcf.SimpleName}/final"
    errorStrategy 'retry'
    maxRetries 2

    time '8h'
    memory '20 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam

    output:
    path "${bed.SimpleName}.beagle.bed", emit: bed
    path "${bed.SimpleName}.beagle.bim", emit: bim
    path "${bed.SimpleName}.beagle.fam", emit: fam
    path "${vcfFile.SimpleName}.beagle.log"

    script:
    """
    plink2 --bfile ${bed.SimpleName} --export vcf bgz --out ${bed.SimpleName}
    java -Xmx18g -jar /usr/bin/beagle.jar \
    gtgl=${bed.SimpleName}.vcf.gz \
    out=${bed.SimpleName}.beagle \
    map=${params.beagle_map} \
    gprobs=true
    plink2 --vcf ${bed.SimpleName}.vcf.gz --make-bed --out ${bed.SimpleName}.beagle
    """
}

//   for dir in ${pop_dirs}
//   do
//     sort -k1,1 -k2,2n \$dir/\$dir.toImputation.bed > \$dir/\$dir.toImputation.sorted.bed
//     bedtools merge 
//     plink2 --bfile \$dir/\$dir.toImputation --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 1000 --make-bed --out \$dir/\$dir.toImputation 
//   done

//   bedtools merge
//   for dir in ${pop_dirs}; do plink2 --bed \$dir/\$dir.toImputation.bed --bim \$dir/\$dir.toImputation.bim --fam \$dir/\$dir.toImputation.fam --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 1000 --make-bed --out \$dir/\$dir.toImputation; done
//   for dir in ${pop_dirs}; do echo "\$dir/\$dir.toImputation.bed \$dir/\$dir.toImputation.bim \$dir/\$dir.toImputation.fam" >> bfile_list.txt; done
//   plink2 --pmerge-list bfile_list.txt --out pop_pc_filtered --make-bed

workflow {
    vcfs = Channel.fromPath("${params.inputDir}/*.vcf.gz")
    fixGTAnnot(vcfs)
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