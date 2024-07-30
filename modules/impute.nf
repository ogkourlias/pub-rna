process sort_bed{
    maxRetries 2
    //errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    time '12h'
    memory '16 GB'
    cpus 1

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    path "${params.cohort_name}.sorted.bed", emit: bed
    path "${params.cohort_name}.sorted.bim", emit: bim
    path "${params.cohort_name}.sorted.fam", emit: fam

    script:
    """
    plink2 --bed ${bed} --bim ${bim} --fam ${fam} --make-bed --output-chr MT --out ${params.cohort_name}.sorted
    """
}

process split_by_chr{
    maxRetries 2
    time '12h'
    memory '16 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam
    val chr

    output:
    path "chr${chr}.bed", emit: bed
    path "chr${chr}.bim", emit: bim
    path "chr${chr}.fam", emit: fam
    val chr, emit: chr

    script:
    """
    plink2 --bfile ${bed.BaseName} --chr ${chr} --make-bed --out chr${chr}
    """
}

process harmonize_hg38{
    maxRetries 2
    time '12h'
    memory '32 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam
    path ref_panel
    path ref_panel_tbi
    val chr


    output:
    path "${bed.SimpleName}.harmonised.bed", emit: bed
    path "${bed.SimpleName}.harmonised.bim", emit: bim
    path "${bed.SimpleName}.harmonised.fam", emit: fam
    val chr, emit: chr

    script:
    """
    java -Xmx25g -jar /usr/bin/GenotypeHarmonizer.jar\
     --input ${bed.baseName}\
     --inputType PLINK_BED\
     --ref ${ref_panel}\
     --refType VCF\
     --update-id\
     --output ${bed.SimpleName}.harmonised
    """
}

process plink_to_vcf{
    maxRetries 2
    time '12h'
    memory '16 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam
    val chr

    output:
    path "${bed.SimpleName}.vcf", emit: vcf
    val chr, emit: chr

    script:
    """
    plink2 --bfile ${bed.BaseName} --recode vcf-iid --chr 1-22 --out ${bed.SimpleName}
    """
}

process vcf_fixref_hg38{
    maxRetries 2
    time '12h'
    memory '16 GB'
    cpus 1

    input:
    path input_vcf
    path ref_fa
    path ref_harm_vcf
    path ref_harm_vcf_tbi
    val chr

    output:
    path "${input_vcf.SimpleName}.fixref_hg38.vcf.gz", emit: vcf
    val chr, emit: chr

    script:
    """
    bgzip -c ${input_vcf} > ${input_vcf}.gz
    bcftools index ${input_vcf}.gz

    bcftools +fixref ${input_vcf}.gz -- -f ${ref_fa} -i ${ref_harm_vcf} | \
    bcftools norm --check-ref x -f ${ref_fa} -Oz -o ${input_vcf.SimpleName}.fixref_hg38.vcf.gz
    """
}

process eagle_prephasing{
    maxRetries 2
    time '12h'
    memory '32 GB'
    cpus 4

    input:
    path vcf
    path genetic_map
    path phasing_dir
    val chr

    output:
    path "${vcf.SimpleName}.phased.vcf.gz", emit: vcf
    val chr, emit: chr

    script:
    """
    bcftools index ${vcf}
    eagle --vcfTarget=${vcf} \
    --vcfRef=${phasing_dir}/chr${chr}.bcf \
    --geneticMapFile=${genetic_map} \
    --chrom=${chr} \
    --outPrefix=${vcf.SimpleName}.phased \
    --numThreads=4
    """
}

process minimac_imputation{
    maxRetries 2
    time '12h'
    memory '32 GB'
    cpus 4

    input:
    path vcf
    path imputation_dir
    val chr

    output:
    path "${vcf.SimpleName}.dose.vcf.gz", emit: vcf

    script:
    """
    minimac4 --refHaps ${imputation_dir}/chr${chr}.m3vcf.gz \
    --rsid \
    --haps ${vcf} \
    --prefix ${vcf.SimpleName} \
    --format GT,DS,GP \
    --noPhoneHome \
    --minRatio 0.0000001 \
    --ChunkLengthMb 300 \
    --cpus 4
    """
}

process filter_r2_maf{
    publishDir "${params.out_dir}/postimpute/", mode: 'move', pattern: "*.vcf.gz", overwrite: true
    maxRetries 2
    time '12h'
    memory '16 GB'
    cpus 1

    input:
    path vcf

    output:
    path "${vcf.SimpleName}.filtered.dose.vcf.gz"
    path "${vcf.SimpleName}.filtered.dose.vcf.gz.tbi"

    script:
    """
    bcftools filter -i 'INFO/R2>=0.3 && INFO/MAF>=0.01' ${vcf} -Oz -o ${vcf.SimpleName}.filtered.dose.vcf.gz
    tabix ${vcf.SimpleName}.filtered.dose.vcf.gz
    """
}

workflow {
    sort_bed(tuple("${params.in_dir}/*.bed", "${params.in_dir}/*.bim", "${params.in_dir}/*.fam"))
    chrs = Channel.from( 1..22 )
    split_by_chr(sort_bed.out.bed, sort_bed.out.bim, sort_bed.out.fam, chrs)
    harmonize_hg38(split_by_chr.out.bed, split_by_chr.out.bim, split_by_chr.out.fam, params.ref_harm_vcf, "${params.ref_harm_vcf}.tbi", split_by_chr.out.chr)
    plink_to_vcf(harmonize_hg38.out.bed, harmonize_hg38.out.bim, harmonize_hg38.out.fam, harmonize_hg38.out.chr)
    vcf_fixref_hg38(plink_to_vcf.out.vcf, params.ref_fa, params.ref_harm_vcf, "${params.ref_harm_vcf}.tbi", plink_to_vcf.out.chr)
    eagle_prephasing(vcf_fixref_hg38.out.vcf, params.phasing_map, params.phasing_dir, vcf_fixref_hg38.out.chr)
    minimac_imputation(eagle_prephasing.out.vcf, params.imputation_dir, eagle_prephasing.out.chr)
    filter_r2_maf(minimac_imputation.out.vcf)
}
