def helpMessage() {
    log.info"""
    =======================================================
     eqtlgenimpute v${workflow.manifest.version}
    =======================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run eQTLGenImpute.nf \
    --qcdata CohortName_hg37_genotyped \
    --cohort_name CohortName_hg38_imputed \
    --genome_build GRCh37 \
    --outdir CohortName \
    --target_ref Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --ref_panel_hg38 30x-GRCh38_NoSamplesSorted \
    --eagle_genetic_map genetic_map_hg38_withX.txt.gz \
    --eagle_phasing_reference /phasing_reference/ \
    --minimac_imputation_reference /imputation_reference/ \
    -profile slurm,singularity \
    -resume

    Mandatory arguments:
      --qcdata                          Path to the folder with input unimputed plink files.
      --cohort_name                     Prefix for the output files.

      --outdir                          The output directory where the results will be saved.
      --target_ref                      Reference genome fasta file for the target genome assembly (e.g. GRCh38).
      --ref_panel_hg38                  Reference panel used for strand fixing and GenotypeHarmonizer after LiftOver (GRCh38).
      --eagle_genetic_map               Eagle genetic map file.
      --eagle_phasing_reference         Phasing reference panel for Eagle (1000 Genomes 30x WGS high coverage).
      --minimac_imputation_reference    Imputation reference panel for Minimac4 in M3VCF format (1000 Genomes 30x WGS high coverage).

    Optional arguments:
      --chain_file                      Chain file to translate genomic coordinates from the source assembly to target assembly (e.g. hg19 --> hg38). hg19-->hg38 & hg18-->hg38 works by default.
      --cohort_build                    The genome build to which the cohort is mapped: hg18, GRCh36, hg19, GRCh37, hg38 or GRCh38 (default is hg37, setting this to hg38 skips crossmapping)

    """.stripIndent()
}

// Define set of accepted genome builds:
def genome_builds_accepted = ['hg18', 'GRCh36', 'hg19', 'GRCh37', 'hg38', 'GRCh38']

// Define input channels
Channel
    .fromFilePairs("${params.qcdata}/outputfolder_gen/gen_data_QCd/*.{bed,bim,fam}", size: -1)
    .ifEmpty {exit 1, "Input genotype files not found!"}
    .map { it.flatten() }
    .set { bfile_ch }

Channel
    .fromPath("${params.qcdata}/outputfolder_exp/exp_data_QCd/exp_data_preprocessed.txt")
    .ifEmpty { exit 1, "Expression matrix not found: ${params.qcdata}/outputfolder_exp/exp_data_QCd/exp_data_preprocessed.txt" }
    .set { exp_mat_ch }

Channel
    .fromPath(params.ref_panel_hg38)
    .map { ref -> [file("${ref}.vcf.gz"), file("${ref}.vcf.gz.tbi")] }
    .into { ref_panel_harmonise_genotypes_hg38; ref_panel_fixref_genotypes_hg38; ref_panel_maf }

Channel
    .fromPath( "${params.eagle_phasing_reference}*" )
    .ifEmpty { exit 1, "Eagle phasing reference not found: ${params.eagle_phasing_reference}" }
    .set { phasing_ref_ch }

Channel
    .fromPath( "${params.minimac_imputation_reference}*" )
    .ifEmpty { exit 1, "Minimac4 imputation reference not found: ${params.minimac_imputation_reference}" }
    .set { imputation_ref_ch }

Channel
    .fromPath(params.eagle_genetic_map)
    .ifEmpty { exit 1, "Eagle genetic map file not found: ${params.eagle_genetic_map}" }
    .set { genetic_map_ch }

Channel
    .fromPath(params.target_ref)
    .ifEmpty { exit 1, "Target reference genome file not found: ${params.target_ref}" }
    .into { target_ref_ch; target_ref_ch2 }

if ((params.genome_build in genome_builds_accepted) == false) {
  exit 1, "[Pipeline error] Genome build $params.genome_build not in accepted genome builds: $genome_builds_accepted \n"
}

params.chain_file="$baseDir/data/GRCh37_to_GRCh38.chain"
if (params.genome_build in ["hg18", "GRCh36"]) {
    chain_file="$baseDir/data/hg18ToHg38.over.chain"
} else {
    chain_file=params.chain_file
}

skip_crossmap = params.genome_build in ["hg38", "GRCh38"]

Channel
    .fromPath(chain_file)
    .ifEmpty { exit 1, "CrossMap.py chain file not found: ${chain_file}" }
    .set { chain_file_ch }

// Header log info
log.info """=======================================================
eqtlgenimpute v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']            = 'eqtlgenimpute'
summary['Pipeline Version']         = workflow.manifest.version
summary['Path to QCd input']        = params.qcdata
summary['Cohort build']             = params.genome_build
summary['Chain file']             = chain_file
summary['Skip crossmap']             = skip_crossmap
summary['Harmonisation ref panel hg38']  = params.ref_panel_hg38
summary['Target reference genome hg38'] = params.target_ref
summary['CrossMap chain file']      = params.chain_file
summary['Eagle genetic map']        = params.eagle_genetic_map
summary['Eagle reference panel']    = params.eagle_phasing_reference
summary['Minimac4 reference panel'] = params.minimac_imputation_reference
summary['Max Memory']               = params.max_memory
summary['Max CPUs']                 = params.max_cpus
summary['Max Time']                 = params.max_time
summary['Cohort name']              = params.cohort_name
summary['Output dir']               = params.outdir
summary['Working dir']              = workflow.workDir
summary['Container Engine']         = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']             = "$HOME"
summary['Current user']             = "$USER"
summary['Current path']             = "$PWD"
summary['Working dir']              = workflow.workDir
summary['Script dir']               = workflow.projectDir
summary['Config Profile']           = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="

process crossmap{

    input:
    set val(study_name), file(study_name_bed), file(study_name_bim), file(study_name_fam) from ( skip_crossmap ? Channel.empty() : bfile_ch)
    file chain_file from chain_file_ch.collect()

    output:
    set val(study_name), file("crossmapped_plink.bed"), file("crossmapped_plink.bim"), file("crossmapped_plink.fam") into crossmapped

    when:
    skip_crossmap == false

    shell:
    //Converts BIM to BED and converts the BED file via CrossMap.
    //Finds excluded SNPs and removes them from the original plink file.
    //Then replaces the BIM with CrossMap's output.
    """
    awk '{print \$1,\$4,\$4+1,\$2,\$5,\$6,\$2 "___" \$5 "___" \$6}' ${study_name}.bim > crossmap_input.bed
    CrossMap.py bed ${chain_file} crossmap_input.bed crossmap_output.bed
    awk '{print \$7}' crossmap_input.bed | sort > input_ids.txt
    awk '{print \$7}' crossmap_output.bed | sort > output_ids.txt
    comm -23 input_ids.txt output_ids.txt | awk '{split(\$0,a,"___"); print a[1]}' > excluded_ids.txt
    plink2 --bfile ${study_name} --exclude excluded_ids.txt --make-bed --output-chr MT --out crossmapped_plink --keep-allele-order
    awk -F'\t' 'BEGIN {OFS=FS} {print \$1,\$4,0,\$2,\$5,\$6}' crossmap_output.bed > crossmapped_plink.bim
    """
}

process sort_bed{

    input:
    set val(study_name), file(study_name_bed), file(study_name_bim), file(study_name_fam) from ( skip_crossmap ? bfile_ch : crossmapped)

    output:
    tuple file("sorted.bed"), file("sorted.bim"), file("sorted.fam") into sorted_genotypes_hg38_ch

    script:
    """
    plink2 --bfile ${study_name_bed.simpleName} --make-bed --output-chr MT --out sorted
    """
}

process split_by_chr{

    publishDir "${params.outdir}/preimpute/", mode: 'copy',
        saveAs: {filename -> if (filename.indexOf(".vcf.gz") > 0) filename else null }

    input:
    tuple file(study_name_bed), file(study_name_bim), file(study_name_fam) from sorted_genotypes_hg38_ch
    each chr from Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

    output:
    tuple val(chr), path("chr_${chr}_sorted.bed"), path("chr_${chr}_sorted.bim"), path("chr_${chr}_sorted.fam") into sorted_split_genotypes_hg38_ch

    script:
    """
    plink2 --bfile ${study_name_bed.baseName} --chr ${chr} --make-bed --out chr_${chr}_sorted
    """
}

process harmonize_hg38{

    input:
    tuple val(chr), path(study_name_bed), path(study_name_bim), path(study_name_fam) from sorted_split_genotypes_hg38_ch
    tuple file(vcf_file), file(vcf_file_index) from ref_panel_harmonise_genotypes_hg38.collect()

    output:
    tuple val(chr), file("harmonised.bed"), file("harmonised.bim"), file("harmonised.fam") into harmonised_split_genotypes_hg38_ch

    script:
    """
    java -Xmx25g -jar /usr/bin/GenotypeHarmonizer.jar\
     --input ${study_name_bed.baseName}\
     --inputType PLINK_BED\
     --ref ${vcf_file.simpleName}\
     --refType VCF\
     --update-id\
     --output harmonised
    """
}

process plink_to_vcf{

    input:
    set val(chr), file(study_name_bed), file(study_name_bim), file(study_name_fam) from harmonised_split_genotypes_hg38_ch

    output:
    tuple val(chr), file("chr${chr}_harmonised_hg38.vcf") into harmonized_split_hg38_vcf_ch

    script:
    """
    plink2 --bfile ${study_name_bed.simpleName} --recode vcf-iid --chr 1-22 --out chr${chr}_harmonised_hg38
    """
}

process vcf_fixref_hg38{

    input:
    tuple val(chr), file(input_vcf) from harmonized_split_hg38_vcf_ch
    file fasta from target_ref_ch2.collect()
    set file(vcf_file), file(vcf_file_index) from ref_panel_fixref_genotypes_hg38.collect()

    output:
    tuple val(chr), file("chr${chr}_fixref_hg38.vcf.gz") into fixed_to_filter_split

    script:
    """
    bgzip -c ${input_vcf} > ${input_vcf}.gz
    bcftools index ${input_vcf}.gz

    bcftools +fixref ${input_vcf}.gz -- -f ${fasta} -i ${vcf_file} | \
    bcftools norm --check-ref x -f ${fasta} -Oz -o chr${chr}_fixref_hg38.vcf.gz
    """
}

// process filter_preimpute_vcf{

//     publishDir "${params.outdir}/preimpute/", mode: 'copy',
//         saveAs: {filename -> if (filename == "filtered.vcf.gz") "${params.cohort_name}_preimpute.vcf.gz" else null }

//     input:
//     tuple val(chr), file(input_vcf) from fixed_to_filter_split

//     output:
//     tuple val(chr), file("filtered.vcf.gz"), file("filtered.vcf.gz.csi") into split_vcf_input
//     file("filtered.vcf.gz") into missingness_input_split

//     script:
//     """
//     #Index
//     bcftools index ${input_vcf}

//     #Add tags
//     bcftools +fill-tags ${input_vcf} -Oz -o tagged.vcf.gz

//     #Filter rare and non-HWE variants and those with abnormal alleles and duplicates
//     bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' tagged.vcf.gz |\
//      bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
//      bcftools filter -e "ALT='.'" |\
//      bcftools norm -d all |\
//      bcftools norm -m+any |\
//      bcftools view -m2 -M2 -Oz -o filtered.vcf.gz

//     #Index the output file
//     bcftools index filtered.vcf.gz
//     """
// }

// process calculate_missingness{

//     publishDir "${params.outdir}/preimpute/", mode: 'copy',
//         saveAs: {filename -> if (filename == "genotypes.imiss") "${params.cohort_name}.imiss" else null }

//     input:
//     file 'filtered_??.vcf.gz' from missingness_input_split.collect()

//     output:
//     file "genotypes.imiss" into missing_individuals

//     script:
//     """
//     bcftools concat filtered_??.vcf.gz -Oz > concat.vcf.gz
//     bcftools index -f concat.vcf.gz
//     vcftools --gzvcf concat.vcf.gz --missing-indv --out genotypes
//     """
// }

process eagle_prephasing{

    input:
    tuple val(chromosome), file(vcf) from fixed_to_filter_split
    file genetic_map from genetic_map_ch.collect()
    file phasing_reference from phasing_ref_ch.collect()

    output:
    tuple val(chromosome), file("chr${chromosome}.phased.vcf.gz") into phased_vcf_cf

    script:
    """
    bcftools index ${vcf}
    eagle --vcfTarget=${vcf} \
    --vcfRef=chr${chromosome}.bcf \
    --geneticMapFile=${genetic_map} \
    --chrom=${chromosome} \
    --outPrefix=chr${chromosome}.phased \
    --numThreads=8
    """
}

// process filter_samples{

//     input:
//     set val(chromosome), file(vcf), file(exp_mat) from phased_vcf_cf2

//     output:
//     tuple val(chromosome), file("chr${chromosome}.phased.samplesfiltered.vcf.gz") into phased_vcf_samples_filtered_cf

//     script:
//     """
//     awk '{print \$1}' ${exp_mat} | awk '(NR>1)' > sample_filter.txt
//     bcftools view -S sample_filter.txt --force-samples ${vcf} -Oz -o chr${chromosome}.phased.samplesfiltered.vcf.gz
//     """
// }

process minimac_imputation{

    input:
    set val(chromosome), file(vcf) from phased_vcf_cf
    file imputation_reference from imputation_ref_ch.collect()

    output:
    tuple val(chromosome), file("chr${chromosome}.dose.vcf.gz") into imputed_vcf_cf

    publishDir "${params.outdir}/postimpute/", mode: 'copy', pattern: "*.vcf.gz", overwrite: true
    script:
    """
    minimac4 --refHaps chr${chromosome}.m3vcf.gz \
    --rsid \
    --haps ${vcf} \
    --prefix chr${chromosome} \
    --format GT,DS,GP \
    --noPhoneHome \
    --ChunkLengthMb 50 \
    --minRatio 0.01
    """
}

// process filter_maf{

//     publishDir "${params.outdir}/postimpute/", mode: 'copy', pattern: "*.filtered.vcf.gz", overwrite: true

//     input:
//     set chromosome, file(vcf) from imputed_vcf_cf

//     output:
//     tuple chromosome, file("chr${chromosome}.filtered.vcf.gz") into imputed_vcf_filtered_cf

//     script:
//     """
//     bcftools +fill-tags ${vcf} -Ou -- -t AF,MAF \
//     | bcftools filter -i 'INFO/MAF[0] > 0.01' -Oz -o chr${chromosome}.filtered.vcf.gz
//     """
// }

// process extract_maf_ref{

//     input:
//     tuple file(vcf_file), file(vcf_file_index) from ref_panel_maf.collect()
    
//     output:
//     file("ref_allele_frequencies.txt") into ref_af

//     script:
//     """
//     bcftools \
//     query -f '%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%AF\\t%AF_EUR\\n' \
//     ${vcf_file} > ref_allele_frequencies.txt
//     """
// }

// process extract_maf_target{

//     input:
//     set val(chromosome), file(vcf) from imputed_vcf_filtered_cf

//     output:
//     file("*_AF.txt") into target_af

//     script:
//     """
//     bcftools \
//     query -f '%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%AF\\t%IMPUTED\\t%TYPED\\t%R2\\n' \
//     ${vcf} > ${chromosome}_AF.txt
//     """
// }

// maf_check_ch = target_af.collect().combine(ref_af)

// process compare_MAF{

//     container = "quay.io/eqtlgen/popassign:v0.6"

//     publishDir "${params.outdir}/", mode: 'copy', pattern: "*.html", overwrite: true

//     input:
//     path files from maf_check_ch
//     val workflow_version from workflow.manifest.version

//     output:
//     path "Report_ImputationQc.html" into report_output_ch

//     script:
//     """
//     # Make report
//     cp -L $baseDir/bin/Report_template.Rmd notebook.Rmd

//     R -e 'library(rmarkdown);rmarkdown::render("notebook.Rmd", "html_document", 
//     output_file = "Report_ImputationQc.html", params=list(workflow_version = "${workflow_version}"))'
//     """

// }

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
