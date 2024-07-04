nextflow.enable.dsl=2

process createCovariateMatrix {
  publishDir "${params.outDir}/1_create_covariate_matrix/", mode: 'copy'
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '16 GB'
  cpus 1

  input:
  path inputDir
  
  output:
  path "${params.cohortName}_covariates.txt.gz", emit: matrix
  path "failed_samples.txt"
  
  script:
  """
  # 1. Create covariate matrix
  create_covariate_matrix.py --input_dir ${params.inputDir} \
  --cohort_name ${params.cohortName}

  # 2. Gzip the matrix
  gzip ${params.cohortName}_covariates.txt

  # 3. Merge covariate and WGS PC dataframe
  # merge_covariates_with_wgs.py ${params.cohortName}_covariates.txt.gz ${params.wgsPCFile} ${params.gteFile} ${params.cohortName}
  
  # 4. Gzip output matrix
  #gzip ${params.cohortName}_covariates_wgs_pcs.txt
  """
}

process createGeneCountsMatrix {
  publishDir "${params.outDir}/3_create_gene_counts_matrix/", mode: 'copy'
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '16 GB'
  cpus 1

  input:
  path inputDir
  
  output:
  path "${params.cohortName}_gene_counts.txt.gz"
  
  script:
  """
  # 1. Create gene count matrix
  python3 ${baseDir}/bin/create_gene_counts_matrix.py --input_dir ${params.inputDir} \
  --cohort_name ${params.cohortName}

  # 2. Gzip the matrix
  gzip ${params.cohortName}_gene_counts.txt
  """
}

process geneCountsTMMNormalization {
  publishDir "${params.outDir}/3_gene_counts_tmm_normalization/", mode: 'copy'
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '16 GB'
  cpus 1

  input:
  path geneCountMatrix
  
  output:
  path "${params.cohortName}_gene_counts-TMM.txt.gz", emit: tmmNormalizedMatrix
  path "${params.cohortName}_gene_counts-CPM.txt.gz"
  
  script:
  """
  # 2. Run TMM normalization
  calculate_TMM.R -c ${geneCountMatrix} \
  -t ./${params.cohortName}_gene_counts-TMM.txt \
  -o ./${params.cohortName}_gene_counts-CPM.txt

  # 3. Gzip output
  gzip *.txt
  """
}

process sampleSelection {
    publishDir "${params.outDir}/4_sample_selection/", mode: 'copy'
    maxRetries 2
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    time '12h'
    memory '16 GB'
    cpus 1

    input:
    path geneCountMatrix
    path covariateMatrix

    output:
    path "${params.cohortName}_gene_counts-TMM.samplesFiltered.txt.gz", emit: geneCountMatrix
    path "${params.cohortName}_covariates_wgs_pcs.samplesFiltered.txt.gz", emit: covariateMatrix

    script:
    """
    # Filter gene count matrix
    python3 ${baseDir}/bin/filter_samples.py ${geneCountMatrix} ${params.sampleList}

    # Filter covariate matrix
    python3 ${baseDir}/bin/filter_samples.py ${covariateMatrix} ${params.sampleList}
    """
}

process normalizeCovariateMatrix {
  publishDir "${params.outDir}/2_normalize_covariate_matrix/", mode: 'copy'
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '16 GB'
  cpus 1

  input:
  path covariateMatrix
  
  output:
  path "${params.cohortName}_covariates_wgs_pcs_normalized.txt.gz"
  
  script:
  """
  # 1. Normalize covariate matrix
  python3 ${baseDir}/bin/normalize_covariate_matrix.py --matrix_path ${covariateMatrix} \
  --cohort_name ${params.cohortName}
  
  # 2. Gzip the matrix
  gzip ${params.cohortName}_covariates_wgs_pcs_normalized.txt
  """
}

process pcaOutlierFiltering {
  publishDir "${params.outDir}/5_gene_counts_pca_outlier_identification/", mode: 'copy'
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '32 GB'
  cpus 1

  input:
  path geneCountMatrix
  
  output:
  path 'pc1_2.txt'
  path 'pc1_2_no_outliers.txt'
  path 'passed_samples.txt'
  path 'outliers.png'
  path "${params.cohortName}_gene_counts-TMM.ProbesWithZeroVarianceRemoved.txt.gz", emit: geneCountMatrix
  
  shell:
  '''
  # 2. Do PCA
  java -Xms20g -jar !{baseDir}/bin/eqtl-mapping-pipeline.jar	\
      --mode normalize \
      --in !{geneCountMatrix} \
      --out $(pwd) \
      --logtransform \
      --qqnorm \
      --adjustPCA \
      --maxnrpcaremoved 0 \
      --stepsizepcaremoval 0

  # 3. Write first PCs out to file  
  zcat *.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors.txt.gz \
  | awk 'BEGIN {OFS="\\t"} NR == 1 {print "Sample","Comp1","Comp2"} NR \
  > 1 {print $1,$2,$3}' > pc1_2.txt

  # 4. Identify outliers based on z-score
  identify_pca_outliers.py --pc_file !{geneCountMatrix} \
  --z_score_threshold !{params.pcaOutlierZScoreThreshold}

  # 5. Remove outliers
  java -Xms20g -jar !{baseDir}/bin/eqtl-mapping-pipeline.jar \
    --mode normalize \
    --sampleInclude passed_samples.txt \
    --out ./ \
    --in !{geneCountMatrix}
  '''
}

process quantileNormalization {
  publishDir "${params.outDir}/6_quantile_normalization/", mode: 'copy'
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '16 GB'
  cpus 1

  input:
  path geneCountMatrix
  
  output:
  path '*.quantile_transform.txt.gz'

  script:
  """
  python3 ${baseDir}/bin/quantile_transform.py ${geneCountMatrix}
  """
}

process geneForceNormal {
  publishDir "${params.outDir}/7_gene_force_normal/", mode: 'copy'
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '16 GB'
  cpus 1

  input:
  path geneCountMatrix
  
  output:
  path '*.geneForceNormal.txt.gz'

  script:
  """
  force_normal_distribution.py ${geneCountMatrix} 1
  """
}

process geneCountsCovariateCorrection {
  publishDir "${params.outDir}/8_covariate_correction/", mode: 'copy'
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '16 GB'
  cpus 1

  input:
  path geneCountMatrix
  path covariateMatrix
  
  
  output:
  path "*.txt"
  path "*.txt.gz"
  path "${params.cohortName}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz", emit: geneCountMatrix

  script:
  """
  # 1. Covariate PCA 
  python3 ${baseDir}/bin/pca.py ${covariateMatrix} covariates

  # 2. Extract and transpose PCs that explain variance greater than a certain threshold
  python3 ${baseDir}/bin/extract_and_transpose_pcs.py --covariate_pcs covariates_PCs.txt \
    --covariates_explained_variance explainedVariance.txt \
    --explained_variance_threshold 0.99 

  # 3. Gzip the covariate matrix
  gzip covariates-pca_PCs-filtered-transpose.txt

  # 5. Regress covariates
  java -jar ${baseDir}/bin/2023-12-08-Regression.jar \
    ${geneCountMatrix} \
    covariates-pca_PCs-filtered-transpose.txt.gz \
    ./${params.cohortName}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz
  """
}

process residualsPCA {
  publishDir "${params.outDir}/9_residuals_pca/", mode: 'copy'
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '16 GB'
  cpus 1

  input:
  path geneCountMatrix

  output:
  path "*.png"
  path "*.txt"
  path "residuals_PCs.txt", emit: residualPCs

  script:
  """
  pca.py ${geneCountMatrix} residuals
  """
}

process covariateAndResidualsCorrection {
  publishDir "${params.outDir}/10_final_covariate_correction/", mode: 'copy'
  maxRetries 2
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  time '12h'
  memory '16 GB'
  cpus 1
  
  input:
  path residualsPCs
  path covariateMatrix
  path forceNormalGeneCountMatrix

  output:
  path "*.gz"

  script:
  """
  # 1. Merge covariate matrix with residuals PCs
  merge_matrices.py ${residualsPCs} ${covariateMatrix} ${params.cohortName}

  # 2. Normalize the matrix
  normalize_covariate_matrix.py --matrix_path ${params.cohortName}_covariates_wgs_pcs_res_pcs.txt --cohort_name ${params.cohortName}

  # 3. Gzip matrix
  gzip ${params.cohortName}_covariates_wgs_pcs_normalized.txt

  # 4. Do a PCA over the merged matrix
  pca.py ${params.cohortName}_covariates_wgs_pcs_normalized.txt.gz covariates

  # 5. Extract and transpose PCs that explain variance greater than a certain threhsold
  python3 ${baseDir}/bin/extract_and_transpose_pcs.py --covariate_pcs covariates_PCs.txt \
    --covariates_explained_variance explainedVariance.txt \
    --explained_variance_threshold 0.99 

  # 6. Gzip the covariate matrix
  gzip covariates-pca_PCs-filtered-transpose.txt

  # 8. Regress covariates and residuals
  java -jar ${baseDir}/bin/2023-12-08-Regression.jar \
    ${forceNormalGeneCountMatrix} \
    covariates-pca_PCs-filtered-transpose.txt.gz \
    ./${params.cohortName}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz
  """
}

workflow {
  createCovariateMatrix(params.inputDir)
  createGeneCountsMatrix(params.inputDir)
  geneCountsTMMNormalization(createGeneCountsMatrix.out)
  // sampleSelection(geneCountsTMMNormalization.out.tmmNormalizedMatrix, createCovariateMatrix.out.matrix)
  normalizeCovariateMatrix(createCovariateMatrix.out.matrix)
  pcaOutlierFiltering(geneCountsTMMNormalization.out.tmmNormalizedMatrix)
  quantileNormalization(pcaOutlierFiltering.out.geneCountMatrix)
  geneForceNormal(quantileNormalization.out)
  geneCountsCovariateCorrection(geneForceNormal.out, normalizeCovariateMatrix.out)
  residualsPCA(geneCountsCovariateCorrection.out.geneCountMatrix)
  covariateAndResidualsCorrection(residualsPCA.out.residualPCs, createCovariateMatrix.out.matrix, geneForceNormal.out)
}