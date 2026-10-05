# ----------------------------------------------------------------------------------
# MASTER'S THESIS SCRIPT:
# Molecular signature analysis for stratification and prediction in Parkinson's disease
# Santos Antequera Fernández - Master's in Bioinformatics - Universidad Europea de Madrid
# ----------------------------------------------------------------------------------

# SCRIPT SUMMARY:
# 1. Environment setup, data download, and cleaning (GSE99039).
# 2. Calculation of biological pathway activity scores (pathway scores).
# 3. Training of Machine Learning models for clinical prediction.
# 4. Visualization and evaluation of model performance.
# 5. Unsupervised clustering analysis to identify patient subtypes.

# ----------------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# SECTION 1: INITIAL SETUP AND DATA LOADING
# -----------------------------------------------------------------------------

## 1.1. Working Directory and R Environment
# -----------------------------------------------------------------------------
# NOTE: This path is local and must be modified by other users.
setwd("Seleccione_su_directorio_de_trabajo")
# To speed up the workflow, we can load a previously saved R environment.
# This allows us to resume the analysis without having to rerun the
# longest and most computationally intensive calculations.
# load("entorno_TFM") # Initial load of unfiltered scores.
# load("entorno_TFM_finalML2") # Load of Machine Learning results.
# load("entorno_TFM_ML_plots") # Load of plotting objects.
load("entorno_TFM_TODO") # Loads the final complete environment.


## 1.2. Package Installation and Loading
# -----------------------------------------------------------------------------

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# We install and load 'pathMED', the main library for calculating
# biological pathway scores.
# BiocManager::install("pathMED")
library(pathMED)
### Documentation and examples:
# browseVignettes("pathMED")



# 1.3. Data Download from GEO using the GEOquery Package
# -----------------------------------------------------------------------------

#if (!require("GEOquery")) {
#  BiocManager::install("GEOquery")
#}

library(GEOquery) # This package allows downloading public data from the NCBI 
# Gene Expression Omnibus (GEO) repository

### Store the data in the gse variable
gse <- getGEO("GSE99039", GSEMatrix = TRUE)

# Check how many objects were downloaded
length(gse)

# Access the object:
gse_data <- gse[[1]]

## Extract the expression matrix and metadata

# Normalized expression matrix (genes x samples) 
exprs_matrix <- exprs(gse_data) # No need to transpose the matrix; it is in the 
# perfect format for the getScores() function

# Phenotypic data (sample metadata)
pheno_data <- pData(gse_data)

# Gene annotation data 
feature_data <- fData(gse_data)




# -----------------------------------------------------------------------------
# SECTION 2: PREPROCESSING AND CLEANING OF THE EXPRESSION MATRIX
# -----------------------------------------------------------------------------
# OBJECTIVE: Transform the expression matrix so that rows represent
# unique genes with standard symbols (e.g., "TP53"), instead of ambiguous
# probe IDs (e.g., "1053_at"). This step is essential for pathway analysis.

# Step 1. # We use annotation data (feature_data) to map probe IDs
# (rownames of exprs_matrix) to official gene symbols.
gene_symbols <- feature_data$`Gene Symbol`[match(rownames(exprs_matrix), feature_data$ID)]

# Step 2. Assign them as row names in the expression matrix
rownames(exprs_matrix) <- gene_symbols

# Step 3. Remove rows with missing or empty gene symbols
exprs_matrix <- exprs_matrix[!is.na(rownames(exprs_matrix)) & rownames(exprs_matrix) != "", ]

# Step 4. Aggregate duplicated genes by symbol
# Multiple probes may map to the same gene. To avoid errors in methods
# that require unique names, we aggregate them by gene symbol.
# Here we use the median as a more robust metric against outliers.

# CRITICAL NOTE: This approach (mean or median) is suitable for exploratory analysis,
# but more precise methods exist, such as:
#  - selecting the probe with the highest variance (most informative),
#  - using primary probe annotations,
#  - or applying weighted models based on probe quality.
# Here we use the median for simplicity and compatibility with functions such as `getScores()`.

# Validation before aggregation
duplicated_genes <- sum(duplicated(rownames(exprs_matrix)))
cat("Genes duplicados a combinar:", duplicated_genes, "\n") 

# Aggregate using median
exprs_matrix <- aggregate(exprs_matrix, by = list(Gene = rownames(exprs_matrix)), FUN = median)
rownames(exprs_matrix) <- exprs_matrix$Gene
exprs_matrix$Gene <- NULL # Remove auxiliary column

# Validation after aggregation (should be 0)
duplicated_genes_after <- sum(duplicated(rownames(exprs_matrix)))
cat("Genes duplicados después de agrupar:", duplicated_genes_after, "\n")





# -----------------------------------------------------------------------------
# SECTION 3: CALCULATION OF PATHWAY ACTIVITY SCORES (MOLECULAR SCORES)
# ----------------------------------------------------------------------------- 
# OBJECTIVE: Transform gene-level expression data into biological pathway-level
# scores. This reduces data dimensionality and allows interpreting results
# within a functional biological context.
# We test a combination of 6 scoring methods and 7 pathway databases
# to identify the most robust biological signals.

# 1. List of geneSets to use
gene_sets <- c("kegg", "reactome", "go_bp", "go_mf", "go_cc", "disgenet", "hpo")

# 2.1. singscore method. Loop to calculate and store results
for (gs in gene_sets) {
  scores <- getScores(exprs_matrix, geneSets = gs, method = "singscore",cores = 10)
  assign(paste0("scores_", gs), scores)  # Creates the variable: scores_kegg, etc.
  cat("\nPrimeros valores de", gs, ":\n")
  print(scores[1:5, 1:5])  # Displays the first values
}

# 2.2. Z-score method. Loop to calculate and store results
for (gs in gene_sets) {
  zscores <- getScores(exprs_matrix, geneSets = gs, method = "Z-score",cores = 10)
  assign(paste0("zscores_", gs), zscores)  # Creates the variable: zscores_kegg, etc.
  cat("\n[zScore] Primeros valores de", gs, ":\n")
  print(zscores[1:5, 1:5])
}

# 2.3. GSVA method. Loop to calculate and store results
for (gs in gene_sets) {
  gsvascores <- getScores(exprs_matrix, geneSets = gs, method = "GSVA",cores = 10)
  assign(paste0("gsvascores_", gs), gsvascores)  # Creates the variable: gsvascores_kegg, etc.
  cat("\n[gsvaScore] Primeros valores de", gs, ":\n")
  print(gsvascores[1:5, 1:5])
}

# 2.4. ssGSEA method. Loop to calculate and store results
for (gs in gene_sets) {
  ssgscores <- getScores(exprs_matrix, geneSets = gs, method = "ssGSEA",cores = 10)
  assign(paste0("ssgscores_", gs), ssgscores)  # Creates the variable: ssgscores_kegg, etc.
  cat("\n[ssGSEAScore] Primeros valores de", gs, ":\n")
  print(ssgscores[1:5, 1:5])
}

# 2.5. Plage method. Loop to calculate and store results
for (gs in gene_sets) {
  plagescores <- getScores(exprs_matrix, geneSets = gs, method = "Plage",cores = 10)
  assign(paste0("plagescores_", gs), plagescores)  # Creates the variable: plagescores_kegg, etc.
  cat("\n[PlageScore] Primeros valores de", gs, ":\n")
  print(plagescores[1:5, 1:5])
}

# 2.6. norm_FGSEA method. Loop to calculate and store results
# First we need to install and load the fgsea package for it to work
# Note: This method requires more time to compute the scores.

# BiocManager::install("fgsea")
library(fgsea)

for (gs in gene_sets) {
  nFscores <- getScores(exprs_matrix, geneSets = gs, method = "norm_FGSEA",cores = 10)
  assign(paste0("nFscores_", gs), nFscores)  # Creates the variable: nFscores_kegg, etc.
  cat("\n[norm_FGSEAScore] Primeros valores de", gs, ":\n")
  print(nFscores[1:5, 1:5])
}



# 3. To avoid repeating the entire process above to inspect the first scores, we do:
# Base names of your objects 
score_names <- c("kegg", "reactome", "go_bp", "go_mf", "go_cc", "disgenet", "hpo")

# Loop to print 5x5 of each scores_*
for (name in score_names) {
  cat("\nPrimeros valores de scores_", name, " (singscore):\n", sep = " ")
  print(get(paste0("scores_", name))[1:5, 1:5])
  
  cat("\nPrimeros valores de scores_", name, " (Z-score):\n", sep = " ")
  print(get(paste0("zscores_", name))[1:5, 1:5])
  
  cat("\nPrimeros valores de scores_", name, " (GSVA):\n", sep = " ")
  print(get(paste0("gsvascores_", name))[1:5, 1:5])
  
  cat("\nPrimeros valores de scores_", name, " (ssGSEA):\n", sep = " ")
  print(get(paste0("ssgscores_", name))[1:5, 1:5])
  
  cat("\nPrimeros valores de scores_", name, " (Plage):\n", sep = " ")
  print(get(paste0("plagescores_", name))[1:5, 1:5])
  
  cat("\nPrimeros valores de scores_", name, " (norm_FGSEA):\n", sep = " ")
  print(get(paste0("nFscores_", name))[1:5, 1:5])
}


#### To annotate identifiers with their corresponding descriptions, we use the ann2term() function
## First for singscore 
# 1. Create scores list
scores_list <- list(
  kegg     = scores_kegg,
  reactome = scores_reactome,
  go_bp    = scores_go_bp,
  go_mf    = scores_go_mf,
  go_cc    = scores_go_cc,
  disgenet = scores_disgenet,
  hpo      = scores_hpo
)
annotations <- lapply(scores_list, ann2term)
for (name in names(annotations)) {
  cat("\nPrimeras filas de las anotaciones de", name, ":\n", sep = " ")
  print(head(annotations[[name]]))
}


## Next for Z-score
zscores_list <- list(
  kegg     = zscores_kegg,
  reactome = zscores_reactome,
  go_bp    = zscores_go_bp,
  go_mf    = zscores_go_mf,
  go_cc    = zscores_go_cc,
  disgenet = zscores_disgenet,
  hpo      = zscores_hpo
)
zannotations <- lapply(zscores_list, ann2term)
for (name in names(zannotations)) {
  cat("\nPrimeras filas de las anotaciones (z-score) de", name, ":\n", sep = " ")
  print(head(zannotations[[name]]))
}


## Next for GSVA
gsvascores_list <- list(
  kegg     = gsvascores_kegg,
  reactome = gsvascores_reactome,
  go_bp    = gsvascores_go_bp,
  go_mf    = gsvascores_go_mf,
  go_cc    = gsvascores_go_cc,
  disgenet = gsvascores_disgenet,
  hpo      = gsvascores_hpo
)
gsvaannotations <- lapply(gsvascores_list, ann2term)
for (name in names(gsvaannotations)) {
  cat("\nPrimeras filas de las anotaciones (GSVA) de", name, ":\n", sep = " ")
  print(head(gsvaannotations[[name]]))
}


## Next for ssGSEA
ssgscores_list <- list(
  kegg     = ssgscores_kegg,
  reactome = ssgscores_reactome,
  go_bp    = ssgscores_go_bp,
  go_mf    = ssgscores_go_mf,
  go_cc    = ssgscores_go_cc,
  disgenet = ssgscores_disgenet,
  hpo      = ssgscores_hpo
)
ssgannotations <- lapply(ssgscores_list, ann2term)
for (name in names(ssgannotations)) {
  cat("\nPrimeras filas de las anotaciones (ssGSEA) de", name, ":\n", sep = " ")
  print(head(ssgannotations[[name]]))
}

## Next for Plage
plagescores_list <- list(
  kegg     = plagescores_kegg,
  reactome = plagescores_reactome,
  go_bp    = plagescores_go_bp,
  go_mf    = plagescores_go_mf,
  go_cc    = plagescores_go_cc,
  disgenet = plagescores_disgenet,
  hpo      = plagescores_hpo
)
plageannotations <- lapply(plagescores_list, ann2term)
for (name in names(plageannotations)) {
  cat("\nPrimeras filas de las anotaciones (Plage) de", name, ":\n", sep = " ")
  print(head(plageannotations[[name]]))
}


## Next for norm_FGSEA
nFscores_list <- list(
  kegg     = nFscores_kegg,
  reactome = nFscores_reactome,
  go_bp    = nFscores_go_bp,
  go_mf    = nFscores_go_mf,
  go_cc    = nFscores_go_cc,
  disgenet = nFscores_disgenet,
  hpo      = nFscores_hpo
)
nFannotations <- lapply(nFscores_list, ann2term)
for (name in names(nFannotations)) {
  cat("\nPrimeras filas de las anotaciones (norm_FGSEA) de", name, ":\n", sep = " ")
  print(head(nFannotations[[name]]))
}




save.image("entorno_TFM") # This environment contains all calculated scores.
# If the following ML part needs to be rerun, this environment can be loaded
# to perform everything cleanly from the start, without extraneous objects or libraries.



# -----------------------------------------------------------------------------
# SECTION 4: REDUCTION OF THE NUMBER OF PATHWAYS FOR MACHINE LEARNING MODEL TRAINING
# -----------------------------------------------------------------------------

## Example of the size of each database
sapply(list(scores_kegg, scores_reactome, scores_go_bp, scores_go_mf, scores_go_cc, scores_hpo, scores_disgenet), nrow)
# [1]   345  2483 12377  4398  1796  8698 10730
# Choosing 750 is fine, but ideally a dynamic filter should be applied; however, due to memory
# limitations, this option is chosen, as otherwise the computation time would be excessive
# (in previous tests calculating 1 or 2 go_bp models took several hours)

## Function to filter by variance (top N)
filter_top_var <- function(score_matrix, topN = 750) {
  n_before <- nrow(score_matrix)
  
  # Cleaning: Remove rows with non-finite values (NA, Inf, -Inf)
  # that could cause errors in variance calculation.
  clean_mat <- score_matrix[complete.cases(score_matrix), ]
  clean_mat <- clean_mat[apply(clean_mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
  n_after <- nrow(clean_mat)
  cat("Eliminadas", n_before - n_after, "features con NA/Inf (de", n_before, ")\n")
  
  # Calculate variance for each pathway (row).
  varianceScores <- apply(clean_mat, 1, function(x) var(x, na.rm = TRUE))
  
  # Sort pathways from highest to lowest variance.
  varianceScores <- sort(varianceScores, decreasing = TRUE)
  
  # Select the names of the 'topN' most variable pathways.
  nFeatures <- min(length(varianceScores), topN)
  selected <- names(varianceScores)[1:nFeatures]
  
  # Return the submatrix with only the selected pathways.
  return(clean_mat[selected, , drop = FALSE])
}


## Complete list of scores, applying the reduction function
scores_all <- list(
  singscore = list(
    kegg     = filter_top_var(scores_kegg, 750),
    reactome = filter_top_var(scores_reactome, 750),
    go_bp    = filter_top_var(scores_go_bp, 750),
    go_mf    = filter_top_var(scores_go_mf, 750),
    go_cc    = filter_top_var(scores_go_cc, 750),
    disgenet = filter_top_var(scores_disgenet, 750),
    hpo      = filter_top_var(scores_hpo, 750)
  ),
  Zscore = list(
    kegg     = filter_top_var(zscores_kegg, 750),
    reactome = filter_top_var(zscores_reactome, 750),
    go_bp    = filter_top_var(zscores_go_bp, 750),
    go_mf    = filter_top_var(zscores_go_mf, 750),
    go_cc    = filter_top_var(zscores_go_cc, 750),
    disgenet = filter_top_var(zscores_disgenet, 750),
    hpo      = filter_top_var(zscores_hpo, 750)
  ),
  GSVA = list(
    kegg     = filter_top_var(gsvascores_kegg, 750),
    reactome = filter_top_var(gsvascores_reactome, 750),
    go_bp    = filter_top_var(gsvascores_go_bp, 750),
    go_mf    = filter_top_var(gsvascores_go_mf, 750),
    go_cc    = filter_top_var(gsvascores_go_cc, 750),
    disgenet = filter_top_var(gsvascores_disgenet, 750),
    hpo      = filter_top_var(gsvascores_hpo, 750)
  ),
  ssGSEA = list(
    kegg     = filter_top_var(ssgscores_kegg, 750),
    reactome = filter_top_var(ssgscores_reactome, 750),
    go_bp    = filter_top_var(ssgscores_go_bp, 750),
    go_mf    = filter_top_var(ssgscores_go_mf, 750),
    go_cc    = filter_top_var(ssgscores_go_cc, 750),
    disgenet = filter_top_var(ssgscores_disgenet, 750),
    hpo      = filter_top_var(ssgscores_hpo, 750)
  ),
  Plage = list(
    kegg     = filter_top_var(plagescores_kegg, 750),
    reactome = filter_top_var(plagescores_reactome, 750),
    go_bp    = filter_top_var(plagescores_go_bp, 750),
    go_mf    = filter_top_var(plagescores_go_mf, 750),
    go_cc    = filter_top_var(plagescores_go_cc, 750),
    disgenet = filter_top_var(plagescores_disgenet, 750),
    hpo      = filter_top_var(plagescores_hpo, 750)
  ),
  norm_FGSEA = list(
    kegg     = filter_top_var(nFscores_kegg, 750),
    reactome = filter_top_var(nFscores_reactome, 750),
    go_bp    = filter_top_var(nFscores_go_bp, 750),
    go_mf    = filter_top_var(nFscores_go_mf, 750),
    go_cc    = filter_top_var(nFscores_go_cc, 750),
    disgenet = filter_top_var(nFscores_disgenet, 750),
    hpo      = filter_top_var(nFscores_hpo, 750)
  )
)





# -----------------------------------------------------------------------------
# SECTION 5: SUPERVISED CLASSIFICATION WITH pathMED
# -----------------------------------------------------------------------------
# OBJECTIVE: Train and evaluate models capable of predicting the condition of a
# subject (IPD or Control) from their pathway activity profiles.

# --------------------------
# 1) Response variable
# --------------------------

# Extract clinical variable of interest from `pheno_data`.
pheno_data$Response <- as.character(pheno_data$`disease label:ch1`)  # key: character for pathMED

# Explicitly define positive ("IPD") and negative ("CONTROL") classes.
pos_class <- "IPD"
neg_class <- "CONTROL"

# Filter metadata table to keep only samples
# belonging to these two groups.
pheno_bin <- pheno_data[pheno_data$Response %in% c(pos_class, neg_class), , drop=FALSE]
# Check class balance.
table(pheno_bin$Response)


# --------------------------
# Exploratory PCA: IPD vs Control
# --------------------------
# Principal Component Analysis (PCA) allows visualizing data structure in 2D.
# It is a first step to assess whether samples from both groups (IPD, Control)
# are separable based on pathway scores.
library(ggplot2)

# --- 1. Create a folder to save all plots ---
# Create directory if it does not exist
if(!dir.exists("ML_Plots")) dir.create("ML_Plots")
pca_output_dir <- "ML_Plots/PCA_Exploratorio"
if (!dir.exists(pca_output_dir)) dir.create(pca_output_dir, recursive = TRUE)

cat("Iniciando la generación de plots PCA para todas las combinaciones...\n")

# --- 2. Nested loop to generate a PCA for each method and database combination.
for (method in names(scores_all)) {
  for (db in names(scores_all[[method]])) {
    
    file_tag <- paste(method, db, sep = "_")
    cat(sprintf("Generando PCA para: %s\n", file_tag))
    
    # a) Select score matrix dynamically
    mat_PCA <- scores_all[[method]][[db]]
    
    # b) Sample alignment: ensure score matrix and metadata
    #    contain exactly the same samples in the same order.
    common_samples_PCA <- intersect(colnames(mat_PCA), rownames(pheno_bin))
    mat_PCA_sub <- mat_PCA[, common_samples_PCA]
    pheno_bin_PCA <- pheno_bin[common_samples_PCA, ]
    
    # Skip if there are insufficient data for PCA
    if(ncol(mat_PCA_sub) < 3) {
      cat(sprintf("  -> Saltando %s: no hay suficientes muestras.\n", file_tag))
      next
    }
    
    # c) PCA calculation. Matrix is transposed (t()) because prcomp expects
    #    samples in rows and features in columns.
    pca <- prcomp(t(mat_PCA_sub), scale. = TRUE)
    
    # d) Create data.frame for ggplot2 visualization.
    pca_df <- data.frame(
      Sample = rownames(pca$x),
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      Group = pheno_bin_PCA$Response
    )
    
    # e) Plot creation, including the percentage of explained variance
    #    by each principal component on the axes.
    pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 3, alpha = 0.8) +
      labs(
        title = paste("PCA (CONTROL vs IPD) -", file_tag),
        x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")
      ) +
      theme_minimal()
    
    # f) Save plot with unique filename
    ggsave(
      filename = file.path(pca_output_dir, paste0("PCA_", file_tag, ".png")),
      plot = pca_plot,
      width = 7, height = 5, dpi = 300
    )
  }
}

cat("\n--- Plots PCA generados y guardados en la carpeta 'PCA_Exploratorio' ---\n")


# --------------------------
# 2) Definition of Machine Learning Models
# --------------------------
modelsList <- methodsML(
  algorithms = c("rf","knn","xgbTree","glm","lda"), 
  outcomeClass = "character",   # mandatory in pathMED
  tuneLength = 10  # Ideally higher, but 10 was chosen due to computational and time constraints            
)

# --------------------------
# 3) Train models
# --------------------------
# This is the core of the Machine Learning analysis.
# Load required libraries for algorithms and parallelization.

library(pathMED)
library(caret)        
library(doParallel)   # parallelization
library(randomForest) # for rf
library(kknn)         # for improved knn
library(xgboost)      # for xgbTree
library(MASS)

# Fixed base seed
base_seed <- 1234

# To run from scratch
results_df <- data.frame(
  ScoreMethod = character(),
  Database = character(),
  BestModel = character(),
  Accuracy = numeric(),
  BalancedAcc = numeric(),
  MCC = numeric(),
  Recall = numeric(),
  Specificity = numeric(),
  Precision = numeric(),
  F1 = numeric(),
  stringsAsFactors = FALSE
)
# Create directory if it does not exist
if(!dir.exists("ML_Plots")) dir.create("ML_Plots")


# Set up parallel computing cluster using 10 CPU cores.
# This substantially accelerates training and cross-validation.
# using the doParallel library
cl <- makeCluster(10)   
registerDoParallel(cl)

# Main loop for model training
for (method in names(scores_all)) {
  for (db in names(scores_all[[method]])) {
    
    # Align samples correctly between input_mat and meta_bin
    common_samples <- intersect(colnames(scores_all[[method]][[db]]), rownames(pheno_bin))
    if(length(common_samples) < 2){
      cat("Saltando", method, db, "- no hay suficientes muestras\n")
      next
    }
    
    # Create final training matrices using only common samples.
    input_mat <- scores_all[[method]][[db]][, common_samples, drop=FALSE]
    meta_bin <- pheno_bin[common_samples, , drop=FALSE]
    
    # Reorder meta_bin according to input_mat for safety
    # This prevents label assignment errors across samples.
    meta_bin <- meta_bin[colnames(input_mat), , drop=FALSE]
    
    # Verify that after filtering, data from both classes are still present
    # (IPD and Control) to train a classification model.
    clase_check <- table(meta_bin$Response)
    if (length(clase_check) < 2) {
      cat("Saltando", method, db, "- solo hay una clase presente después de alinear\n")
    } else {
      cat("Alineación correcta para", method, db, "-", 
          clase_check[pos_class], pos_class, "y", clase_check[neg_class], neg_class, "\n")
    }
    
    
    # Set unique seed for each combination.
    # This ensures that upon re-running the script, cross-validation results
    # (which involve randomness in data partitioning) will be identical.
    seed_iter <- base_seed + match(method, names(scores_all)) * 100 + match(db, names(scores_all[[method]]))
    set.seed(seed_iter)
    
    # Confirm seed and combination
    cat("Entrenando:", method, "-", db, "con semilla:", seed_iter, "\n")
    
    # This is the main function that trains and evaluates defined models.
    # Uses nested cross-validation to obtain a robust and reliable
    # estimate of model performance on unseen data.
    trainedModel <- trainModel(
      inputData = input_mat,
      metadata = meta_bin,
      var2predict = "Response",
      positiveClass = pos_class,   # "IPD"
      models = modelsList,
      Koutter = 5,
      Kinner = 3,
      repeatsCV = 3
    )
    
    # Across all tested models (rf, knn, etc.), pathMED selects the best one.
    # Extract its performance metrics.
    stats <- trainedModel$stats
    best  <- trainedModel$model$method
    
    # Store metrics
    results_df <- rbind(
      results_df,
      data.frame(
        ScoreMethod = method,
        Database = db,
        BestModel = best,
        Accuracy = stats["accuracy", best],
        BalancedAcc = stats["balacc", best],
        MCC = (stats["mcc", best] + 1) / 2,  # MCC normalized to [0,1]
        Recall = stats["recall", best],
        Specificity = stats["specificity", best],
        Precision = stats["precision", best],
        F1 = stats["fscore", best]
      )
    )
    
    # Save full 'trainedModel' object. It is a large object containing
    # the final model, predictions, etc. Needed later for plots.
    save(trainedModel, file = paste0("ML_Plots/trainedModel_", method, "_", db, ".RData"))
    
    # The 'trainedModel' object can consume significant RAM, so remove it
    # after finishing each combination.
    rm(trainedModel)
    gc()
    
    print(results_df)
    
    # Save results table after each iteration
    save(results_df, file="results_todo_ML2.RData")
  }
}

# Once loop finishes, stop parallel cluster to release CPU cores.
stopCluster(cl)


# Save results table after each iteration to a .csv file
write.csv(results_df, "results_todo_ML2.csv", row.names = FALSE)



# Estimated time for methods with kegg as database: 2 minutes (value for initial runs)
# Estimated time for methods filtered to 750 scores: 5 minutes (value for initial runs)

save.image("entorno_TFM_finalML2") # This environment contains scores 
# and Machine Learning results with pathMED






# =============================================================================
# SECTION 6: VISUALIZATION AND ANALYSIS OF MACHINE LEARNING RESULTS
# =============================================================================
# OBJECTIVE: Once models are trained, evaluate and visualize performance
# to draw conclusions. In this section, generate a series of plots comparing
# different approaches (scoring methods and databases) and analyze the behavior
# of each trained model in detail.


# -----------------------------------------------------------------------------
# 1. Metrics Heatmap by ScoreMethod and Database 
# -----------------------------------------------------------------------------
# This heatmap provides an overview of performance across all tested combinations. 
library(reshape2)
library(ggplot2)

metrics <- c("Accuracy","BalancedAcc","MCC","Recall","Specificity","Precision","F1")

# The `melt` function transforms `results_df` from "wide" format
# (one column per metric) to "long" format (one column for metric name
# and one for its value). This format is required by ggplot for `facet_wrap`.
df_melt <- melt(results_df, id.vars = c("ScoreMethod","Database"), measure.vars = metrics)

# Create plot with ggplot2
p_heat <- ggplot(df_melt, aes(x = Database, y = ScoreMethod, fill = value)) +
  geom_tile() +
  facet_wrap(~variable) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Heatmap de métricas por ScoreMethod y Database", fill="Valor")

ggsave("ML_Plots/Heatmap_metrics.png", plot=p_heat, width=10, height=6, dpi=300)


# -----------------------------------------------------------------------------
# 2. ROC and PR curves for best models (using subsample.preds)
# -----------------------------------------------------------------------------
# These curves evaluate model discriminative power.
# Generated using "Out-of-Fold" (OOF) predictions, which are predictions made
# on unseen data during cross-validation training folds. This provides a realistic
# and reliable estimate of generalization performance.

library(pROC)
library(PRROC)

# Create output directory if it does not exist
if (!dir.exists("ML_Plots/ROC_PR")) dir.create("ML_Plots/ROC_PR", recursive = TRUE)

# Explicitly define classes
pos_class <- "IPD"
neg_class <- "CONTROL"

for (i in 1:nrow(results_df)) {
  method <- results_df$ScoreMethod[i]
  db     <- results_df$Database[i]
  
  # Load `trainedModel` object saved during training.
  load(paste0("ML_Plots/trainedModel_", method, "_", db, ".RData"))  
  
  # Extract OOF ("Out-of-Fold") predictions.
  oof <- trainedModel$subsample.preds
  if (is.null(oof)) { 
    cat("Sin subsample.preds para", method, db, "\n") 
    next 
  }
  
  # Check required columns (IPD, CONTROL, and obs)
  if (!all(c(pos_class, neg_class, "obs") %in% colnames(oof))) {
    cat("Esperaba columnas", pos_class, ",", neg_class, "y obs en", method, db, "\n")
    cat("Columnas disponibles:", paste(colnames(oof), collapse = ", "), "\n")
    next
  }
  
  # Prepare data: 'obs' represents ground truth (actual patient status)
  # and 'prob' is the probability assigned to the positive class (IPD).
  obs  <- factor(oof$obs, levels = c(neg_class, pos_class))
  prob <- oof[[pos_class]]  # positive class probability (IPD)
  
  # Verify both classes are present
  if (length(unique(obs)) < 2) {
    cat("Saltando ROC/PR para", method, db, "- solo una clase presente\n")
    next
  }
  
  # ========================
  # ROC Curve.
  # Shows tradeoff between True Positive Rate (sensitivity) and
  # False Positive Rate (1 - specificity). A perfect model yields
  # an area under the curve (AUC) of 1.
  # ========================
  roc_obj <- roc(
    response = obs,
    predictor = prob,
    levels = c(neg_class, pos_class),
    direction = "<"
  )
  
  png(paste0("ML_Plots/ROC_PR/ROC_", method, "_", db, ".png"), width = 900, height = 700)
  plot(roc_obj, main = paste("ROC (OOF):", method, "-", db), col = "#1C6DD0", lwd = 3)
  legend("bottomright", legend = sprintf("AUC = %.3f", auc(roc_obj)), col = "#1C6DD0", lwd = 3)
  dev.off()
  
  # ========================
  # PR Curve.
  # Informative for imbalanced classes. Shows tradeoff between
  # Precision (fraction of predicted positives that are true positives)
  # and Recall (fraction of actual positives detected).
  # ========================
  pr_obj <- pr.curve(
    scores.class0 = prob[obs == pos_class],
    scores.class1 = prob[obs == neg_class], # Note: PRROC inverts class order
    curve = TRUE
  )
  
  png(paste0("ML_Plots/ROC_PR/PR_", method, "_", db, ".png"), width = 900, height = 700)
  plot(pr_obj$curve[, 1], pr_obj$curve[, 2],
       type = "l", lwd = 3, col = "#E91E63",
       xlab = "Recall", ylab = "Precision",
       main = paste("PR (OOF):", method, "-", db))
  legend("bottomright", legend = sprintf("AUPRC = %.3f", pr_obj$auc.integral), 
         col = "#E91E63", lwd = 3)
  dev.off()
  
  cat("ROC y PR (OOF) generados para", method, db, "\n")
}



# -----------------------------------------------------------------------------
# 3. Confusion Matrices 
# -----------------------------------------------------------------------------
# OBJECTIVE: Analyze model classification errors in detail. The confusion matrix
# reports how many patients were classified correctly (True Positives/Negatives)
# versus incorrectly (False Positives/Negatives). Again, OOF predictions
# are used for realistic performance estimation.

library(caret)

# Output folder for plots if not existing
if (!dir.exists("ML_Plots/Confusion_matrix")) dir.create("ML_Plots/Confusion_matrix", recursive = TRUE)

# Explicitly define classes (performed previously)
pos_class <- "IPD"
neg_class <- "CONTROL"

for (i in 1:nrow(results_df)) {
  
  method_imp <- results_df$ScoreMethod[i]
  db_imp     <- results_df$Database[i]
  
  # === 1) Load trained model ===
  model_file <- paste0("ML_Plots/trainedModel_", method_imp, "_", db_imp, ".RData")
  if (!file.exists(model_file)) {
    cat("Modelo no encontrado para", method_imp, db_imp, "\n")
    next
  }
  load(model_file)  # Restores 'trainedModel' object into environment.
  
  # === 2) Extract OOF (subsample) predictions ===
  oof <- trainedModel$subsample.preds
  if (is.null(oof) || !all(c(pos_class, neg_class, "obs") %in% colnames(oof))) {
    cat("subsample.preds no disponible o sin columnas esperadas (", 
        pos_class, ",", neg_class, ", obs ) para", method_imp, db_imp, "\n")
    next
  }
  
  # === 3) Ground truth and predicted classes from OOF probabilities ===
  # The model yields probabilities (e.g., 0.8 probability of IPD).
  # For the confusion matrix, assign categorical prediction ('IPD' or 'CONTROL').
  # Use standard threshold of 0.5: if probability >= 50%, predict 'IPD'.
  obs  <- factor(oof$obs, levels = c(neg_class, pos_class))
  prob <- oof[[pos_class]]
  pred <- factor(ifelse(prob >= 0.5, pos_class, neg_class), levels = c(neg_class, pos_class))
  
  # Check both classes
  if (length(unique(obs)) < 2) {
    cat("Saltando Confusion Matrix para", method_imp, db_imp, "- solo hay una clase en OOF\n")
    next
  }
  
  # === 4) Confusion matrix (OOF) ===
  cm <- caret::confusionMatrix(pred, obs, positive = pos_class)
  
  # Print to console
  cat("\nMatriz de confusión (OOF) para", method_imp, db_imp, ":\n")
  print(cm)
  
  # === 5) Save plot (fourfoldplot) ===
  png(paste0("ML_Plots/Confusion_matrix/ConfMatrix_", method_imp, "_", db_imp, ".png"),
      width = 800, height = 600)
  # `fourfoldplot` provides a visual representation of the matrix. The area of each
  # quadrant is proportional to sample count in that cell.
  fourfoldplot(cm$table,
               color = c("#FF9999", "#99CCFF"),
               conf.level = 0,
               margin = 1,
               main = paste("Confusion Matrix (OOF):", method_imp, "-", db_imp))
  dev.off()
  
  cat("Matriz de confusión generada para", method_imp, db_imp, "\n")
}

# -----------------------------------------------------------------------------
# 4. Calibration Plots
# -----------------------------------------------------------------------------
# OBJECTIVE: Evaluate whether model predicted probabilities are reliable.
# In a well-calibrated model, if it predicts 80% probability across a group of samples,
# approximately 80% of those samples should truly belong to the positive class.
# This is critical for assessing prediction confidence.

library(dplyr)
library(ggplot2)

if (!dir.exists("ML_Plots/Calibration_plot")) dir.create("ML_Plots/Calibration_plot", recursive = TRUE)

# Explicitly define classes (performed previously)
pos_class <- "IPD"
neg_class <- "CONTROL"

for (i in 1:nrow(results_df)) {
  method <- results_df$ScoreMethod[i]
  db     <- results_df$Database[i]
  
  # Load trained model 
  model_file <- paste0("ML_Plots/trainedModel_", method, "_", db, ".RData")
  if (!file.exists(model_file)) {
    cat("Modelo no encontrado para", method, db, "\n")
    next
  }
  load(model_file)
  
  # Extract OOF predictions.
  oof <- trainedModel$subsample.preds
  if (is.null(oof) || !all(c(pos_class, neg_class, "obs") %in% colnames(oof))) {
    cat("Sin OOF (subsample.preds) con columnas", pos_class, "/", neg_class, "/obs para", method, db, "\n")
    if (!is.null(oof)) cat("Columnas disponibles:", paste(colnames(oof), collapse = ", "), "\n")
    next
  }
  
  # 1) Prepare data for analysis ===
  # Ground truth in numeric format (0 for CONTROL, 1 for IPD)
  # to calculate means and calibration error.
  obs_fac <- factor(oof$obs, levels = c(neg_class, pos_class))
  if (length(unique(obs_fac)) < 2) {
    cat("Saltando Calibration para", method, db, "- solo una clase presente en OOF\n")
    next
  }
  prob    <- oof[[pos_class]]                 # probability of IPD (positive class)
  obs_num <- as.integer(obs_fac == pos_class) # 1=IPD, 0=CONTROL
  
  # 2) Brier score (overall calibration)
  # Overall metric measuring calibration error: mean squared difference
  # between predicted probabilities and observed binary outcomes (0 or 1).
  # Values closer to 0 reflect superior calibration.
  brier <- mean((prob - obs_num)^2, na.rm = TRUE)
  
  # 3) Calibration curve by probability bins (deciles by default)
  #    Safeguard if few unique values exist
  unique_probs <- sum(!is.na(unique(prob)))
  if (unique_probs < 3) {
    # With minimal variation, use fixed bins
    cuts <- seq(0, 1, by = 0.2)
  } else {
    cuts <- quantile(prob, probs = seq(0, 1, by = 0.1), na.rm = TRUE, names = FALSE)
    # Avoid duplicate boundary cuts from ties
    if (any(duplicated(cuts))) {
      cuts <- unique(cuts)
      # Ensure inclusion of 0 and 1
      if (min(cuts) > 0) cuts <- c(0, cuts)
      if (max(cuts) < 1) cuts <- c(cuts, 1)
    }
  }
  # Ensure at least 3 cuts to form bins
  if (length(cuts) < 3) cuts <- unique(sort(c(0, cuts, 1)))
  
  # Group predictions into probability bins (e.g., 0.1 to 0.2).
  # For each bin, compute mean predicted probability and observed positive proportion.
  df_cal <- data.frame(prob = prob, obs = obs_num) %>%
    mutate(bin = cut(prob, breaks = cuts, include.lowest = TRUE, right = TRUE)) %>%
    group_by(bin, .drop = TRUE) %>%
    summarise(
      pred_mean = mean(prob, na.rm = TRUE), # Mean predicted probability in bin
      obs_rate  = mean(obs,  na.rm = TRUE), # Observed positive rate in bin
      n         = dplyr::n(),
      .groups   = "drop"
    ) %>%
    filter(is.finite(pred_mean), is.finite(obs_rate))
  
  # 4) Plot and save
  # A perfectly calibrated model aligns along the 45-degree dashed diagonal
  # (where predicted probability equals observed proportion).
  p <- ggplot(df_cal, aes(x = pred_mean, y = obs_rate)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point() +
    geom_line() +
    labs(
      title = paste("Calibration (OOF):", method, "-", db),
      subtitle = sprintf("Brier score = %.4f", brier),
      x = paste0("Predicted probability of ", pos_class, " (bin mean)"),
      y = paste0("Observed rate of ", pos_class)
    ) +
    theme_minimal(base_size = 12)
  
  ggsave(filename = paste0("ML_Plots/Calibration_plot/Calibration_", method, "_", db, ".png"),
         plot = p, width = 8, height = 6, dpi = 150)
  
  cat("Calibration (OOF) para", method, db, "- Brier:", sprintf("%.4f", brier), "\n")
}


# -----------------------------------------------------------------------------
# 5. Feature Importance + ann2term in plots
#    Exports ALL.csv, top20.csv, and top5 plot
#    (with progress logging and progress bar for KNN)
# -----------------------------------------------------------------------------
# OBJECTIVE: Identify which pathways (features) contribute most significantly
# to distinguishing IPD patients from Controls.
# Essential for biological interpretability: elucidates not only THAT the model works,
# but HOW it arrives at predictions.

# Load required libraries
library(caret)
library(pathMED)
library(ggplot2)
library(dplyr)
library(stringr)

#--- Logging function ---------------------------------------------------------
VERBOSE <- TRUE
logf <- function(fmt, ...) {
  if (!VERBOSE) return(invisible(NULL))
  cat(format(Sys.time(), "%H:%M:%S"), "-", sprintf(fmt, ...), "\n")
  flush.console()
}

# Function to normalize pathway IDs (e.g., "GO.123" -> "GO:123")
normalize_ids <- function(x){
  x <- sub("^X(?=\\d)", "", x, perl=TRUE)
  x <- sub("^GO\\.", "GO:", x); x <- sub("^HP\\.", "HP:", x)
  x <- sub("^R[._-]?HSA[._-]", "R-HSA-", x)
  sub("[._](\\d{3,7})$", "-\\1", x, perl=TRUE)
}

# Function to retrieve annotations for a given method and database combination.
get_ann <- function(score_method, gene_set){
  prefix <- c(singscore="scores", Zscore="zscores", GSVA="gsvascores",
              ssGSEA="ssgscores", Plage="plagescores", norm_FGSEA="nFscores")[score_method]
  if (is.na(prefix)) return(NULL)
  nm <- paste0(prefix, "_", gene_set)
  mat <- tryCatch(get(nm, envir=.GlobalEnv), error=function(e) NULL)
  if (is.null(mat)) return(NULL)
  tryCatch(ann2term(mat), error=function(e)
    data.frame(ID=rownames(mat), term=NA_character_, stringsAsFactors=FALSE))
}

# SPECIAL helper to compute feature importance in KNN models.
# KNN models do not possess intrinsic importance metrics.
# This function implements permutation importance:
# 1. Measures model accuracy on original data.
# 2. For each pathway, randomly shuffles feature values across samples.
# 3. Measures resulting performance drop.
# 4. Performance degradation indicates pathway importance.
imp_knn <- function(model, nsim=5, seed=123){
  if (is.null(model$trainingData)) stop("train$trainingData no disponible.")
  set.seed(seed)
  df <- model$trainingData
  y  <- if (is.factor(df$.outcome)) df$.outcome else factor(df$.outcome)
  X  <- df[, setdiff(names(df), intersect(names(df), c(".outcome",".weights"))), drop=FALSE]
  base_acc <- mean(predict(model, X) == y)
  
  p_total <- ncol(X) * nsim
  done <- 0L
  pb <- NULL
  if (interactive()) pb <- txtProgressBar(min=0, max=p_total, style=3)
  on.exit({ if (!is.null(pb)) close(pb) }, add=TRUE)
  
  drop_acc <- numeric(ncol(X))
  names(drop_acc) <- colnames(X)
  
  for (j in seq_along(colnames(X))) {
    vals <- numeric(nsim)
    for (s in seq_len(nsim)) {
      Xp <- X; Xp[[j]] <- sample(Xp[[j]])
      vals[s] <- base_acc - mean(predict(model, Xp) == y)
      done <- done + 1L
      if (!is.null(pb)) setTxtProgressBar(pb, done)
    }
    m <- mean(vals, na.rm=TRUE)
    drop_acc[j] <- ifelse(m < 0, 0, m)
  }
  data.frame(Term=names(drop_acc), Overall=as.numeric(drop_acc), check.names=FALSE)
}

#--- Define output directories -----------------------------------------------
dir_plots <- "ML_Plots/Importance/Plots"
dir_csv   <- "ML_Plots/Importance/Anotations"
dir.create(dir_plots, recursive=TRUE, showWarnings=FALSE)
dir.create(dir_csv,   recursive=TRUE, showWarnings=FALSE)

#--- MAIN LOOP: Calculate and visualize feature importance per model ---
for(i in seq_len(nrow(results_df))){
  sm <- results_df$ScoreMethod[i]
  gs <- results_df$Database[i]
  alg <- results_df$BestModel[i]
  tag <- paste(alg, sm, gs, sep="_") # Unique tag for files.
  f   <- file.path("ML_Plots", paste0("trainedModel_", sm, "_", gs, ".RData"))
  
  logf("(%d/%d) Iniciando %s", i, nrow(results_df), tag)
  
  # Check model file exists before attempting to load.
  if(!file.exists(f)){
    logf("AVISO: Modelo no encontrado, saltando: %s", f)
    next # If not found, skip to next iteration.
  }
  
  # Measure load time in R and print execution duration in seconds.
  t_load <- proc.time()[3]; load(f); t_load <- proc.time()[3]-t_load
  logf("Modelo cargado en %.2fs", t_load)
  
  mdl <- trainedModel$model
  
  # --- 1. Importance Calculation ---
  t_imp <- proc.time()[3]
  # Determine which pathways (features) were most influential in model predictions.
  # Method differs by algorithm family.
  if (isTRUE(grepl("^knn$", mdl$method, ignore.case=TRUE))) {
    # For KNN, apply permutation feature importance: shuffle feature values
    # and quantify performance degradation.
    nfeat <- if (!is.null(mdl$finalModel$xNames)) length(mdl$finalModel$xNames) else NA_integer_
    logf("Cálculo importancia por permutación (KNN, nsim=5, p=%s)", ifelse(is.na(nfeat), "?", nfeat))
    imp <- tryCatch({
      df <- imp_knn(mdl, nsim=5, seed=123)
      tibble(Term=df$Term, Overall=df$Overall)
    }, error=function(e){ logf("ERROR imp_knn: %s", e$message); tibble(Term=character(), Overall=numeric()) })
  } else {
    # For tree-based models (Random Forest, XGBoost), use standard caret `varImp`.
    logf("Cálculo varImp (%s)", mdl$method)
    imp <- tryCatch({
      vi <- varImp(mdl, scale=FALSE)
      df <- if (inherits(vi, "varImp.train")) vi$importance else vi
      if (!"Overall" %in% colnames(df)) colnames(df)[1] <- "Overall"
      tibble(Term=rownames(df), Overall=df$Overall)
    }, error=function(e){ logf("ERROR varImp: %s", e$message); tibble(Term=character(), Overall=numeric()) })
  }
  t_imp <- proc.time()[3]-t_imp
  logf("Importancia calculada en %.2fs (n=%d filas)", t_imp, nrow(imp))
  if (nrow(imp) == 0){ logf("AVISO: Sin importancia válida para %s. Saltando.", tag); next }
  
  # --- 2. Pathway Term Annotation ---
  # Map pathway accession IDs (e.g., "hsa04110") to descriptive biological terms (e.g., "Cell cycle").
  t_ann <- proc.time()[3]
  ann <- get_ann(sm, gs)
  logf("Anotación ann2term: %d términos", if (is.null(ann)) 0L else nrow(ann))
  
  # Merge importance data with annotation metadata.
  imp <- imp %>%
    arrange(desc(Overall)) %>% # Sort pathways from most to least important.
    mutate(ID_norm = normalize_ids(Term)) %>% # Standardize IDs for alignment.
    left_join({ if(is.null(ann)) tibble(ID=character(), term=character()) else as_tibble(ann) },
              by = c("ID_norm"="ID")) %>% # Merge tables.
    mutate(TermPlot = ifelse(is.na(term) | term=="", Term, term)) %>% # Fallback to ID if term unavailable.
    select(Term = TermPlot, Importance = Overall, ID = ID_norm) # Select and rename columns.
  
  t_ann <- proc.time()[3]-t_ann
  logf("Anotado en %.2fs", t_ann)
  
  # --- 3. EXPORT CSVs AND PLOT ---
  all_csv_path <- file.path(dir_csv, paste0(tag, "_all_importances.csv"))
  
  # Save FULL feature importance list to CSV.
  # 'write.csv2' uses semicolon as delimiter.
  write.csv2(imp, all_csv_path, row.names = FALSE)
  logf("CSV con TODAS las importancias guardado en: %s", all_csv_path)
  
  top20_csv_path <- file.path(dir_csv, paste0(tag, "_top20_importances.csv"))
  
  # Export top 20 features to CSV.
  write.csv2(imp %>% slice_head(n = 20), top20_csv_path, row.names = FALSE)
  logf("CSV con el TOP 20 guardado en: %s", top20_csv_path)
  
  # Barplot of top 5 most important pathways.
  top5 <- imp %>% slice_head(n = 5)
  if (nrow(top5) > 0) {
    # `str_wrap` wraps long labels for formatting.
    top5$Term <- stringr::str_wrap(top5$Term, width = 45)
    # Order plot bars by descending importance.
    top5$Term <- factor(top5$Term, levels = rev(top5$Term))
    
    p <- ggplot(top5, aes(x = Term, y = Importance)) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      coord_flip() + # Horizontal orientation for label legibility.
      labs(
        title = paste("Top 5 Importancia:", tag),
        x = "Término",
        y = "Importancia (Overall)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.y  = element_text(size = 12),
        axis.text.x  = element_text(size = 12),
        axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin  = margin(15, 20, 15, 20)
      )
    
    png_path <- file.path(dir_plots, paste0("importance_plot_", tag, ".png"))
    ggsave(png_path, p, width = 10, height = 7, dpi = 300, bg = "white")
    logf("Plot TOP 5 guardado: %s", png_path)
  } else {
    logf("AVISO: No hay datos para generar el plot de top 5.")
  }
  logf("--- Completado %s ---", tag)
}

# --- Final message ---
logf("========== PROCESO FINALIZADO ==========")
cat("Resultados guardados en las siguientes carpetas:\n")
cat("-> Plots (.png):", normalizePath(dir_plots), "\n")
cat("-> Anotaciones (.csv):", normalizePath(dir_csv), "\n")



# -----------------------------------------------------------------------------
# 6. Heatmaps (PNG) with IPD on the left and CONTROL on the right
# -----------------------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)
library(pathMED)
library(grid)

# General Heatmap Configuration
dir.create("ML_Plots/Heatmaps", showWarnings = FALSE, recursive = TRUE)
DO_ZSCORE <- TRUE
TOP_N     <- 100
PNG_W     <- 2200; PNG_H <- 2000; PNG_RES <- 200

# Class Definitions (performed also in prior visualizations)
pos_class <- "IPD"
neg_class <- "CONTROL"

# Cleans pathway IDs for consistency.
norm_ids <- function(ids){
  x <- ids
  x <- sub("^X(?=\\d)", "", x, perl=TRUE)
  x <- sub("^GO\\.", "GO:", x)
  x <- sub("^HP\\.", "HP:", x)
  x <- sub("^R[._-]?HSA[._-]?", "R-HSA-", x)
  x <- sub("[._](\\d{3,7})$", "-\\1", x, perl=TRUE)
  make.unique(x)
}
# Computes Z-score across rows (pathways). Standardizes features,
# highlighting relative deviation patterns rather than raw scale.
row_z <- function(m){
  m <- as.matrix(m)
  mu <- rowMeans(m, na.rm=TRUE)
  sdv <- apply(m, 1, sd, na.rm=TRUE); sdv[sdv==0 | is.na(sdv)] <- 1
  sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
}
# Selects 'n' rows with highest variance.
top_var <- function(m, n=100){
  if(nrow(m)<=n) return(m)
  v <- apply(m, 1, var, na.rm=TRUE)
  m[order(v, decreasing=TRUE)[seq_len(n)], , drop=FALSE]
}
# Orders samples placing IPD first, followed by CONTROL.
order_cols_pos_neg <- function(y_chr, pos=pos_class, neg=neg_class){
  if (all(c(pos,neg) %in% unique(y_chr))) c(which(y_chr==pos), which(y_chr==neg)) else order(y_chr)
}
# Annotates pathway IDs with descriptive terms using ann2term.
annotate_all <- function(row_ids, values_for_ann){
  m <- as.matrix(setNames(values_for_ann, row_ids))
  ann <- tryCatch(ann2term(m), error=function(e) NULL)
  if (is.null(ann)) return(row_ids)
  mm <- match(row_ids, rownames(ann))
  out <- row_ids
  ok <- !is.na(mm) & !is.na(ann$term[mm]) & nzchar(ann$term[mm])
  out[ok] <- ann$term[mm][ok]
  out
}

# Main Loop to Generate Heatmaps ---
files <- list.files("ML_Plots", pattern="^trainedModel_.*\\.RData$", full.names=TRUE)
for (f in files){
  load(f)  # -> trainedModel
  
  # Extract training data used for this specific model.
  td <- if(!is.null(trainedModel$model$trainingData)) trainedModel$model$trainingData else trainedModel$trainingData
  if (is.null(td) || !(".outcome" %in% names(td))) {
    cat("Sin trainingData utilizable en:", f, "- se salta\n")
    rm(trainedModel); gc(); next
  }
  
  # Separate outcome labels (y) from features (X).
  y <- as.character(td$.outcome)
  X <- td[, setdiff(names(td), ".outcome"), drop=FALSE]  # samples x features
  # Transpose matrix to standard format: pathways in rows, samples in columns.
  M <- t(as.matrix(X))                                   # rows=pathways, cols=samples
  rownames(M) <- norm_ids(rownames(M))
  
  # Order columns (samples): IPD first, then CONTROL.
  ord <- order_cols_pos_neg(y)   # IPD -> CONTROL
  M   <- M[, ord, drop=FALSE]
  y   <- y[ord]
  
  # Apply row-wise Z-score scaling (if enabled).
  Mz <- if (DO_ZSCORE) row_z(M) else as.matrix(M)
  
  # ann2term labels
  rownames(Mz) <- make.unique(annotate_all(rownames(Mz), rowMeans(Mz, na.rm=TRUE)))
  
  # Outcome annotation bar
  lvl <- if (all(c(pos_class,neg_class) %in% unique(y))) c(pos_class,neg_class) else sort(unique(y))
  y_fac <- factor(y, levels = lvl)
  
  # Create column annotation bar identifying sample groups.
  ha <- HeatmapAnnotation(
    Outcome = y_fac,
    col = list(Outcome = structure(circlize::rand_color(length(levels(y_fac))), names=levels(y_fac)))
  )
  col_fun <- circlize::colorRamp2(c(-2,0,2), c("#2b6cb0","#f7fafc","#c53030"))
  model_name <- gsub("^ML_Plots/trainedModel_|\\.RData$", "", f)
  
  # FULL
  png(file.path("ML_Plots/Heatmaps", paste0(model_name, "_heatmap_FULL.png")),
      width=PNG_W, height=PNG_H, res=PNG_RES)
  draw(Heatmap(Mz, name=if(DO_ZSCORE) "z-score" else "score",
               col=col_fun, top_annotation=ha,
               show_row_names=FALSE, show_column_names=FALSE,
               column_title=paste0(model_name, " (FULL)"),
               use_raster=TRUE, raster_quality=1,
               cluster_columns=FALSE))   # <-- maintains IPD|CONTROL order fixed
  dev.off()
  
  # TOP-N
  Mz_top <- top_var(Mz, TOP_N)
  png(file.path("ML_Plots/Heatmaps", paste0(model_name, "_heatmap_TOP", TOP_N, ".png")),
      width=PNG_W, height=PNG_H, res=PNG_RES)
  draw(Heatmap(Mz_top, name=if(DO_ZSCORE) "z-score" else "score",
               col=col_fun, top_annotation=ha,
               show_row_names=TRUE, row_names_gp=grid::gpar(fontsize=6.5),
               show_column_names=FALSE,
               column_title=paste0(model_name, " (Top-", TOP_N, " var)"),
               use_raster=TRUE, raster_quality=1,
               cluster_columns=FALSE))   # <-- maintains IPD|CONTROL order fixed
  dev.off()
  
  cat("Heatmaps PNG generados (IPD|CONTROL separados) para:", model_name, "\n")
  rm(trainedModel); gc()
}

# Save environment containing visualization objects.
save.image("entorno_TFM_ML_plots") 

# -----------------------------------------------------------------------------
# 7. Stratification across ALL combinations
#
# - Computes ALL pairwise cluster comparisons (C1-C2, C1-C3, etc.).
# - Appends Max_RCSI and Max_Prop_IPD columns to the final summary table.
# - Generates ALL plots and files from the original pipeline.
# - Displays real-time progress bar.
# - Generates dedicated subfolder per combination.
# - Uses doParallel and foreach for speedup.
# -----------------------------------------------------------------------------
# OBJECTIVE: Switch to an UNSUPERVISED learning framework to determine
# whether patient samples naturally group into distinct subtypes (clusters)
# based purely on pathway activity profiles, blind to clinical disease labels.
# This can reveal hidden biological heterogeneity within patient cohorts.
#
# STRATEGY: Use the M3C (Monte Carlo Consensus Clustering) algorithm, a robust
# approach to establish optimal cluster number (K) and assign samples.
# Evaluate ALL score/database combinations in parallel.

# --- 0. LOAD REQUIRED LIBRARIES ---
suppressPackageStartupMessages({
  library(M3C); library(ggplot2); library(ComplexHeatmap); library(circlize)
  library(cluster); library(pathMED); library(stringr); library(foreach)
  library(doParallel); library(dplyr); library(progressr); library(parallel)
})

# --- 1. GENERAL CONFIGURATION ---
pos_class <- "IPD"; neg_class <- "CONTROL"
output_dir <- "ML_Plots/Clustering"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# M3C algorithm parameters. Higher values increase stability at computational expense.
M3C_MAX_K <- 6; M3C_ITERS <- 50; M3C_REPS_REAL <- 100; M3C_REPS_REF <- 100

# --- 2. GENERATE TASK LIST ---
# `expand.grid` produces a table of all method and database combinations.
combinations <- expand.grid(score_method=names(scores_all), gene_set=names(scores_all[[1]]))
cat(sprintf("Se van a procesar %d combinaciones en total.\n", nrow(combinations)))

# --- 3. CLUSTERING ANALYSIS FUNCTION ---
# Define helper functions outside loop
# Normalize pathway IDs (e.g., "GO.123" -> "GO:123")
normalize_ids <- function(x){
  x <- sub("^X(?=\\d)", "", x, perl=TRUE); x <- sub("^GO\\.", "GO:", x)
  x <- sub("^HP\\.", "HP:", x); x <- sub("^R[._-]?HSA[._-]", "R-HSA-", x)
  sub("[._](\\d{3,7})$", "-\\1", x, perl=TRUE)
}
# Row-wise Z-score calculation
row_z <- function(m) {
  mu <- rowMeans(m, na.rm=TRUE)
  sdv <- apply(m, 1, sd, na.rm=TRUE); sdv[sdv==0 | is.na(sdv)] <- 1
  sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
}
# Map pathway IDs to terms
map_ids_to_terms <- function(mat_with_ids, ids){
  ann <- tryCatch(ann2term(mat_with_ids), error=function(e) NULL)
  if (is.null(ann)) return(make.unique(ids))
  if (!"ID" %in% names(ann)) ann$ID <- rownames(ann)
  if (!"term" %in% names(ann)) return(make.unique(ids))
  ann$ID_norm <- normalize_ids(ann$ID)
  term_map <- setNames(ann$term, ann$ID_norm)
  terms <- term_map[normalize_ids(ids)]
  terms[is.na(terms) | !nzchar(terms)] <- ids[is.na(terms) | !nzchar(terms)]
  make.unique(terms)
}

# Main parallelized clustering function
run_clustering_analysis <- function(score_method, gene_set) {
  
  file_tag <- paste(score_method, gene_set, sep = "_")
  combination_dir <- file.path(output_dir, file_tag)
  if (!dir.exists(combination_dir)) dir.create(combination_dir, recursive = TRUE)
  
  # `tryCatch` ensures workflow continuation if an individual combination fails.
  tryCatch({
    mat <- as.matrix(scores_all[[score_method]][[gene_set]])
    if(nrow(mat) < 2 || ncol(mat) < 2) {
      return(data.frame(Combination=file_tag, Status="Skipped", Optimal_K=NA, Max_RCSI=NA, Silhouette_Avg=NA, Separation_PValue=NA, Max_Prop_IPD=NA))
    }
    
    # 1. Run M3C and extract key metrics
    set.seed(123)
    res <- M3C(mydata=mat, maxK=M3C_MAX_K, iters=M3C_ITERS, repsreal=M3C_REPS_REAL, repsref=M3C_REPS_REF, seed=123, removeplots=TRUE)
    df <- res$scores # Contains metrics such as RCSI to identify optimal K.
    clusters <- res$assignments # Cluster assignments per sample.
    max_rcsi_value <- round(max(df$RCSI, na.rm = TRUE), 4) # Highest RCSI score
    # (Relative Cluster Stability Index) obtained across evaluated K values.
    
    # 2. Align data and calculate quality / clinical relevance metrics
    outcome <- as.character(pheno_bin$Response); names(outcome) <- rownames(pheno_bin)
    common <- intersect(colnames(mat), names(outcome))
    mat_sub <- mat[, common, drop=FALSE]
    if (is.null(names(clusters))) names(clusters) <- colnames(mat)
    clusters_sub <- clusters[common]
    outcome_sub  <- outcome[common]
    
    # a) Clinical Relevance: Test whether unsupervised clusters associate
    #    with true IPD/Control diagnoses using statistical testing.
    tab <- table(Cluster = clusters_sub, Outcome = outcome_sub)
    test_res <- if (any(tab < 5)) fisher.test(tab) else chisq.test(tab)
    prop_table <- prop.table(tab, 1)
    max_prop_ipd <- if (pos_class %in% colnames(prop_table)) round(max(prop_table[, pos_class]), 4) else 0
    
    # b) Clustering Quality: Evaluate cluster separation
    #    using mean Silhouette width (values near 1 indicate well-separated clusters).
    dmat <- dist(t(mat_sub), method="euclidean")
    sil <- silhouette(as.integer(factor(clusters_sub)), dmat)
    sil_mean <- mean(sil[, "sil_width"], na.rm=TRUE)
    
    # 3. Generate summary plots and files (RCSI, Heatmap, PCA, etc.)
    
    # RCSI plot and assignments
    png(file.path(combination_dir, paste0("M3C_RCSI_", file_tag, ".png")), width=1200, height=1000, res=150)
    plot(df$K, df$RCSI, type="b", pch=19, xlab="K", ylab="RCSI", main=paste("Selección de K (", file_tag, ")"))
    arrows(df$K, df$RCSI - df$RCSI_SE, df$K, df$RCSI + df$RCSI_SE, angle=90, code=3, length=0.05)
    dev.off()
    write.csv2(clusters, file.path(combination_dir, paste0("M3C_clusters_", file_tag, ".csv")))
    
    # Contingency table and test results
    write.csv2(as.data.frame(tab), file.path(combination_dir, paste0("M3C_", file_tag, "_contingency.csv")))
    capture.output({
      cat("Tabla de contingencia:\n"); print(tab)
      cat("\nTest Estadístico:\n"); print(test_res)
    }, file = file.path(combination_dir, paste0("M3C_", file_tag, "_test.txt")))
    
    # Heatmap
    ord <- order(clusters_sub)
    Mz  <- row_z(mat_sub[, ord, drop = FALSE])
    ha  <- HeatmapAnnotation(Cluster = factor(clusters_sub[ord]))
    col_fun <- colorRamp2(c(-2,0,2), c("#2b6cb0","#f7fafc","#c53030"))
    png(file.path(combination_dir, paste0("M3C_heatmap_", file_tag, ".png")), width=2000, height=1600, res=200)
    draw(Heatmap(Mz, name="z", col=col_fun, top_annotation=ha, show_row_names=FALSE, show_column_names=FALSE, cluster_columns=FALSE, column_title=paste0(file_tag, " — Heatmap")))
    dev.off()
    
    # PCA
    pca <- prcomp(t(mat_sub), scale. = TRUE)
    pc_df <- data.frame(Sample=colnames(mat_sub), PC1=pca$x[,1], PC2=pca$x[,2], Cluster=factor(clusters_sub), Outcome=factor(outcome_sub))
    expl <- summary(pca)$importance[2, 1:2] * 100
    p1 <- ggplot(pc_df, aes(PC1, PC2, color=Cluster)) + geom_point(size=2.2) + theme_minimal(base_size=12) + labs(title="PCA por Cluster", x=sprintf("PC1 (%.1f%%)", expl[1]), y=sprintf("PC2 (%.1f%%)", expl[2]))
    p2 <- ggplot(pc_df, aes(PC1, PC2, color=Outcome)) + geom_point(size=2.2) + theme_minimal(base_size=12) + labs(title="PCA por Outcome", x=sprintf("PC1 (%.1f%%)", expl[1]), y=sprintf("PC2 (%.1f%%)", expl[2]))
    ggsave(file.path(combination_dir, paste0("M3C_PCA_byCluster_", file_tag, ".png")), p1, width=7, height=5, dpi=300)
    ggsave(file.path(combination_dir, paste0("M3C_PCA_byOutcome_", file_tag, ".png")), p2, width=7, height=5, dpi=300)
    
    # Composition barplot
    p_bar <- ggplot(as.data.frame(prop.table(tab, 1)), aes(x = Cluster, y = Freq, fill = Outcome)) + geom_col(width=0.7, color="grey30") + theme_minimal(base_size=12) + scale_y_continuous(labels=scales::percent_format()) + labs(title="Composición por cluster", y="Proporción", x="Cluster")
    ggsave(file.path(combination_dir, paste0("M3C_composition_barplot_", file_tag, ".png")), p_bar, width=6.5, height=4.8, dpi=300)
    
    # Silhouette
    write.csv2(as.data.frame(sil[, 1:3]), file.path(combination_dir, paste0("M3C_silhouette_", file_tag, ".csv")))
    png(file.path(combination_dir, paste0("M3C_silhouette_plot_", file_tag, ".png")), width=1800, height=900, res=200)
    plot(sil, main=sprintf("Silhouette (media = %.3f)", sil_mean), col=2:(length(unique(clusters_sub))+1), border=NA)
    dev.off()
    
    # --- 4. TOP DIFFERENTIALLY ACTIVE PATHWAYS ACROSS ALL PAIRS ---
    # For each cluster pair (e.g., C1 vs C2), test for differentially active pathways.
    grpC <- factor(clusters_sub)
    if (nlevels(grpC) >= 2) {
      
      cluster_pairs <- combn(levels(grpC), 2, simplify = FALSE)
      
      # Apply pairwise Wilcoxon tests across all cluster comparisons.
      all_diff_results <- lapply(cluster_pairs, function(pair) {
        c1 <- pair[1]; c2 <- pair[2]
        idx1 <- which(grpC == c1); idx2 <- which(grpC == c2)
        
        # Perform Wilcoxon rank-sum test per pathway feature.
        pvals <- apply(mat_sub, 1, function(x) wilcox.test(x[idx1], x[idx2])$p.value)
        padj  <- p.adjust(pvals, method="fdr")
        mean1 <- rowMeans(mat_sub[, idx1, drop=FALSE], na.rm=TRUE)
        mean2 <- rowMeans(mat_sub[, idx2, drop=FALSE], na.rm=TRUE)
        
        # Build results data.frame (means, p-values, direction).
        res_pair <- data.frame(
          Comparison = paste0("C", c1, "_vs_C", c2), Pathway_ID = rownames(mat_sub),
          Term = map_ids_to_terms(mat_sub, rownames(mat_sub)),
          Mean_C1 = mean1, Mean_C2 = mean2, Diff = mean1 - mean2,
          pval = pvals, padj = padj,
          Direction = ifelse(mean1 > mean2, paste0("↑ C", c1), paste0("↑ C", c2))
        )
        # Dynamically assign cluster-specific mean column names.
        names(res_pair)[names(res_pair) == "Mean_C1"] <- paste0("Mean_C", c1)
        names(res_pair)[names(res_pair) == "Mean_C2"] <- paste0("Mean_C", c2)
        return(res_pair)
      })
      # Combine all pairwise results into single table.
      res_clust_all <- dplyr::bind_rows(all_diff_results)
      # Sort by comparison and adjusted p-value.
      res_clust_all <- res_clust_all %>% arrange(Comparison, padj)
      
      write.csv2(res_clust_all, file.path(combination_dir, paste0("M3C_topPathways_AllPairs_", file_tag, ".csv")), row.names=FALSE)
      
      # Generate dotplot per cluster comparison pair.
      for (pair in cluster_pairs) {
        c1 <- pair[1]; c2 <- pair[2]
        comp_name <- paste0("C", c1, "_vs_C", c2)
        
        # Filter top 10 most significant pathways.
        top10 <- res_clust_all %>% filter(Comparison == comp_name) %>% head(10)
        
        if (nrow(top10) >= 2) {
          # Dotplot visualizing top 10 differential pathways.
          p_dot <- ggplot(top10, aes(x=Diff, y=reorder(str_wrap(Term, 50), Diff))) +
            geom_vline(xintercept=0, linetype=2, linewidth=0.4) +
            geom_point(aes(size=-log10(padj), color=Diff)) +
            scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
            scale_size_continuous(name="-log10(FDR)") +
            labs(title=paste("Top-10 Rutas:", comp_name), x=paste0("Diferencia (C",c1," - C",c2,")"), y="Ruta") +
            theme_minimal(base_size=14) + theme(plot.title=element_text(hjust=0.5, face="bold"))
          
          ggsave(file.path(combination_dir, paste0("M3C_topPathways_dotplot_", comp_name, "_", file_tag, ".png")), p_dot, width=9, height=6, dpi=300)
        }
      }
    }
    
    # 5. Return summary metrics row
    return(data.frame(
      Combination=file_tag, Status="Completed", Optimal_K=length(unique(clusters)),
      Max_RCSI=max_rcsi_value, Silhouette_Avg=round(sil_mean, 4),
      Separation_PValue=test_res$p.value, Max_Prop_IPD=max_prop_ipd
    ))
  }, error = function(e) {
    return(data.frame(Combination=file_tag, Status=paste("Error:", e$message), Optimal_K=NA, Max_RCSI=NA, Silhouette_Avg=NA, Separation_PValue=NA, Max_Prop_IPD=NA))
  })
}

# --- 4. PARALLEL EXECUTION WITH PROGRESS TRACKING ---
# Configure parallel cluster using N-1 CPU cores.
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("\nIniciando procesamiento paralelo con %d núcleos...\n", num_cores))
with_progress({
  # Parallel loop via `foreach` executing `run_clustering_analysis`.
  # `.combine = 'rbind'` binds metric summary rows.
  p <- progressor(steps = nrow(combinations))
  summary_results <- foreach(i = 1:nrow(combinations), .combine = 'rbind', .packages = c("M3C", "ggplot2", "ComplexHeatmap", "circlize", "cluster", "pathMED", "stringr", "dplyr", "parallel")) %dopar% {
    # Update progress bar per iteration.
    p(sprintf("Procesando %s_%s", combinations$score_method[i], combinations$gene_set[i])) 
    run_clustering_analysis(combinations$score_method[i], combinations$gene_set[i])
  }
})
# Stop parallel cluster to release compute resources.
stopCluster(cl)
cat("\n--- PROCESO PARALELO FINALIZADO ---\n")

# --- 5. FINAL RESULTS EXPORT ---
if (!is.null(summary_results) && nrow(summary_results) > 0) {
  # Sort final results by Silhouette score descending.
  summary_results <- summary_results %>% arrange(desc(Silhouette_Avg))
  print("Tabla Resumen de Resultados de Clustering:")
  print(summary_results)
  write.csv2(summary_results, file.path(output_dir, "CLUSTERING_SUMMARY_RESULTS.csv"), row.names = FALSE)
  cat(sprintf("\nTabla resumen guardada en: %s\n", file.path(output_dir, "CLUSTERING_SUMMARY_RESULTS.csv")))
} else {
  cat("No se generaron resultados. Revisa los mensajes de error.\n")
}

# Save complete workspace containing all final objects and results.
save.image("entorno_TFM_TODO")




sessionInfo()
