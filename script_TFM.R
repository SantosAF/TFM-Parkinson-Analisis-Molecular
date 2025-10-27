# ----------------------------------------------------------------------------------
# SCRIPT TFM:
# Análisis de firmas moleculares para la estratificación y la predicción en la enfermedad de Parkinson
# Santos Antequera Fernández - Máster en Bioinformática - Universidad Europea de Madrid
# ----------------------------------------------------------------------------------

# RESUMEN DEL SCRIPT:
# 1. Preparación del entorno, descarga y limpieza de datos (GSE99039).
# 2. Cálculo de puntuaciones de actividad de rutas biológicas (pathway scores).
# 3. Entrenamiento de modelos de Machine Learning para predicción clínica.
# 4. Visualización y evaluación del rendimiento de los modelos.
# 5. Análisis de clustering no supervisado para identificar subtipos de pacientes.

# ----------------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# SECCIÓN 1: CONFIGURACIÓN INICIAL Y CARGA DE DATOS
# -----------------------------------------------------------------------------

## 1.1. Directorio de Trabajo y Entorno de R
# -----------------------------------------------------------------------------
# NOTA: Esta ruta es local y deberá ser modificada por otros usuarios.
setwd("C:/Users/Usuario/Desktop/Máster Bioinformática/Módulo 10 TFM/Contenido")
# Para agilizar el trabajo, podemos cargar un entorno de R previamente guardado.
# Esto nos permite reanudar el análisis sin tener que ejecutar de nuevo los
# cálculos más largos y pesados.
# load("entorno_TFM") # Carga inicial de scores sin filtrar.
# load("entorno_TFM_finalML2") # Carga de los resultados del Machine Learning.
# load("entorno_TFM_ML_plots") # Carga de los objetos para las gráficas.
load("entorno_TFM_TODO") # Carga el entorno completo final.


## 1.2. Instalación y Carga de Librerías
# -----------------------------------------------------------------------------

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# Instalamos y cargamos 'pathMED', la librería principal para el cálculo
# de puntuaciones de rutas biológicas.
# BiocManager::install("pathMED")
library(pathMED)
### Documentación y ejemplos:
# browseVignettes("pathMED")



# 1.3. Descarga de Datos desde GEO usando el paquete GEOquery
# -----------------------------------------------------------------------------

#if (!require("GEOquery")) {
#  BiocManager::install("GEOquery")
#}

library(GEOquery) # Este paquete nos permite descargar datos públicos desde el repositorio 
#Gene Expression Omnibus (GEO) de NCBI

### Guardamos los datos en la variable gse
gse <- getGEO("GSE99039", GSEMatrix = TRUE)

# Ver cuántos objetos se descargaron
length(gse)

# Accedemos al objeto:
gse_data <- gse[[1]]

## Extraemos la matriz de expresión y metadatos

# Matriz de expresión normalizada (genes x muestras) 
exprs_matrix <- exprs(gse_data) # No necesitamos transponer la matriz, es el 
# formato perfecto para la función getScores()

# Datos fenotípicos (información de las muestras)
pheno_data <- pData(gse_data)

# Datos de anotación de genes 
feature_data <- fData(gse_data)




# -----------------------------------------------------------------------------
# SECCIÓN 2: PREPROCESAMIENTO Y LIMPIEZA DE LA MATRIZ DE EXPRESIÓN
# -----------------------------------------------------------------------------
# OBJETIVO: Transformar la matriz de expresión para que las filas representen
# genes únicos con símbolos estándar (ej. "TP53"), en lugar de IDs de sondas
# ambiguos (ej. "1053_at"). Este paso es fundamental para el análisis de rutas.

# Paso 1. # Usamos los datos de anotación (feature_data) para traducir los IDs de las sondas
# (rownames de exprs_matrix) a símbolos de gen oficiales.
gene_symbols <- feature_data$`Gene Symbol`[match(rownames(exprs_matrix), feature_data$ID)]

# Paso 2. Asignarlos como nombres de fila en la matriz de expresión
rownames(exprs_matrix) <- gene_symbols

# Paso 3. Eliminar las filas sin nombre de gen o vacías
exprs_matrix <- exprs_matrix[!is.na(rownames(exprs_matrix)) & rownames(exprs_matrix) != "", ]

# Paso 4. Agrupar genes duplicados por símbolo
# Varias sondas pueden mapear al mismo gen. Para evitar errores en métodos
# que requieren nombres únicos, los agregamos por símbolo de gen.
# Aquí usamos la mediana como medida más robusta frente a valores atípicos.

# NOTA CRÍTICA: Esta aproximación (media o mediana) es válida para exploración,
# pero existen métodos más precisos, como:
#  - seleccionar la sonda con mayor varianza (más informativa),
#  - usar anotaciones de sondas principales,
#  - o aplicar modelos ponderados según calidad de sonda.
# Aquí usamos la mediana por simplicidad y compatibilidad con funciones como `getScores()`.

# Validación antes de agrupar
duplicated_genes <- sum(duplicated(rownames(exprs_matrix)))
cat("Genes duplicados a combinar:", duplicated_genes, "\n") 

# Agrupar usando la mediana
exprs_matrix <- aggregate(exprs_matrix, by = list(Gene = rownames(exprs_matrix)), FUN = median)
rownames(exprs_matrix) <- exprs_matrix$Gene
exprs_matrix$Gene <- NULL # Eliminamos la columna auxiliar

# Validación después de agrupar (debería salir 0)
duplicated_genes_after <- sum(duplicated(rownames(exprs_matrix)))
cat("Genes duplicados después de agrupar:", duplicated_genes_after, "\n")





# -----------------------------------------------------------------------------
# SECCIÓN 3: CÁLCULO DE PUNTUACIONES DE ACTIVIDAD DE RUTAS (SCORES MOLECUALRES)
# ----------------------------------------------------------------------------- 
# OBJETIVO: Transformar los datos de expresión a nivel de gen en puntuaciones
# a nivel de ruta biológica. Esto reduce la dimensionalidad de los datos y
# permite interpretar los resultados en un contexto biológico funcional.
# Probamos una combinación de 6 métodos de puntuación y 7 bases de datos de
# rutas para identificar las señales biológicas más robustas.

# 1. Lista de geneSets que quieres usar
gene_sets <- c("kegg", "reactome", "go_bp", "go_mf", "go_cc", "disgenet", "hpo")

# 2.1. Método singscore. Bucle para calcular y guardar los resultados
for (gs in gene_sets) {
  scores <- getScores(exprs_matrix, geneSets = gs, method = "singscore",cores = 10)
  assign(paste0("scores_", gs), scores)  # Crea la variable: scores_kegg, etc.
  cat("\nPrimeros valores de", gs, ":\n")
  print(scores[1:5, 1:5])  # Muestra los primeros valores
}

# 2.2. Método Z-score. Bucle para calcular y guardar los resultados
for (gs in gene_sets) {
  zscores <- getScores(exprs_matrix, geneSets = gs, method = "Z-score",cores = 10)
  assign(paste0("zscores_", gs), zscores)  # Crea la variable: zscores_kegg, etc.
  cat("\n[zScore] Primeros valores de", gs, ":\n")
  print(zscores[1:5, 1:5])
}

# 2.3. Método GSVA. Bucle para calcular y guardar los resultados
for (gs in gene_sets) {
  gsvascores <- getScores(exprs_matrix, geneSets = gs, method = "GSVA",cores = 10)
  assign(paste0("gsvascores_", gs), gsvascores)  # Crea la variable: gsvascores_kegg, etc.
  cat("\n[gsvaScore] Primeros valores de", gs, ":\n")
  print(gsvascores[1:5, 1:5])
}

# 2.4. Método ssGSEA. Bucle para calcular y guardar los resultados
for (gs in gene_sets) {
  ssgscores <- getScores(exprs_matrix, geneSets = gs, method = "ssGSEA",cores = 10)
  assign(paste0("ssgscores_", gs), ssgscores)  # Crea la variable: ssgscores_kegg, etc.
  cat("\n[ssGSEAScore] Primeros valores de", gs, ":\n")
  print(ssgscores[1:5, 1:5])
}

# 2.5. Método Plage. Bucle para calcular y guardar los resultados
for (gs in gene_sets) {
  plagescores <- getScores(exprs_matrix, geneSets = gs, method = "Plage",cores = 10)
  assign(paste0("plagescores_", gs), plagescores)  # Crea la variable: plagescores_kegg, etc.
  cat("\n[PlageScore] Primeros valores de", gs, ":\n")
  print(plagescores[1:5, 1:5])
}

# 2.6. Método norm_FGSEA. Bucle para calcular y guardar los resultados
# Primero necesitamos instalar y cargar el paquete fgsea para que funcione
# Nota: Este método necesita más tiempo para obtener los scores.

# BiocManager::install("fgsea")
library(fgsea)

for (gs in gene_sets) {
  nFscores <- getScores(exprs_matrix, geneSets = gs, method = "norm_FGSEA",cores = 10)
  assign(paste0("nFscores_", gs), nFscores)  # Crea la variable: nFscores_kegg, etc.
  cat("\n[norm_FGSEAScore] Primeros valores de", gs, ":\n")
  print(nFscores[1:5, 1:5])
}



# 3. Para no repetir todo el proceso anterior si queremos ver los primeros scores, hacemos lo siguiente:
# Nombres base de tus objetos 
score_names <- c("kegg", "reactome", "go_bp", "go_mf", "go_cc", "disgenet", "hpo")

# Bucle para imprimir 5x5 de cada scores_*
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


#### Para anotar los identificadores con sus descripciones correspondientes, usamos la función ann2term()
## Primero para singscore 
# 1. Crear lista de scores
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


## Después para Z-score
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


## Después para GSVA
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


## Después para ssGSEA
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

## Después para Plage
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


## Después para norm_FGSEA
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




save.image("entorno_TFM") # Este entorno contiene todos los scores calculados.
# En caso de tener que repetir la siguiente parte de ML se cargaría este entorno
# para poder hacer todo de manera correcta desde el inicio, sin objetos o librerías
# sobrantes



# -----------------------------------------------------------------------------
# SECCIÓN 4: REDUCCIÓN DE NÚMERO DE RUTAS PARA EL ENTRENAMIENTO DE MODELOS DE MACHINE LEARNING
# -----------------------------------------------------------------------------

## Ejemplo de lo que ocupa cada base de datos
sapply(list(scores_kegg, scores_reactome, scores_go_bp, scores_go_mf, scores_go_cc, scores_hpo, scores_disgenet), nrow)
# [1]   345  2483 12377  4398  1796  8698 10730
# Escoger 750 está bien pero lo suyo sería hacer un filtro dinámico, sin embargo, debido a las limitaciones
# de memoria se opta por esta opción, ya que si no el tiempo de calculo sería demasiado largo
# (en anteriores pruebas para calcular 1 o dos modelos de go_bp se tardaron varias horas)

## Función para filtrar por varianza (top N)
filter_top_var <- function(score_matrix, topN = 750) {
  n_before <- nrow(score_matrix)
  
  # Limpieza. Eliminar filas con valores no finitos (NA, Inf, -Inf)
  # que podrían causar errores en el cálculo de la varianza.
  clean_mat <- score_matrix[complete.cases(score_matrix), ]
  clean_mat <- clean_mat[apply(clean_mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
  n_after <- nrow(clean_mat)
  cat("Eliminadas", n_before - n_after, "features con NA/Inf (de", n_before, ")\n")
  
  # Calcular la varianza para cada ruta (fila).
  varianceScores <- apply(clean_mat, 1, function(x) var(x, na.rm = TRUE))
  
  # Ordenar las rutas de mayor a menor varianza.
  varianceScores <- sort(varianceScores, decreasing = TRUE)
  
  # Seleccionar los nombres de las 'topN' rutas más variables.
  nFeatures <- min(length(varianceScores), topN)
  selected <- names(varianceScores)[1:nFeatures]
  
  # Devolver la submatriz con solo las rutas seleccionadas.
  return(clean_mat[selected, , drop = FALSE])
}


## Lista completa de scores, aplicando la función de reducción
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
# SECCIÓN 5: CLASIFICACIÓN SUPERVISADA CON pathMED
# -----------------------------------------------------------------------------
# OBJETIVO: Entrenar y evaluar modelos capaces de predecir la condición de un
# sujeto (IPD o Control) a partir de sus perfiles de actividad de rutas.

# --------------------------
# 1) Variable respuesta
# --------------------------

# Extraemos la variable clínica de interés del `pheno_data`.
pheno_data$Response <- as.character(pheno_data$`disease label:ch1`)  # clave: character para pathMED

# Definimos explícitamente la clase positiva ("IPD") y negativa ("CONTROL").
pos_class <- "IPD"
neg_class <- "CONTROL"

# Filtramos la tabla de metadatos para quedarnos solo con las muestras
# que pertenecen a estos dos grupos.
pheno_bin <- pheno_data[pheno_data$Response %in% c(pos_class, neg_class), , drop=FALSE]
# Comprobamos el balance de clases.
table(pheno_bin$Response)


# --------------------------
# PCA exploratorio: IPD vs Control
# --------------------------
# El Análisis de Componentes Principales (PCA) nos permite visualizar la estructura
# de los datos en 2D. Es un primer paso para evaluar si las muestras de los dos
# grupos (IPD, Control) son separables basándose en los scores de las rutas.
library(ggplot2)

# --- 1. Crear una carpeta para guardar todos los plots ---
# Crear carpeta si no existe
if(!dir.exists("ML_Plots")) dir.create("ML_Plots")
pca_output_dir <- "ML_Plots/PCA_Exploratorio"
if (!dir.exists(pca_output_dir)) dir.create(pca_output_dir, recursive = TRUE)

cat("Iniciando la generación de plots PCA para todas las combinaciones...\n")

# --- 2. Bucle anidado para generar un PCA para cada combinación de método y base de datos.
for (method in names(scores_all)) {
  for (db in names(scores_all[[method]])) {
    
    file_tag <- paste(method, db, sep = "_")
    cat(sprintf("Generando PCA para: %s\n", file_tag))
    
    # a) Escoger la matriz de scores dinámicamente
    mat_PCA <- scores_all[[method]][[db]]
    
    # b) Alineación de muestras: nos aseguramos de que la matriz de scores y
    #    los metadatos contienen exactamente las mismas muestras en el mismo orden.
    common_samples_PCA <- intersect(colnames(mat_PCA), rownames(pheno_bin))
    mat_PCA_sub <- mat_PCA[, common_samples_PCA]
    pheno_bin_PCA <- pheno_bin[common_samples_PCA, ]
    
    # Salta si no hay suficientes datos para el PCA
    if(ncol(mat_PCA_sub) < 3) {
      cat(sprintf("  -> Saltando %s: no hay suficientes muestras.\n", file_tag))
      next
    }
    
    # c) Cálculo del PCA. Se transpone la matriz (t()) porque prcomp espera
    #    muestras en las filas y características en las columnas.
    pca <- prcomp(t(mat_PCA_sub), scale. = TRUE)
    
    # d) Creación de un data.frame para la visualización con ggplot2.
    pca_df <- data.frame(
      Sample = rownames(pca$x),
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      Group = pheno_bin_PCA$Response
    )
    
    # e) Creación del gráfico, incluyendo el porcentaje de varianza explicada
    #    por cada componente principal en los ejes.
    pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 3, alpha = 0.8) +
      labs(
        title = paste("PCA (CONTROL vs IPD) -", file_tag),
        x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")
      ) +
      theme_minimal()
    
    # f) Guardar el plot con un nombre de archivo único
    ggsave(
      filename = file.path(pca_output_dir, paste0("PCA_", file_tag, ".png")),
      plot = pca_plot,
      width = 7, height = 5, dpi = 300
    )
  }
}

cat("\n--- Plots PCA generados y guardados en la carpeta 'PCA_Exploratorio' ---\n")


# --------------------------
# 2) Definición de los Modelos de Machine Learning
# --------------------------
modelsList <- methodsML(
  algorithms = c("rf","knn","xgbTree","glm","lda"), 
  outcomeClass = "character",   # obligatorio en pathMED
  tuneLength = 10  # Lo ideal sería poner más pero por las limitaciones computacionales y de tiempo se escogio 10            
)

# --------------------------
# 3) Entrenar modelos
# --------------------------
# Este es el núcleo del análisis de Machine Learning.
# Cargamos las librerías necesarias para los algoritmos y la paralelización.

library(pathMED)
library(caret)        
library(doParallel)   # paralelización
library(randomForest) # para rf
library(kknn)         # para knn mejorado
library(xgboost)      # para xgbTree
library(MASS)

# Base seed fija
base_seed <- 1234

# Para hacerlo desde el principio
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
# Crear carpeta si no existe
if(!dir.exists("ML_Plots")) dir.create("ML_Plots")


# Configuramos un clúster de computación paralela para usar 10 núcleos de la CPU.
# Esto acelera drásticamente el proceso de entrenamiento y validación cruzada.
# con la librería doParallel
cl <- makeCluster(10)   
registerDoParallel(cl)

# Bucle principal para el entrenamiento de los modelos
for (method in names(scores_all)) {
  for (db in names(scores_all[[method]])) {
    
    # Alinear muestras correctamente entre input_mat y meta_bin
    common_samples <- intersect(colnames(scores_all[[method]][[db]]), rownames(pheno_bin))
    if(length(common_samples) < 2){
      cat("Saltando", method, db, "- no hay suficientes muestras\n")
      next
    }
    
    # Creamos las matrices finales para el entrenamiento, usando solo las muestras comunes.
    input_mat <- scores_all[[method]][[db]][, common_samples, drop=FALSE]
    meta_bin <- pheno_bin[common_samples, , drop=FALSE]
    
    # Reordenar meta_bin según input_mat por seguridad
    # Esto evita errores de asignación de etiquetas a las muestras.
    meta_bin <- meta_bin[colnames(input_mat), , drop=FALSE]
    
    # Comprobamos que, tras el filtrado, todavía tenemos datos de las dos clases
    # (IPD y Control) para poder entrenar un modelo de clasificación.
    clase_check <- table(meta_bin$Response)
    if (length(clase_check) < 2) {
      cat("Saltando", method, db, "- solo hay una clase presente después de alinear\n")
    } else {
      cat("Alineación correcta para", method, db, "-", 
          clase_check[pos_class], pos_class, "y", clase_check[neg_class], neg_class, "\n")
    }
    
    
    # Establecemos una semilla única para cada combinación.
    # Esto garantiza que si volvemos a ejecutar el script, los resultados de la
    # validación cruzada (que tiene un componente aleatorio en la división de datos)
    # serán exactamente los mismos.
    seed_iter <- base_seed + match(method, names(scores_all)) * 100 + match(db, names(scores_all[[method]]))
    set.seed(seed_iter)
    
    # Confirmación de semilla y combinación
    cat("Entrenando:", method, "-", db, "con semilla:", seed_iter, "\n")
    
    # Esta es la función principal que entrena y evalúa los modelos definidos.
    # Utiliza una validación cruzada anidada para obtener una estimación robusta
    # y fiable del rendimiento del modelo en datos nuevos
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
    
    # De todos los modelos probados (rf, knn, etc.), pathMED elige el mejor.
    # Extraemos sus métricas de rendimiento.
    stats <- trainedModel$stats
    best  <- trainedModel$model$method
    
    #Guardar métricas
    results_df <- rbind(
      results_df,
      data.frame(
        ScoreMethod = method,
        Database = db,
        BestModel = best,
        Accuracy = stats["accuracy", best],
        BalancedAcc = stats["balacc", best],
        MCC = (stats["mcc", best] + 1) / 2,  # MCC normalizado a [0,1]
        Recall = stats["recall", best],
        Specificity = stats["specificity", best],
        Precision = stats["precision", best],
        F1 = stats["fscore", best]
      )
    )
    
    # Guardamos el objeto 'trainedModel' completo. Es un objeto grande que contiene
    # el modelo final, predicciones, etc. Lo necesitaremos más tarde para los gráficos.
    save(trainedModel, file = paste0("ML_Plots/trainedModel_", method, "_", db, ".RData"))
    
    # El objeto 'trainedModel' puede ocupar mucha memoria RAM, por ello lo eliminamos
    # cada vez que termina con una de las combinaciones
    rm(trainedModel)
    gc()
    
    print(results_df)
    
    # Guardar la tabla de resultados después de cada iteración
    save(results_df, file="results_todo_ML2.RData")
  }
}

# Una vez terminado el bucle, detenemos el clúster de computación paralela
# para liberar los núcleos de la CPU.
stopCluster(cl)


# Guardar la tabla de resultados después de cada iteración en un archivo .csv
write.csv(results_df, "results_todo_ML2.csv", row.names = FALSE)



# Tiempo estimado con aquellos metodos con kegg como base de datos: 2 minutos (valor para los primeros casos)
# Tiempo estimado para los que se les aplicó el filtro de 750 scores: 5 minutos (valor para los primeros casos)

save.image("entorno_TFM_finalML2") # Este entorno contiene los scores  
# y los resultados de Machine Learning con pathMED






# =============================================================================
# SECCIÓN 6: VISUALIZACIÓN Y ANÁLISIS DE RESULTADOS DE MACHINE LEARNING
# =============================================================================
# OBJETIVO: Una vez entrenados los modelos, necesitamos interpretar y visualizar
# su rendimiento para poder sacar conclusiones. En esta sección generaremos
# una serie de gráficos que nos ayudarán a comparar los diferentes enfoques
# (métodos de score y bases de datos) y a entender en detalle el comportamiento
# de cada modelo entrenado.


# -----------------------------------------------------------------------------
# 1. Heatmap de métricas por ScoreMethod y Database 
# -----------------------------------------------------------------------------
# Este heatmap nos proporciona una vista de pájaro del rendimiento de todas la combinaciones 
# probadas. 
library(reshape2)
library(ggplot2)

metrics <- c("Accuracy","BalancedAcc","MCC","Recall","Specificity","Precision","F1")

# La función `melt` transforma nuestra tabla `results_df` de un formato "ancho"
# (una columna por cada métrica) a un formato "largo" (una columna para el
# nombre de la métrica y otra para su valor). Este formato es el que necesita
# ggplot para crear los paneles con `facet_wrap`.
df_melt <- melt(results_df, id.vars = c("ScoreMethod","Database"), measure.vars = metrics)

# Creamos el gráfico con ggplot2
p_heat <- ggplot(df_melt, aes(x = Database, y = ScoreMethod, fill = value)) +
  geom_tile() +
  facet_wrap(~variable) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Heatmap de métricas por ScoreMethod y Database", fill="Valor")

ggsave("ML_Plots/Heatmap_metrics.png", plot=p_heat, width=10, height=6, dpi=300)


# -----------------------------------------------------------------------------
# 2. ROC y PR de los mejores modelos (con subsample.preds)
# -----------------------------------------------------------------------------
# Estas curvas evalúan la capacidad de discriminación de cada modelo.
# Se generan usando las predicciones "Out-of-Fold" (OOF), que son las
# predicciones que el modelo hace sobre datos que no ha visto durante su
# entrenamiento en la validación cruzada. Esto nos da una estimación muy
# realista y fiable de cómo se comportará el modelo con datos nuevos.

library(pROC)
library(PRROC)

# Crear carpeta de salida si no existe
if (!dir.exists("ML_Plots/ROC_PR")) dir.create("ML_Plots/ROC_PR", recursive = TRUE)

# Definir clases explícitamente
pos_class <- "IPD"
neg_class <- "CONTROL"

for (i in 1:nrow(results_df)) {
  method <- results_df$ScoreMethod[i]
  db     <- results_df$Database[i]
  
  # Cargamos el objeto `trainedModel` que guardamos durante el entrenamiento.
  load(paste0("ML_Plots/trainedModel_", method, "_", db, ".RData"))  
  
  # Extraemos las predicciones OOF ("Out-of-Fold").
  oof <- trainedModel$subsample.preds
  if (is.null(oof)) { 
    cat("Sin subsample.preds para", method, db, "\n") 
    next 
  }
  
  # Comprobar columnas necesarias (IPD, CONTROL y obs)
  if (!all(c(pos_class, neg_class, "obs") %in% colnames(oof))) {
    cat("Esperaba columnas", pos_class, ",", neg_class, "y obs en", method, db, "\n")
    cat("Columnas disponibles:", paste(colnames(oof), collapse = ", "), "\n")
    next
  }
  
  # Preparamos los datos: 'obs' es la verdad (lo que realmente era el paciente)
  # y 'prob' es la probabilidad que el modelo asignó a que fuera de la clase positiva (IPD).
  obs  <- factor(oof$obs, levels = c(neg_class, pos_class))
  prob <- oof[[pos_class]]  # probabilidad de la clase positiva (IPD)
  
  # Comprobar que hay ambas clases
  if (length(unique(obs)) < 2) {
    cat("Saltando ROC/PR para", method, db, "- solo una clase presente\n")
    next
  }
  
  # ========================
  # Curva ROC.
  # Muestra el balance entre la Tasa de Verdaderos Positivos (sensibilidad) y la
  # Tasa de Falsos Positivos (1 - especificidad). Un modelo perfecto tendría un
  # área bajo la curva (AUC) de 1.
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
  # Curva PR.
  # Es especialmente informativa cuando las clases están desbalanceadas. Muestra
  # el balance entre la Precisión (qué porcentaje de positivos predichos son correctos)
  # y el Recall (qué porcentaje de los positivos reales se encontraron).
  # ========================
  pr_obj <- pr.curve(
    scores.class0 = prob[obs == pos_class],
    scores.class1 = prob[obs == neg_class], # Nota: PRROC invierte las clases
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
# 3. Matrices de confusión 
# -----------------------------------------------------------------------------
# OBJETIVO: Analizar en detalle los errores del modelo. La matriz de confusión
# es una tabla que nos dice exactamente cuántos pacientes fueron clasificados
# correctamente (Verdaderos Positivos/Negativos) y cuántos incorrectamente
# (Falsos Positivos/Negativos). De nuevo, usamos las predicciones OOF para
# obtener una visión realista del rendimiento.

library(caret)

# Carpeta para guardar plots si no existe
if (!dir.exists("ML_Plots/Confusion_matrix")) dir.create("ML_Plots/Confusion_matrix", recursive = TRUE)

# Definir las clases explícitamente (ya hecho anteriormente)
pos_class <- "IPD"
neg_class <- "CONTROL"

for (i in 1:nrow(results_df)) {
  
  method_imp <- results_df$ScoreMethod[i]
  db_imp     <- results_df$Database[i]
  
  # === 1) Cargar modelo entrenado ===
  model_file <- paste0("ML_Plots/trainedModel_", method_imp, "_", db_imp, ".RData")
  if (!file.exists(model_file)) {
    cat("Modelo no encontrado para", method_imp, db_imp, "\n")
    next
  }
  load(model_file)  # Esto restaura el objeto 'trainedModel' en el entorno.
  
  # === 2) Extraemos las predicciones OOF (subsample) ===
  oof <- trainedModel$subsample.preds
  if (is.null(oof) || !all(c(pos_class, neg_class, "obs") %in% colnames(oof))) {
    cat("subsample.preds no disponible o sin columnas esperadas (", 
        pos_class, ",", neg_class, ", obs ) para", method_imp, db_imp, "\n")
    next
  }
  
  # === 3) Verdades y clases predichas desde prob OOF ===
  # El modelo nos da una probabilidad (ej. 0.8 de ser IPD). Para la matriz
  # de confusión, necesitamos una decisión final ('IPD' o 'CONTROL').
  # Usamos un umbral estándar de 0.5: si la probabilidad es >= 50%, predecimos 'IPD'.
  obs  <- factor(oof$obs, levels = c(neg_class, pos_class))
  prob <- oof[[pos_class]]
  pred <- factor(ifelse(prob >= 0.5, pos_class, neg_class), levels = c(neg_class, pos_class))
  
  # Comprobar ambas clases
  if (length(unique(obs)) < 2) {
    cat("Saltando Confusion Matrix para", method_imp, db_imp, "- solo hay una clase en OOF\n")
    next
  }
  
  # === 4) Matriz de confusión (OOF) ===
  cm <- caret::confusionMatrix(pred, obs, positive = pos_class)
  
  # Mostrar en consola
  cat("\nMatriz de confusión (OOF) para", method_imp, db_imp, ":\n")
  print(cm)
  
  # === 5) Guardar plot (fourfoldplot) ===
  png(paste0("ML_Plots/Confusion_matrix/ConfMatrix_", method_imp, "_", db_imp, ".png"),
      width = 800, height = 600)
  # `fourfoldplot` es una forma gráfica de ver la matriz. El área de cada
  # cuadrante es proporcional al número de muestras en esa categoría.
  fourfoldplot(cm$table,
               color = c("#FF9999", "#99CCFF"),
               conf.level = 0,
               margin = 1,
               main = paste("Confusion Matrix (OOF):", method_imp, "-", db_imp))
  dev.off()
  
  cat("Matriz de confusión generada para", method_imp, db_imp, "\n")
}

# -----------------------------------------------------------------------------
# 4. Gráficos de Calibración
# -----------------------------------------------------------------------------
# OBJETIVO: Evaluar si las probabilidades que predice el modelo son fiables.
# Un modelo bien calibrado es aquel en el que si predice un 80% de probabilidad
# para un grupo de muestras, aproximadamente el 80% de esas muestras
# pertenecen realmente a la clase positiva. Esto es crucial para la confianza
# en las predicciones del modelo.

library(dplyr)
library(ggplot2)

if (!dir.exists("ML_Plots/Calibration_plot")) dir.create("ML_Plots/Calibration_plot", recursive = TRUE)

# Definir clases explícitamente (realizado anteriormente)
pos_class <- "IPD"
neg_class <- "CONTROL"

for (i in 1:nrow(results_df)) {
  method <- results_df$ScoreMethod[i]
  db     <- results_df$Database[i]
  
  # Cargar modelo entrenado 
  model_file <- paste0("ML_Plots/trainedModel_", method, "_", db, ".RData")
  if (!file.exists(model_file)) {
    cat("Modelo no encontrado para", method, db, "\n")
    next
  }
  load(model_file)
  
  # Extraemos las predicciones OOF.
  oof <- trainedModel$subsample.preds
  if (is.null(oof) || !all(c(pos_class, neg_class, "obs") %in% colnames(oof))) {
    cat("Sin OOF (subsample.preds) con columnas", pos_class, "/", neg_class, "/obs para", method, db, "\n")
    if (!is.null(oof)) cat("Columnas disponibles:", paste(colnames(oof), collapse = ", "), "\n")
    next
  }
  
  # 1) Preparar los datos para el análisis ===
  # Necesitamos la verdad en formato numérico (0 para CONTROL, 1 para IPD)
  # para poder calcular medias y errores.
  obs_fac <- factor(oof$obs, levels = c(neg_class, pos_class))
  if (length(unique(obs_fac)) < 2) {
    cat("Saltando Calibration para", method, db, "- solo una clase presente en OOF\n")
    next
  }
  prob    <- oof[[pos_class]]                # prob de IPD (clase positiva)
  obs_num <- as.integer(obs_fac == pos_class) # 1=IPD, 0=CONTROL
  
  # 2) Brier score (calibración global)
  # Es una métrica global que mide la calibración. Es el error cuadrático medio
  # entre las probabilidades predichas y los valores reales (0 o 1).
  # Un valor más cercano a 0 indica una mejor calibración.
  brier <- mean((prob - obs_num)^2, na.rm = TRUE)
  
  # 3) Curva de calibración por bins (deciles por defecto)
  #    Salvaguarda si hay pocos valores únicos
  unique_probs <- sum(!is.na(unique(prob)))
  if (unique_probs < 3) {
    # Con muy poca variación, usar bins fijos
    cuts <- seq(0, 1, by = 0.2)
  } else {
    cuts <- quantile(prob, probs = seq(0, 1, by = 0.1), na.rm = TRUE, names = FALSE)
    # Evitar límites duplicados por empates
    if (any(duplicated(cuts))) {
      cuts <- unique(cuts)
      # garantizar incluye 0 y 1
      if (min(cuts) > 0) cuts <- c(0, cuts)
      if (max(cuts) < 1) cuts <- c(cuts, 1)
    }
  }
  # Si por algún motivo quedaron <3 cortes, asegurar al menos 3 para formar bins
  if (length(cuts) < 3) cuts <- unique(sort(c(0, cuts, 1)))
  
  # Agrupamos las predicciones en "cajones" (bins) de probabilidad (ej. todas
  # las predicciones entre 0.1 y 0.2). Para cada cajón, calculamos la
  # probabilidad media predicha y la proporción real de casos positivos.
  df_cal <- data.frame(prob = prob, obs = obs_num) %>%
    mutate(bin = cut(prob, breaks = cuts, include.lowest = TRUE, right = TRUE)) %>%
    group_by(bin, .drop = TRUE) %>%
    summarise(
      pred_mean = mean(prob, na.rm = TRUE), # Probabilidad media predicha en el bin
      obs_rate  = mean(obs,  na.rm = TRUE), # Proporción real de positivos en el bin
      n         = dplyr::n(),
      .groups   = "drop"
    ) %>%
    filter(is.finite(pred_mean), is.finite(obs_rate))
  
  # 4) Plot y guardado
  # Un modelo perfectamente calibrado tendría todos sus puntos sobre la línea
  # diagonal (donde la probabilidad predicha es igual a la proporción observada).
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
# 5. Importancia de variables + ann2term en plots
#   Exporta ALL.csv, top20.csv y plot top5
#   (con logs de progreso y barra para KNN)
# -----------------------------------------------------------------------------
# OBJETIVO: Descubrir qué rutas (características) son las más importantes para
# que cada modelo pueda distinguir entre pacientes IPD y Controles.
# Este paso es fundamental para la interpretabilidad biológica: no solo
# queremos saber si el modelo funciona, sino también CÓMO funciona.

# Cargar librerías necesarias
library(caret)
library(pathMED)
library(ggplot2)
library(dplyr)
library(stringr)

#--- Función de logging ---------------------------------------------------------
VERBOSE <- TRUE
logf <- function(fmt, ...) {
  if (!VERBOSE) return(invisible(NULL))
  cat(format(Sys.time(), "%H:%M:%S"), "-", sprintf(fmt, ...), "\n")
  flush.console()
}

# Función para normalizar los IDs de las rutas (ej. "GO.123" -> "GO:123")
normalize_ids <- function(x){
  x <- sub("^X(?=\\d)", "", x, perl=TRUE)
  x <- sub("^GO\\.", "GO:", x); x <- sub("^HP\\.", "HP:", x)
  x <- sub("^R[._-]?HSA[._-]", "R-HSA-", x)
  sub("[._](\\d{3,7})$", "-\\1", x, perl=TRUE)
}

# Función para obtener las anotaciones correspondientes a una combinación
# de método y base de datos.
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

# Helper ESPECIAL para calcular la importancia en modelos KNN.
# Los modelos como KNN no tienen una medida de importancia intrínseca.
# Esta función implementa la "importancia por permutación":
# 1. Mide el rendimiento del modelo con los datos originales.
# 2. Para cada ruta, desordena aleatoriamente sus valores.
# 3. Vuelve a medir el rendimiento.
# 4. La caída en el rendimiento indica cuán importante era esa ruta.
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

#--- Definir directorios de salida -----------------------------------------------
dir_plots <- "ML_Plots/Importance/Plots"
dir_csv   <- "ML_Plots/Importance/Anotations"
dir.create(dir_plots, recursive=TRUE, showWarnings=FALSE)
dir.create(dir_csv,   recursive=TRUE, showWarnings=FALSE)

#--- BUCLE PRINCIPAL: Calcular y visualizar importancia para cada modelo ---
for(i in seq_len(nrow(results_df))){
  sm <- results_df$ScoreMethod[i]
  gs <- results_df$Database[i]
  alg <- results_df$BestModel[i]
  tag <- paste(alg, sm, gs, sep="_") # Etiqueta única para los archivos.
  f   <- file.path("ML_Plots", paste0("trainedModel_", sm, "_", gs, ".RData"))
  
  logf("(%d/%d) Iniciando %s", i, nrow(results_df), tag)
  
  # Comprobamos si el archivo del modelo existe antes de intentar cargarlo.
  if(!file.exists(f)){
    logf("AVISO: Modelo no encontrado, saltando: %s", f)
    next # Si no existe, saltamos a la siguiente iteración.
  }
  
  # Esta parte mide el tiempo que tarda en cargar un archivo 
  # en R y luego imprime un mensaje en la consola mostrando ese tiempo en segundos.
  t_load <- proc.time()[3]; load(f); t_load <- proc.time()[3]-t_load
  logf("Modelo cargado en %.2fs", t_load)
  
  mdl <- trainedModel$model
  
  # --- 1. Cálculo de Importancia ---
  t_imp <- proc.time()[3]
  # Determinamos qué rutas (características) fueron más influyentes para las predicciones del modelo.
  # El método varía según el algoritmo.
  if (isTRUE(grepl("^knn$", mdl$method, ignore.case=TRUE))) {
    # Para KNN, usamos la "importancia por permutación": desordenamos cada ruta
    # y medimos cuánto empeora el modelo. Una gran caída en el rendimiento
    # significa que la ruta era muy importante.
    nfeat <- if (!is.null(mdl$finalModel$xNames)) length(mdl$finalModel$xNames) else NA_integer_
    logf("Cálculo importancia por permutación (KNN, nsim=5, p=%s)", ifelse(is.na(nfeat), "?", nfeat))
    imp <- tryCatch({
      df <- imp_knn(mdl, nsim=5, seed=123)
      tibble(Term=df$Term, Overall=df$Overall)
    }, error=function(e){ logf("ERROR imp_knn: %s", e$message); tibble(Term=character(), Overall=numeric()) })
  } else {
    # Para modelos como Random Forest o XGBoost, usamos la función estándar
    # `varImp` de caret, que es mucho más eficiente.
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
  
  # --- 2. Anotación de términos ---
  # Los resultados de importancia nos dan IDs de rutas (ej. "hsa04110").
  # Ahora los traducimos a nombres descriptivos (ej. "Cell cycle").
  t_ann <- proc.time()[3]
  ann <- get_ann(sm, gs)
  logf("Anotación ann2term: %d términos", if (is.null(ann)) 0L else nrow(ann))
  
  # Unimos la tabla de importancia con la tabla de anotaciones.
  imp <- imp %>%
    arrange(desc(Overall)) %>% # Ordenamos las rutas de más a menos importantes.
    mutate(ID_norm = normalize_ids(Term)) %>% # Limpiamos los IDs para que coincidan.
    left_join({ if(is.null(ann)) tibble(ID=character(), term=character()) else as_tibble(ann) },
              by = c("ID_norm"="ID")) %>% # Unimos las dos tablas.
    mutate(TermPlot = ifelse(is.na(term) | term=="", Term, term)) %>% # Si no hay nombre, usamos el ID.
    select(Term = TermPlot, Importance = Overall, ID = ID_norm) # Seleccionamos y renombramos columnas.
  
  t_ann <- proc.time()[3]-t_ann
  logf("Anotado en %.2fs", t_ann)
  
  # --- 3. EXPORTACIÓN DE CSVs Y PLOT ---
  all_csv_path <- file.path(dir_csv, paste0(tag, "_all_importances.csv"))
  
  # Guardamos la lista COMPLETA de importancias en un archivo CSV.
  # 'write.csv2' usa punto y coma como separador, ideal para abrir en Excel en Europa.
  write.csv2(imp, all_csv_path, row.names = FALSE)
  logf("CSV con TODAS las importancias guardado en: %s", all_csv_path)
  
  top20_csv_path <- file.path(dir_csv, paste0(tag, "_top20_importances.csv"))
  
  # Hacemos lo mismo para el archivo del top 20.
  write.csv2(imp %>% slice_head(n = 20), top20_csv_path, row.names = FALSE)
  logf("CSV con el TOP 20 guardado en: %s", top20_csv_path)
  
  # Creamos un gráfico de barras con las 5 rutas más importantes.
  top5 <- imp %>% slice_head(n = 5)
  if (nrow(top5) > 0) {
    # `str_wrap` ajusta los nombres largos para que quepan bien en el gráfico.
    top5$Term <- stringr::str_wrap(top5$Term, width = 45)
    # Este truco ordena las barras en el gráfico de mayor a menor importancia.
    top5$Term <- factor(top5$Term, levels = rev(top5$Term))
    
    p <- ggplot(top5, aes(x = Term, y = Importance)) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      coord_flip() + # Pone las barras en horizontal para que los nombres se lean mejor.
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

# --- Mensaje final ---
logf("========== PROCESO FINALIZADO ==========")
cat("Resultados guardados en las siguientes carpetas:\n")
cat("-> Plots (.png):", normalizePath(dir_plots), "\n")
cat("-> Anotaciones (.csv):", normalizePath(dir_csv), "\n")



# -----------------------------------------------------------------------------
# 6. Heatmaps (PNG) con IPD a la izquierda y CONTROL a la derecha
# -----------------------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)
library(pathMED)
library(grid)

# Configuración General del Heatmap
dir.create("ML_Plots/Heatmaps", showWarnings = FALSE, recursive = TRUE)
DO_ZSCORE <- TRUE
TOP_N     <- 100
PNG_W     <- 2200; PNG_H <- 2000; PNG_RES <- 200

# Definición de Clases (realizado tambien en anteriores visualizaciones)
pos_class <- "IPD"
neg_class <- "CONTROL"

# Limpia los IDs de las rutas para que sean consistentes.
norm_ids <- function(ids){
  x <- ids
  x <- sub("^X(?=\\d)", "", x, perl=TRUE)
  x <- sub("^GO\\.", "GO:", x)
  x <- sub("^HP\\.", "HP:", x)
  x <- sub("^R[._-]?HSA[._-]?", "R-HSA-", x)
  x <- sub("[._](\\d{3,7})$", "-\\1", x, perl=TRUE)
  make.unique(x)
}
# Calcula el Z-score para cada fila (ruta). Esto estandariza los datos,
# resaltando los patrones relativos (qué está más alto/bajo de lo normal para ESA ruta)
# en lugar de los valores absolutos.
row_z <- function(m){
  m <- as.matrix(m)
  mu <- rowMeans(m, na.rm=TRUE)
  sdv <- apply(m, 1, sd, na.rm=TRUE); sdv[sdv==0 | is.na(sdv)] <- 1
  sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
}
# Selecciona las 'n' filas con la mayor varianza.
top_var <- function(m, n=100){
  if(nrow(m)<=n) return(m)
  v <- apply(m, 1, var, na.rm=TRUE)
  m[order(v, decreasing=TRUE)[seq_len(n)], , drop=FALSE]
}
# Ordena las muestras para que todas las IPD aparezcan primero
# y luego todas las de CONTROL. Esto es esencial para la comparación visual.
order_cols_pos_neg <- function(y_chr, pos=pos_class, neg=neg_class){
  if (all(c(pos,neg) %in% unique(y_chr))) c(which(y_chr==pos), which(y_chr==neg)) else order(y_chr)
}
# Traduce los IDs de las rutas a nombres descriptivos usando ann2term.
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

# Bucle Principal para Generar Heatmaps ---
files <- list.files("ML_Plots", pattern="^trainedModel_.*\\.RData$", full.names=TRUE)
for (f in files){
  load(f)  # -> trainedModel
  
  # Extraemos los datos exactos que se usaron para entrenar este modelo.
  td <- if(!is.null(trainedModel$model$trainingData)) trainedModel$model$trainingData else trainedModel$trainingData
  if (is.null(td) || !(".outcome" %in% names(td))) {
    cat("Sin trainingData utilizable en:", f, "- se salta\n")
    rm(trainedModel); gc(); next
  }
  
  # Separamos las etiquetas (y) de las características (X).
  y <- as.character(td$.outcome)
  X <- td[, setdiff(names(td), ".outcome"), drop=FALSE]  # muestras x features
  # Transponemos la matriz para tener el formato estándar: rutas en filas, muestras en columnas.
  M <- t(as.matrix(X))                                  # filas=rutas, cols=muestras
  rownames(M) <- norm_ids(rownames(M))
  
  # Aplicamos la función para ordenar las columnas (muestras): primero IPD, luego CONTROL.
  ord <- order_cols_pos_neg(y)   # IPD -> CONTROL
  M   <- M[, ord, drop=FALSE]
  y   <- y[ord]
  
  # Aplicamos el escalado Z-score a las filas (si está activado).
  Mz <- if (DO_ZSCORE) row_z(M) else as.matrix(M)
  
  # Etiquetas ann2term
  rownames(Mz) <- make.unique(annotate_all(rownames(Mz), rowMeans(Mz, na.rm=TRUE)))
  
  # Barra Outcome
  lvl <- if (all(c(pos_class,neg_class) %in% unique(y))) c(pos_class,neg_class) else sort(unique(y))
  y_fac <- factor(y, levels = lvl)
  
  # Creamos la barra de color superior que identificará a qué grupo pertenece cada muestra.
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
               cluster_columns=FALSE))   # <-- mantiene IPD|CONTROL fijos
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
               cluster_columns=FALSE))   # <-- mantiene IPD|CONTROL fijos
  dev.off()
  
  cat("Heatmaps PNG generados (IPD|CONTROL separados) para:", model_name, "\n")
  rm(trainedModel); gc()
}

# Guardamos el entorno con los objetos usados en esta sección de visualización.
save.image("entorno_TFM_ML_plots") 

# -----------------------------------------------------------------------------
# 7.Estratificación para TODAS las combinaciones
#
# - Realiza TODAS las comparaciones pairwise entre clusters (C1-C2, C1-C3, etc.).
# - Columnas Max_RCSI y Max_Prop_IPD a la tabla resumen final.
# - Genera TODOS los plots y archivos del script original.
# - Muestra una barra de progreso en tiempo real.
# - Crea una subcarpeta para cada combinación.
# - Utiliza doParallel y foreach para acelerar el proceso.
# -----------------------------------------------------------------------------
# OBJETIVO: En las secciones anteriores, usamos Machine Learning SUPERVISADO para
# predecir las etiquetas que ya conocíamos (IPD vs. Control). Ahora, cambiamos a
# un enfoque NO SUPERVISADO. El objetivo es descubrir si los pacientes se agrupan
# de forma natural en subtipos (clusters) basándonos únicamente en sus perfiles
# de actividad de rutas, sin usar las etiquetas de enfermedad. Esto podría
# revelar una heterogeneidad biológica desconocida dentro del grupo de pacientes.
#
# ESTRATEGIA: Usaremos el algoritmo M3C (Monte Carlo Consensus Clustering), que es
# un método muy robusto para encontrar el número óptimo de clusters y asignar
# cada muestra a uno de ellos. Repetiremos este análisis para TODAS las
# combinaciones de score/base de datos y lo haremos en paralelo para ahorrar tiempo.

# --- 0. CARGAR LIBRERÍAS NECESARIAS ---
suppressPackageStartupMessages({
  library(M3C); library(ggplot2); library(ComplexHeatmap); library(circlize)
  library(cluster); library(pathMED); library(stringr); library(foreach)
  library(doParallel); library(dplyr); library(progressr); library(parallel)
})

# --- 1. CONFIGURACIÓN GENERAL ---
pos_class <- "IPD"; neg_class <- "CONTROL"
output_dir <- "ML_Plots/Clustering"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# Parámetros para el algoritmo M3C. Valores más altos dan más robustez pero tardan más.
M3C_MAX_K <- 6; M3C_ITERS <- 50; M3C_REPS_REAL <- 100; M3C_REPS_REF <- 100

# --- 2. GENERAR LA LISTA DE TAREAS ---
# `expand.grid` crea una tabla con todas las combinaciones posibles de métodos de
# score y bases de datos. Esta será nuestra "lista de tareas" para el bucle paralelo.
combinations <- expand.grid(score_method=names(scores_all), gene_set=names(scores_all[[1]]))
cat(sprintf("Se van a procesar %d combinaciones en total.\n", nrow(combinations)))

# --- 3. LA GRAN FUNCIÓN DE ANÁLISIS ---
# Definimos las funciones helper fuera del bucle
# Función para normalizar los IDs de las rutas (ej. "GO.123" -> "GO:123")
normalize_ids <- function(x){
  x <- sub("^X(?=\\d)", "", x, perl=TRUE); x <- sub("^GO\\.", "GO:", x)
  x <- sub("^HP\\.", "HP:", x); x <- sub("^R[._-]?HSA[._-]", "R-HSA-", x)
  sub("[._](\\d{3,7})$", "-\\1", x, perl=TRUE)
}
# Cálculo de Z-score
row_z <- function(m) {
  mu <- rowMeans(m, na.rm=TRUE)
  sdv <- apply(m, 1, sd, na.rm=TRUE); sdv[sdv==0 | is.na(sdv)] <- 1
  sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
}
# Anotación de rutas
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

# Esta es la función principal que se ejecutará en paralelo
run_clustering_analysis <- function(score_method, gene_set) {
  
  file_tag <- paste(score_method, gene_set, sep = "_")
  combination_dir <- file.path(output_dir, file_tag)
  if (!dir.exists(combination_dir)) dir.create(combination_dir, recursive = TRUE)
  
  # `tryCatch` es un mecanismo de seguridad: si algo falla dentro de este bloque,
  # el script no se detendrá. En su lugar, registrará el error y continuará.
  tryCatch({
    mat <- as.matrix(scores_all[[score_method]][[gene_set]])
    if(nrow(mat) < 2 || ncol(mat) < 2) {
      return(data.frame(Combination=file_tag, Status="Skipped", Optimal_K=NA, Max_RCSI=NA, Silhouette_Avg=NA, Separation_PValue=NA, Max_Prop_IPD=NA))
    }
    
    # 1. Ejecutar M3C y extraer resultados clave
    set.seed(123)
    res <- M3C(mydata=mat, maxK=M3C_MAX_K, iters=M3C_ITERS, repsreal=M3C_REPS_REAL, repsref=M3C_REPS_REF, seed=123, removeplots=TRUE)
    df <- res$scores # Contiene métricas como el RCSI para elegir el mejor K.
    clusters <- res$assignments # Nos dice a qué cluster pertenece cada muestra.
    max_rcsi_value <- round(max(df$RCSI, na.rm = TRUE), 4) # busca la puntuación RCSI 
    # (Índice de Estabilidad Relativa del Clúster) más alta obtenida en el análisis de clustering.
    
    # 2. Alinear datos y calcular Métricas de Calidad y Relevancia Clínica
    outcome <- as.character(pheno_bin$Response); names(outcome) <- rownames(pheno_bin)
    common <- intersect(colnames(mat), names(outcome))
    mat_sub <- mat[, common, drop=FALSE]
    if (is.null(names(clusters))) names(clusters) <- colnames(mat)
    clusters_sub <- clusters[common]
    outcome_sub  <- outcome[common]
    
    # a) Relevancia Clínica: ¿Los clusters que encontramos "a ciegas" se corresponden
    #    con las etiquetas reales de IPD/Control? Lo medimos con un test estadístico.
    tab <- table(Cluster = clusters_sub, Outcome = outcome_sub)
    test_res <- if (any(tab < 5)) fisher.test(tab) else chisq.test(tab)
    prop_table <- prop.table(tab, 1)
    max_prop_ipd <- if (pos_class %in% colnames(prop_table)) round(max(prop_table[, pos_class]), 4) else 0
    
    # b) Calidad del Clustering: ¿Qué tan bien separados están los clusters?
    #    Lo medimos con el Coeficiente de Silueta (Silhouette). Un valor cercano a 1
    #    indica clusters muy bien definidos y separados.
    dmat <- dist(t(mat_sub), method="euclidean")
    sil <- silhouette(as.integer(factor(clusters_sub)), dmat)
    sil_mean <- mean(sil[, "sil_width"], na.rm=TRUE)
    
    # 3. Generar plots y archivos generales (RCSI, Heatmap, PCA, etc.)
    
    # RCSI plot y asignaciones
    png(file.path(combination_dir, paste0("M3C_RCSI_", file_tag, ".png")), width=1200, height=1000, res=150)
    plot(df$K, df$RCSI, type="b", pch=19, xlab="K", ylab="RCSI", main=paste("Selección de K (", file_tag, ")"))
    arrows(df$K, df$RCSI - df$RCSI_SE, df$K, df$RCSI + df$RCSI_SE, angle=90, code=3, length=0.05)
    dev.off()
    write.csv2(clusters, file.path(combination_dir, paste0("M3C_clusters_", file_tag, ".csv")))
    
    # Tabla de contingencia y test
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
    
    # Barplot de composición
    p_bar <- ggplot(as.data.frame(prop.table(tab, 1)), aes(x = Cluster, y = Freq, fill = Outcome)) + geom_col(width=0.7, color="grey30") + theme_minimal(base_size=12) + scale_y_continuous(labels=scales::percent_format()) + labs(title="Composición por cluster", y="Proporción", x="Cluster")
    ggsave(file.path(combination_dir, paste0("M3C_composition_barplot_", file_tag, ".png")), p_bar, width=6.5, height=4.8, dpi=300)
    
    # Silhouette
    write.csv2(as.data.frame(sil[, 1:3]), file.path(combination_dir, paste0("M3C_silhouette_", file_tag, ".csv")))
    png(file.path(combination_dir, paste0("M3C_silhouette_plot_", file_tag, ".png")), width=1800, height=900, res=200)
    plot(sil, main=sprintf("Silhouette (media = %.3f)", sil_mean), col=2:(length(unique(clusters_sub))+1), border=NA)
    dev.off()
    
    # --- 4. TOP RUTAS DIFERENCIALES PARA TODOS LOS PARES ---
    # Este es un paso clave de interpretación. Para cada par de clusters (ej. C1 vs C2),
    # buscamos las rutas que tienen una actividad significativamente diferente.
    grpC <- factor(clusters_sub)
    if (nlevels(grpC) >= 2) {
      
      cluster_pairs <- combn(levels(grpC), 2, simplify = FALSE)
      
      # Usamos `lapply` para realizar la comparación para todos los pares.
      all_diff_results <- lapply(cluster_pairs, function(pair) {
        c1 <- pair[1]; c2 <- pair[2]
        idx1 <- which(grpC == c1); idx2 <- which(grpC == c2)
        
        # Para cada ruta, realizamos un test de Wilcoxon para ver si hay diferencias.
        pvals <- apply(mat_sub, 1, function(x) wilcox.test(x[idx1], x[idx2])$p.value)
        padj  <- p.adjust(pvals, method="fdr")
        mean1 <- rowMeans(mat_sub[, idx1, drop=FALSE], na.rm=TRUE)
        mean2 <- rowMeans(mat_sub[, idx2, drop=FALSE], na.rm=TRUE)
        
        # Crea una tabla con los resultados de la comparación (medias, p-valor, etc.).
        res_pair <- data.frame(
          Comparison = paste0("C", c1, "_vs_C", c2), Pathway_ID = rownames(mat_sub),
          Term = map_ids_to_terms(mat_sub, rownames(mat_sub)),
          Mean_C1 = mean1, Mean_C2 = mean2, Diff = mean1 - mean2,
          pval = pvals, padj = padj,
          Direction = ifelse(mean1 > mean2, paste0("↑ C", c1), paste0("↑ C", c2))
        )
        # Nombra las columnas de media de forma dinámica (ej. "Mean_C1", "Mean_C2").
        names(res_pair)[names(res_pair) == "Mean_C1"] <- paste0("Mean_C", c1)
        names(res_pair)[names(res_pair) == "Mean_C2"] <- paste0("Mean_C", c2)
        return(res_pair)
      })
      # Une los resultados de todas las comparaciones en una única tabla.
      res_clust_all <- dplyr::bind_rows(all_diff_results)
      # Ordena la tabla por comparación y p-valor ajustado.
      res_clust_all <- res_clust_all %>% arrange(Comparison, padj)
      
      write.csv2(res_clust_all, file.path(combination_dir, paste0("M3C_topPathways_AllPairs_", file_tag, ".csv")), row.names=FALSE)
      
      # Bucle para crear un gráfico para cada par de clústeres.
      for (pair in cluster_pairs) {
        c1 <- pair[1]; c2 <- pair[2]
        comp_name <- paste0("C", c1, "_vs_C", c2)
        
        # Selecciona las 10 rutas más significativas para la comparación actual.
        top10 <- res_clust_all %>% filter(Comparison == comp_name) %>% head(10)
        
        if (nrow(top10) >= 2) {
          # Crea un dotplot que visualiza las 10 rutas principales.
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
    
    # 5. Devolver el resumen final
    # Devuelve una única fila con el resumen de métricas para esta combinación.
    return(data.frame(
      Combination=file_tag, Status="Completed", Optimal_K=length(unique(clusters)),
      Max_RCSI=max_rcsi_value, Silhouette_Avg=round(sil_mean, 4),
      Separation_PValue=test_res$p.value, Max_Prop_IPD=max_prop_ipd
    ))
  }, error = function(e) {
    return(data.frame(Combination=file_tag, Status=paste("Error:", e$message), Optimal_K=NA, Max_RCSI=NA, Silhouette_Avg=NA, Separation_PValue=NA, Max_Prop_IPD=NA))
  })
}

# --- 4. EJECUCIÓN PARALELA CON BARRA DE PROGRESO ---
# Configura el clúster para usar N-1 núcleos de la CPU, acelerando el proceso.
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("\nIniciando procesamiento paralelo con %d núcleos...\n", num_cores))
with_progress({
  # Bucle paralelo `foreach`: ejecuta `run_clustering_analysis` para cada combinación.
  # `.combine = 'rbind'` une las filas de resumen devueltas por cada ejecución.
  p <- progressor(steps = nrow(combinations))
  summary_results <- foreach(i = 1:nrow(combinations), .combine = 'rbind', .packages = c("M3C", "ggplot2", "ComplexHeatmap", "circlize", "cluster", "pathMED", "stringr", "dplyr", "parallel")) %dopar% {
    # Actualiza la barra de progreso en cada iteración.
    p(sprintf("Procesando %s_%s", combinations$score_method[i], combinations$gene_set[i])) 
    run_clustering_analysis(combinations$score_method[i], combinations$gene_set[i])
  }
})
# Detiene el clúster paralelo para liberar los recursos.
stopCluster(cl)
cat("\n--- PROCESO PARALELO FINALIZADO ---\n")

# --- 5. RESULTADOS FINALES ---
if (!is.null(summary_results) && nrow(summary_results) > 0) {
  # Ordena los resultados finales por el Silhouette Score (de mejor a peor).
  summary_results <- summary_results %>% arrange(desc(Silhouette_Avg))
  print("Tabla Resumen de Resultados de Clustering:")
  print(summary_results)
  write.csv2(summary_results, file.path(output_dir, "CLUSTERING_SUMMARY_RESULTS.csv"), row.names = FALSE)
  cat(sprintf("\nTabla resumen guardada en: %s\n", file.path(output_dir, "CLUSTERING_SUMMARY_RESULTS.csv")))
} else {
  cat("No se generaron resultados. Revisa los mensajes de error.\n")
}

# Guarda el entorno completo con todos los objetos y resultados finales.
save.image("entorno_TFM_TODO")




sessionInfo()
