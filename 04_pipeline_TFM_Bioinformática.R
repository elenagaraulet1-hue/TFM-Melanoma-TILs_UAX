### IMPACTO DE LA ABUNDANCIA DE SUBPOBLACIONES DE TILS EN LA SUPERVIVENCIA DEL MELANOMA
## Autor: Elena Garaulet Belda

#------------------------------------------------------------------------------
########  BLOQUE 0: Establecemos el entorno de trabajo.  ----------------------
if (!require("here")) install.packages("here")
library(here)

# Creamos toda la estructura de directorios de una sola vez de forma segura
directorios <- c("data/raw", "data/processed", "results/plots", "results/tables", "logs", "scripts")
for (dir in directorios) {
  if(!dir.exists(here(dir))) {
    dir.create(here(dir), recursive = TRUE)
  }
}
message("Entorno de trabajo y carpetas preparadas con éxito.")

#------------------------------------------------------------------------------
########  BLOQUE 1: Descarga de datos originales  ------------------------------
# Objetivo: Descargar datos de TCGA-SKCM y guardarlos en data/raw
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(here)

# 1. Configuración de Log
if (!dir.exists(here("logs"))) dir.create(here("logs"), recursive = TRUE)
# nombre: 
log_file <- here("logs", paste0("log_descarga_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")) 
con_log <- file(log_file, open = "wt") 

sink(con_log, append = TRUE, split = TRUE)
sink(con_log, append = TRUE, split = FALSE, type = "message") 
message("--- LOG INICIADO: ", Sys.time(), " ---")


# 2. Ejecución con manejo de errores (tryCatch)
datos_ok <- tryCatch({
  message("Paso 1: Iniciando consulta a GDC (TCGA-SKCM) - Filtrando muestras metastásicas ..." )
  
  query <- GDCquery(
    project = "TCGA-SKCM",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Metastatic"
  )
  
  GDCdownload(query)
  data <- GDCprepare(query)
  
  # Guardar en raw
  if (!dir.exists(here("data/raw"))) dir.create(here("data/raw"), recursive = TRUE)
  saveRDS(data, here("data/raw/tcga_skcm_metastatic.rds"))
  
  message("ÉXITO: Datos guardados en data/raw/tcga_skcm_metastatic.rds")
data
  
}, error = function(e) {
  message("ERROR CRÍTICO: ", e$message)
  return(NULL)
})

# 3. Resumen estadístico (validación de descarga de datos)
if (!is.null(datos_ok)) {
  # Accedemos a los metadatos clínicos 
  col_data <- colData(datos_ok)
  
  resumen_datos <- data.frame(
    Muestra = c("Total muestras descargadas (Metastásicas)"),
    Cantidad = c(ncol(datos_ok))
  )
  
  if (!dir.exists(here("results/tables"))) dir.create(here("results/tables"), recursive = TRUE)
  write.csv(resumen_datos, here("results/tables/resumen_datos_filtrados.csv"), row.names = FALSE)
  
  message("Resumen de datos generado en results/tables/resumen_datos_filtrados.csv")
  
  # Mostrar en consola
  print(resumen_datos)
}

#------------------------------------------------------------------------------
########  BLOQUE 2: Limpieza y normalización de los datos. --------------------
library(SummarizedExperiment)
library(edgeR)
library(here)
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager") 
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# 1. Configuración de Log (cargamos la lógica que ya creamos)

message("--- INICIO DE LIMPIEZA Y NORMALIZACIÓN: ", Sys.time(), " ---")

# 2. Cargar datos brutos
data <- readRDS(here("data/raw/tcga_skcm_metastatic.rds"))
  # Extraemos la matriz de conteos
counts <- assay(data) 
message("Dimensiones originales: ", nrow(counts), " genes, ", ncol(counts), " muestras.")


# 3. Filtrado: Eliminar genes con baja expresión
# Umbral: genes con 10 counts en el 20% de las muestras
counts_cpm <- edgeR::cpm(counts) 
keep <- rowSums(counts_cpm > 1) >= (0.2 * ncol(counts))
counts_filtrados <- counts[keep, ]

# Guardamos el número de genes: 
num_genes <- nrow(counts_filtrados)
message("Genes filtrados. Genes restantes: ", num_genes)
print(paste("Número de genes filtrados:", num_genes))

# 4. Normalización (TPM - esencial para deconvolución)
# Usaremos la normalización CPM log-transformada para análisis estadísticos
dge <- DGEList(counts = counts_filtrados)
dge <- calcNormFactors(dge, method = "TMM")
matriz_norm <- cpm(dge, log = FALSE) # Normalización TPM/CPM

# 5. conversión de IDs Ensembl a Symbol. 
rownames(matriz_norm) <- gsub("\\..*", "", rownames(matriz_norm))
ids <- mapIds(org.Hs.eg.db, keys = rownames(matriz_norm), column = "SYMBOL", 
              keytype = "ENSEMBL", multiVals = "first")
  # Limpieza de datos: 
matriz_final <- matriz_norm[!is.na(ids), ]
rownames(matriz_final) <- ids[!is.na(ids)]
matriz_final <- matriz_final[!duplicated(rownames(matriz_final)), ]

# 6. Guardar
if (!dir.exists(here("data/processed"))) dir.create(here("data/processed"), recursive = TRUE)
saveRDS(matriz_final, here("data/processed/tcga_skcm_normalizado.rds"))

message("ÉXITO: Matriz normalizada guardada. Genes finales:", nrow(matriz_final))
print(paste("Número de genes en matriz normalizada:", nrow(matriz_final)))

message("Estadísticas de la matriz final:")
print(summary(colSums(matriz_final))) 


#------------------------------------------------------------------------------
########  BLOQUE 3: Deconvolución CIBERSORT  ----------------------------------
#if (!require("devtools")) install.packages("devtools")
if (!require("immunedeconv")) devtools::install_github("omnideconv/immunedeconv")
library(here)
library(tidyr)
library(tibble)
library(here)

message("--- Descarga de archivos...: ", Sys.time(), " ---")
if (!dir.exists(here("scripts"))) dir.create(here("scripts"), recursive = TRUE)
if (!dir.exists(here("data/raw"))) dir.create(here("data/raw"), recursive = TRUE)

# 1. URLs de los archivos
url_cibersort <- "https://raw.githubusercontent.com/SiYangming/CIBERSORT-DATA/main/CIBERSORT.R"
url_lm22 <- "https://raw.githubusercontent.com/SiYangming/CIBERSORT-DATA/main/LM22.txt"

# 2. Descarga de archivos 
if (!file.exists(here("scripts/CIBERSORT.R"))) {
  message("Descargando CIBERSORT.R...") 
  download.file(url_cibersort, here("scripts/CIBERSORT.R")) 
} else {
    message("CIBERSORT.R ya existe, saltando descarga.")
}

if (!file.exists(here("data/raw/LM22.txt"))) {
  message("Descargando LM22.txt...")
  download.file(url_lm22, here("data/raw/LM22.txt"))
} else {
  message("LM22.txt ya existe, saltando descarga.")
}

# 3. Instalación y carga. 
paquetes <- c("e1071", "preprocessCore")
if (length(setdiff(paquetes, installed.packages()[, "Package"]))) {
  install.packages("e1071")
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("preprocessCore")
}

source(here("scripts/CIBERSORT.R"))


message("--- Ejecutando CIBERSORT...: ", Sys.time(), " ---")

# 1. Preparar archivo de mezcla: 
if (!dir.exists(here("data/processed"))) dir.create(here("data/processed"), recursive = TRUE)

cibersort_input <- cbind(GeneSymbol = rownames(matriz_final), as.data.frame(matriz_final))
write.table(cibersort_input, here("data/processed/mixture_to_cibersort.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# 2. Resultados CIBERSORT
resultados_ciber <- CIBERSORT(here("data/raw/LM22.txt"), here("data/processed/mixture_to_cibersort.txt"), 
                              perm = 100, QN = FALSE)

rownames(resultados_ciber) <- gsub("\\.", "-", rownames(resultados_ciber))

# 3. Guardar en la estructura de carpetas definida
if (!dir.exists(here("results/tables"))) dir.create(here("results/tables"), recursive = TRUE)
write.csv(resultados_ciber, here("results/tables/Resultados_TILs_SKCM.csv"))

# Filtrar muestras con p-valor significativo
muestras_validas <- resultados_ciber[resultados_ciber[, "P-value"] < 0.05, ]
message("Muestras con deconvolución significativa (p < 0.05): ", nrow(muestras_validas), " de ", nrow(resultados_ciber))

# Guardar solo las muestras de calidad
write.csv(muestras_validas, here("results/tables/Resultados_TILs_Calidad_Alta.csv"))


message("CIBERSORT finalizado. Muestras procesadas: ", nrow(resultados_ciber)) 
message("Muestras con deconvolución significativa (p < 0.05): ", nrow(muestras_validas), " de ", nrow(resultados_ciber))
print(paste("CIBERSORT finalizado. Muestras procesadas: ", nrow(resultados_ciber)))



#------------------------------------------------------------------------------
########  BLOQUE 3.1: HETEROGENEIDAD CELULAR  ----------------------------------

message("Paisaje inmunitario: Heterogeneidad celular")

# 1. Visualización de Heterogeneidad Inmune: 
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# 2. Preparación de datos
df_long_bar <- as.data.frame(resultados_ciber) %>%
  dplyr::select(-c(`P-value`, Correlation, RMSE)) %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion")

# 3. Orden por CD8+: validación de la existencia de la columna 
col_cd8_name <- grep("CD8", colnames(resultados_ciber), value = TRUE)[1]

if (is.na(col_cd8_name)) {
  message("ADVERTENCIA: No se encontró columna 'CD8'. Ordenando por la primera columna disponible.")
  col_cd8_name <- colnames(resultados_ciber)[1]
}

orden_muestras <- resultados_ciber[order(resultados_ciber[, col_cd8_name], decreasing = TRUE), ]
df_long_bar$Sample <- factor(df_long_bar$Sample, levels = rownames(orden_muestras))

# 4. Generación del gráfico:
p_bar <- ggplot(df_long_bar, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(ncol(resultados_ciber))) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +
  labs(title = paste("Heterogeneidad Inmune (Ordenado por", col_cd8_name, ")"), 
       x = "Pacientes (TCGA)", y = "Proporción")

# 5. Guardado
if (!dir.exists(here("results/plots"))) dir.create(here("results/plots"), recursive = TRUE)
ggsave(here("results/plots/1.heterogeneidad_inmune_barplot.pdf"), plot = p_bar, width = 12, height = 6)
message("Barplot de Heterogeneidad inmune guardado en results/plots")




#------------------------------------------------------------------------------
########  BLOQUE 4: Deconvolución quanTIseq - Validación -----------------------
library(immunedeconv)
library(tidyr)    
library(tibble)   
library(dplyr)   
library(ggplot2)
library(ggpubr)

message("EJECUTANDO VALIDACIÓN CON QUANTISEQ...")

# 2. Deconvolución
res_quantiseq_raw <- deconvolute(matriz_final, "quantiseq")

# 3. Transformar formato para análisis comparativo
# Eliminamos columna cell_type: 
df_temp <- res_quantiseq_raw %>% column_to_rownames("cell_type")
# Transponemos para tener el mismo formato
res_q_df <- as.data.frame(t(df_temp))
# Limpieza: 
rownames(res_q_df) <- gsub("\\.", "-", rownames(res_q_df))
print(colnames(res_q_df))

# 4. Guardar resultados de la validación
write.csv(res_q_df, here("results/tables/Resultados_Validacion_quantiseq_SKCM.csv"))

message("Validación con quanTIseq completada. Muestras obtenidas: ", nrow(res_q_df))
print(paste("quanTIseq finalizado. Muestras procesadas: ", nrow(res_q_df)))



message("--- CORRELACIÓN ENTRE ALGORITMOS ---")
# 1. Cruzar muestras comunes
muestras_comunes <- intersect(rownames(resultados_ciber), rownames(res_q_df))

# 2. Extraer CD8 (buscamos nombres similares en ambos)
col_cd8_ciber <- "T cells CD8" 
col_cd8_quant <- grep("CD8", colnames(res_q_df), ignore.case = TRUE, value = TRUE)[1]

# Verifica qué encontró
message("Columna CD8 en CIBERSORT: ", col_cd8_ciber)
message("Columna CD8 en quanTIseq: ", col_cd8_quant)
if (is.na(col_cd8_ciber) || is.na(col_cd8_quant)) {
  stop("¡Error! No se encontró la subpoblación CD8 en alguno de los métodos. Revisa colnames().")
}

# 3. Creación del dataframe 
df_validacion <- data.frame(
  Sample = muestras_comunes,
  CIBERSORT = as.numeric(resultados_ciber[muestras_comunes, col_cd8_ciber]),
  quanTIseq = as.numeric(res_q_df[muestras_comunes, col_cd8_quant])
)
df_validacion[df_validacion < 0] <- 0


# 3. Gráfico de Correlación de Validación
p_val <- ggplot(df_validacion, aes(x = CIBERSORT, y = quanTIseq)) +
  geom_point(color = "#88A0A8", alpha = 0.6) +
  geom_smooth(method = "lm", color = "#D4A373", fill = "#E9EDC9") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(title = "Validación Cruzada: CIBERSORTx vs quanTIseq",
       subtitle = paste("Correlación en células T CD8+ (", col_cd8_ciber, " vs ", col_cd8_quant, ")"),
       x = "Proporción CIBERSORTx", y = "Proporción quanTIseq")

ggsave(here("results/plots/2.validacion_cruzada_cd8.pdf"), plot = p_val, width = 7, height = 6)
message("Gráfico de validación guardado en results/plots/validacion_cruzada_cd8.pdf")






#------------------------------------------------------------------------------
########  BLOQUE 5: ANÁLISIS DE SUPERVIVENCIA ---------------------------------
# Objetivo: ¿Tener más células T CD8+ o Macrófagos M2 hace que el paciente sobreviva más tiempo?
library(survival)
library(survminer)
library(here)
library(dplyr)
library(SummarizedExperiment)

message("--- INICIO DE ANÁLISIS DE SUPERVIVENCIA ---")
# 1. Asegurar que tenemos los datos clínicos listos
# Extraer los datos clínicos a un dataframe
data_raw <- readRDS(here("data/raw/tcga_skcm_metastatic.rds"))
clinico_df <- as.data.frame(colData(data_raw))
clinico_df$sample <- rownames(clinico_df)

# 2. Cargar resultados CIBERSORT
res_ciber <- read.csv(here("results/tables/Resultados_TILs_SKCM.csv"), row.names = 1)
res_ciber$sample <- rownames(res_ciber)

# 3. Merge: 
df_final <- merge(res_ciber, clinico_df, by = "sample")

# 4. Continuar con la construcción de variables...
df_final$death_time <- as.numeric(df_final$days_to_death)
df_final$follow_up_time <- as.numeric(df_final$days_to_last_follow_up)
df_final$time <- ifelse(is.na(df_final$death_time), df_final$follow_up_time, df_final$death_time)
df_final$status <- ifelse(df_final$vital_status == "Dead", 1, 0)

df_final <- df_final[!is.na(df_final$time) & !is.na(df_final$status), ]
df_final <- df_final[df_final$time > 0, ]

message("Pacientes disponibles tras el merge: ", nrow(df_final))

# 5. Estandarización y agrupación de variables: 
df_final$age_scaled <- scale(as.numeric(df_final$age_at_diagnosis))
df_final$CD8_scaled <- scale(df_final$`T.cells.CD8`)
df_final$M1_scaled <- scale(df_final$Macrophages.M1)
df_final$M2_scaled <- scale(df_final$Macrophages.M2)

# Crear el Ratio M1/M2 (añadimos 0.001 para evitar división por cero)
df_final$ratio_M1_M2 <- df_final$Macrophages.M1 / (df_final$Macrophages.M2 + 0.001)
df_final$ratio_scaled <- scale(df_final$ratio_M1_M2)

df_final$ajcc_group <- as.character(df_final$ajcc_pathologic_stage)
df_final$ajcc_group[df_final$ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
df_final$ajcc_group[df_final$ajcc_pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC")] <- "Stage II"
df_final$ajcc_group[df_final$ajcc_pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"
df_final$ajcc_group[df_final$ajcc_pathologic_stage %in% c("Stage IV")] <- "Stage IV"
df_final$ajcc_group[df_final$ajcc_pathologic_stage %in% c("Stage 0")] <- "Stage I"

df_final <- df_final[!is.na(df_final$ajcc_group), ]
df_final$ajcc_group <- factor(df_final$ajcc_group, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))

message("Muestras finales para el modelo de Cox: ", nrow(df_final))
print(table(df_final$ajcc_group, useNA = "always"))

# 6. Modelo Cox CD8+: 
formula_cox_cd8 <- Surv(time, status) ~ CD8_scaled + age_scaled + ajcc_group
modelo_cox_cd8 <- coxph(formula_cox_cd8, data = df_final)
forest_cd8 <- ggforest(modelo_cox_cd8, data = df_final, main = "Hazard Ratio: Infiltración CD8+")
# Ver resumen
summary(modelo_cox_cd8)

# 7. Modelo de Cox M1: 
formula_cox_m1 <- Surv(time, status) ~ M1_scaled + age_scaled + ajcc_group
modelo_cox_m1 <- coxph(formula_cox_m1, data = df_final)
forest_m1 <- ggforest(modelo_cox_m1, data = df_final, main = "Hazard Ratio: Infiltración Macrófagos M1")
# Ver resumen
summary(modelo_cox_m1)

# 8. Modelo Cox M2: 
formula_cox_m2 <- Surv(time, status) ~ M2_scaled + age_scaled + ajcc_group
modelo_cox_m2 <- coxph(formula_cox_m2, data = df_final)
forest_m2 <- ggforest(modelo_cox_m2, data = df_final, main = "Hazard Ratio: Infiltración Macrófagos M2")
# Ver resumen
summary(modelo_cox_m2)


# Asegurar carpeta y guardar
ggsave(here("results/plots/3.forest_plot_cd8.pdf"), forest_cd8, width = 10, height = 6)
ggsave(here("results/plots/4.forest_plot_m1.pdf"), forest_m1, width = 10, height = 6)
ggsave(here("results/plots/5.forest_plot_m2.pdf"), forest_m2, width = 10, height = 6)

message("Análisis de supervivencia completado. Archivos guardados en results/plots/")




#------------------------------------------------------------------------------
########  BLOQUE 6: CURVAS DE KAPLAN-MEIER  -----------------------------------
library(survminer)

message("--- CURVA DE KAPLAN-MEIER: CD8+ ---")

# 1. Estratificación basada en la mediana: Alta/Baja infiltración 
mediana_cd8 <- median(df_final$`T.cells.CD8`, na.rm = TRUE)
df_final$CD8_group <- ifelse(df_final$`T.cells.CD8` > mediana_cd8, "Alto CD8", "Bajo CD8")
df_final$CD8_group <- factor(df_final$CD8_group, levels = c("Bajo CD8", "Alto CD8"))

# 2. Ajustar el modelo de Kaplan-Meier
fit_cd8 <- survfit(Surv(time, status) ~ CD8_group, data = df_final)

# 3. Graficar
km_cd8 <- ggsurvplot(
  fit_cd8, 
  data = df_final,
  pval = TRUE,             # P-valor del log-rank test
  conf.int = TRUE,         # Intervalos de confianza
  risk.table = TRUE,       # Tabla de pacientes en riesgo
  palette = c("#F8AFA8", "#89CFF0"),
  ggtheme = theme_minimal(base_size = 14),
  title = "Supervivencia según infiltración de T CD8+",
  xlab = "Días",
  ylab = "Probabilidad de supervivencia",
  legend.title = "Grupo CD8",
  legend.labs = c("Bajo", "Alto"), 
  font.main = c(16, "bold"),
  risk.table.fontsize = 4
)

# 4. Guardar
here("results/plots/...")
pdf('6.kaplan_meier_cd8.pdf', width= 10, height=8)
print(km_cd8, newpage= FALSE)
dev.off()


library(gridExtra)
archivo_pdf <- here("results/plots/6.kaplan_meier_cd8.pdf")
pdf(archivo_pdf, width = 9, height = 7)
arrange_ggsurvplots(list(km_cd8), print = TRUE)
dev.off()
message("Curva Kaplan-Meier guardada en results/plots/")


message("--- CURVA DE KAPLAN-MEIER: ESTADIOS TUMORALES ---")
# 1. Crear una variable binaria para los estadios
df_final$stage_binary <- ifelse(df_final$ajcc_group %in% c("Stage I", "Stage II"), 
                                "Temprano (I-II)", "Avanzado (III-IV)")
df_final$stage_binary <- factor(df_final$stage_binary, levels = c("Temprano (I-II)", "Avanzado (III-IV)"))

# 2. Ajustar el modelo de Kaplan-Meier para estadios
fit_km_stage <- survfit(Surv(time, status) ~ stage_binary, data = df_final)

# 3. Graficar 
km_plot_stage <- ggsurvplot (
  fit_km_stage, 
  data = df_final,
  pval = TRUE,              # P-valor para ver si la diferencia es significativa
  conf.int = TRUE,          # Intervalos de confianza
  risk.table = TRUE,        # Tabla de pacientes en riesgo
  palette = c("#F7D794", "#94B49F"), 
  ggtheme = theme_minimal(base_size = 14),
  title = "Supervivencia según Estadio Clínico",
  xlab = "Días",
  ylab = "Probabilidad de supervivencia",
  legend.title = "Estadio",
  legend.labs = c("Temprano (I-II)", "Avanzado (III-IV)"), 
  font.main = c(16, "bold")
)

# 4. Guardar
pdf("results/plots/7.kaplan_meier_stage.pdf", width= 10, height=8)
print(km_plot_stage, newpage= FALSE)
dev.off()
message("Gráfico de estadios guardado en results/plots/")


message("--- CURVA DE KAPLAN-MEIER: MACRÓFAGOS M1 ---")
# 1. Estratificación
df_final$Macrophages.M1 <- as.numeric(as.character(df_final$Macrophages.M1))
mediana_M1 <- median(df_final$Macrophages.M1, na.rm = TRUE)
df_final$M1_group <- ifelse(df_final$Macrophages.M1 > mediana_M1, "Alto", "Bajo")
df_final$M1_group <- factor(df_final$M1_group, levels = c("Bajo", "Alto"))

# 2. Ajustar modelo
fit_m1 <- survfit(Surv(time, status) ~ M1_group, data = df_final)

# 3. Graficar
km_m1 <- ggsurvplot(
  fit_m1, data = df_final, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
  palette = c("#F7D794", "#778beb"), # Usamos verde para M1 (color "saludable/bueno")
  title = "Supervivencia según infiltración de Macrófagos M1",
  xlab = "Días", ylab = "Probabilidad de supervivencia",
  legend.title = "Nivel M1", legend.labs = c("Bajo", "Alto")
)

print(km_m1)

# 4. Guardar
pdf("results/plots/8.kaplan_meier_m1.pdf", width= 10, height=8)
print(km_m1, newpage= FALSE)
dev.off()

message("--- CURVA DE KAPLAN-MEIER: MACRÓFAGOS M2 ---")
# 1. Estratificación basada en la mediana: Alta/Baja infiltración 
df_final$Macrophages.M2 <- as.numeric(as.character(df_final$Macrophages.M2))
mediana_M2 <- median(df_final$`Macrophages.M2`, na.rm = TRUE)
df_final$M2_group <- ifelse(df_final$`Macrophages.M2` > mediana_M2, "Alto", "Bajo")
df_final$M2_group <- factor(df_final$M2_group, levels = c("Bajo", "Alto"))

# 2. Ajustar el modelo de Kaplan-Meier
fit_m2 <- survfit(Surv(time, status) ~ M2_group, data = df_final)

# 3. Graficar
km_m2 <- ggsurvplot(
  fit_m2, 
  data = df_final,
  pval = TRUE,             # P-valor del log-rank test
  conf.int = TRUE,         # Intervalos de confianza
  risk.table = TRUE,       # Tabla de pacientes en riesgo
  palette = c("#F7D794", "#778beb"),
  ggtheme = theme_minimal(base_size = 14),
  title = "Supervivencia según infiltración de Macrófagos",
  xlab = "Días",
  ylab = "Probabilidad de supervivencia",
  legend.title = "Nivel M2",
  legend.labs = c("Bajo", "Alto"), 
  font.main = c(16, "bold"),
  risk.table.fontsize = 4
)

print(km_m2)


# 4. Guardar
pdf("results/plots/9.kaplan_meier_m2.pdf", width= 10, height=8)
print(km_m2, newpage= FALSE)
dev.off()
message("Curva Kaplan-Meier guardada en results/plots/")



message("--- CURVA DE KAPLAN-MEIER: RATIO M1/M2 ---")
mediana_ratio <- median(df_final$ratio_M1_M2, na.rm = TRUE)
df_final$ratio_group <- ifelse(df_final$ratio_M1_M2 > mediana_ratio, "Ratio Alto (M1)", "Ratio Bajo (M2)")
df_final$ratio_group <- factor(df_final$ratio_group, levels = c("Ratio Bajo (M2)", "Ratio Alto (M1)"))

fit_ratio <- survfit(Surv(time, status) ~ ratio_group, data = df_final)
km_ratio <- ggsurvplot(fit_ratio, data = df_final, pval = TRUE, conf.int = TRUE,risk.table = TRUE,
                       palette = c("#E76F51", "#2A9D8F"), 
                       title = "Supervivencia según Ratio M1/M2",
                       legend.labs = c("Dominio M2 (Bajo)", "Dominio M1 (Alto)"),
                       font.main = c(16, "bold"),
                       risk.table.fontsize = 4
)

print(km_ratio)

# 4. Guardar
pdf("results/plots/10.kaplan_meier_ratio.pdf", width= 10, height=8)
print(km_ratio, newpage= FALSE)
dev.off()
message("Curva Kaplan-Meier guardada en results/plots/")


#------------------------------------------------------------------------------
########  BLOQUE 7: ANÁLISIS MULTIVARIANTE  ---------------------------
message("--- INICIO DE ANÁLISIS MULTIVARIANTE ---")
library(survival)
library(survminer)
library(car)

# 1. Anotamos de nuevo las variables: 
df_final$age_scaled <- scale(as.numeric(df_final$age_at_diagnosis))
df_final$CD8_scaled <- scale(df_final$`T.cells.CD8`)
df_final$M1_scaled  <- scale(df_final$Macrophages.M1)
df_final$M2_scaled  <- scale(df_final$Macrophages.M2)
df_final$ratio_scaled <- scale(df_final$ratio_M1_M2)

# 2. Modelo multivariante: 
formula_multi <- Surv(time, status) ~ CD8_scaled + M1_scaled + M2_scaled + age_scaled + ajcc_group
modelo_multi <- coxph(formula_multi, data = df_final)

# Ver resumen: 
summary(modelo_multi)


# 3. Modelo multivariante con el ratio M1/M2. 
formula_multi_ratio <- Surv(time, status) ~ CD8_scaled + ratio_scaled + age_scaled + ajcc_group
modelo_multi_ratio <- coxph(formula_multi_ratio, data = df_final)
# Ver resumen: 
summary(modelo_multi_ratio)

# 4. Validación: Proporcionalidad de riesgos y Multicolinealidad
message("Validación del modelo (Test Schoenfeld y VIF):")
print(cox.zph(modelo_multi))
print(vif(modelo_multi))


# 5. Forest Plot: 
forest_multi <- ggforest(modelo_multi, data = df_final, main = "HR: Modelo Completo")
forest_ratio <- ggforest(modelo_multi_ratio, data = df_final, main = "HR: Modelo con Ratio M1/M2")


# Guardado: 
ggsave(here("results/plots/11.forest_plot_multivariante.pdf"), plot = forest_multi, width = 10, height = 6)
ggsave(here("results/plots/12.forest_plot_multivariante_ratio.pdf"), plot = forest_ratio, width = 10, height = 6)
message("Análisis multivariante completado y guardado.")




#------------------------------------------------------------------------------
########  BLOQUE 8: CORRELACIÓN E HETEROGENEIDAD INMUNE ----------------------
message("--- INICIO COCORRELACIÓN E HETEROGENEIDAD INMUNE ---")
library(ggplot2)
library(ggpubr)
library(corrplot)

# 1. Cálculo de Spearman: 
res_cor <- cor.test(df_final$`T.cells.CD8`, 
                    df_final$Macrophages.M2, 
                    method = "spearman", 
                    exact = FALSE)

# 2. Resultados: 
message("Rho de Spearman: ", round(res_cor$estimate, 3))
message("P-value: ", format.pval(res_cor$p.value, digits = 3))

# 3. Gráfico CD8+ y M2: 
cor_plot <- ggplot(df_final, aes(x = `T.cells.CD8`, y = Macrophages.M2)) +
  geom_point(aes(color = ajcc_group), alpha = 0.6) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(title = "Correlación entre Células T CD8+ y Macrófagos M2",
       subtitle = paste("Spearman Rho =", round(res_cor$estimate, 2), 
                        "| p-value =", format.pval(res_cor$p.value, digits = 2)),
       x = "Infiltración Células T CD8+",
       y = "Infiltración Macrófagos M2",
       color = "Estadio") +
  theme_minimal()

# Guardar
ggsave(here("results/plots/13.correlacion_CD8_M2.pdf"), plot = cor_plot, width = 8, height = 6)
print(cor_plot)

# 4. Matriz de Correlación para ver asociaciones entre todas las poblaciones
vars_inmunes <- df_final[, c("T.cells.CD8", "Macrophages.M1", "Macrophages.M2", "ratio_M1_M2")]
cor_matrix <- cor(vars_inmunes, method = "spearman", use = "pairwise.complete.obs")

png(here("results/plots/14.matriz_correlacion_inmune.png"), width = 800, height = 800)
corrplot(cor_matrix, method = "color", type = "upper", addCoef.col = "black", tl.col = "black", title = "Matriz de correlación inmune")
dev.off()



message("Gráfico de correlación guardado y valores calculados.")



#------------------------------------------------------------------------------
########  BLOQUE 9: MODELO PREDICTIVO: SCORE DE RIESGO  ---------------------------
message("--- INICIO DEL MODELO PREDICTIVO ---")
library(pROC)

# 1. Score de riesgo: 
df_final$risk_score <- predict(modelo_multi_ratio, type = "lp")

# 2.Estratificación: 
mediana_score <- median(df_final$risk_score, na.rm = TRUE)
df_final$risk_group <- factor(ifelse(df_final$risk_score > mediana_score, "Alto Riesgo", "Bajo Riesgo"), 
                              levels = c("Bajo Riesgo", "Alto Riesgo"))

 # Ajustar el modelo de Kaplan-Meier: 
fit_risk <- survfit(Surv(time, status) ~ risk_group, data = df_final)

# 3. Gráfico de curva Kaplan-Meier: 
km_risk <- ggsurvplot (
  fit_risk, 
  data = df_final,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#C3B1E1", "#E9C4A6"),
  ggtheme = theme_minimal(base_size = 14) ,
  title = "Supervivencia según Score de Riesgo (Ratio M1/M2)",
  subtitle = "Basado en CD8, Ratio M1/M2, Edad y Estadio",
  xlab = "Días",
  ylab = "Probabilidad de supervivencia",
  legend.title = "Grupo",
  legend.labs = c("Bajo Riesgo", "Alto Riesgo"), 
  risk.table.labs = c("Bajo Riesgo", "Alto Riesgo"), 
  font.main = c(16, "bold")
) 

 
# Guardar
pdf("results/plots/15.kaplan_meier_score_riesgo.pdf", width= 10, height=8)
print(km_risk, newpage= FALSE)
dev.off()
message("Curva Kaplan-Meier guardada en results/plots/")


# 4. Análisis de la curva ROC
  # Crear el objeto ROC: comparación entre "risk_score" y "status" (vital_status)
roc_obj <- roc(df_final$status, df_final$risk_score)

  # 2. Calcular el AUC
auc_val <- auc(roc_obj)
print(paste("El AUC de tu modelo es:", round(auc_val, 3)))

# 3. Graficar la Curva ROC
roc_plot <- ggroc(roc_obj, color = "#F8AFA8", size = 1.5) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray") +
  theme_minimal(base_size = 14) +
  labs(title = "Capacidad Predictiva del Modelo (ROC)",
       subtitle = paste("AUC =:", round(auc_val, 3)))

# 4. Guardar
ggsave(here("results/plots/16.curva_roc_modelo.pdf"), plot = roc_plot, width = 7, height = 6)
message("Curva ROC guardada. ¡Tu modelo tiene un poder predictivo de ", round(auc_val, 2), "!")


#------------------------------------------------------------------------------
########  BLOQUE 10: DESCRIPCIÓN DE LA COHORTE ANALIZADA  ----------------------
library(gt)
library(dplyr)

message("--- TABLA DE FLUJO DE MUESTRAS  ---")
# 1. Extracción de valores 
val_metastaticas <- ncol(data) # Muestras descargadas 
val_genes_brutos <- nrow(counts) # Genes originales
val_genes_filtrados <- num_genes # Genes tras filtrado 
val_genes_finales <- nrow(matriz_final) # Genes tras limpiar IDs
val_finales_analisis <- nrow(df_final) # Muestras finales (tras merge y limpieza)
val_excluidas <- val_metastaticas - val_finales_analisis
val_eventos <- sum(df_final$status == 1) # Cuenta los "1" en la columna status

# 2. Crear el dataframe 
datos_flujo <- data.frame(
  Etapa = c(
    "Muestras Metastásicas obtenidas (TCGA-SKCM)", 
    "Genes en la matriz bruta de conteos",
    "Genes tras el filtrado: CPM >1 en ≥20 % muestras",
    "Genes en la matriz normalizada final (TMM-CPM)",
    "Excluidas (datos clínicos faltantes o tiempo ≤ 0)",
    "Muestras finales para análisis",
    "Eventos observados (fallecimientos)"
  ),
  N = c(val_metastaticas, val_genes_brutos, val_genes_filtrados, 
        val_genes_finales, val_excluidas, val_finales_analisis, val_eventos)
)

# 3. Generar la tabla profesional con gt
tabla_flujo_final <- datos_flujo %>%
  gt() %>%
  tab_header(
    title = "Flujo del Proceso Analítico",
    subtitle = "Recuento automatizado de muestras y genes"
  ) %>%
  cols_label(
    Etapa = "Etapa del flujo analítico",
    N = "N"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  fmt_number(columns = N, decimals = 0, sep_mark = ".") %>%
  tab_options(table.width = pct(90))

# 4. Guardar
if (!dir.exists(here("results/tables"))) dir.create(here("results/tables"), recursive = TRUE)
gtsave(tabla_flujo_final, here("results/tables/Tabla_Flujo_Analitico_Completa.docx"))

message("Tabla de flujo generada dinámicamente con N=", val_finales_analisis)


message("--- PACIENTES INCLUIDOS EN ANÁLISIS DE SUPERVIVENCIA ---")
library(gtsummary)
library(gt)


# 1. Preparación limpia
df_tabla <- df_final %>%
  mutate(
    vital_status = factor(vital_status, levels = c("Alive", "Dead"), labels = c("Vivo", "Fallecido")),
    `T.cells.CD8` = `T.cells.CD8` * 100,
    Macrophages.M1 = Macrophages.M1 * 100,
    Macrophages.M2 = Macrophages.M2 * 100
  ) %>%
  select(age_at_diagnosis, gender, ajcc_group, `T.cells.CD8`, Macrophages.M1, Macrophages.M2, vital_status)

# 2. Creación de la tabla
tabla1 <- df_tabla %>%
  tbl_summary(
    by = vital_status,
    label = list(
      age_at_diagnosis ~ "Edad al diagnóstico",
      ajcc_group       ~ "Estadio AJCC",
      `T.cells.CD8`    ~ "Infiltración T CD8+ (%)",
      Macrophages.M1   ~ "Macrófagos M1 (%)",
      Macrophages.M2   ~ "Macrófagos M2 (%)"
    ),
    statistic = all_continuous() ~ "{mean} (± {sd})",
    missing = "no" 
  ) %>%
  add_p() %>%
  bold_labels()

# 3. Guardado
tabla1_gt <- as_gt(tabla1)
gtsave(tabla1_gt, here("results/tables/Descripción de Cohorte.docx"))
message("Tabla Descripción de Cohorte generada correctamente con N=78.")

# Tabla 2: Estadios (Distribución)
tabla_estadios <- df_final %>%
  group_by(ajcc_group) %>%
  summarise(n = n(), Porcentaje = round(n()/nrow(df_final)*100, 1)) %>%
  gt() %>%
  tab_header(title = "Distribución de Pacientes por Estadio")

gtsave(tabla_estadios, here("results/tables/Tabla_Estadios.docx"))



# Resumen estadístico de resultados cibersort:
# 1. Cargar la librería necesaria
library(dplyr)

# 2. Cargar los resultados de CIBERSORT
res_ciber <- read.csv(here("results/tables/Resultados_TILs_SKCM.csv"), row.names = 1)

# 3. Calcular medias y desviación estándar para las columnas de interés
# Nota: Asegúrate de que los nombres coincidan exactamente con las columnas del CSV
resumen_stats <- res_ciber %>%
  summarise(
    Media_CD8 = mean(`T.cells.CD8`, na.rm = TRUE),
    DE_CD8 = sd(`T.cells.CD8`, na.rm = TRUE),
    
    Media_M1 = mean(Macrophages.M1, na.rm = TRUE),
    DE_M1 = sd(Macrophages.M1, na.rm = TRUE),
    
    Media_M2 = mean(Macrophages.M2, na.rm = TRUE),
    DE_M2 = sd(Macrophages.M2, na.rm = TRUE)
  )

# 4. Ver resultados
print(resumen_stats)

# Opcional: Guardar este resumen en un archivo de texto para tu TFM
write.csv(resumen_stats, here("results/tables/Resumen_Estadistico_TILs.csv"), row.names = FALSE)









# IMPORTANTE: Cerramos el log_run
close_log <- function() {
  sink(type = "message")
  sink(type = "output")
  close(con_log)
  message("Log cerrado correctamente en: ", log_file)
}

