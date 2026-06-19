# TFM-Melanoma-TILs_UAX
# Impacto de la Abundancia de Subpoblaciones de TILs en la Supervivencia del Melanoma

Este repositorio contiene el flujo de trabajo computacional y estadístico desarrollado para mi Trabajo de Fin de Máster (TFM).

## 📊 Estructura del Proyecto
* `scripts/`: Scripts ordenados en R para la descarga (TCGAbiolinks), deconvolución (CIBERSORTx/quanTIseq) y modelos de Cox.
* `results/`: Tablas de caracterización clínica y gráficos de supervivencia generados.

## 🛠️ Requisitos y Dependencias
Para replicar este análisis, es necesario disponer de R (v4.x) y las siguientes librerías:
* `TCGAbiolinks`
* `CIBERSORTx` (Requiere token/matriz LM22 oficial)
* `survival` y `survminer`
* `here`, `dplyr`, `gt`

## 🚀 Instrucciones de Uso
El proyecto está estructurado utilizando la librería `here`. Basta con clonar este repositorio, mantener la estructura de directorios y ejecutar los scripts en orden numérico. Los datos se descargarán automáticamente desde el GDC Data Portal.
