# =============================================================================
# Análisis de Expresión Diferencial - GSE240033
#Comparación:
# Osteosarcoma Tissue vs Skin
# =============================================================================

# 1. LIBRERIAS 
library(tidyverse)
library(data.table)
library(pheatmap)
library(ggplot2)
library(limma)
library(EnhancedVolcano)
library(GEOquery)
library(corrplot)

# 2. CARGA Y EXPLORACIÓN INICIAL DE LOS DATOS
expr_RNA_69 <- read.delim("GSE240033_VIGOR.txt.gz", check.names = FALSE)
expr_RNA_69 <- as.data.frame(expr_RNA_69)

# Exploración inicial
head(expr_RNA_69[, 1:5])
colnames(expr_RNA_69)[1:10]
dim(expr_RNA_69)

# Verificar si son datos crudos (conteos enteros)
mat_check <- as.matrix(expr_RNA_69[, -1])   # excluir columna de genes
any(mat_check %% 1 != 0)                    # FALSE = son conteos enteros
sum(mat_check %% 1 != 0)
range(mat_check)
head(sort(unique(as.vector(mat_check))), 20)

# 3. PREPARACIÓN DE LA MATRIZ DE EXPRESIÓN 
gene_names <- as.character(expr_RNA_69[["Name"]]) #guardar nombres de genes

expr_RNA_69 <- expr_RNA_69[, colnames(expr_RNA_69) != "Name"] #eliminar columna Name
expr_RNA_69 <- as.data.frame(
  lapply(expr_RNA_69, function(x) as.numeric(as.character(x)))) #forzar valor numérico en todas las col.
  
gene_names_unique <- make.unique(gene_names, sep = ".") #reasignar/eliminar duplicados

rownames(expr_RNA_69) <- gene_names_unique

# Verificación
stopifnot(head(rownames(expr_RNA_69), 3) == head(gene_names_unique, 3))
cat("Rownames asignados correctamente:",head(rownames(expr_RNA_69), 5), "\n") 
head(rownames(expr_RNA_69))
nrow(expr_RNA_69)
summary(as.vector(as.matrix(expr_RNA_69)))

# Visualización distribución cruda
hist(as.numeric(as.matrix(expr_RNA_69)), breaks = 100,
     main = "Distribución cruda", xlab = "Counts")
boxplot(log2(expr_RNA_69 + 1), outline = FALSE,
        main = "Boxplot log2(counts+1)", las = 2, cex.axis = 0.6)

# 4. METADATOS (GEO)
options(download.file.method = "auto") #mejora la conexión
gseRNA69     <- getGEO("GSE240033", GSEMatrix = TRUE)
gseRNA69     <- gseRNA69[[1]]
metadataRNA69 <- pData(gseRNA69)

head(metadataRNA69)
colnames(metadataRNA69)
table(metadataRNA69$`tissue:ch1`)

# 5. cORRECCIONES: ALINEACIÓN SAMPLE IDs 
metadataRNA69$sampleID <- sub("_S[0-9]+$", "", metadataRNA69$title)

# Diagnóstico de coincidencia
sum(metadataRNA69$sampleID %in% colnames(expr_RNA_69))
setdiff(metadataRNA69$sampleID, colnames(expr_RNA_69))
setdiff(colnames(expr_RNA_69), metadataRNA69$sampleID)

# CORRECCIÓN: renombrar las tres muestras duplicadas que quedaron
metadataRNA69$sampleID[
  metadataRNA69$sampleID %in%
    c("MN03_bone_pre", "MN05_bone_post", "MN11_bone_pre")
] <- c("MN03_bone_pre_2", "MN05_bone_post_2", "MN11_bone_pre_2")

# Verificar que ahora coincidan las 69 muestras
sum(metadataRNA69$sampleID %in% colnames(expr_RNA_69))

# 6. FILTRAR LÍNEAS CELULARES
# CORRECCIÓN: "celllLine" → "cellLine"
metadataRNA69$`tissue:ch1` <- gsub(
  "celllLine", "cellLine",
  metadataRNA69$`tissue:ch1`)

metadata_filtrado <- metadataRNA69[
  metadataRNA69$`tissue:ch1` != "cellLine", ]

table(metadata_filtrado$`tissue:ch1`)

# 7. MATRIZ FILTRADA DE EXPRESIÓN (alineada con metadata)
expr_filtrado <- expr_RNA_69[
  , match(metadata_filtrado$sampleID, colnames(expr_RNA_69))]

# Verificación de alineación
all(colnames(expr_filtrado) == metadata_filtrado$sampleID)

dim(expr_filtrado)
dim(metadata_filtrado)

# 8. SUBCONJUNTO TUMOR vs SKIN 
metadata_TS <- metadata_filtrado[
  metadata_filtrado$`tissue:ch1` %in% c("osteosarcomaTissue", "skin"), ]

table(metadata_TS$`tissue:ch1`)

expr_TS <- expr_filtrado[, metadata_TS$sampleID]
dim(expr_TS)

# 9. TRANSFORMACIÓN Y FILTRADO DE GENES 
expr_log <- log2(expr_TS + 1)

# Filtrar genes con expresión media > 1 (en escala log2)
keep      <- rowMeans(expr_log) > 1
expr_filt <- expr_log[keep, ]

cat("Genes antes del filtrado:", nrow(expr_log), "\n")
cat("Genes después del filtrado:", nrow(expr_filt), "\n")
summary(rowMeans(expr_log))

# 10. PREPROCESAMIENTO Y CONTROL DE CALIDAD 
dev.new(width = 16, height = 8) #Abrir otra ventada
grupo      <- metadata_TS$`tissue:ch1`
col_tejido <- ifelse(grupo == "osteosarcomaTissue", "#E74C3C", "#3498DB")

ann_colors_qc <- list(
  Tejido = c(osteosarcomaTissue = "#E74C3C", skin = "#3498DB"))

# ETAPA 1 — DATOS SIN PROCESAR (expr_TS: conteos en escala original)
##Boxplot
boxplot(log2(expr_TS + 1),
        las      = 2,
        cex.axis = 0.45,
        outline  = FALSE,
        col      = col_tejido,
        main     = "Datos preliminares — Distribución por muestra/GSE240033",
        ylab     = "log2(counts + 1) [solo visualización]")

legend("topright",
       legend = c("Osteosarcoma", "Skin"),
       fill   = c("#E74C3C", "#3498DB"),
       bty    = "n", cex = 0.8)

##Gráfico de densidad
expr_TS_log <- log2(expr_TS + 1)  
n_muestras  <- ncol(expr_TS_log)
cols_dens   <- colorRampPalette(c("#E74C3C", "#3498DB"))(n_muestras)

plot(density(expr_TS_log[, 1], na.rm = TRUE),
     main  = "Datos preliminares — Densidad de expresión\nGSE240033",
     xlab  = "log2(counts + 1)",
     ylab  = "Densidad",
     col   = cols_dens[1],
     lwd   = 1.2,
     ylim  = c(0, 0.5))
for (i in 2:n_muestras) {
  lines(density(expr_TS_log[, i], na.rm = TRUE),
        col = cols_dens[i], lwd = 1.2)
}
legend("topright",
       legend = c("Osteosarcoma", "Skin"),
       col    = c("#E74C3C", "#3498DB"),
       lwd    = 2, bty = "n", cex = 0.8)

##Corrplot 
pdf(
  "correlacion_GSE240033.pdf",
  width = 14,
  height = 14
)

corrplot.mixed(
  cor_preliminar,
  
  lower = "number",
  upper = "circle",
  
  lower.col = "black",
  upper.col = colorRampPalette(
    c("white", "#1A3A6B")
  )(200),
  
  number.cex = 0.40,
  
  tl.cex = 0.60,
  tl.col = "black",
  tl.srt = 45,
  tl.pos = "lt",
  
  cl.cex = 0.80,
  
  diag = "n",
  
  mar = c(0, 0, 3, 0)
)

title(
  "Correlación entre muestras\nGSE240033",
  cex.main = 1.3
)

dev.off()


# ETAPA 2 — DATOS PROCESADOS (expr_filt: transformados + filtrados)
##Boxplot 
boxplot(expr_filt,
        las      = 2,
        cex.axis = 0.45,
        outline  = FALSE,
        col      = col_tejido,
        main     = "DESPUÉS del procesamiento — Distribución por muestra\nGSE240033",
        ylab     = "log2(counts + 1)")

legend("topright",
       legend = c("Osteosarcoma", "Skin"),
       fill   = c("#E74C3C", "#3498DB"),
       bty    = "n", cex = 0.8)

##Gráfico de densidad 
plot(density(expr_filt[, 1], na.rm = TRUE),
     main  = "DESPUÉS del procesamiento — Densidad de expresión\nGSE240033",
     xlab  = "log2(counts + 1)",
     ylab  = "Densidad",
     col   = cols_dens[1],
     lwd   = 1.2,
     ylim  = c(0, 0.5))
for (i in 2:ncol(expr_filt)) {
  lines(density(expr_filt[, i], na.rm = TRUE),
        col = cols_dens[i], lwd = 1.2)}
legend("topright",
       legend = c("Osteosarcoma", "Skin"),
       col    = c("#E74C3C", "#3498DB"),
       lwd    = 2, bty = "n", cex = 0.8)

##Corplot 
cor_procesado <- cor(expr_filt, use = "pairwise.complete.obs")

corrplot.mixed(cor_procesado,
               lower         = "number",
               upper         = "circle",
               lower.col     = "black",
               upper.col     = colorRampPalette(c("white", "#1A3A6B"))(200),
               number.cex    = 0.45,
               tl.cex        = 0.45,
               tl.col        = "red",
               tl.pos        = "lt",
               mar           = c(0, 0, 2, 0))
title("Procesado — Correlación entre muestras\nGSE240033",cex.main = 0.9)

# ETAPA 3 — Visualizaciones de estructura global (solo datos procesados)
##PCA (variabilidad global y separación entre grupos) 
pca     <- prcomp(t(expr_filt), scale. = TRUE)
var_exp <- pca$sdev^2 / sum(pca$sdev^2)

plot(100 * var_exp[1:10],
     type = "b",
     xlab = "Componente Principal",
     ylab = "% Varianza explicada",
     main = "Scree Plot — GSE240033",
     pch  = 19, col = "steelblue")

pca_df <- data.frame(
  Sample = colnames(expr_filt),
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2],
  Grupo  = grupo)

ggplot(pca_df, aes(PC1, PC2, color = Grupo, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 2.5) +
  scale_color_manual(
    values = c(osteosarcomaTissue = "#E74C3C", skin = "#3498DB"),
    labels = c("Osteosarcoma", "Skin")
  ) +
  labs(
    title = "PCA — Osteosarcoma vs Skin (GSE240033)",
    x     = paste0("PC1 (", round(100 * var_exp[1], 1), "%)"),
    y     = paste0("PC2 (", round(100 * var_exp[2], 1), "%)"),
    color = "Tejido"
  ) +
  theme_bw()

##Agrupamiento jerárquico (detección de muestras atípicas)
d  <- dist(t(expr_filt))
hc <- hclust(d, method = "complete")
plot(hc,
     main   = "Agrupamiento jerárquico — GSE240033",
     xlab   = "",
     sub    = "Distancia euclídea, enlace completo",
     cex    = 0.6,
     labels = paste0(colnames(expr_filt), "\n(", grupo, ")"))

##Heatmap de correlación 
ann_col_qc <- data.frame(
  Tejido    = grupo,
  row.names = colnames(expr_filt))

pheatmap(cor_procesado,
         annotation_col    = ann_col_qc,
         annotation_colors = ann_colors_qc,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = FALSE,
         show_colnames = FALSE,
         main          = "Heatmap de correlación — GSE240033 (procesado)",
         fontsize      = 8)

# 11. EXPRESIÓN DIFERENCIAL — LIMMA 
# Referencia = skin → logFC positivo = mayor expresión en osteosarcoma
grupoTS <- factor(
  metadata_TS$`tissue:ch1`,
  levels = c("skin", "osteosarcomaTissue"))

table(grupoTS)

design    <- model.matrix(~ grupoTS)
colnames(design)   # "(Intercept)" "grupoTSosteosarcomaTissue"

fit       <- lmFit(expr_filt, design)
fit       <- eBayes(fit)

coef_name      <- colnames(design)[2]

##Resultados
resultados_DEG <- topTable(
  fit,
  coef          = coef_name,
  number        = Inf,
  adjust.method = "BH",
  sort.by       = "P")

head(resultados_DEG, 20)
summary(resultados_DEG$adj.P.Val)

DEG_sig <- subset(
  resultados_DEG,
  adj.P.Val < 0.05 & abs(logFC) > 1)

cat("Total genes analizados:", nrow(resultados_DEG), "\n")
cat("DEGs significativos (FDR<0.05, |logFC|>1):", nrow(DEG_sig), "\n")
cat("  Sobreexpresados en osteosarcoma (logFC > 1):",sum(DEG_sig$logFC > 1), "\n")
cat("  Sobreexpresados en skin (logFC < -1):",sum(DEG_sig$logFC < -1), "\n")

# Visualización - Volcano plot
p_volcano <- EnhancedVolcano(resultados_DEG,
  lab         = rownames(resultados_DEG),
  x           = "logFC",
  y           = "adj.P.Val",
  pCutoff     = 0.05,
  FCcutoff    = 1,
  title       = "Osteosarcoma vs Skin",
  subtitle    = "GSE240033 · limma · BH · |logFC| > 1",
  xlab        = "log2 Fold Change (positivo = ↑ en Osteosarcoma)",
  col         = c("grey60", "steelblue", "orange", "firebrick"),
  colAlpha    = 0.6,
  pointSize   = 2
)

ggsave(
  "volcano_GSE240033_Osteosarcoma_vs_Skin.png",
  plot   = p_volcano,
  width  = 12,
  height = 10,
  units  = "in",
  dpi    = 300,
  bg     = "white")

##Heatmap Top DEGs 
n_top <- min(50, nrow(DEG_sig))

if (n_top >= 2) {
  top_genes <- rownames(head(DEG_sig[order(DEG_sig$adj.P.Val), ], n_top))
  
  ann_col_deg <- data.frame(
    Tejido    = grupoTS,
    row.names = colnames(expr_filt)
  )
  
  pheatmap(
    as.matrix(expr_filt[top_genes, ]),
    annotation_col    = ann_col_deg,
    annotation_colors = ann_colors_qc,
    scale             = "row",
    show_colnames     = FALSE,
    fontsize_row      = 6,
    cluster_cols      = TRUE,
    cluster_rows      = TRUE,
    main              = paste0("Top ", n_top,
                               " DEGs — Osteosarcoma vs Skin\nGSE240033")
  )
} else {
  cat("Menos de 2 DEGs significativos\n")
}

# 12. EXPORTAR RESULTADOS 
write.csv(resultados_DEG,"GSE240033_DEG_osteosarcoma_vs_skin_todos.csv",row.names = TRUE)
write.csv(DEG_sig,"GSE240033_DEG_osteosarcoma_vs_skin_sig.csv",row.names = TRUE)

