# =============================================================================
# Análisis de Expresión Diferencial — GSE239948
# Comparaciones:
#   a) Ubicación tumoral: Femur, Humerus, Radius, Tibia, Ulna (10 contrastes) - No se utilizo
#   b) Tamaño de raza: Giant vs Large
# =============================================================================

# 1. LIBRERÍAS 
library(tidyverse)
library(data.table)
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(corrplot)
library(EnhancedVolcano)

# 2. CARGA Y EXPLORACIÓN INICIAL DE LOS DATOS 
expr43 <- fread("GSE239948_CCOGC.txt.gz")
expr43 <- as.data.frame(expr43)

# Exploración inicial
dim(expr43)
head(expr43[, 1:5])
str(expr43[, 1:5])

# Verificar si son datos crudos (conteos enteros)
mat_check43 <- as.matrix(expr43[, -1])   # excluir columna de genes
any(mat_check43 %% 1 != 0)                    # FALSE = son conteos enteros
sum(mat_check43 %% 1 != 0)
range(mat_check43)
head(sort(unique(as.vector(mat_check43))), 20)

# 3. PREPARACIÓN DE LA MATRIZ DE EXPRESIÓN
colnames(expr43)[1] <- "Gene" #renombrar columna de genes

# resolver/eliminar duplicados conservando la fila de mayor expresión total
expr_unique43 <- expr43 %>%
  mutate(total_expression = rowSums(across(-Gene), na.rm = TRUE)) %>%
  group_by(Gene) %>%
  slice_max(total_expression, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Genes antes de deduplicación:", nrow(expr43), "\n")
cat("Genes tras deduplicación:", nrow(expr_unique43), "\n")

# construir matriz con genes como rownames
expr_unique43   <- as.data.frame(expr_unique43)
rownames(expr_unique43) <- expr_unique43$Gene

expr_matrix43 <- expr_unique43 %>%
  select(-Gene, -total_expression) %>%
  as.matrix()

mode(expr_matrix43) <- "numeric"

# Verificación
cat("NA en la matriz:", sum(is.na(expr_matrix43)), "\n")
cat("Rownames asignados:", head(rownames(expr_matrix43), 5), "\n")
cat("Dimensiones:", dim(expr_matrix43), "\n")
summary(as.numeric(expr_matrix43))

# Visualización distribución cruda
hist(as.numeric(as.matrix(expr43)), breaks = 100,
     main = "Distribución cruda", xlab = "Counts")

# 4. METADATOS (GEO)
options(download.file.method = "auto") #mejora la conexión
gse43      <- getGEO("GSE239948", GSEMatrix = TRUE)
gse43      <- gse43[[1]]
metadata43 <- pData(gse43)

colnames(metadata43)
head(metadata43[, c("title", "geo_accession")])

# Limpiar columna de ubicación tumoral
metadata43$location <- gsub(
  "tumor location: ", "",
  metadata43$`tumor location:ch1`)

print(table(metadata43$location))
print(table(metadata43$`breed:ch1`))

# 5. ALINEACIÓN DE SAMPLE IDs 
rownames(metadata43) <- metadata43$title

# Diagnóstico de coincidencia
cat("Muestras en matriz:", ncol(expr_matrix43), "\n")
cat("Muestras en metadata:", nrow(metadata43), "\n")
cat("Coincidencias:", sum(colnames(expr_matrix43) %in% rownames(metadata43)), "\n")

setdiff(colnames(expr_matrix43), rownames(metadata43))   # en matriz pero no en metadatos
setdiff(rownames(metadata43), colnames(expr_matrix43))   # en metadatos pero no en matriz

# Alinear metadata con el orden de columnas de la matriz
metadata43 <- metadata43[
  match(colnames(expr_matrix43), rownames(metadata43)), ]

(all(colnames(expr_matrix43) == rownames(metadata43)))

# ── 6. MATRIZ FILTRADA DE EXPRESIÓN (alineada con metadata)
# Eliminar genes con varianza = 0 (constantes en todas las muestras)
var_global <- apply(expr_matrix43, 1, var)
expr_matrix43 <- expr_matrix43[var_global != 0, ]

cat("Genes tras eliminar varianza=0:", nrow(expr_matrix43), "\n")
(sum(apply(expr_matrix43, 1, var) == 0) == 0)

# 8. SUBCONJUNTOS DE DATOS 
# (a).Subconjunto por ubicación tumoral
# Excluir "Other site" por ser una categoría heterogénea no interpretable
metadata_loc <- metadata43[
  metadata43$location %in%
    c("Femur", "Humerus", "Radius", "Tibia", "Ulna"), ]

metadata_loc$location <- droplevels(factor(metadata_loc$location))

print(table(metadata_loc$location))

expr_loc <- expr_matrix43[, rownames(metadata_loc)]
(all(colnames(expr_loc) == rownames(metadata_loc)))
cat("Dimensiones subconjunto ubicación:", dim(expr_loc), "\n")

# (b).Subconjunto por raza (Medianas-grandes y Gigantes)
metadata43$size_group <- case_when(
  metadata43$`breed:ch1` %in% c(
    "Great Dane", "Great Pyrenees", "Irish Wolfhound"
  ) ~ "Giant",
  
  metadata43$`breed:ch1` %in% c(
    "Mixed Breed", "Greyhound", "Labrador Retriever",
    "Golden Retriever", "Belgian Sheepdog", "Doberman Pinscher",
    "Old English Sheepdog", "Australian Cattle Dog",
    "Staffordshire Bull Terrier", "Weimaraner", "Rhodesian Ridgeback"
  ) ~ "Large",
  
  TRUE ~ NA_character_
)

# Razas sin clasificar
cat("Razas sin grupo asignado:\n")
print(unique(metadata43$`breed:ch1`[is.na(metadata43$size_group)]))

metadata_size <- metadata43[!is.na(metadata43$size_group), ]
metadata_size$size_group <- factor(metadata_size$size_group,
                                   levels = c("Large", "Giant"))

print(table(metadata_size$size_group))

expr_size <- expr_matrix43[, rownames(metadata_size)]
(all(colnames(expr_size) == rownames(metadata_size)))
cat("Dimensiones subconjunto raza:", dim(expr_size), "\n")

# 9. TRANSFORMACIÓN Y FILTRADO DE GENES 
# Los datos son conteos crudos (enteros) → transformación log2(x+1)
# para estabilizar varianza y reducir la asimetría positiva

##Transformación y filtrado — UBICACIÓN 
expr_loc_log  <- log2(expr_loc + 1)

# Conservar top 5000 genes más variables 
vars_loc      <- apply(expr_loc_log, 1, var)
top5000_loc   <- names(sort(vars_loc, decreasing = TRUE)[1:5000])
expr_loc_filt <- expr_loc_log[top5000_loc, ]

cat("Genes tras filtrado ubicación (top 5000):", nrow(expr_loc_filt), "\n")
summary(as.numeric(expr_loc_filt))

##Transformación y filtrado — RAZA 
expr_size_log  <- log2(expr_size + 1)

vars_size      <- apply(expr_size_log, 1, var)
top5000_size   <- names(sort(vars_size, decreasing = TRUE)[1:5000])
expr_size_filt <- expr_size_log[top5000_size, ]

cat("Genes tras filtrado raza (top 5000):", nrow(expr_size_filt), "\n")

# 10. PREPROCESAMIENTO Y CONTROL DE CALIDAD
# Colores por grupo de ubicación para visualizaciones
n_loc      <- length(levels(metadata_loc$location))
col_loc    <- setNames(
  RColorBrewer::brewer.pal(max(n_loc, 3), "Set1")[1:n_loc],
  levels(metadata_loc$location))
col_muestra_loc <- col_loc[as.character(metadata_loc$location)]

#Colores raza
grupo_raza      <- metadata_size$size_group
col_raza <- ifelse(grupo_raza == "Giant", "#E74C3C", "#3498DB")

ann_colors_Raza <- list(
  Raza = c(Giant = "#E74C3C", Large = "#3498DB"))

# ETAPA 1 — DATOS SIN PROCESAR 

##Boxplot — Ubicación 
par(mar = c(9, 5, 4, 10), xpd = TRUE)
boxplot(expr_loc_log,
        las      = 2,
        cex.axis = 0.45,
        outline  = FALSE,
        col      = col_muestra_loc,
        main     = paste0("Datos preliminares — Distribución por muestra",
                          "\nGSE239948 (Ubicación tumoral)"),
        ylab     = "log2(counts + 1)")
usr <- par("usr")
legend(x      = usr[2] * 1.01,
       y      = usr[4],
       legend = names(col_loc),
       fill   = col_loc,
       bty    = "n", cex = 0.8,
       title  = "Ubicación")
par(mar = c(5, 4, 4, 2), xpd = FALSE)
dev.off()

##Boxplot — Razas 
par(mar = c(9, 5, 4, 10), xpd = TRUE)
boxplot(expr_size_log,
        las      = 2,
        cex.axis = 0.45,
        outline  = FALSE,
        col      = col_raza,
        main     = paste0("Datos preliminares — Distribución por muestra",
                          "\nGSE239948 (Razas)"),
        ylab     = "log2(counts + 1)")
usr <- par("usr")
legend(x      = usr[2] * 1.01,
       y      = usr[4],
       legend = names(col_loc),
       fill   = col_loc,
       bty    = "n", cex = 0.8,
       title  = "Razas")
par(mar = c(5, 4, 4, 2), xpd = FALSE)
dev.off()

##Gráfico de densidad — Ubicación 
n_muestras_loc <- ncol(expr_loc_log)
cols_dens_loc  <- col_muestra_loc

plot(density(expr_loc_log[, 1], na.rm = TRUE),
     main  = paste0("Datos preliminares — Densidad de expresión",
                    "\nGSE239948 (Ubicación tumoral)"),
     xlab  = "log2(counts + 1)",
     ylab  = "Densidad",
     col   = cols_dens_loc[1],
     lwd   = 1.2,
     ylim  = c(0, 0.4))
for (i in 2:n_muestras_loc) {
  lines(density(expr_loc_log[, i], na.rm = TRUE),
        col = cols_dens_loc[i], lwd = 1.2)
}
legend("topright",
       legend = names(col_loc),
       col    = col_loc,
       lwd    = 2, bty = "n", cex = 0.8,
       title  = "Ubicación")

##Gráfico de densidad — Raza 
n_muestras_raza <- ncol(expr_size_log)
cols_dens_raza  <- col_raza

plot(density(expr_size_log[, 1], na.rm = TRUE),
     main  = paste0("Datos preliminares — Densidad de expresión",
                    "\nGSE239948 (Raza)"),
     xlab  = "log2(counts + 1)",
     ylab  = "Densidad",
     col   = col_raza[1],
     lwd   = 1.2,
     ylim  = c(0, 0.4))
for (i in 2:n_muestras_raza) {
  lines(density(expr_size_log[, i], na.rm = TRUE),
        col = col_raza[i], lwd = 1.2)
}
legend("topright",
       legend = names(col_raza),
       col    = col_raza,
       lwd    = 2, bty = "n", cex = 0.8,
       title  = "Raza")

##Corrplot — Ubicación 
cor_raza_antes <- cor(expr_size_log, use = "pairwise.complete.obs")

corrplot.mixed(cor_raza_antes,
               lower      = "number",
               upper      = "circle",
               lower.col  = "black",
               upper.col  = colorRampPalette(c("white", "#1A3A6B"))(200),
               number.cex = 0.4,
               tl.cex     = 0.4,
               tl.col     = "red",
               tl.pos     = "lt")
title(paste0("Datos preliminares — Correlación entre muestras",
             "\nGSE239948 (Raza)"))
dev.off()

# =============================================================================
# ETAPA B — DATOS PROCESADOS (expr_loc_filt: top 5000 genes más variables)
# =============================================================================

# --- B.1 Boxplot DESPUÉS del procesamiento — Ubicación -----------------------
png("boxplot_despues_ubicacion_GSE239948.png",
    width = 2400, height = 800, res = 150)
par(mar = c(9, 5, 4, 10), xpd = TRUE)
boxplot(expr_loc_filt,
        las      = 2,
        cex.axis = 0.45,
        outline  = FALSE,
        col      = col_muestra_loc,
        main     = paste0("DESPUÉS del procesamiento — Distribución por muestra",
                          "\nGSE239948 (Ubicación tumoral)"),
        ylab     = "log2(counts + 1)")
usr <- par("usr")
legend(x      = usr[2] * 1.01,
       y      = usr[4],
       legend = names(col_loc),
       fill   = col_loc,
       bty    = "n", cex = 0.8,
       title  = "Ubicación")
par(mar = c(5, 4, 4, 2), xpd = FALSE)
dev.off()

# --- B.2 Densidad DESPUÉS del procesamiento — Ubicación ----------------------
plot(density(expr_loc_filt[, 1], na.rm = TRUE),
     main  = paste0("DESPUÉS del procesamiento — Densidad de expresión",
                    "\nGSE239948 (Ubicación tumoral)"),
     xlab  = "log2(counts + 1)",
     ylab  = "Densidad",
     col   = cols_dens_loc[1],
     lwd   = 1.2,
     ylim  = c(0, 0.4))
for (i in 2:n_muestras_loc) {
  lines(density(expr_loc_filt[, i], na.rm = TRUE),
        col = cols_dens_loc[i], lwd = 1.2)
}
legend("topright",
       legend = names(col_loc),
       col    = col_loc,
       lwd    = 2, bty = "n", cex = 0.8,
       title  = "Ubicación")

# --- B.3 Corrplot DESPUÉS del procesamiento — Ubicación ----------------------
cor_loc_despues <- cor(expr_loc_filt, use = "pairwise.complete.obs")

png("corrplot_despues_ubicacion_GSE239948.png",
    width = 1800, height = 1800, res = 150)
corrplot.mixed(cor_loc_despues,
               lower      = "number",
               upper      = "circle",
               lower.col  = "black",
               upper.col  = colorRampPalette(c("white", "#1A3A6B"))(200),
               number.cex = 0.4,
               tl.cex     = 0.4,
               tl.col     = "red",
               tl.pos     = "lt")
title(paste0("DESPUÉS del procesamiento — Correlación entre muestras",
             "\nGSE239948 (Ubicación tumoral)"))
dev.off()

# =============================================================================
# ETAPA C — VISUALIZACIONES DE ESTRUCTURA GLOBAL (datos procesados)
# =============================================================================

# --- C.1 PCA — Ubicación tumoral ---------------------------------------------
top500_loc <- names(sort(apply(expr_loc_filt, 1, var),
                         decreasing = TRUE)[1:500])
pca_loc    <- prcomp(t(expr_loc_filt[top500_loc, ]), scale. = TRUE)
var_loc    <- round(100 * pca_loc$sdev^2 / sum(pca_loc$sdev^2), 1)

pca_df_loc <- data.frame(
  PC1      = pca_loc$x[, 1],
  PC2      = pca_loc$x[, 2],
  Location = metadata_loc$location
)

ggplot(pca_df_loc, aes(PC1, PC2, color = Location)) +
  geom_point(size = 3) +
  scale_color_manual(values = col_loc) +
  labs(
    title = "PCA — Ubicación tumoral (GSE239948)",
    x     = paste0("PC1 (", var_loc[1], "%)"),
    y     = paste0("PC2 (", var_loc[2], "%)")
  ) +
  theme_bw()

# --- C.2 Heatmap de distancias entre muestras --------------------------------
sample_dist <- dist(t(expr_loc_filt))

pheatmap(as.matrix(sample_dist),
         annotation_col = data.frame(
           Ubicacion  = metadata_loc$location,
           row.names  = colnames(expr_loc_filt)
         ),
         show_rownames = FALSE,
         show_colnames = FALSE,
         main          = "Distancias entre muestras — GSE239948 (Ubicación)")

# --- C.3 PCA — Raza ----------------------------------------------------------
top500_size <- names(sort(apply(expr_size_filt, 1, var),
                          decreasing = TRUE)[1:500])
pca_size    <- prcomp(t(expr_size_filt[top500_size, ]), scale. = TRUE)
var_size    <- round(100 * pca_size$sdev^2 / sum(pca_size$sdev^2), 1)

col_size <- c(Large = "#3498DB", Giant = "#E74C3C")

pca_df_size <- data.frame(
  PC1   = pca_size$x[, 1],
  PC2   = pca_size$x[, 2],
  Raza  = metadata_size$size_group
)

ggplot(pca_df_size, aes(PC1, PC2, color = Raza)) +
  geom_point(size = 3) +
  scale_color_manual(values = col_size) +
  labs(
    title = "PCA — Tamaño de raza (GSE239948)",
    x     = paste0("PC1 (", var_size[1], "%)"),
    y     = paste0("PC2 (", var_size[2], "%)")
  ) +
  theme_bw()

# 11. EXPRESIÓN DIFERENCIAL — LIMMA 
## ANÁLISIS A — UBICACIÓN TUMORAL (10 contrastes entre 5 localizaciones)
design_loc <- model.matrix(~ 0 + location, data = metadata_loc)
colnames(design_loc) <- levels(metadata_loc$location)

fit_loc  <- lmFit(expr_loc_filt, design_loc)

contrast_loc <- makeContrasts(
  Femur_vs_Tibia    = Femur - Tibia,
  Femur_vs_Humerus  = Femur - Humerus,
  Femur_vs_Radius   = Femur - Radius,
  Femur_vs_Ulna     = Femur - Ulna,
  Humerus_vs_Radius = Humerus - Radius,
  Humerus_vs_Tibia  = Humerus - Tibia,
  Humerus_vs_Ulna   = Humerus - Ulna,
  Radius_vs_Tibia   = Radius - Tibia,
  Radius_vs_Ulna    = Radius - Ulna,
  Tibia_vs_Ulna     = Tibia - Ulna,
  levels = design_loc
)

fit2_loc <- contrasts.fit(fit_loc, contrast_loc)
fit2_loc <- eBayes(fit2_loc)

# Función auxiliar para extraer resultados y exportar por contraste
extraer_DEG_localización <- function(fit, coef_name, archivo) {
  res <- topTable(fit,
                  coef          = coef_name,
                  number        = Inf,
                  adjust.method = "BH",
                  sort.by       = "P")
  deg <- subset(res, adj.P.Val < 0.05 & abs(logFC) > 1)
  cat(coef_name, "— DEGs significativos:", nrow(deg), "\n")
  write.csv(data.frame(Gene = rownames(res), res),
            archivo, row.names = FALSE)
  return(res)
}

res_Femur_Tibia    <- extraer_DEG_localización(fit2_loc, "Femur_vs_Tibia",
                                  "GSE239948_Femur_vs_Tibia.csv")
res_Femur_Humerus  <- extraer_DEG_localización(fit2_loc, "Femur_vs_Humerus",
                                  "GSE239948_Femur_vs_Humerus.csv")
res_Femur_Radius   <- extraer_DEG_localización(fit2_loc, "Femur_vs_Radius",
                                  "GSE239948_Femur_vs_Radius.csv")
res_Femur_Ulna     <- extraer_DEG_localización(fit2_loc, "Femur_vs_Ulna",
                                  "GSE239948_Femur_vs_Ulna.csv")
res_Humerus_Radius <- extraer_DEG_localización(fit2_loc, "Humerus_vs_Radius",
                                  "GSE239948_Humerus_vs_Radius.csv")
res_Humerus_Tibia  <- extraer_DEG_localización(fit2_loc, "Humerus_vs_Tibia",
                                  "GSE239948_Humerus_vs_Tibia.csv")
res_Humerus_Ulna   <- extraer_DEG_localización(fit2_loc, "Humerus_vs_Ulna",
                                  "GSE239948_Humerus_vs_Ulna.csv")
res_Radius_Tibia   <- extraer_DEG_localización(fit2_loc, "Radius_vs_Tibia",
                                  "GSE239948_Radius_vs_Tibia.csv")
res_Radius_Ulna    <- extraer_DEG_localización(fit2_loc, "Radius_vs_Ulna",
                                  "GSE239948_Radius_vs_Ulna.csv")
res_Tibia_Ulna     <- extraer_DEG_localización(fit2_loc, "Tibia_vs_Ulna",
                                  "GSE239948_Tibia_vs_Ulna.csv")

# Volcano plot — contraste representativo (Humero vs Tibia)
EnhancedVolcano(res_Humerus_Tibia,
                lab      = rownames(res_Humerus_Tibia),
                x        = "logFC",
                y        = "adj.P.Val",
                pCutoff  = 0.05,
                FCcutoff = 1,
                title    = "Humero vs Tibia",
                subtitle = "GSE239948 · limma · BH · |logFC|>1")

# Tabla resumen de todos los contrastes de ubicación
tabla_resumen_loc <- data.frame(
  Gene                   = rownames(res_Femur_Tibia),
  Femur_vs_Tibia_logFC   = res_Femur_Tibia$logFC,
  Femur_vs_Tibia_adjP    = res_Femur_Tibia$adj.P.Val,
  Femur_vs_Humerus_logFC = res_Femur_Humerus$logFC,
  Femur_vs_Humerus_adjP  = res_Femur_Humerus$adj.P.Val,
  Femur_vs_Radius_logFC  = res_Femur_Radius$logFC,
  Femur_vs_Radius_adjP   = res_Femur_Radius$adj.P.Val,
  Femur_vs_Ulna_logFC    = res_Femur_Ulna$logFC,
  Femur_vs_Ulna_adjP     = res_Femur_Ulna$adj.P.Val,
  Humerus_vs_Radius_logFC = res_Humerus_Radius$logFC,
  Humerus_vs_Radius_adjP  = res_Humerus_Radius$adj.P.Val,
  Humerus_vs_Tibia_logFC  = res_Humerus_Tibia$logFC,
  Humerus_vs_Tibia_adjP   = res_Humerus_Tibia$adj.P.Val,
  Humerus_vs_Ulna_logFC   = res_Humerus_Ulna$logFC,
  Humerus_vs_Ulna_adjP    = res_Humerus_Ulna$adj.P.Val,
  Radius_vs_Tibia_logFC   = res_Radius_Tibia$logFC,
  Radius_vs_Tibia_adjP    = res_Radius_Tibia$adj.P.Val,
  Radius_vs_Ulna_logFC    = res_Radius_Ulna$logFC,
  Radius_vs_Ulna_adjP     = res_Radius_Ulna$adj.P.Val,
  Tibia_vs_Ulna_logFC     = res_Tibia_Ulna$logFC,
  Tibia_vs_Ulna_adjP      = res_Tibia_Ulna$adj.P.Val
)

write.csv(tabla_resumen_loc,"GSE239948_tabla_resumen_ubicacion.csv",row.names = FALSE)

# ANÁLISIS B — TAMAÑO DE RAZA: Giant vs Large
design_size <- model.matrix(~ 0 + size_group, data = metadata_size)
colnames(design_size) <- levels(metadata_size$size_group)

fit_size <- lmFit(expr_size_filt, design_size)

contrast_size <- makeContrasts(
  Giant_vs_Large = Giant - Large,
  levels         = design_size)

fit2_size <- contrasts.fit(fit_size, contrast_size)
fit2_size <- eBayes(fit2_size)

res_Giant_Large <- topTable(
  fit2_size,
  coef          = "Giant_vs_Large",
  number        = Inf,
  adjust.method = "BH",
  sort.by       = "P"
)

deg_Giant_Large <- subset(res_Giant_Large,
                          adj.P.Val < 0.05 & abs(logFC) > 1)

cat("Total genes analizados:", nrow(res_Giant_Large), "\n")
cat("DEGs significativos (FDR<0.05, |logFC|>1):", nrow(deg_Giant_Large), "\n")
cat("  Sobreexpresados en Giant (logFC > 1):",sum(deg_Giant_Large$logFC > 1), "\n")
cat("  Sobreexpresados en Large (logFC < -1):",sum(deg_Giant_Large$logFC < -1), "\n")

head(res_Giant_Large, 10)

write.csv(data.frame(Gene = rownames(res_Giant_Large), res_Giant_Large),"GSE239948_Giant_vs_Large.csv",row.names = FALSE)

EnhancedVolcano(res_Giant_Large,
                lab      = rownames(res_Giant_Large),
                x        = "logFC",
                y        = "adj.P.Val",
                pCutoff  = 0.05,
                FCcutoff = 1,
                title    = "Giant vs Large",
                subtitle = "GSE239948 · limma · BH · |logFC|>1",
                xlab     = "log2 FC (positivo = ↑ en Giant)")
