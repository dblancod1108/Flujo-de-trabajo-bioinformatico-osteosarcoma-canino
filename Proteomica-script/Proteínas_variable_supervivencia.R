# =============================================================================
# Análisis de Expresión Diferencial Proteómica - Etiqueta supervivencia
# =============================================================================

# 1. LIBRERIAS
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(writexl)

# 2. CARGA Y EXPLORACIÓN INICIAL DE LOS DATOS
ruta_proteinas <- "C:/Users/danie/OneDrive/Documentos/Trabajo de Grado Maestría Bioinformática/Proteínasxmuestras_proteómica_280626.xlsx"
ruta_metadatos <- "C:/Users/danie/OneDrive/Documentos/Trabajo de Grado Maestría Bioinformática/Metadatos_proteómica.xlsx"

prot_raw <- read_excel(ruta_proteinas, sheet = 3, col_names = TRUE,
                       na = c("", "NA", "N/A", "#N/A"))
colnames(prot_raw)[1] <- "Accesion"

prot_mat <- prot_raw %>%
  column_to_rownames("Accesion") %>%
  as.matrix()
class(prot_mat) <- "numeric"

cat("Dimensiones de la matriz cruda:", nrow(prot_mat), "proteínas ×",
    ncol(prot_mat), "muestras\n")

meta_raw <- read_excel(ruta_metadatos, sheet = 1, col_names = TRUE,
                       na = c("", "NA", "-", "N/A"))
colnames(meta_raw)[1] <- "Variable"

meta <- as.data.frame(meta_raw) %>%
  column_to_rownames("Variable") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Muestra")

cat("Metadatos cargados:", nrow(meta), "muestras ×", ncol(meta), "variables\n")

# Alinear metadatos
meta <- meta[match(colnames(prot_mat), meta$Muestra), ]
rownames(meta) <- meta$Muestra
stopifnot("Desalineación muestra-metadatos" = all(colnames(prot_mat) == meta$Muestra))

##SUBSET -Muestra 1##
col_surv <- "Supervivencia después del diagnóstico (días)"
print(meta[, c("Muestra", col_surv)])

print(meta$Muestra[is.na(meta[[col_surv]])])

# Subset
muestras_validas <- intersect(
  colnames(prot_mat),
  meta$Muestra[!is.na(meta[[col_surv]])])

cat("\nMuestras incluidas en el análisis:", length(muestras_validas),
    "->", paste(muestras_validas, collapse = ", "), "\n")

# Aplicar subset a matriz y metadatos
prot_mat_surv <- prot_mat[, muestras_validas]
meta_surv     <- meta[match(muestras_validas, meta$Muestra), ]
rownames(meta_surv) <- meta_surv$Muestra

stopifnot("Desalineación tras subset" = all(colnames(prot_mat_surv) == meta_surv$Muestra))
print(table(meta_surv[[col_surv]]))

# 3. TRANSFORMACIÓN Y NORMALIZACIÓN
# Verificar que efectivamente no hay NA en la matriz cargada
n_na <- sum(is.na(prot_mat))
cat("\nNúmero de NA en la matriz cargada:", n_na, "\n")

# Transformación log2
prot_log2 <- log2(prot_mat + 0.0001)
prot_log2[is.infinite(prot_log2)] <- NA   # por si algún 0 exacto quedara

# Verificar que log2 no generó NA nuevos (no debería con la matriz ya filtrada)
cat("NA tras log2:", sum(is.na(prot_log2)), "(deben ser 0 si la matriz ya estaba limpia)\n")

# Normalización por mediana
medianas       <- apply(prot_log2, 2, median, na.rm = TRUE)
mediana_global <- median(prot_log2, na.rm = TRUE)
prot_norm      <- sweep(prot_log2, 2, medianas, "-") + mediana_global

cat("\nMedianas post-normalización (deben ser todas iguales):\n")
print(round(apply(prot_norm, 2, median, na.rm = TRUE), 3))

# 4. CONTROL DE CALIDAD VISUAL
# Boxplot 
prot_long <- as.data.frame(prot_norm) %>%
  rownames_to_column("Proteina") %>%
  pivot_longer(-Proteina, names_to = "Muestra", values_to = "log2emPAI") %>%
  left_join(meta_surv %>% select(Muestra, all_of(col_surv)), by = "Muestra")

# Renombrar columna para ggplot (evita problemas con caracteres especiales)
colnames(prot_long)[colnames(prot_long) == col_surv] <- "Supervivencia"

p_box <- ggplot(prot_long, aes(x = Muestra, y = log2emPAI, fill = Supervivencia)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  scale_fill_manual(values = c("Larga" = "#417505", "Corta" = "#D0021B")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold")) +
  labs(title    = "Distribución log2-emPAI por muestra",
       subtitle = paste0("Supervivencia | n = ", ncol(prot_norm),
                         " muestras | post-normalización"),
       x = NULL, y = "log2(emPAI)", fill = "Supervivencia")

print(p_box)

# PCA 
pca_res <- prcomp(t(prot_norm), center = TRUE, scale. = FALSE)
var_exp <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

pca_df <- as.data.frame(pca_res$x[, 1:2]) %>%
  rownames_to_column("Muestra") %>%
  left_join(meta_surv %>% select(Muestra, all_of(col_surv)), by = "Muestra")
colnames(pca_df)[colnames(pca_df) == col_surv] <- "Supervivencia"

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                            color = Supervivencia,
                            label = Muestra)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(size = 3.5, fontface = "bold") +
  scale_color_manual(values = c("Larga" = "#417505", "Corta" = "#D0021B")) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold")) +
  labs(title    = "PCA — Osteosarcoma canino",
       subtitle = "Coloreado por supervivencia",
       x = paste0("PC1 (", var_exp[1], "% varianza)"),
       y = paste0("PC2 (", var_exp[2], "% varianza)"),
       color = "Supervivencia")

if (dev.cur() != 1) dev.off()
print(p_pca)

# 5. ANÁLISIS DIFERENCIAL CON LIMMA
# Factor de grupo: Adulto = referencia, Geriatrico = caso
surv_factor <- factor(meta_surv[[col_surv]], levels = c("Larga", "Corta"))

print(table(surv_factor))
print(data.frame(Muestra = meta_surv$Muestra, Supervivencia = surv_factor))

# Matriz de diseño
design_surv <- model.matrix(~ surv_factor)
colnames(design_surv) <- c("Larga", "Corta_vs_Larga")
rownames(design_surv) <- meta_surv$Muestra

print(design_surv)

# Ajuste del modelo
# Nota: con matriz sin NA, lmFit corre de forma más estable y
# los grados de libertad de eBayes son exactos (no hay pesos implícitos)
fit_surv    <- lmFit(prot_norm, design_surv)
fit_surv_eb <- eBayes(fit_surv, trend = TRUE)# trend = TRUE: modela la relación media-varianza, recomendado para emPAI

Resultados_supervivencia_proteómica <- topTable(
  fit_surv_eb,
  coef          = "Corta_vs_Larga",
  number        = Inf,
  adjust.method = "BH",
  sort.by       = "P"
) %>%
  rownames_to_column("Accesion") %>%
  mutate(
    Regulacion = case_when(
      adj.P.Val < 0.05 & logFC >  1 ~ "Up en Corta",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down en Corta",
      TRUE                           ~ "No significativo"
    ))

print(Resultados_supervivencia_proteómica %>%
        filter(Regulacion != "No significativo") %>%
        select(Accesion, logFC, adj.P.Val, Regulacion))


Resultados_superv_proteómica_anotados <- Resultados_supervivencia_proteómica %>%
  left_join(
    prot_supervivencia_anotada %>%
      select(Accesion, Gene_Symbol),
    by = "Accesion")

###Volcano plot###
# Etiquetas: top 15 por FDR entre las significativas
top_labels <- Resultados_superv_proteómica_anotados %>%
  filter(
    Regulacion != "No significativo",
    !is.na(Gene_Symbol),
    Gene_Symbol != ""
  ) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 15)


# Si hay menos de 5 proteínas significativas anotadas,
# completar con las de menor FDR
if (nrow(top_labels) < 5) {
  
  top_labels <- Resultados_superv_proteómica_anotados %>%
    filter(
      !is.na(Gene_Symbol),
      Gene_Symbol != ""
    ) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 15)
}

# Volcano plot
p_volcano_surv <- ggplot(
  Resultados_superv_proteómica_anotados,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Regulacion
  )
) +
  
  # Puntos
  geom_point(
    alpha = 0.7,
    size = 2.5
  ) +
  
  # Etiquetas
  geom_text_repel(
    data = top_labels,
    aes(label = Gene_Symbol),
    size = 3,
    fontface = "bold",
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3,
    force = 2,
    segment.color = "grey40",
    segment.size = 0.3
  ) +
  
  # Umbral FDR
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.5
  ) +
  
  # Umbral de logFC
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.5
  ) +
  
  # Colores
  scale_color_manual(
    values = c(
      "Up en Corta"      = "#D0021B",
      "Down en Corta"    = "#417505",
      "No significativo" = "grey70"
    )
  ) +
  
  # Tema
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10,
      color = "grey40"
    ),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  # Títulos
  labs(
    title = "Volcano Plot — Supervivencia Corta vs Larga",
    subtitle = paste0(
      "Osteosarcoma canino | FDR (BH) | n = 7 muestras | ",
      sum(
        Resultados_superv_proteómica_anotados$Regulacion != 
          "No significativo"
      ),
      " proteínas significativas"
    ),
    x = "log2(Fold Change) [Corta / Larga]",
    y = "-log10(FDR ajustado)"
  )
 

print(p_volcano_surv)
write.csv(Resultados_supervivencia_proteómica,"Tabla_resultados_supervivencia.csv")
