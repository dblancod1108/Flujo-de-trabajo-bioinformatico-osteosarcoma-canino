# =============================================================================
# Análisis de Expresión Diferencial Proteómica - Etiqueta edad
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

prot_raw <- read_excel(ruta_proteinas, sheet = 4, col_names = TRUE,
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

# prot_norm es la matriz final que entra a limma

# 4. CONTROL DE CALIDAD VISUAL

# Boxplot de distribuciones por muestra

prot_long <- as.data.frame(prot_norm) %>%
  rownames_to_column("Proteina") %>%
  pivot_longer(-Proteina, names_to = "Muestra", values_to = "log2emPAI") %>%
  left_join(meta %>% select(Muestra,`Edad (años)`), by = "Muestra")

p_box <- ggplot(prot_long, aes(x = Muestra, y = log2emPAI, fill = `Edad (años)`)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  scale_fill_manual(values = c("Adulto" = "#F5A623", "Geriatrico" = "#50B83C")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold")) +
  labs(title    = "Distribución log2-emPAI por muestra",
       subtitle = "Coloreado por estadio de edad | post-normalización",
       x = NULL, y = "log2(emPAI)", fill = "Estadio de edad")

print(p_box)
ggsave("01_boxplot_muestras_edad.pdf", plot = p_box, width = 10, height = 6)


# PCA 
pca_res <- prcomp(t(prot_norm), center = TRUE, scale. = FALSE)
var_exp <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

pca_df <- as.data.frame(pca_res$x[, 1:2]) %>%
  rownames_to_column("Muestra") %>%
  left_join(meta %>% select(Muestra, `Edad (años)`), by = "Muestra")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                            color = `Edad (años)`,
                            label = Muestra)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(size = 3.5, fontface = "bold") +
  scale_color_manual(values = c("Adulto"    = "#F5A623",   # ← solo UNA escala
                                "Geriatrico" = "#50B83C")) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold")) +
  labs(
    title    = "PCA — Osteosarcoma canino",
    subtitle = "Coloreado por estadio de edad",
    x        = paste0("PC1 (", var_exp[1], "% varianza)"),
    y        = paste0("PC2 (", var_exp[2], "% varianza)"),
    color    = "Estadio de edad"
  )

# Resetear dispositivo gráfico antes de imprimir
if (dev.cur() != 1) dev.off()

print(p_pca)
ggsave("02_PCA_edad.pdf", plot = p_pca, width = 8, height = 7)

###ANÁLISIS DIFERENCIAL CON LIMMA###
# Factor de grupo: Adulto = referencia, Geriatrico = caso
edad_factor <- factor(meta$`Edad (años)`, levels = c("Adulto", "Geriatrico"))

print(table(edad_factor))
print(data.frame(Muestra = meta$Muestra, Edad = edad_factor))

# Matriz de diseño
design_edad <- model.matrix(~ edad_factor)
colnames(design_edad) <- c("Adulto", "Geriatrico_vs_Adulto")
rownames(design_edad) <- meta$Muestra

print(design_edad)

# Ajuste del modelo
# Nota: con matriz sin NA, lmFit corre de forma más estable y
# los grados de libertad de eBayes son exactos (no hay pesos implícitos)
fit_edad    <- lmFit(prot_norm, design_edad)
fit_edad2 <- eBayes(fit_edad, trend = TRUE)# trend = TRUE: modela la relación media-varianza, recomendado para emPAI

Resultados_edad_proteómica <- topTable(
  fit_edad2,
  coef          = "Geriatrico_vs_Adulto",
  number        = Inf,
  adjust.method = "BH",
  sort.by       = "P"
) %>%
  rownames_to_column("Accesion") %>%
  mutate(
    Regulacion = case_when(
      adj.P.Val < 0.05 & logFC >  1 ~ "Up en Geriatrico",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down en Geriatrico",
      TRUE                           ~ "No significativo"
    ))

Resultados_edad_proteómica_anotados <- Resultados_edad_proteómica %>%
  left_join(
    prot_edad_anotada %>%
      select(Accesion, Gene_Symbol),
    by = "Accesion")

# VOLCANO PLOT 
# Etiquetas: top 15 por FDR entre las significativas
top_labels <- Resultados_edad_proteómica_anotados %>%
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
  
  top_labels <- Resultados_edad_proteómica_anotados %>%
    filter(
      !is.na(Gene_Symbol),
      Gene_Symbol != ""
    ) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 15)
}

# Volcano plot
p_volcano <- ggplot(
  Resultados_edad_proteómica_anotados,
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
      "Up en Geriatrico"   = "#50B83C",
      "Down en Geriatrico" = "#F5A623",
      "No significativo"   = "grey70"
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
    title = "Volcano Plot — Geriátrico vs Adulto",
    subtitle = paste0(
      "Osteosarcoma canino | FDR (BH) | n = ",
      nrow(meta),
      " muestras | ",
      sum(
        Resultados_edad_proteómica_anotados$Regulacion !=
          "No significativo"
      ),
      " proteínas significativas"
    ),
    x = "log2(Fold Change) [Geriátrico / Adulto]",
    y = "-log10(FDR ajustado)"
  )

print(p_volcano)
ggsave("04_volcano_Edad.pdf", plot = p_volcano, width = 9, height = 8)
ggsave("04_volcano_Edad.png", plot = p_volcano, width = 9, height = 8, dpi = 300)

print(Resultados_edad_proteómica %>%
        filter(Regulacion != "No significativo") %>%
        select(Accesion, logFC, adj.P.Val, Regulacion))

write.csv(Resultados_edad_proteómica,"Tabla_resultados_edad.csv")
