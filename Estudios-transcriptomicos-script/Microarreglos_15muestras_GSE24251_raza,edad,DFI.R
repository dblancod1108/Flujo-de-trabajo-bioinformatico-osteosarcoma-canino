# =============================================================================
# Análisis de Expresión Diferencial — GSE24251
# Comparaciones:
#   a) Raza: Gigantes vs Medianos/grandes               
#   b) Edad: Geriátrico vs. Adulto
# =============================================================================

# 1. LIBRERÍAS
library(GEOquery)
library(openxlsx)
library(dplyr)
library(corrplot)
library(pheatmap)
library(limma)

# 2. CARGA Y EXPLORACIÓN INICIAL DE LOS DATOS
options(download.file.method = "libcurl")
gsetmicro15 <- getGEO("GSE24251", GSEMatrix = TRUE, getGPL = TRUE)
gsetmicro15 <- gsetmicro15[[1]]
class(gsetmicro15)

# 3. PREPARACIÓN DE LA MATRIZ DE EXPRESIÓN
exprmicro15  <- exprs(gsetmicro15)
genesmicro15 <- fData(gsetmicro15)
metadatamicro15 <- pData(gsetmicro15)

dim(exprmicro15)
summary(as.vector(as.matrix(exprmicro15)))

# Inspeccionar rango para determinar si los datos requieren transformación
range(exprmicro15, na.rm = TRUE)
max(exprmicro15, na.rm = TRUE)

# GSE24251 contiene intensidades que requieren transformación log2.
# Se utiliza un offset de 1 para evitar log2(0).
exprmicro15 <- log2(exprmicro15 + 1)
# Verificación 
summary(as.vector(as.matrix(exprmicro15)))
range(exprmicro15, na.rm = TRUE)

# 4. PREPROCESAMIENTO Y CONTROL DE CALIDAD
boxplot(
  exprmicro15,
  las       = 2,
  cex.axis  = 0.6,
  outline   = FALSE,
  main      = "Distribución de expresión por muestra — GSE24251",
  ylab      = "Intensidad de expresión",
  xlab      = "Muestras")

abline(h   = median(exprmicro15, na.rm = TRUE),lty = 2,col = "red")

plot(
  density(exprmicro15[, 1], na.rm = TRUE),
  main = "Distribución de densidad de expresión — GSE24251",
  xlab = "Intensidad de expresión",
  ylab = "Densidad",
  col  = "gray40",
  lwd  = 2)

for (i in 2:ncol(exprmicro15)) {
  lines(
    density(exprmicro15[, i], na.rm = TRUE),
    lwd = 1.5)}

legend(
  "topright",
  legend = colnames(exprmicro15),
  lwd    = 2,
  cex    = 0.6,
  bty    = "n")

cor_micro15 <- cor(
  exprmicro15,
  use = "pairwise.complete.obs",
  method = "pearson")

corrplot(
  cor_micro15,
  method = "color",
  order = "hclust",
  tl.cex = 0.7,
  tl.col = "black",
  main = "Correlación entre muestras — GSE24251")

pca15 <- prcomp(
  t(exprmicro15),
  scale. = TRUE)

var_exp15 <- round(
  100 * pca15$sdev^2 /
    sum(pca15$sdev^2),
  1)

# Crear tabla para graficar
pca_df15 <- data.frame(
  Muestra = rownames(pca15$x),
  PC1 = pca15$x[, 1],
  PC2 = pca15$x[, 2]
)

plot(
  pca_df15$PC1,
  pca_df15$PC2,
  pch = 19,
  xlab = paste0("PC1 (", var_exp15[1], "%)"),
  ylab = paste0("PC2 (", var_exp15[2], "%)"),
  main = "PCA de muestras — GSE24251"
)

text(
  pca_df15$PC1,
  pca_df15$PC2,
  labels = pca_df15$Muestra,
  pos = 3,
  cex = 0.7
)

pheatmap(
  cor_micro15,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Heatmap de correlación entre muestras — GSE24251",
  fontsize = 8)

# 5. ANOTACIÓN Y FILTRADO DE SONDAS
tabla_expmicro15 <- as.data.frame(exprmicro15)
tabla_expmicro15$GeneSymbol <- genesmicro15$`Gene Symbol`
tabla_expmicro15$EntrezID <- genesmicro15$ENTREZ_GENE_ID

dim(tabla_expmicro15)

# Eliminar sondas sin gen
tabla_expmicrof15 <-
  tabla_expmicro15 %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != "",
    GeneSymbol != "0",
    !grepl(
      "^[0-9]+$",
      GeneSymbol))

cat("Sondas con anotación génica:",nrow(tabla_expmicrof15),"\n")

# Limpiar símbolos
tabla_expmicrof15$GeneSymbol <- trimws(tabla_expmicrof15$GeneSymbol)

tabla_expmicrof15$GeneSymbol <-
  sapply(
    strsplit(
      tabla_expmicrof15$GeneSymbol,
      " /// "),`[`,1)

# Verificación
sum(grepl(" ",tabla_expmicrof15$GeneSymbol))

# Eliminar símbolos puramente numéricos
tabla_expmicrof15 <- tabla_expmicrof15 %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != "",
    GeneSymbol != "0",
    !grepl("^[0-9]+$", GeneSymbol))

# Verificación
sum(grepl(" ",tabla_expmicrof15$GeneSymbol))
dim(tabla_expmicrof15)

#Identificar valores faltantes
na_total <- sum(is.na(exprmicro15))
cat("Número total de valores NA:",na_total,"\n")

#Filtrado de baja expresión
sample_cols15 <- colnames(exprmicro15)

# Para datos log2, se puede utilizar un umbral de intensidad.
# Se conserva una sonda si presenta intensidad > 5
# en al menos el 50% de las muestras.

keep_expr15 <- rowSums(
  tabla_expmicrof15[, sample_cols15] > 5,
  na.rm = TRUE
) >= ceiling(length(sample_cols15) * 0.50)


tabla_expmicrof15 <- tabla_expmicrof15[keep_expr15,]

cat("Sondas después del filtrado de baja expresión:",nrow(tabla_expmicrof15),"\n")

# 6. SUMARIZACIÓN (SONDAS A GENES)
tabla_expmicrof15$varianza <- apply(
  tabla_expmicrof15[, sample_cols15],
  1,var,na.rm = TRUE)

tabla_filtradamicro15 <- tabla_expmicrof15 %>%
  group_by(GeneSymbol) %>%
  slice_max(
    order_by = varianza,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(-varianza)


# Verificación
cat("Genes únicos después de sumarización:",nrow(tabla_filtradamicro15),"\n")
any(duplicated(tabla_filtradamicro15$GeneSymbol))

# 7. METADATOS Y CLASIFICACIÓN DE LAS VARIABLES
#Extracción de raza
raza_col <- grep(
  "breed|genotype|race",
  colnames(metadatamicro15),
  ignore.case = TRUE,
  value = TRUE)

print(raza_col)

table(metadatamicro15[[raza_col[1]]])

metadatos_limpiomicro15 <- metadatamicro15 %>%
  transmute(
    Titulo = as.character(title),
    
    GEO_ID = as.character(geo_accession),
    
    Raza = trimws(
      as.character(
        .data[[raza_col[1]]]
      )
    ),
    
    Edad_años = as.numeric(
      metadatamicro15[
        ,
        "age (years):ch1"
      ]
    ),
    
    Grupo_DFI = case_when(
      
      grepl(
        "DFI<100",
        title,
        ignore.case = TRUE
      ) ~ "DFI_corto",
      
      grepl(
        "DFI>300",
        title,
        ignore.case = TRUE
      ) ~ "DFI_largo",
      
      TRUE ~ NA_character_))

# Clasificar por tamaño de raza
metadatos_limpiomicro15$Grupo_Raza <- case_when(
  metadatos_limpiomicro15$Raza %in%
    c(
      "Greyhound",
      "Labrador Ret.",
      "Doberman",
      "Golden Ret."
    ) ~ "Mediana_Grande",
  
  metadatos_limpiomicro15$Raza %in%
    c(
      "Rottweiler",
      "Mix"
    ) ~ "Gigante",
  
  TRUE ~ NA_character_)

table(metadatos_limpiomicro15$Grupo_Raza, useNA = "always")

# Clasificar por edad
metadatos_limpiomicro15$Grupo_Tamanio <-
  case_when(
    
    metadatos_limpiomicro15$Grupo_Raza ==
      "Gigante" ~ "Giant",
    
    metadatos_limpiomicro15$Grupo_Raza ==
      "Mediana_Grande" ~ "Large",
    
    TRUE ~ NA_character_
  )

metadatos_limpiomicro15$Etapa_Vida <- case_when(
  # Gigantes
  metadatos_limpiomicro15$Grupo_Tamanio ==
    "Giant" &
    metadatos_limpiomicro15$Edad_años >= 8 ~
    "Geriatrico",
  
  metadatos_limpiomicro15$Grupo_Tamanio ==
    "Giant" &
    metadatos_limpiomicro15$Edad_años < 8 ~
    "Adulto_Senior",
  
  
  # Medianas-grandes
  metadatos_limpiomicro15$Grupo_Tamanio ==
    "Large" &
    metadatos_limpiomicro15$Edad_años >= 10 ~
    "Geriatrico",
  
  metadatos_limpiomicro15$Grupo_Tamanio ==
    "Large" &
    metadatos_limpiomicro15$Edad_años < 10 ~
    "Adulto_Senior",
  
  TRUE ~ NA_character_)

table(metadatos_limpiomicro15$Etapa_Vida, useNA = "always")

# Resumen 
metadatos_limpiomicro15 %>%
  filter(
    !is.na(Etapa_Vida)
  ) %>%
  group_by(
    Grupo_Raza,
    Etapa_Vida
  ) %>%
  summarise(
    n = n(),
    edad_media = mean(
      Edad_años,
      na.rm = TRUE
    ),
    edad_min = min(
      Edad_años,
      na.rm = TRUE
    ),
    edad_max = max(
      Edad_años,
      na.rm = TRUE
    ),
    .groups = "drop")

# Alinear metadatos con muestras
rownames(
  metadatos_limpiomicro15
) <- metadatos_limpiomicro15$GEO_ID


stopifnot(all(colnames(exprmicro15) == rownames(metadatos_limpiomicro15)))

# Renombrar muestras
sum_expmicro15 <- colnames(
  exprmicro15)

nuevos_nombres15 <- paste0(
  "Muestra_",
  seq_along(sum_expmicro15))

# Renombrar expresión
colnames(
  tabla_filtradamicro15
)[
  colnames(tabla_filtradamicro15)
  %in% sum_expmicro15
] <- nuevos_nombres15

# Renombrar metadatos
rownames(
  metadatos_limpiomicro15
) <- nuevos_nombres15

# Exportar expresión y metadatos
tabla_finalmicro15 <- tabla_filtradamicro15 %>%
  select(
    GeneSymbol,
    EntrezID,
    starts_with("Muestra_"))

write.xlsx(tabla_finalmicro15,file = "GSE24251_genesxmuestras_final.xlsx",rowNames = FALSE)

metadatos_finalmicro15 <-
  as.data.frame(
    t(metadatos_limpiomicro15))
colnames(metadatos_finalmicro15) <- nuevos_nombres15

write.xlsx(metadatos_finalmicro15,file = "GSE24251_metadata_organizado.xlsx",rowNames = TRUE)

# 8. CONTROL DE CALIDAD ENTRE MUESTRAS
#Colores
col_raza15 <- case_when(
  
  metadatos_limpiomicro15$Grupo_Raza ==
    "Gigante" ~ "#E67E22",
  
  metadatos_limpiomicro15$Grupo_Raza ==
    "Mediana_Grande" ~ "#8E44AD",
  
  TRUE ~ "grey60")

boxplot(
  expr_genes15,
  las = 2,
  cex.axis = 0.6,
  outline = FALSE,
  col = col_raza15,
  main = "Distribución de expresión por muestra — GSE24251",
  ylab = "Intensidad (log2)")

legend(
  "topright",
  legend = c(
    "Gigante",
    "Mediana_Grande"
  ),
  fill = c(
    "#E67E22",
    "#8E44AD"
  ),
  bty = "n"
)

plot(
  density(
    expr_genes15[, 1],
    na.rm = TRUE
  ),
  main = "Distribución de expresión — GSE24251",
  xlab = "Intensidad (log2)",
  col = col_raza15[1]
)


for (i in 2:ncol(expr_genes15)) {
  
  lines(
    density(
      expr_genes15[, i],
      na.rm = TRUE
    ),
    col = col_raza15[i]
  )
}


legend(
  "topright",
  legend = c(
    "Gigante",
    "Mediana_Grande"
  ),
  col = c(
    "#E67E22",
    "#8E44AD"
  ),
  lwd = 2,
  bty = "n"
)

cor_micro15 <- cor(
  expr_genes15,
  use = "pairwise.complete.obs")


corrplot(
  cor_micro15,
  method = "color",
  order = "hclust",
  tl.cex = 0.7,
  main = "Correlación entre muestras — GSE24251")


pca15 <- prcomp(
  t(expr_genes15),
  scale. = TRUE)


var_exp15 <- round(
  100 *
    pca15$sdev^2 /
    sum(pca15$sdev^2),
  1
)


pca_df15 <- data.frame(
  
  PC1 = pca15$x[, 1],
  
  PC2 = pca15$x[, 2],
  
  Grupo_Raza =
    metadatos_limpiomicro15$Grupo_Raza,
  
  Raza =
    metadatos_limpiomicro15$Raza,
  
  Grupo_DFI =
    metadatos_limpiomicro15$Grupo_DFI,
  
  Etapa_Vida =
    metadatos_limpiomicro15$Etapa_Vida)

ggplot(
  pca_df15,
  aes(
    PC1,
    PC2,
    color = Grupo_Raza,
    label = Raza
  )
) +
  geom_point(size = 4) +
  geom_text(
    vjust = -0.7,
    size = 2.8
  ) +
  scale_color_manual(
    values = c(
      Gigante = "#E67E22",
      Mediana_Grande = "#8E44AD"
    ),
    na.value = "grey60"
  ) +
  labs(
    title = "PCA — GSE24251",
    x = paste0(
      "PC1 (",
      var_exp15[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      var_exp15[2],
      "%)"
    )
  ) +
  theme_bw()


annotation_col15 <- data.frame(
  
  Grupo_Raza =
    metadatos_limpiomicro15$Grupo_Raza,
  
  Grupo_DFI =
    metadatos_limpiomicro15$Grupo_DFI,
  
  Etapa_Vida =
    metadatos_limpiomicro15$Etapa_Vida,
  
  row.names =
    colnames(expr_genes15)
)


pheatmap(
  cor_micro15,
  annotation_col = annotation_col15,
  clustering_distance_rows =
    "euclidean",
  clustering_distance_cols =
    "euclidean",
  main =
    "Heatmap de correlación — GSE24251",
  fontsize = 8)

# 9. EXPRESIÓN DIFERENCIAL — LIMMA
# Preparar matriz de expresión por genes
expr_genes15 <- tabla_finalmicro15 %>%
  select(
    -GeneSymbol,
    -EntrezID
  ) %>%
  as.matrix()

rownames(expr_genes15) <- tabla_finalmicro15$GeneSymbol

# Verificación
stopifnot(all(colnames(expr_genes15) == rownames(metadatos_limpiomicro15)))

##ANÁLISIS A — RAZA
grupo_raza15 <- factor(
  metadatos_limpiomicro15$Grupo_Raza,
  levels = c(
    "Mediana_Grande",
    "Gigante"))

table(grupo_raza15, useNA = "always")

valid_raza15 <-!is.na(grupo_raza15)
expr_raza15 <- expr_genes15[, valid_raza15]

# Factor
grupo_raza15_valid <-
  droplevels(
    grupo_raza15[valid_raza15])

meta_raza15 <-
  metadatos_limpiomicro15[
    valid_raza15,]

stopifnot(all(colnames(expr_raza15) == rownames(meta_raza15)))

#Diseño
# Mediana_Grande = referencia
# logFC positivo = mayor expresión/ganancia en Gigante

design_raza15 <- model.matrix(
  ~ grupo_raza15_valid)
colnames(design_raza15)
design_raza15

#Ajustar modelo
fit_raza15 <- lmFit(
  expr_raza15,
  design_raza15)

fit_raza15 <- eBayes(fit_raza15)

#Resultados
resultados_raza15 <- topTable(
  fit_raza15,
  coef = "grupo_raza15_validGigante",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)


DEG_raza15 <- resultados_raza15 %>%
  filter(adj.P.Val < 0.05,abs(logFC) > 1)

# Volcano
EnhancedVolcano(
  resultados_raza15,
  lab = rownames(resultados_raza15),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "GSE24251 — Gigante vs. Mediana_Grande",
  subtitle = "limma · BH · |logFC| > 1",
  xlab =
    "log2 FC (positivo = mayor ganancia en Gigante)")

##ANÁLISIS B — EDAD
grupo_edad15 <- factor(
  metadatos_limpiomicro15$Etapa_Vida,
  levels = c(
    "Adulto_Senior",
    "Geriatrico"))

table(grupo_edad15,useNA = "always")

#Filtrar muestras
valid_edad15 <-!is.na(grupo_edad15)
expr_edad15 <- expr_genes15[, valid_edad15]

#Factos
grupo_edad15_valid <-
  droplevels(grupo_edad15[valid_edad15])

meta_edad15 <- metadatos_limpiomicro15[valid_edad15,]

stopifnot(all(colnames(expr_edad15) == rownames(meta_edad15)))

#Diseño
# Adulto_Senior = referencia
# logFC positivo = mayor expresión/ganancia en Geriátrico

design_edad15 <- model.matrix(~ grupo_edad15_valid)
colnames(design_edad15)

#Ajustar modelo
fit_edad15 <- lmFit(
  expr_edad15,design_edad15)

fit_edad15 <- eBayes(fit_edad15)

#Resultados
resultados_edad15 <- topTable(
  fit_edad15,
  coef = "grupo_edad15_validGeriatrico",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)


DEG_edad15 <- resultados_edad15 %>%
  filter(adj.P.Val < 0.05,abs(logFC) > 1)

# Volcano
EnhancedVolcano(
  resultados_edad15,
  lab = rownames(resultados_edad15),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "GSE24251 — Geriátrico vs. Adulto_Senior",
  subtitle = "limma · BH · |logFC| > 1",
  xlab =
    "log2 FC (positivo = mayor ganancia en Geriátrico)")

# ANÁLISIS SECUNDARIO — DFI CORTO VS. DFI LARGO
grupo_dfi15 <- factor(
  metadatos_limpiomicro15$Grupo_DFI,
  levels = c(
    "DFI_largo",
    "DFI_corto"))

table(grupo_dfi15,useNA = "always")

#Filtrar muestras
valid_dfi15 <-!is.na(grupo_dfi15)
expr_dfi15 <- expr_genes15[, valid_dfi15]

#Factor
grupo_dfi15_valid <-
  droplevels(
    grupo_dfi15[valid_dfi15])


meta_dfi15 <-
  metadatos_limpiomicro15[
    valid_dfi15,]

stopifnot(all(colnames(expr_dfi15) == rownames(meta_dfi15)))

#Diseño
# DFI_largo = referencia
# logFC positivo = mayor expresión/ganancia en DFI_corto

design_dfi15 <- model.matrix(
  ~ grupo_dfi15_valid)

#Ajustar modelo
fit_dfi15 <- lmFit(
  expr_dfi15,
  design_dfi15)

fit_dfi15 <- eBayes(fit_dfi15)

#Resultados
resultados_dfi15 <- topTable(
  fit_dfi15,
  coef = "grupo_dfi15_validDFI_corto",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P")

DEG_dfi15 <- resultados_dfi15 %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) > 1)

#Volcano
EnhancedVolcano(
  resultados_dfi15,
  lab = rownames(resultados_dfi15),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "GSE24251 — DFI corto vs. DFI largo",
  subtitle = "Análisis secundario · limma · BH",
  xlab =
    "log2 FC (positivo = mayor ganancia en DFI corto)")
