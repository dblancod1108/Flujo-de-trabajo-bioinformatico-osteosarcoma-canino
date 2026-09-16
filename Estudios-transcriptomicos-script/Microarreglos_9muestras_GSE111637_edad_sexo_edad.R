# =============================================================================
# Análisis de Expresión Diferencial — GSE111637
# Comparaciones:
#   a) Raza: Gigantes vs Medianos/grandes               
#   b) Sexo: Hembra vs. Macho                                          
#   c) Edad: Geriátrico vs. Adulto
# =============================================================================

# 1. LIBRERÍAS
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(openxlsx)
library(corrplot)

# 2. CARGA Y EXPLORACIÓN INICIAL DE LOS DATOS 
options(download.file.method = "libcurl")
gset111 <- getGEO("GSE111637", GSEMatrix = TRUE, getGPL = TRUE)
gset111 <- gset111[[1]]
class(gset111)

# 3. PREPARACIÓN DE LA MATRIZ DE EXPRESIÓN
expr111  <- exprs(gset111)
genes111 <- fData(gset111)    # anotación de sondas (plataforma Agilent)
meta111  <- pData(gset111)    # metadatos de muestras

dim(expr111)                  # ~180880 sondas × 9 muestras
colnames(expr111)
head(rownames(expr111))
summary(as.vector(expr111))
range(expr111, na.rm = TRUE)

# 4. PREPROCESAMIENTO Y CONTROL DE CALIDAD
boxplot(expr111,
        las      = 2,
        cex.axis = 0.7,
        outline  = FALSE,
        col      = "steelblue",
        main     = "Distribución log ratio por muestra — GSE111637",
        ylab     = "log2(tumor/normal)",
        xlab     = "")
abline(h = 0, lty = 2, col = "red")   # línea en 0 = sin cambio de copias

plot(density(expr111[, 1], na.rm = TRUE),
     main = "Distribución de densidad de expresión — GSE111637",
     xlab = "log2(tumor/normal)",
     ylab = "Densidad",
     lwd  = 2,
     col  = "gray40")

# Añadir las demás muestras
for(i in 2:ncol(expr111)){
  lines(density(expr111[, i], na.rm = TRUE),
        lwd = 1.5,
        col = "gray70")}

abline(v = 0,lty = 2,col = "red")

legend("topright",
       legend = c("Muestras", "Referencia (log2 ratio = 0)"),
       col    = c("gray70", "red"),
       lty    = c(1, 2),
       lwd    = c(1.5, 1),
       bty    = "n")

cor_111 <- cor(expr111,
               use = "pairwise.complete.obs",
               method = "pearson")

corrplot(cor_111,
         method = "color",
         order = "hclust",
         type = "upper",
         tl.cex = 0.8,
         tl.col = "black",
         addCoef.col = "black",
         number.cex = 0.5,
         main = "Correlación entre muestras — GSE111637")

pheatmap(cor_111,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Heatmap de correlación entre muestras — GSE111637",
         fontsize = 9,
         border_color = NA)

# 5. ANOTACIÓN Y FILTRADO DE SONDAS
print(colnames(genes111))

gene_col <- grep("gene.symbol|GENE_SYMBOL|SystematicName|GeneName",
  colnames(genes111),
  ignore.case = TRUE,
  value = TRUE)

gene_col

tabla_expr111 <- as.data.frame(expr111)
tabla_expr111$GeneSymbol <- genes111[[gene_col[1]]]

sample_cols <- colnames(expr111)

tabla_expr111f <- tabla_expr111 %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != "",
    GeneSymbol != "0",
    !grepl("^[0-9]+$", GeneSymbol))

cat("Sondas con anotación génica:",nrow(tabla_expr111f),"\n")

# Filtrar sondas con demasiados valores faltantes
na_prop <- rowMeans(is.na(tabla_expr111f[, sample_cols]))

tabla_expr111f <- tabla_expr111f[na_prop < 0.5,]

cat("Sondas tras filtrado de NAs:",nrow(tabla_expr111f),"\n")

#Limpieza de símbolos 
tabla_expr111f$GeneSymbol <- trimws(tabla_expr111f$GeneSymbol)

tabla_expr111f$GeneSymbol <- sapply(
  strsplit(
    tabla_expr111f$GeneSymbol,
    " /// |;|,"
  ),`[`,1)

tabla_expr111f$GeneSymbol <- trimws(tabla_expr111f$GeneSymbol)

# 6. SUMARIZACIÓN (SONDAS A GENES)
tabla_expr111f$varianza <- apply(
  tabla_expr111f[, sample_cols],
  1,var,na.rm = TRUE)

tabla_filtrada111 <- tabla_expr111f %>%
  group_by(GeneSymbol) %>%
  slice_max(
    order_by = varianza,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(-varianza)

cat("Genes únicos tras sumarización:",nrow(tabla_filtrada111),"\n")

# 7. METADATOS Y CLASIFICACIÓN DE RAZAS
colnames(meta111)

# Extraer información de raza y características
metadatos_limpios <- meta111 %>%
  transmute(
    GEO_ID   = geo_accession,
    Titulo   = title,
    Tejido   = as.character(`tissue:ch1`),
    Raza     = as.character(`breed:ch1`),
    Edad     = as.character(`age:ch1`),
    Sexo     = as.character(`gender:ch1`))

# Ver información original
table(metadatos_limpios$Raza, useNA = "always")
table(metadatos_limpios$Sexo, useNA = "always")
table(metadatos_limpios$Edad, useNA = "always")

# Clasificar según los grupos definidos
metadatos_limpios$Grupo_Raza <- case_when(
  metadatos_limpios$Raza %in% c(
    "Greyhound",
    "Mixed Breed (Labrador Retriever/X)",
    "Boxer"
  ) ~ "Mediana_Grande",
  
  metadatos_limpios$Raza %in% c(
    "Rottweiler",
    "Rottweiler X (Canine Mixed Breed)",
    "Mixed Breed (German Shepherd/Collie/X)",
    "Bull Mastiff",
    "Great Pyrenees"
  ) ~ "Gigante",
  
  TRUE ~ NA_character_
)

table(metadatos_limpios$Grupo_Raza, useNA = "always")

# Alinear rownames con colnames de la matriz
rownames(metadatos_limpios) <- metadatos_limpios$GEO_ID
stopifnot(all(colnames(expr111) == rownames(metadatos_limpios)))

# 8. CONTROL DE CALIDAD ENTRE MUESTRAS
expr_genes111 <- as.matrix(tabla_filtrada111[, sample_cols])
rownames(expr_genes111) <- tabla_filtrada111$GeneSymbol

# Colores por grupo de raza
col_grupos111 <- ifelse(
  metadatos_limpios$Grupo_Raza == "Gigante",
  "#E74C3C",       # rojo = Gigante
  ifelse(
    metadatos_limpios$Grupo_Raza == "Mediana_Grande",
    "#3498DB",     # azul = Mediana_Grande
    "grey60"       # gris = sin clasificar
  )
)

boxplot(expr_genes111,
        las      = 2,
        cex.axis = 0.7,
        outline  = FALSE,
        col      = col_grupos111,
        main     = "Distribución log ratio (genes) por muestra",
        ylab     = "log2(tumor/normal)")
abline(h = 0, lty = 2, col = "black")
legend("topright",
       legend = c("Gigante", "Mediana_Grande"),
       fill   = c("#E74C3C", "#3498DB"),
       bty    = "n", cex = 0.8)

cor_111 <- cor(expr_genes111, use = "pairwise.complete.obs")
corrplot(cor_111,
         method  = "color",
         order   = "hclust",
         tl.cex  = 0.8,
         main    = "Correlación entre muestras — GSE111637")

pca_111 <- prcomp(
  t(expr_genes111),
  scale. = TRUE)

var_exp <- round(
  100 * pca_111$sdev^2 /
    sum(pca_111$sdev^2),1)

plot(
  pca_111$x[, 1],
  pca_111$x[, 2],
  xlab = paste0("PC1 (", var_exp[1], "%)"),
  ylab = paste0("PC2 (", var_exp[2], "%)"),
  main = "PCA — GSE111637",
  pch = 19)

text(
  pca_111$x[, 1],
  pca_111$x[, 2],
  labels = colnames(expr_genes111),
  pos = 3,
  cex = 0.7)

pheatmap(cor_111,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main     = "Heatmap correlación muestras — GSE111637",
         fontsize = 9,
         annotation_col = data.frame(
           Grupo = metadatos_limpios$Grupo_Raza,
           row.names = colnames(expr_genes111)
         ))

# 9. EXPRESIÓN DIFERENCIAL - LIMMA
#Preparar matriz de expresión
expr_matrix111 <- expr_genes111

stopifnot(all(colnames(expr_matrix111) == rownames(metadatos_limpios)))

##ANÁLISIS A — RAZA
table(metadatos_limpios$Grupo_Raza,useNA = "always")
# Filtrar solo muestras con grupo de raza asignado
valid_raza <- !is.na(metadatos_limpios$Grupo_Raza)
expr_raza111 <- expr_matrix111[, valid_raza]
meta_raza111 <- metadatos_limpios[valid_raza,]

stopifnot(all(colnames(expr_raza111) == rownames(meta_raza111)))

# Factor
grupo_raza111 <- factor(
  meta_raza111$Grupo_Raza,
  levels = c(
    "Gigante",
    "Mediana_Grande"))

table(grupo_raza111)

#Diseño
design_raza111 <- model.matrix(~0 + grupo_raza111)
colnames(design_raza111) <- c(
  "Gigante",
  "Mediana_Grande")
rownames(design_raza111) <-rownames(meta_raza111)

design_raza111

#Contraste
contraste_raza111 <- makeContrasts(
  Gigante_vs_MedianaGrande =
    Gigante - Mediana_Grande,
  levels = design_raza111)

#Ajuste Modelo
fit_raza111 <- lmFit(expr_raza111,design_raza111)

fit2_raza111 <- contrasts.fit(fit_raza111,contraste_raza111)
fit2_raza111 <- eBayes(fit2_raza111)

#Resultados
resultados_raza111 <- topTable(
  fit2_raza111,
  coef = "Gigante_vs_MedianaGrande",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P")

head(resultados_raza111,20)
summary(resultados_raza111$adj.P.Val)

deg_raza111 <- resultados_raza111 %>%
  filter(adj.P.Val < 0.05,abs(logFC) > 0.5)

#Volcano
EnhancedVolcano(
  resultados_raza111,
  lab = rownames(resultados_raza111),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 0.5,
  title = "Gigante vs Mediana_Grande",
  subtitle = "GSE111637 · limma · BH · |logFC| > 0.5",
  xlab = "log2 FC (Gigante / Mediana_Grande)")

##ANÁLISIS B — SEXO
table(metadatos_limpios$Sexo,useNA = "always")

metadatos_limpios$Grupo_Sexo <- case_when(
  metadatos_limpios$Sexo == "FS" ~ "Hembra",
  metadatos_limpios$Sexo == "MC" ~ "Macho",
  TRUE ~ NA_character_)

table(metadatos_limpios$Grupo_Sexo,useNA = "always")

#Filtrar muestras
valid_sexo111 <- !is.na(metadatos_limpios$Grupo_Sexo)
expr_sexo111 <- expr_matrix111[,valid_sexo111]
meta_sexo111 <- metadatos_limpios[valid_sexo111,]

stopifnot(all(colnames(expr_sexo111) ==rownames(meta_sexo111)))

#Factor
grupo_sexo111 <- factor(
  meta_sexo111$Grupo_Sexo,
  levels = c(
    "Hembra",
    "Macho"))

table(grupo_sexo111)

#Diseño
design_sexo111 <- model.matrix(~0 + grupo_sexo111)
colnames(design_sexo111) <- c("Hembra","Macho")
rownames(design_sexo111) <-rownames(meta_sexo111)

design_sexo111

#Contraste
contraste_sexo111 <- makeContrasts(
  Hembra_vs_Macho =
    Hembra - Macho,
  levels = design_sexo111)

#Ajuste modelo
fit_sexo111 <- lmFit(expr_sexo111,design_sexo111)

fit2_sexo111 <- contrasts.fit(fit_sexo111,contraste_sexo111)
fit2_sexo111 <- eBayes(fit2_sexo111)

#Resultados
resultados_sexo111 <- topTable(
  fit2_sexo111,
  coef = "Hembra_vs_Macho",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P")

deg_sexo111 <- resultados_sexo111 %>%
  filter(adj.P.Val < 0.05,abs(logFC) > 0.5)

#Volcano
EnhancedVolcano(
  resultados_sexo111,
  lab = rownames(resultados_sexo111),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 0.5,
  title = "Hembra vs Macho",
  subtitle = "GSE111637 · limma · BH · |logFC| > 0.5",
  xlab = "log2 FC (Hembra / Macho)")

##ANÁLISIS C - EDAD
metadatos_limpios$Edad_num <- as.numeric(
  gsub("[^0-9.]","",
    metadatos_limpios$Edad))

data.frame(Original = metadatos_limpios$Edad,Numerico = metadatos_limpios$Edad_num)

#Clasificación dependiente de raza
metadatos_limpios$Grupo_Edad <- case_when(
  
  metadatos_limpios$Grupo_Raza == "Mediana_Grande" &
    metadatos_limpios$Edad_num >= 10 ~
    "Geriatrico",
  
  metadatos_limpios$Grupo_Raza == "Mediana_Grande" &
    metadatos_limpios$Edad_num < 10 ~
    "Adulto_senior",
  
  metadatos_limpios$Grupo_Raza == "Gigante" &
    metadatos_limpios$Edad_num >= 8 ~
    "Geriatrico",
  
  metadatos_limpios$Grupo_Raza == "Gigante" &
    metadatos_limpios$Edad_num < 8 ~
    "Adulto_senior",
  
  TRUE ~ NA_character_)

table(metadatos_limpios$Grupo_Edad,useNA = "always")

metadatos_limpios %>%filter(!is.na(Grupo_Edad)) %>%
  group_by(Grupo_Raza,Grupo_Edad) %>%
  summarise(
    n = n(),
    edad_media = mean(
      Edad_num,
      na.rm = TRUE),
    edad_min = min(
      Edad_num,
      na.rm = TRUE),
    edad_max = max(
      Edad_num,
      na.rm = TRUE),
    .groups = "drop")

#Filtrar muestras para edad
valid_edad111 <- !is.na(metadatos_limpios$Grupo_Edad)
expr_edad111 <- expr_matrix111[,valid_edad111]
meta_edad111 <- metadatos_limpios[valid_edad111,]

stopifnot(all(colnames(expr_edad111) ==rownames(meta_edad111)))

#Diseño y grupo
grupo_edad111 <- factor(
  meta_edad111$Grupo_Edad,
  levels = c("Geriatrico","Adulto_senior"))

table(grupo_edad111)

design_edad111 <- model.matrix(~0 + grupo_edad111)
colnames(design_edad111) <- c("Geriatrico","Adulto_senior")
rownames(design_edad111) <-rownames(meta_edad111)

design_edad111

#Contraste
contraste_edad111 <- makeContrasts(
  Geriatrico_vs_Adulto =
    Geriatrico - Adulto_senior,
  levels = design_edad111)

#Ajuste de modelo
fit_edad111 <- lmFit(expr_edad111,design_edad111)

fit2_edad111 <- contrasts.fit(fit_edad111,contraste_edad111)
fit2_edad111 <- eBayes(fit2_edad111)

#Resultados
resultados_edad111 <- topTable(
  fit2_edad111,
  coef = "Geriatrico_vs_Adulto",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P")

head(resultados_edad111,20)
summary(resultados_edad111$adj.P.Val)

deg_edad111 <- resultados_edad111 %>%filter(adj.P.Val < 0.05,abs(logFC) > 0.5)

#Volcano
EnhancedVolcano(
  resultados_edad111,
  lab = rownames(resultados_edad111),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 0.5,
  title = "Geriátrico vs Adulto_senior",
  subtitle = "GSE111637 · umbral de edad dependiente de raza",
  xlab = "log2 FC (Geriátrico / Adulto_senior)")

