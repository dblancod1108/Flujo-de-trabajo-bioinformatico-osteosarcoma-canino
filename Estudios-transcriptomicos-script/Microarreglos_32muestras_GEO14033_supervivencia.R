# =============================================================================
# Análisis de Expresión Diferencial — GSE14033
# Comparación:
#   a) Supervivencia: Corta vs. Larga                                                           
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
gsetmicro <- getGEO("GSE14033", GSEMatrix = TRUE, getGPL = TRUE)
gsetmicro <- gsetmicro[[1]]
class(gsetmicro)

# 3. PREPARACIÓN DE LA MATRIZ DE EXPRESIÓN
exprmicro32  <- exprs(gsetmicro)      # matriz sondas x muestras (log2-ratios)
genesmicro32 <- fData(gsetmicro)      # anotacion de sondas

dim(exprmicro32)
summary(as.vector(as.matrix(exprmicro32)))
colnames(genesmicro32)

# 4. PREPROCESAMIENTO Y CONTROL DE CALIDAD PRELIMINAR
cat("Matriz inicial:",nrow(exprmicro32), "sondas ×",ncol(exprmicro32), "muestras\n")

# Valores faltantes
cat("Valores NA:",sum(is.na(exprmicro32)),"\n")

# Valores infinitos
cat("Valores infinitos:",sum(is.infinite(exprmicro32)),"\n")

# Varianza por sonda
varianza_inicial32 <- apply(exprmicro32,1,var,na.rm = TRUE)

cat("Sondas con varianza cero:",sum(varianza_inicial32 == 0, na.rm = TRUE),"\n")

boxplot(
  exprmicro32,
  las = 2,
  cex.axis = 0.6,
  outline = FALSE,
  main = "Distribución por muestras — GSE14033",
  ylab = "log2 ratio (M-value)")

abline(h = 0,col = "blue",lty = 2)

plot(
  density(
    exprmicro32[, 1],
    na.rm = TRUE
  ),
  main = "Distribución de expresión — GSE14033",
  xlab = "log2 ratio (M-value)",
  col = "gray40"
)

for (i in 2:ncol(exprmicro32)) {
  
  lines(
    density(
      exprmicro32[, i],
      na.rm = TRUE
    ),
    col = "gray70"
  )}

abline(v = 0,col = "blue",lty = 2)

cor_micro32 <- cor(
  exprmicro32,
  use = "pairwise.complete.obs")

corrplot(
  cor_micro32,
  method = "color",
  order = "hclust",
  tl.cex = 0.5,
  main = "Correlación entre muestras — GSE14033")

keep_pca32 <- is.finite(varianza_inicial32) &varianza_inicial32 > 0

expr_pca32 <- exprmicro32[keep_pca32,]

# Eliminar sondas con valores NA
keep_complete32 <- complete.cases(expr_pca32)

expr_pca32 <- expr_pca32[keep_complete32,]

cat("Sondas utilizadas para PCA:",nrow(expr_pca32),"\n")

pca32 <- prcomp(
  t(expr_pca32),
  center = TRUE,
  scale. = TRUE)

var_exp32 <- round(100 * pca32$sdev^2 /sum(pca32$sdev^2),1)

dev.new(width = 12, height = 10)
plot(
  pca32$x[, 1],
  pca32$x[, 2],
  xlab = paste0(
    "PC1 (",
    var_exp32[1],
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    var_exp32[2],
    "%)"
  ),
  main = "PCA de muestras — GSE14033",
  pch = 19,
  col = ifelse(
    grepl(
      "SHORT",
      colnames(exprmicro32),
      ignore.case = TRUE
    ),
    "#E74C3C",
    "#3498DB"))

text(
  pca32$x[, 1],
  pca32$x[, 2],
  labels = colnames(exprmicro32),
  cex = 0.6,
  pos = 3
)

legend(
  "topleft",
  legend = c(
    "Corta",
    "Larga"
  ),
  col = c(
    "#E74C3C",
    "#3498DB"
  ),
  pch = 19,
  bty = "n")

# 5. ANOTACIÓN Y FILTRADO DE SONDAS
tabla_expmicro32 <- as.data.frame(exprmicro32)

tabla_expmicro32$GeneSymbol <- genesmicro32$GENE_SYMBOL
tabla_expmicro32$GB_ACC     <- genesmicro32$GB_ACC   # GenBank como ID alternativo

dim(tabla_expmicro32)

# Eliminar sondas sin gen
tabla_expmicrof32 <- tabla_expmicro32 %>%
  filter(!is.na(GeneSymbol),
         GeneSymbol != "")

dim(tabla_expmicrof32)

tabla_expmicrof32$GeneSymbol <- trimws(tabla_expmicrof32$GeneSymbol)# Eliminar columna GeneSymbol 

# Si hay multiples simbolos separados por " " espacio (GPL5117), quedarse con el primero
tabla_expmicrof32$GeneSymbol <- sapply(
  strsplit(tabla_expmicrof32$GeneSymbol, " "), `[`, 1)

# Inspeccionar 
problematicos <- tabla_expmicrof32$GeneSymbol[grepl(" ", tabla_expmicrof32$GeneSymbol)]
head(problematicos, 20)
head(nchar(problematicos), 10)# exactos

tabla_expmicrof32 <- tabla_expmicro32 %>%
  filter(!is.na(GeneSymbol),
         GeneSymbol != "",
         GeneSymbol != "0",
         !grepl("^[0-9]+$", GeneSymbol))  # elimina entradas puramente numéricas

# Verificar: no debe quedar ningún GeneSymbol con espacio interno
sum(grepl(" ", tabla_expmicrof32$GeneSymbol))  
head(sort(tabla_expmicrof32$GeneSymbol), 20)   # revisar visualmente

# FILTRADO DE SONDAS DE BAJA SEÑAL
#  NOTA: En arrays de 2 colores los valores son log2-ratios y PUEDEN ser
#  negativos. El umbral media > 5 aplica en arrays de 1 canal (Affymetrix).
#  Aqui se usa un criterio de varianza mínima: se eliminan sondas en el
#  percentil 10 inferior de varianza (sondas sin señal biologica util).

sum_expmicro32 <- colnames(exprmicro32)   # nombres originales de muestras

varianza_ini <- apply(tabla_expmicrof32[, sum_expmicro32], 1, var, na.rm = TRUE)
umbral_var   <- quantile(varianza_ini, 0.10, na.rm = TRUE)
keep         <- varianza_ini > umbral_var
tabla_expmicrof32 <- tabla_expmicrof32[keep, ]

cat("Sondas tras filtrado por varianza (>P10):", nrow(tabla_expmicrof32), "\n")

# 6. SUMARIZACIÓN (SONDAS A GENES)
# Criterio: sonda con MAYOR varianza por gen (mayor información biológica)
tabla_expmicrof32$varianza <- apply(
  tabla_expmicrof32[, sum_expmicro32], 1, var, na.rm = TRUE)

tabla_filtradamicro32 <- tabla_expmicrof32 %>%
  group_by(GeneSymbol) %>%
  slice_max(order_by = varianza, n = 1, with_ties = FALSE) %>%
  ungroup()

tabla_filtradamicro32 <- tabla_filtradamicro32 %>% select(-varianza)# Limpieza de columna auxiliar

any(duplicated(tabla_filtradamicro32$GeneSymbol))  #FALSE 

cat("Genes unicos tras sumarización:", nrow(tabla_filtradamicro32), "\n")

# 7. RENOMBRAR MUESTRAR Y EXPORTAR EXPRESIÓN 
nuevos_nombres32 <- paste0("Muestra_", seq_along(sum_expmicro32))
nombres_map32    <- setNames(nuevos_nombres32, sum_expmicro32)

colnames(tabla_filtradamicro32)[colnames(tabla_filtradamicro32) 
        %in% sum_expmicro32] <- nuevos_nombres32

tabla_finalmicro32 <- tabla_filtradamicro32 %>%
  select(GeneSymbol, GB_ACC, starts_with("Muestra_"))

write.xlsx(tabla_finalmicro32,file= "GSE14033_genesxmuestras_final.xlsx",rowNames = FALSE)

# 8. METADATOS
metadatamicro32 <- pData(gsetmicro)
colnames(metadatamicro32)

metadatos_limpiomicro32 <- metadatamicro32 %>%
  transmute(
    Titulo         = as.character(title),
    GEO_ID         = as.character(geo_accession),
    Grupo_Survival = ifelse(grepl("SHORT", title, ignore.case = TRUE),
                            "Short", "Long"),
    # Extraer dias de supervivencia de characteristics_ch1
    # Ejemplo: "Short_36 Days" -> 36
    Dias_Survival  = as.numeric(
      gsub("[^0-9]", "",
           gsub("[Dd]ays?", "",
                as.character(characteristics_ch1)))),
    Canal_Cy5      = as.character(source_name_ch1),  # muestra tumoral
    Canal_Cy3      = as.character(source_name_ch2)   # referencia pool
  )

# Validación
metadatos_limpiomicro32 <- metadatos_limpiomicro32[colnames(exprmicro32), ]
all(colnames(exprs(gsetmicro)) == rownames(metadatos_limpiomicro32))  

# Renombrar filas igual que las muestras renombradas
rownames(metadatos_limpiomicro32) <- nuevos_nombres32

# Exportar metadatos transpuestos (variables x muestras)
metadatos_finalmicro32 <- as.data.frame(t(metadatos_limpiomicro32))
colnames(metadatos_finalmicro32) <- nuevos_nombres32

write.xlsx(metadatos_finalmicro, file = "GSE14033_metadata_organizado.xlsx",rowNames = TRUE)

summary(metadatos_limpiomicro32)
table(metadatos_limpiomicro32$Grupo_Survival)

# 9. CONTROL DE CALIDAD 
boxplot(Dias_Survival ~ Grupo_Survival,
        data    = metadatos_limpiomicro32,
        col     = c("#3498DB", "#E74C3C"),
        main    = "Supervivencia por grupo — GSE14033",
        ylab    = "Dias de supervivencia",
        xlab    = "Grupo",
        outline = FALSE)

stripchart(Dias_Survival ~ Grupo_Survival,
           data     = metadatos_limpiomicro32,
           method   = "jitter",
           pch      = 19, cex = 0.8, col = "black",
           add      = TRUE, vertical = TRUE)


# 10. EXPRESIÓN DIFERENCIAL - LIMMA
#Preparar matriz de expresión
expr_matrix32 <- tabla_finalmicro32 %>%
  select(-GB_ACC, -GeneSymbol) %>%
  as.matrix()

rownames(expr_matrix32) <- tabla_finalmicro32$GeneSymbol
class(expr_matrix32)

all(colnames(expr_matrix32) == rownames(metadatos_limpiomicro32))  # validación
table(metadatos_limpiomicro32$Grupo_Survival)

#Factor
grupo_survival32 <- trimws(as.character(metadatos_limpiomicro32$Grupo_Survival))

grupo_survival32 <- factor(grupo_survival32,levels = c("Long", "Short"))

table(grupo_survival32, useNA = "always")

#Seleccionar muestras validas
valid_surv32 <-!is.na(grupo_survival32)
expr_surv32 <-expr_matrix32[,valid_surv32]
grupo_surv32 <-grupo_survival32[valid_surv32]

meta_surv32 <-metadatos_limpiomicro32[valid_surv32,]

cat("Muestras utilizadas:",ncol(expr_surv32),"\n")
cat("Supervivencia corta:",sum(grupo_surv32 == "Short"),"\n")
cat("Supervivencia larga:",sum(grupo_surv32 == "Long"),"\n")

# Matriz de diseño
surv_design32 <- model.matrix(~ grupo_surv32)
colnames(surv_design32)
surv_design32

# Ajustar modelo
fit_surv32 <- lmFit(expr_surv32, surv_design32)
fit_surv32 <- eBayes(fit_surv32)

#Resultados
resultados_surv32 <- topTable(fit_surv32,
                            coef          = "grupo_surv32Short",
                            number        = Inf,
                            adjust.method = "BH")

DEG_surv32 <- resultados_surv32 %>%
  dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 1)


# Volcano plot
volcanoplot(fit_surv32,
            coef      = "grupo_surv32Short",
            highlight = 15,
            names     = rownames(expr_surv32),
            main      = "Volcano: Short vs Long survivors\nGSE14033",
            xlab      = "log2 Fold Change (Short / Long)")
abline(v = c(-1, 1),      lty = 2, col = "gray50")
abline(h = -log10(0.05),  lty = 2, col = "red")
