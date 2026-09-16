# =============================================================================
# Análisis de Expresión Diferencial — GSE180303
# Comparaciones:
#   a) Raza: Gigantes vs Medianos/grandes               
#   b) Sexo: Hembra vs. Macho                                          
#   c) Supervivencia: Corta vs. Larga                                                           
#   d) Edad: Geriátrico vs. Adulto
# =============================================================================

# 1. LIBRERÍAS
library(GEOquery)
library(openxlsx)
library(dplyr)
library(corrplot)
library(pheatmap)
library(survival)
library(limma)

# 2. CARGA Y EXPLORACIÓN INICIAL DE LOS DATOS 
options(download.file.method = "auto")
options(download.file.method = "libcurl")#En caso que no cargue el gset/Evita errores de conexión

gsetmicro <- getGEO("GSE180303", GSEMatrix = TRUE)
gsetmicro <- gsetmicro[[1]]   # Convertir lista a ExpressionSet
class(gsetmicro)

# 3. PREPARACIÓN DE LA MATRIZ DE EXPRESIÓN
exprmicro <- exprs(gsetmicro) # matriz de expresión (sondas × muestras)
genesmicro <- fData(gsetmicro) # anotación sondas

dim(exprmicro)
summary(as.vector(as.matrix(exprmicro)))

# 4. PREPROCESAMIENTO Y CONTROL DE CALIDAD
boxplot(exprmicro,
        las = 2,
        cex.axis = 0.6,
        main = "Distribución por muestras")

plot(density(exprmicro[,1]),
     main="Distribución de expresión",
     xlab="Expresión")
for(i in 2:ncol(exprmicro)){
  lines(density(exprmicro[,i]))}

cor_micro <- cor(exprmicro)
corrplot(cor_micro,
         method="color",
         order="hclust",
         tl.cex=0.5,
         main="Correlación entre muestras")

pca <- prcomp(t(exprmicro), scale. = TRUE)
plot(pca$x[,1],pca$x[,2],
     xlab="PC1",ylab="PC2",
     main="PCA de muestras")
text(pca$x[,1],pca$x[,2],
     labels=colnames(exprmicro),
     cex=0.6,pos=3)

dev.new(width = 16, height = 8) #Abrir otra ventada

pheatmap(cor_micro,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         main="Heatmap correlación muestras")

# 5. ANOTACIÓN Y FILTRADO DE SONDAS
tabla_expmicro <- as.data.frame(exprmicro) #Convertir a dataframe

tabla_expmicro$GeneSymbol <- genesmicro$`Gene Symbol`
tabla_expmicro$EntrezID <- genesmicro$ENTREZ_GENE_ID

#Eliminar sondas sin anotación
tabla_expmicrof <- tabla_expmicro %>%
  filter(!is.na(GeneSymbol),
         GeneSymbol != "")
dim(tabla_expmicro)
dim(tabla_expmicrof)

#Eliminar espacios
tabla_expmicrof$GeneSymbol <- trimws(tabla_expmicrof$GeneSymbol)

#Asociar cada sonda a un solo gen 
tabla_expmicrof$GeneSymbol <- sapply(
  strsplit(tabla_expmicrof$GeneSymbol, " /// "),`[`,1)

#Filtrado de baja expresión 
sum_expmicro   <- colnames(exprmicro)   # nombres originales de muestras
keep <- rowMeans(tabla_expmicrof[, sum_expmicro]) > 5
tabla_expmicrof <- tabla_expmicrof[keep, ]

cat("Sondas tras filtrado de baja expresión:", nrow(tabla_expmicrof), "\n")

# 6. SUMARIZACIÓN (SONDAS A GENES)
#Calcular varianza (mayor)
tabla_expmicrof$varianza <- apply(
  tabla_expmicrof[, sum_expmicro],
  1,var)

#Seleccionar sonda con mayor varianza
tabla_filtradamicro <- tabla_expmicrof %>%
  group_by(GeneSymbol) %>%
  slice_max(order_by = varianza, n = 1) %>%
  ungroup()

#Limpieza final
tabla_filtradamicro <- tabla_filtradamicro %>%select(-varianza)#eliminar columna auxiliar (varianza)
any(duplicated(tabla_filtradamicro$GeneSymbol))#TRUE si un valor ya aparecio antes

# 7. RENOMBRAR MUESTRAR Y EXPORTAR EXPRESIÓN 
nuevos_nombres <- paste0("Muestra_", seq_along(sum_expmicro))
nombres_map    <- setNames(nuevos_nombres, sum_expmicro)  # mapa nombre_original → Muestra_N

colnames(tabla_filtradamicro)[
  colnames(tabla_filtradamicro) %in% sum_expmicro
] <- nuevos_nombres

tabla_finalmicro <- tabla_filtradamicro %>%
  select(GeneSymbol, EntrezID, starts_with("Muestra_"))

write.xlsx(tabla_finalmicro,
           file = "GSE180303_genesxmuestras_final.xlsx",
           rowNames = FALSE)

# 8. METADATOS
metadatamicro <- pData(gsetmicro)
colnames(metadatamicro)

metadatos_limpiomicro <- metadatamicro %>%
  transmute(
    Sexo                       = as.character(`Sex:ch1`),
    Raza                       = as.character(`breed:ch1`),
    Esterilizacion             = as.character(`neutered status:ch1`),
    Edad_años                  = as.numeric(`age (years):ch1`),
    Peso_kg                    = as.numeric(`body weight (kgs):ch1`),
    N_T.carboplatino           = as.numeric(`number of carboplatin treatments:ch1`),
    N_T.doxorrubicina          = as.numeric(`number of doxorubicin treatments:ch1`),
    Tiempo_progresion_dias     = as.numeric(`time to progression (days):ch1`),
    Supervivencia_general_dias = as.numeric(`overall survival (days):ch1`),
    Conteo_monocitos           = as.numeric(`monocyte counts (1000 cells/μl):ch1`))

metadatos_limpiomicro <- metadatos_limpiomicro[colnames(exprmicro), ]#Alinear metadatos con expresión
all(colnames(exprs(gsetmicro)) == rownames(metadatos_limpiomicro)) #Verificación

rownames(metadatos_limpiomicro) <- nuevos_nombres#Renombrar filas igual que muestras

metadatos_finalmicro <- as.data.frame(t(metadatos_limpiomicro)) #Transponer metadatos (variablesxmuestras)
colnames(metadatos_finalmicro) <- nuevos_nombres

write.xlsx(metadatos_finalmicro,file = "GSE180303_metadata_organizado.xlsx",rowNames = TRUE)

summary(metadatos_limpiomicro)

# 9. CONTROL DE CALIDAD DE VARIABLES CLÍNICAS
boxplot(metadatos_limpiomicro[,c("Edad_años","Peso_kg","N_T.carboplatino","N° T.doxorrubicina",
                                 "Tiempo de progresión(días)","Supervivencia general",
                                 "Conteo de monocitos" )],
        las = 2,
        cex.axis = 0.6,
        outline = FALSE)


# 10. EXPRESIÓN DIFERENCIAL - LIMMA
#Preparar matriz de expresión
expr_matrix <- tabla_finalmicro %>%
  select(-EntrezID,-GeneSymbol) %>% as.matrix() #Quitar columna de genes y EntrezID

rownames(expr_matrix) <- tabla_finalmicro$GeneSymbol#Asignar el nombre de los genes
class(expr_matrix)

all(colnames(expr_matrix) == rownames(metadatos_limpiomicro))#Validación

##ANÁLISIS A — RAZA
table(metadatos_limpiomicro$Raza)

metadatos_limpiomicro$Grupo_Raza <- case_when(
  metadatos_limpiomicro$Raza %in% c(
    "Great Dane", "Irish Wolfhound", "Saint Bernard",
    "Newfoundland", "Great Pyrenees", "Anatolian Shep",
    "Akbash", "Rottweiler", "Mix breed", "Mix", "Malamute"
  ) ~ "Gigante",
  
  metadatos_limpiomicro$Raza %in% c(
    "Border Collie", "Australian Heeler", "Airedale", "Elkhound",
    "Standard Poodle", "Labrador Ret", "Lab", "Golden Retriever",
    "Golden Ret", "Belgian Tervuren", "Australian Shepherd",
    "Doberman", "Greyhound", "Vizsla", "German Shep",
    "Bernese Mtn Dog", "Akita"
  ) ~ "Medianas_Grande",
  
  TRUE ~ NA_character_
)

table(metadatos_limpiomicro$Grupo_Raza) #Verificar grupos
summary(metadatos_limpiomicro$Grupo_Raza)

###Expresión diferencial 
expr_genes <- tabla_filtradamicro %>%
  select(starts_with("Muestra_")) %>% as.matrix() #Matriz expresión - Raza

rownames(expr_genes) <- tabla_filtradamicro$GeneSymbol

dim(expr_genes)

(colnames(expr_genes) == rownames(metadatos_limpiomicro)) #Verificación

#Filtrar muestras con raza clasificada
idx <- !is.na(metadatos_limpiomicro$Grupo_Raza)
expr_race <- expr_genes[, idx]
meta_race <- metadatos_limpiomicro[idx, ] 

#Factor/grupo
grupo_raza <- factor(meta_race$Grupo_Raza,
                     levels = c("Gigante","Medianas_Grande"))
table(grupo_raza)

#Diseño
design_raza <- model.matrix(~ grupo_raza)
colnames(design_raza) <-
  c( "Gigante","Medianas_Grandes_vs_Gigante")
rownames(design_raza) <- rownames(meta_race)

design_raza

#Ajuste del modelo
fit_raza <- lmFit(expr_race, design_raza)
fit_raza <- eBayes(fit_raza)

#Resultados
resultados_raza <- topTable(fit_raza,
                            coef = "Medianas_Grandes_vs_Gigante",
                            number = Inf,
                            adjust.method = "BH")

head(resultados_raza)

DEG_raza <- resultados_raza %>%
  filter(adj.P.Val < 0.05 &
           abs(logFC) > 1)

#Volcano Plot
volcanoplot(fit_raza,
            coef      = "Medianas_Grandes_vs_Gigante",
            highlight = 10,
            names     = rownames(expr_race),
            main      = "Volcano: Medianas_Grande vs Gigante")

##ANÁLISIS B — SEXO
table(metadatos_limpiomicro$Sexo)

grupo_sexo <- trimws(as.character(metadatos_limpiomicro$Sexo))
grupo_sexo <- factor(grupo_sexo, levels = c("Male", "Female"))# Convertir a factor y definir referencia

table(grupo_sexo)#Verificar los grupos
length(grupo_sexo)

valid_sexo  <- !is.na(grupo_sexo)
expr_sexo   <- expr_matrix[, valid_sexo]
grupo_sexo  <- grupo_sexo[valid_sexo]
meta_sexo   <- metadatos_limpiomicro[valid_sexo, ]

#Matriz de diseño y modelo
sex_design <- model.matrix(~ grupo_sexo)
colnames(sex_design) <- c("Macho", "Hembra_vs_Macho")
rownames(sex_design) <- rownames(meta_sexo)

sex_design

#Ajustar modelo
fit_sexo <- lmFit(expr_sexo, sex_design)
fit_sexo <- eBayes(fit_sexo)

#Resultados
resultados_sexo <- topTable(
  fit_sexo,coef = "Hembra_vs_Macho",
  number = Inf,adjust.method = "BH",sort.by = "P")

DEG_sexo <- resultados_sexo %>%filter(adj.P.Val < 0.05,
    abs(logFC) > 1)

#Volcano Plot
volcanoplot(fit_sexo,
            coef= "Hembra_vs_Macho",
            highlight=10,
            names=rownames(expr_sexo),
            main = "Volcano: Hembra vs Macho")

##ANÁLISIS C - SUPERVIVENCIA
#Distribución (Umbral 180 días)
summary(metadatos_limpiomicro$Supervivencia_general_dias)
hist(metadatos_limpiomicro$Supervivencia_general_dias,
     breaks = 20,
     main   = "Distribución supervivencia (días)",
     xlab   = "Días",
     col    = "steelblue")
abline(v = 180, col = "red", lwd = 2, lty = 2)
legend("topright",
       legend = "Umbral 180 días (6 meses)",
       col    = "red", lty = 2, lwd = 2)

sum(is.na(metadatos_limpiomicro$Supervivencia_general_dias))

metadatos_limpiomicro$Supervivencia_general_dias

summary(metadatos_limpiomicro$Supervivencia_general_dias)
head(metadatos_limpiomicro$Supervivencia_general_dias,20)
sum(is.na(metadatos_limpiomicro$Supervivencia_general_dias))
str(metadatos_limpiomicro$Supervivencia_general_dias)

##Crear grupo
UMBRAL_DIAS <- 180

metadatos_limpiomicro$Grupo_Supervivencia <- case_when(
  metadatos_limpiomicro$Supervivencia_general_dias < UMBRAL_DIAS ~ "Corto_plazo",
  metadatos_limpiomicro$Supervivencia_general_dias >= UMBRAL_DIAS ~ "Largo_plazo",
  TRUE ~ NA_character_
)

table(metadatos_limpiomicro$Grupo_Supervivencia, useNA = "always")
sum(!is.na(metadatos_limpiomicro$Grupo_Supervivencia))

dim(meta_surv)
length(grupo_surv)
dim(expr_surv)
table(grupo_surv, useNA = "always")


metadatos_limpiomicro %>%
  group_by(Grupo_Supervivencia) %>%
  summarise(
    n           = n(),
    media_dias  = mean(Supervivencia_general_dias, na.rm = TRUE),
    min_dias    = min(Supervivencia_general_dias,  na.rm = TRUE),
    max_dias    = max(Supervivencia_general_dias,  na.rm = TRUE)
  )

boxplot(
  Supervivencia_general_dias ~ Grupo_Supervivencia,
  data    = metadatos_limpiomicro,
  col     = c("salmon", "steelblue"),
  ylab    = "Días",
  main    = "Supervivencia por grupo",
  outline = FALSE
)
abline(h = UMBRAL_DIAS, lty = 2, col = "red")


valid_surv <- !is.na(metadatos_limpiomicro$Grupo_Supervivencia)

expr_surv  <- expr_matrix[, valid_surv]          # submatriz de expresión
meta_surv  <- metadatos_limpiomicro[valid_surv, ] # submetadata

cat("Muestras en análisis de supervivencia:", ncol(expr_surv), "\n")
cat("  Corto plazo (<180 días):",
    sum(meta_surv$Grupo_Supervivencia == "Corto_plazo"), "\n")
cat("  Largo plazo (≥180 días):",
    sum(meta_surv$Grupo_Supervivencia == "Largo_plazo"), "\n")

all(colnames(expr_surv) == rownames(meta_surv))# Verificar 

#Diseño y factor
grupo_surv <- factor(
  meta_surv$Grupo_Supervivencia,
  levels = c("Largo_plazo","Corto_plazo")   # Largo_plazo = referencia
)

design_surv <- model.matrix(
  ~ grupo_surv)

colnames(design_surv) <-
  c(
    "Largo_plazo",
    "Corto_plazo_vs_Largo_plazo")

rownames(design_surv) <-
  rownames(meta_surv)

table(design_surv)   

fit_surv  <- lmFit(expr_surv, design_surv)
fit_surv  <- eBayes(fit_surv)

# Resultados 
resultados_surv <- topTable(
  fit_surv,
  coef = "Corto_plazo_vs_Largo_plazo",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P")

DEG_surv <- resultados_surv %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) > 1)

#Volcano Plot
volcanoplot(fit_surv,
            coef= "Corto_plazo_vs_Largo_plazo",
            highlight=10,
            names=rownames(expr_surv),
            main = "Volcano: Supervivencia: Corto plazo vs Largo plazo")

##ANÁLISIS D - EDAD

class(metadatos_limpiomicro$Edad_años)   #"numeric"
summary(metadatos_limpiomicro$Edad_años)
sum(is.na(metadatos_limpiomicro$Edad_años))

#Clasificar edad según umbral dependiente de la raza 
table(metadatos_limpiomicro$Grupo_Raza, useNA = "always")

metadatos_limpiomicro$Grupo_Edad <- case_when(
  metadatos_limpiomicro$Grupo_Raza == "Medianas_Grande" &
    metadatos_limpiomicro$Edad_años >= 10 ~ "Geriatrico",
  metadatos_limpiomicro$Grupo_Raza == "Medianas_Grande" &
    metadatos_limpiomicro$Edad_años <  10 ~ "Adulto_senior",
  
  metadatos_limpiomicro$Grupo_Raza == "Gigante" &
    metadatos_limpiomicro$Edad_años >= 8  ~ "Geriatrico",
  metadatos_limpiomicro$Grupo_Raza == "Gigante" &
    metadatos_limpiomicro$Edad_años <  8  ~ "Adulto_senior",
  
  TRUE ~ NA_character_ 
)

table(metadatos_limpiomicro$Grupo_Edad, useNA = "always")

# Resumen de verificación: edad real por grupo y raza
metadatos_limpiomicro %>%
  filter(!is.na(Grupo_Edad)) %>%
  group_by(Grupo_Raza, Grupo_Edad) %>%
  summarise(
    n          = n(),
    edad_media = mean(Edad_años),
    edad_min   = min(Edad_años),
    edad_max   = max(Edad_años),
    .groups    = "drop"
  ) %>%
  print()

# Visualización: distribución de edad por raza con líneas de umbral
boxplot(
  Edad_años ~ Grupo_Raza,
  data    = metadatos_limpiomicro,
  col     = c("#E74C3C", "#3498DB"),
  ylab    = "Edad (años)",
  main    = "Edad por grupo de raza — GSE180303",
  outline = FALSE)

# Líneas de umbral correspondientes a cada grupo (posición x aproximada)
segments(x0 = 0.6, x1 = 1.4, y0 = 8,  col = "red",  lty = 2, lwd = 2)  # Gigante
segments(x0 = 1.6, x1 = 2.4, y0 = 10, col = "red",  lty = 2, lwd = 2)  # Medianas_Grande
legend("topright",
       legend = "Umbral Geriátrico (8 años Gigante / 10 años Mediana)",
       col = "red", lty = 2, lwd = 2, bty = "n", cex = 0.7)

##Expresión diferencial
# Fitrar muestras validas
valid_edad <- !is.na(metadatos_limpiomicro$Grupo_Edad)
expr_edad  <- expr_matrix[, valid_edad]
meta_edad  <- metadatos_limpiomicro[valid_edad, ]

# Verificar alineación
stopifnot(all(colnames(expr_edad) == rownames(meta_edad)))

# Grupo y Diseño
grupo_edad <- factor(
  meta_edad$Grupo_Edad,
  levels = c(
    "Adulto_senior",
    "Geriatrico"
  ))

table(grupo_edad)

design_edad <- model.matrix(
  ~ grupo_edad)
colnames(design_edad) <-
  c(
    "Adulto_senior",
    "Geriatrico_vs_Adulto_senior")
rownames(design_edad) <-rownames(meta_edad)

design_edad

fit_edad <- lmFit(expr_edad, design_edad)
fit_edad <- eBayes(fit_edad)

# Resultados
resultados_edad <- topTable(
  fit_edad,
  coef = "Geriatrico_vs_Adulto_senior",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P")

DEG_edad <- resultados_edad %>%
  filter(adj.P.Val < 0.05,abs(logFC) > 1)

# Volcano plot 
volcanoplot(fit_edad,
            coef      = "Geriatrico_vs_Adulto_senior",
            highlight = 10,
            names     = rownames(expr_edad),
            main      = "Volcano: Geriátrico vs Adulto_senior")
abline(v = c(-1, 1),      lty = 2, col = "gray50")
abline(h = -log10(0.05),  lty = 2, col = "red")









