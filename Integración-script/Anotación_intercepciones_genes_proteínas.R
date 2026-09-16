# =============================================================================
#  ANOTACIÓN E INTEGRACIÓN DE CAPAS PROTEÓMICA Y TRANSCRIPTÓMICA
# =============================================================================

library(readxl)
library(writexl)
library(httr)
library(jsonlite)
library(ggplot2)
library(ggrepel)
library(ggVennDiagram)
library(patchwork)
library(biomaRt)
library(dplyr)

# Conectar a Ensembl canino
mart_perro <- useEnsembl(
  biomart = "genes",
  dataset = "clfamiliaris_gene_ensembl")

mart_perro
listFilters(mart_perro)
listAttributes(mart_perro)

#Validación identificadores BioMart canino
"uniprot_gn_id" %in% listFilters(mart_perro)$name
c(
  "uniprot_gn_id",
  "external_gene_name",
  "ensembl_gene_id",
  "description"
) %in% listAttributes(mart_perro)$name

#Construcción del conjunto de accession
todos_accessions <- unique(c(
  Resultados_sexo_proteómica$Accesion,
  Resultados_edad_proteómica$Accesion,
  Resultados_supervivencia_proteómica$Accesion))

#Extracción del mapeo UniProt → Gene Symbol
mapeo_ensembl_raw <- getBM(
  attributes = c(
    "uniprot_gn_id",
    "external_gene_name",
    "ensembl_gene_id",
    "description"
  ),
  filters = "uniprot_gn_id",
  values  = todos_accessions,
  mart    = mart_perro)

cat("Total accessions únicos a mapear:", length(todos_accessions), "\n")

# FUENTE 1 — biomaRt / Ensembl

# Intentar también con uniprotswissprot para las proteínas que no mapearon
sin_mapeo_ensembl <- setdiff(todos_accessions, mapeo_ensembl_raw$uniprot_gn_id)
cat("Accessions mapeados en primer intento:", nrow(mapeo_ensembl_raw), "\n")
cat("Sin mapeo, intentando filtro alternativo:", length(sin_mapeo_ensembl), "\n")

if (length(sin_mapeo_ensembl) > 0) {
  
  mapeo_ensembl_alt <- getBM(
    attributes = c(
      "uniprotswissprot",
      "external_gene_name",
      "ensembl_gene_id",
      "description"
    ),
    filters = "uniprotswissprot",
    values = sin_mapeo_ensembl,
    mart = mart_perro
  ) %>%
    mutate(
      uniprotswissprot = as.character(uniprotswissprot),
      external_gene_name = as.character(external_gene_name),
      ensembl_gene_id = as.character(ensembl_gene_id),
      description = as.character(description)
    ) %>%
    rename(
      uniprot_gn_id = uniprotswissprot
    )
  
  # Asegurar que ambas tablas tengan tipos compatibles
  mapeo_ensembl_raw <- mapeo_ensembl_raw %>%
    mutate(
      uniprot_gn_id = as.character(uniprot_gn_id),
      external_gene_name = as.character(external_gene_name),
      ensembl_gene_id = as.character(ensembl_gene_id),
      description = as.character(description)
    )
  
  mapeo_ensembl_raw <- bind_rows(
    mapeo_ensembl_raw,
    mapeo_ensembl_alt
  )}

# Limpiar: un Gene Symbol por Accession (priorizar el más informativo)
mapeo_ensembl <- mapeo_ensembl_raw %>%
  filter(!is.na(external_gene_name) & external_gene_name != "") %>%
  group_by(uniprot_gn_id) %>%
  slice(1) %>%
  ungroup() %>%
  rename(
    Accesion        = uniprot_gn_id,
    Gene_Ensembl    = external_gene_name,
    Ensembl_ID      = ensembl_gene_id,
    Descripcion_Ensembl = description
  )

cat("Accessions con Gene Symbol Ensembl:", nrow(mapeo_ensembl), "\n")
cat("Accessions sin cobertura Ensembl:",
    length(setdiff(todos_accessions, mapeo_ensembl$Accesion)), "\n")

# FUENTE 2 — UniProt REST API (para los TrEMBL sin cobertura en Ensembl)
# Función corregida con simplifyVector = FALSE
get_uniprot_annotation <- function(accession) {
  url  <- paste0("https://rest.uniprot.org/uniprotkb/", accession, ".json")
  resp <- tryCatch(GET(url), error = function(e) NULL)
  
  if (is.null(resp) || status_code(resp) != 200) {
    return(data.frame(
      Accesion        = accession,
      Gene_UniProt    = NA_character_,
      Ensembl_API     = NA_character_,
      Descripcion_API = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  
  # simplifyVector = FALSE evita el error de longitud 0
  datos <- fromJSON(
    content(resp, "text", encoding = "UTF-8"),
    simplifyVector = FALSE
  )
  
  # Gene Symbol
  gene <- tryCatch({
    val <- datos$genes[[1]]$geneName$value
    if (length(val) == 0 || is.null(val)) NA_character_
    else as.character(val[1])
  }, error = function(e) NA_character_)
  
  # Ensembl ID desde cross-references
  ensembl <- tryCatch({
    refs <- datos$uniProtKBCrossReferences
    ensembl_entries <- Filter(function(x) {
      !is.null(x$database) && x$database == "Ensembl"
    }, refs)
    if (length(ensembl_entries) == 0) NA_character_
    else as.character(ensembl_entries[[1]]$id[1])
  }, error = function(e) NA_character_)
  
  # Descripción (recommended name o submitted name como respaldo)
  desc <- tryCatch({
    val <- datos$proteinDescription$recommendedName$fullName$value
    if (length(val) == 0 || is.null(val)) {
      val2 <- datos$proteinDescription$submittedName[[1]]$fullName$value
      if (length(val2) == 0 || is.null(val2)) NA_character_
      else as.character(val2[1])
    } else as.character(val[1])
  }, error = function(e) NA_character_)
  
  data.frame(
    Accesion        = accession,
    Gene_UniProt    = gene,
    Ensembl_API     = ensembl,
    Descripcion_API = desc,
    stringsAsFactors = FALSE
  )
}

# Consultar solo los que NO mapearon con Ensembl (evita consultas innecesarias)
accessions_para_api <- setdiff(todos_accessions, mapeo_ensembl$Accesion)

#Consulta masiva accession pendientes (Vía API)
if (length(accessions_para_api) > 0) {
  mapeo_api <- bind_rows(
    lapply(seq_along(accessions_para_api), function(i) {
      acc <- accessions_para_api[i]
      if (i %% 20 == 0) cat("  Procesados:", i, "/",
                            length(accessions_para_api), "\n")
      Sys.sleep(0.15)   # pausa para no saturar la API
      get_uniprot_annotation(acc)
    })
  )
  cat("Respuestas recibidas de UniProt API:", nrow(mapeo_api), "\n")
  cat("Con Gene Symbol vía API:",
      sum(!is.na(mapeo_api$Gene_UniProt)), "\n")
} else {
  cat("Todos los accessions mapearon con Ensembl — API no necesaria.\n")
  mapeo_api <- data.frame(
    Accesion        = character(),
    Gene_UniProt    = character(),
    Ensembl_API     = character(),
    Descripcion_API = character()
  )
}

# CONSOLIDACIÓN: Ensembl + API con priorización explícita
# Todos los accessions como base
base_accessions <- data.frame(Accesion = todos_accessions,
                              stringsAsFactors = FALSE)

mapeo_consolidado <- base_accessions %>%
  
  # Unir resultados de Ensembl
  left_join(
    mapeo_ensembl %>%
      select(Accesion, Gene_Ensembl, Ensembl_ID, Descripcion_Ensembl),
    by = "Accesion"
  ) %>%
  
  # Unir resultados de UniProt API
  left_join(
    mapeo_api %>%
      select(Accesion, Gene_UniProt, Ensembl_API, Descripcion_API),
    by = "Accesion"
  ) %>%
  
  # Reglas de priorización
  mutate(
    
    # Gene Symbol: Ensembl primero, API como respaldo
    Gene_Symbol = case_when(
      !is.na(Gene_Ensembl) & Gene_Ensembl != "" ~ Gene_Ensembl,
      !is.na(Gene_UniProt) & Gene_UniProt != "" ~ Gene_UniProt,
      TRUE ~ NA_character_
    ),
    
    # Ensembl ID: biomaRt primero, API como respaldo
    Ensembl_ID_final = case_when(
      !is.na(Ensembl_ID)  ~ Ensembl_ID,
      !is.na(Ensembl_API) ~ Ensembl_API,
      TRUE ~ NA_character_
    ),
    
    # Descripción: Ensembl primero, API como respaldo
    Descripcion_final = case_when(
      !is.na(Descripcion_Ensembl) ~ Descripcion_Ensembl,
      !is.na(Descripcion_API)     ~ Descripcion_API,
      TRUE ~ NA_character_
    ),
    
    # Fuente del Gene Symbol (para trazabilidad)
    Fuente_Gene = case_when(
      !is.na(Gene_Ensembl) & Gene_Ensembl != "" ~ "Ensembl/biomaRt",
      !is.na(Gene_UniProt) & Gene_UniProt != "" ~ "UniProt API",
      TRUE ~ "Sin anotación"
    )
  ) %>%
  
  select(Accesion, Gene_Symbol, Ensembl_ID_final, Descripcion_final,
         Gene_Ensembl, Gene_UniProt, Fuente_Gene,
         Descripcion_Ensembl, Descripcion_API)

# RESUMEN Y EXPORTACIÓN
cat("Total accessions procesados:", nrow(mapeo_consolidado), "\n")
cat("Con Gene Symbol (cualquier fuente):",
    sum(!is.na(mapeo_consolidado$Gene_Symbol)), "\n")
cat("  → Fuente Ensembl/biomaRt:",
    sum(mapeo_consolidado$Fuente_Gene == "Ensembl/biomaRt", na.rm = TRUE), "\n")
cat("  → Fuente UniProt API:    ",
    sum(mapeo_consolidado$Fuente_Gene == "UniProt API",     na.rm = TRUE), "\n")
cat("Sin anotación en ninguna fuente:",
    sum(mapeo_consolidado$Fuente_Gene == "Sin anotación",   na.rm = TRUE), "\n")
cat("\nCobertura total:",
    round(100 * mean(!is.na(mapeo_consolidado$Gene_Symbol)), 1), "%\n")

write_xlsx(list(
  Consolidado_completo = mapeo_consolidado,
  Solo_con_anotacion   = mapeo_consolidado %>%
    filter(!is.na(Gene_Symbol)),
  Sin_anotacion        = mapeo_consolidado %>%
    filter(is.na(Gene_Symbol)) %>%
    select(Accesion, Fuente_Gene)
),"mapeo_UniProt_GeneSymbol_consolidado.xlsx")


#### Guardar el objeto limpio(solo proteínas con Gene Symbol) para integración transcriptómica-proteómica
mapeo_limpio <- mapeo_consolidado %>%
  filter(!is.na(Gene_Symbol)) %>%
  select(Accesion, Gene_Symbol, Ensembl_ID_final,
         Descripcion_final, Fuente_Gene)

# VERIFICACIÓN PREVIA 
# Verificar que mapeo_limpio existe y tiene las columnas necesarias
stopifnot(
  "mapeo_limpio no está en el entorno — ejecuta el Paso 2 primero" =
    exists("mapeo_limpio"),
  "mapeo_limpio no tiene la columna Accesion" =
    "Accesion" %in% colnames(mapeo_limpio),
  "mapeo_limpio no tiene la columna Gene_Symbol" =
    "Gene_Symbol" %in% colnames(mapeo_limpio)
)


cat("✓ Resultados Sexo disponibles:",nrow(Resultados_sexo_proteómica),  "proteínas\n")
cat("✓ Resultados Edad disponibles:",nrow(Resultados_edad_proteómica),  "proteínas\n")
cat("✓ Resultados Supervivencia disponibles:",nrow(Resultados_supervivencia_proteómica),  "proteínas\n")

# FUNCIÓN CENTRAL DE ANOTACIÓN
# Recibe una tabla DEA y el mapeo, devuelve la tabla enriquecida con
# Gene Symbol, Ensembl ID, descripción funcional y fuente de anotación.
# Incluye diagnóstico detallado de cobertura.

anotar_tabla_dea <- function(tabla_dea,
                             mapeo,
                             nombre_comparacion) {
  
  cat("\n──────────────────────────────────────────\n")
  cat("Anotando:", nombre_comparacion, "\n")
  cat("──────────────────────────────────────────\n")
  
  # Verificar que la columna Accesion existe en la tabla DEA
  if (!"Accesion" %in% colnames(tabla_dea)) {
    stop(paste("La tabla de", nombre_comparacion,
               "no tiene columna 'Accesion'. Columnas disponibles:",
               paste(colnames(tabla_dea), collapse = ", ")))
  }
  
  n_total <- nrow(tabla_dea)
  
  # ── Join con el mapeo ──────────────────────────────────
  tabla_anotada <- tabla_dea %>%
    left_join(
      mapeo %>%
        select(Accesion, Gene_Symbol, Ensembl_ID_final,
               Descripcion_final, Fuente_Gene),
      by = "Accesion"
    ) %>%
    # Reorganizar columnas: identificadores al frente, luego estadísticos
    select(
      Accesion,
      Gene_Symbol,
      Ensembl_ID_final,
      Descripcion_final,
      Fuente_Gene,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B,
      Regulacion,
      everything()   # cualquier columna adicional que ya existiera
    )
  
  # ── Diagnóstico de cobertura ───────────────────────────
  n_con_symbol   <- sum(!is.na(tabla_anotada$Gene_Symbol))
  n_sin_symbol   <- sum( is.na(tabla_anotada$Gene_Symbol))
  pct_cobertura  <- round(100 * n_con_symbol / n_total, 1)
  
  cat("Proteínas totales:          ", n_total, "\n")
  cat("Con Gene Symbol:            ", n_con_symbol,
      "(", pct_cobertura, "% )\n")
  cat("Sin Gene Symbol:            ", n_sin_symbol, "\n")
  
  # Desglose por fuente
  if ("Fuente_Gene" %in% colnames(tabla_anotada)) {
    cat("\nDesglose por fuente de anotación:\n")
    print(table(tabla_anotada$Fuente_Gene, useNA = "always"))
  }
  
  # Cobertura específica entre las proteínas significativas
  sig_fdr <- tabla_anotada %>%
    filter(Regulacion != "No significativo")
  sig_nom <- tabla_anotada %>%
    filter(!is.na(Regulacion) &
             grepl("nominal", Regulacion, ignore.case = TRUE))
  
  if (nrow(sig_fdr) > 0) {
    pct_sig <- round(100 * mean(!is.na(sig_fdr$Gene_Symbol)), 1)
    cat("\nProteínas significativas (FDR):", nrow(sig_fdr),
        "| Con Gene Symbol:", pct_sig, "%\n")
  }
  
  # Advertencia si la cobertura es baja
  if (pct_cobertura < 60) {
    cat("\n⚠ ADVERTENCIA: cobertura de anotación < 60%.\n")
    cat("  Considera revisar si los accessions están en formato\n")
    cat("  UniProt estándar o si hay accessions obsoletos.\n")
  }
  
  return(tabla_anotada)
}

# APLICAR ANOTACIÓN A LAS TRES COMPARACIONES
prot_sexo_anotada <- anotar_tabla_dea(
  Resultados_sexo_proteómica,
  mapeo_limpio,
  "Sexo — Hembra vs Macho")

prot_edad_anotada <- anotar_tabla_dea(
  Resultados_edad_proteómica,
  mapeo_limpio,
  "Edad — Geriátrico vs Adulto")

prot_supervivencia_anotada <- anotar_tabla_dea(
  Resultados_supervivencia_proteómica,
  mapeo_limpio,
  "Supervivencia — Corta vs Larga")

# RESUMEN COMPARATIVO DE COBERTURA 
resumen_cobertura <- data.frame(
  Comparacion = c("Sexo", "Edad", "Supervivencia"),
  N_proteinas = c(
    nrow(prot_sexo_anotada),
    nrow(prot_edad_anotada),
    nrow(prot_supervivencia_anotada)
  ),
  Con_Gene_Symbol = c(
    sum(!is.na(prot_sexo_anotada$Gene_Symbol)),
    sum(!is.na(prot_edad_anotada$Gene_Symbol)),
    sum(!is.na(prot_supervivencia_anotada$Gene_Symbol))
  ),
  Cobertura_pct = c(
    round(100 * mean(!is.na(prot_sexo_anotada$Gene_Symbol)), 1),
    round(100 * mean(!is.na(prot_edad_anotada$Gene_Symbol)), 1),
    round(100 * mean(!is.na(prot_supervivencia_anotada$Gene_Symbol)), 1)
  ),
  Fuente_Ensembl = c(
    sum(prot_sexo_anotada$Fuente_Gene == "Ensembl/biomaRt", na.rm = TRUE),
    sum(prot_edad_anotada$Fuente_Gene == "Ensembl/biomaRt", na.rm = TRUE),
    sum(prot_supervivencia_anotada$Fuente_Gene == "Ensembl/biomaRt", na.rm = TRUE)
  ),
  Fuente_UniProt_API = c(
    sum(prot_sexo_anotada$Fuente_Gene == "UniProt API", na.rm = TRUE),
    sum(prot_edad_anotada$Fuente_Gene == "UniProt API", na.rm = TRUE),
    sum(prot_supervivencia_anotada$Fuente_Gene == "UniProt API", na.rm = TRUE)
  )
)

print(resumen_cobertura)

# TABLAS DE PROTEÍNAS SIN ANOTACIÓN 
sin_anotacion_sexo <- prot_sexo_anotada %>%
  filter(is.na(Gene_Symbol)) %>%
  select(Accesion, logFC, adj.P.Val, Regulacion)

sin_anotacion_edad <- prot_edad_anotada %>%
  filter(is.na(Gene_Symbol)) %>%
  select(Accesion, logFC, adj.P.Val, Regulacion)

sin_anotacion_supervivencia <- prot_supervivencia_anotada %>%
  filter(is.na(Gene_Symbol)) %>%
  select(Accesion, logFC, adj.P.Val, Regulacion)

cat("  Sexo:          ", nrow(sin_anotacion_sexo), "\n")
cat("  Edad:          ", nrow(sin_anotacion_edad), "\n")
cat("  Supervivencia: ", nrow(sin_anotacion_supervivencia), "\n")

# Verificar si alguna proteína significativa quedó sin anotar
sig_sin_anot_edad <- prot_edad_anotada %>%
  filter(Regulacion != "No significativo" & is.na(Gene_Symbol))

if (nrow(sig_sin_anot_edad) > 0) {
  cat("\n⚠ Proteínas SIGNIFICATIVAS sin Gene Symbol en Edad:\n")
  print(sig_sin_anot_edad %>% select(Accesion, logFC, adj.P.Val))
  cat("  Estas proteínas no podrán participar en la integración.\n")
  cat("  Considera buscarlas manualmente en uniprot.org\n")}

# TABLAS FILTRADAS PARA LA INTEGRACIÓN
prot_sexo_para_integracion <- prot_sexo_anotada %>%
  filter(!is.na(Gene_Symbol))

prot_edad_para_integracion <- prot_edad_anotada %>%
  filter(!is.na(Gene_Symbol))

prot_surv_para_integracion <- prot_supervivencia_anotada %>%
  filter(!is.na(Gene_Symbol))

#CONJUNTOS DE REFERENCIA PROTEÓMICOS PARA INTEGRACIÓN
cat("Sexo:          ", nrow(prot_sexo_para_integracion),
    "proteínas con Gene Symbol\n")
cat("Edad:          ", nrow(prot_edad_para_integracion),
    "proteínas con Gene Symbol\n")
cat("Supervivencia: ", nrow(prot_surv_para_integracion),
    "proteínas con Gene Symbol\n")

# EXPORTACIÓN
# Tablas anotadas completas (con y sin Gene Symbol)
write_xlsx(
  list(
    Sexo_completa          = prot_sexo_anotada,
    Edad_completa          = prot_edad_anotada,
    Supervivencia_completa = prot_surv_anotada
  ),
  "proteinas_anotadas_completas.xlsx")

# Tablas filtradas para integración (solo con Gene Symbol)
write_xlsx(
  list(
    Sexo_integracion          = prot_sexo_para_integracion,
    Edad_integracion          = prot_edad_para_integracion,
    Supervivencia_integracion = prot_surv_para_integracion
  ),
  "proteinas_para_integracion.xlsx")

# Tabla de proteínas sin anotación (para revisión manual)
write_xlsx(
  list(
    Sexo_sin_anotacion          = sin_anotacion_sexo,
    Edad_sin_anotacion          = sin_anotacion_edad,
    Supervivencia_sin_anotacion = sin_anotacion_surv
  ),
  "proteinas_sin_Gene_Symbol.xlsx")

# CARGUE Y ESTANDARIZACIÓN DE LAS TABLAS TOP1.000 TRANSCRIPTÓMICAS
# CARGAR LAS TRES TABLAS
ruta_top1000 <- "C:/Users/danie/OneDrive/Documentos/Trabajo de Grado Maestría Bioinformática/Versiones finales/Tablas_TOP1000.xlsx"

top1000_sexo <- read_excel(ruta_top1000,  sheet = 1,
                           col_names = TRUE,
                           na = c("", "NA", "N/A"))
top1000_surv <- read_excel(ruta_top1000,  sheet = 2,
                           col_names = TRUE,
                           na = c("", "NA", "N/A"))
top1000_edad <- read_excel(ruta_top1000,  sheet = 3,
                           col_names = TRUE,
                           na = c("", "NA", "N/A"))

cat("Sexo:          ", nrow(top1000_sexo), "genes\n")
cat("Edad:          ", nrow(top1000_edad), "genes\n")
cat("Supervivencia: ", nrow(top1000_surv), "genes\n")

# VERIFICAR ESTRUCTURA DE COLUMNAS
cols_requeridas <- c("Gene_Symbol", "Score_total",
                     "Direccion_consenso", "logFC_ref")

verificar_columnas <- function(tabla, nombre) {
  cols_faltantes <- setdiff(cols_requeridas, colnames(tabla))
  if (length(cols_faltantes) > 0) {
    cat("⚠ ADVERTENCIA en", nombre, "— columnas faltantes:\n")
    cat("  ", paste(cols_faltantes, collapse = ", "), "\n")
    cat("  Columnas disponibles:", paste(colnames(tabla), collapse = ", "), "\n")
  } else {
    cat("✓", nombre, "— estructura correcta\n")}}


verificar_columnas(top1000_sexo, "TOP1000 Sexo")
verificar_columnas(top1000_edad, "TOP1000 Edad")
verificar_columnas(top1000_surv, "TOP1000 Supervivencia")

# ESTANDARIZACIÓN Y LIMPIEZA
# Garantiza que Gene_Symbol sea texto limpio y sin duplicados,
# que logFC_ref sea numérico, y que Direccion_consenso sea "Up" o "Down".

estandarizar_top1000 <- function(tabla, nombre_variable) {
  
  cat("\n── Estandarizando:", nombre_variable, "──\n")
  
  tabla_limpia <- tabla %>%
    
    # Limpiar Gene_Symbol: sin espacios, sin NAs
    filter(!is.na(Gene_Symbol) & Gene_Symbol != "") %>%
    mutate(
      Gene_Symbol = trimws(Gene_Symbol),  # eliminar espacios al inicio/fin
      
      # Asegurar que logFC_ref sea numérico
      logFC_ref = as.numeric(logFC_ref),
      
      # Estandarizar Direccion_consenso a "Up" / "Down" exactamente
      # (cubre variantes como "UP", "up", "Arriba", "+", etc.)
      Direccion_consenso = case_when(
        toupper(trimws(Direccion_consenso)) %in%
          c("UP", "ARRIBA", "SOBREEXPRESADO", "SOBRE", "+", "1") ~ "Up",
        toupper(trimws(Direccion_consenso)) %in%
          c("DOWN", "ABAJO", "SUBEXPRESADO", "SUB", "-", "-1")   ~ "Down",
        TRUE ~ Direccion_consenso  # conservar si ya está bien
      ),
      
      # Añadir columna de variable clínica para trazabilidad
      Variable_clinica = nombre_variable
    ) %>%
    
    # Eliminar duplicados de Gene_Symbol (conservar el de mayor score)
    arrange(desc(Score_total)) %>%
    distinct(Gene_Symbol, .keep_all = TRUE) %>%
    
    # Mantener solo las columnas relevantes para el Paso 5
    # (más cualquier columna extra que quieras conservar)
    select(
      Gene_Symbol,
      Score_total,
      Direccion_consenso,
      logFC_ref,
      Variable_clinica,
      any_of(c("Estudio_ref", "N_estudios_detectado",
               "N_estudios_misma_dir"))
    )
  
  # Diagnóstico
  cat("Genes tras limpieza:       ", nrow(tabla_limpia), "\n")
  cat("Con logFC_ref numérico:    ",
      sum(!is.na(tabla_limpia$logFC_ref)), "\n")
  cat("Dirección Up:              ",
      sum(tabla_limpia$Direccion_consenso == "Up",   na.rm = TRUE), "\n")
  cat("Dirección Down:            ",
      sum(tabla_limpia$Direccion_consenso == "Down", na.rm = TRUE), "\n")
  
  # Advertencia si hay menos de 1000 genes tras la limpieza
  if (nrow(tabla_limpia) < 1000) {
    cat("⚠ Nota: quedan", nrow(tabla_limpia),
        "genes (< 1000 tras eliminar duplicados/NAs)\n")
  }
  
  return(tabla_limpia)}

top1000_sexo_std <- estandarizar_top1000(top1000_sexo, "Sexo")
top1000_edad_std <- estandarizar_top1000(top1000_edad, "Edad")
top1000_surv_std <- estandarizar_top1000(top1000_surv, "Supervivencia")

# SOLAPAMIENTO ENTRE ETIQUETAS/VARIABLES (análisis exploratorio)
# Cuántos genes del TOP1000 son compartidos entre variables clínicas.
# etiquetas tienen evidencia más amplia de desregulación en OSA.

genes_sexo <- top1000_sexo_std$Gene_Symbol
genes_edad <- top1000_edad_std$Gene_Symbol
genes_surv <- top1000_surv_std$Gene_Symbol

cat("Sexo ∩ Edad:              ",
    length(intersect(genes_sexo, genes_edad)), "genes\n")
cat("Sexo ∩ Supervivencia:     ",
    length(intersect(genes_sexo, genes_surv)), "genes\n")
cat("Edad ∩ Supervivencia:     ",
    length(intersect(genes_edad, genes_surv)), "genes\n")
cat("Sexo ∩ Edad ∩ Supervivencia:",
    length(Reduce(intersect, list(genes_sexo, genes_edad, genes_surv))),
    "genes\n")

# Genes en las tres etiquetas simultáneamente (los más robustos)
genes_triple <- Reduce(intersect, list(genes_sexo, genes_edad, genes_surv))
if (length(genes_triple) > 0) {
  cat("\nGenes en TOP1000 de las TRES etiquetas:\n")
  print(genes_triple)}

# Exportar tablas estandarizadas para revisión
write_xlsx(
  list(
    TOP1000_Sexo          = top1000_sexo_std,
    TOP1000_Edad          = top1000_edad_std,
    TOP1000_Supervivencia = top1000_surv_std
  ),
  "TOP1000_estandarizados.xlsx")

# Organizar en lista para intersección y análisis de concordancia direccional
top1000_std <- list(
  Sexo          = top1000_sexo_std,
  Edad          = top1000_edad_std,
  Supervivencia = top1000_surv_std)

cat("  top1000_std$Sexo          —", nrow(top1000_std$Sexo), "genes\n")
cat("  top1000_std$Edad          —", nrow(top1000_std$Edad), "genes\n")
cat("  top1000_std$Supervivencia —", nrow(top1000_std$Supervivencia), "genes\n")

# INTERSECCIÓN Y ANÁLISIS DE CONCORDANCIA DIRECCIONAL
# VERIFICACIONES PREVIAS
stopifnot(
  "top1000_std no existe" = exists("top1000_std"),
  "top1000_std debe contener Sexo, Edad y Supervivencia" =
    all(c("Sexo", "Edad", "Supervivencia") %in% names(top1000_std))
)

for (obj in c(
  "prot_sexo_para_integracion",
  "prot_edad_para_integracion",
  "prot_surv_para_integracion"
)) {
  
  if (!exists(obj)) {
    stop(
      paste(
        obj,
        "no existe. Ejecuta primero la sección de anotación proteómica."
      )
    )
  }
}


# TABLAS PROTEÓMICA
tablas_prot <- list(
  Sexo          = prot_sexo_para_integracion,
  Edad          = prot_edad_para_integracion,
  Supervivencia = prot_surv_para_integracion
)

# UNIVERSO PROTEÓMICO
# Corresponde a todos los genes con Gene Symbol identificados en la proteómica.
genes_universo_prot <- unique(
  unlist(
    lapply(
      tablas_prot,
      function(x) x$Gene_Symbol
    ) ))

genes_universo_prot <- genes_universo_prot[
  !is.na(genes_universo_prot) &
    genes_universo_prot != ""]

cat("\nGenes únicos identificados en proteómica:",
    length(genes_universo_prot), "\n")

# UNIVERSO TRANSCRIPTÓMICO
genes_universo_trans <- unique(
  c(
    top1000_sexo_std$Gene_Symbol,
    top1000_edad_std$Gene_Symbol,
    top1000_surv_std$Gene_Symbol))

if (!exists("genes_universo_trans")) {
  
  stop(
    paste0(
      "\nNO SE ENCONTRÓ 'genes_universo_trans'.\n\n",
      "Para realizar correctamente la prueba hipergeométrica ",
      "debes proporcionar el conjunto de genes evaluados en ",
      "la transcriptómica.\n\n",
      "IMPORTANTE: NO uses solamente los TOP1000 como universo."
    ))}


genes_universo_trans <- unique(
  genes_universo_trans)

genes_universo_trans <- genes_universo_trans[
  !is.na(genes_universo_trans) &
    genes_universo_trans != ""]

cat("Genes únicos evaluables en transcriptómica:",
  length(genes_universo_trans),"\n")

# UNIVERSO COMÚN
# Solo se consideran genes que podían ser observados en ambas capas ómicas.
universo_comun <- intersect(
  genes_universo_trans,
  genes_universo_prot)

universo_comun <- unique(universo_comun)

cat("Genes del universo común transcriptómica-proteómica:",
  length(universo_comun),"\n")


if (length(universo_comun) == 0) {
  
  stop(
    "El universo común está vacío. ",
    "Revisa que los Gene Symbol de transcriptómica y proteómica ",
    "estén en el mismo formato."
  )
}

# FUNCIÓN CENTRAL DE INTEGRACIÓN
integrar_variable <- function(
    top1000_var,
    prot_var,
    nombre_variable,
    universo_trans,
    universo_prot) {
  
  cat("\n")
  cat("============================================================\n")
  cat("INTEGRACIÓN:", nombre_variable, "\n")
  cat("============================================================\n")
  
  
  # ----------------------------------------------------------
  # 6.1 GENES DEL TOP1000
  # ----------------------------------------------------------
  
  genes_trans_top1000 <- unique(
    top1000_var$Gene_Symbol
  )
  
  genes_trans_top1000 <- genes_trans_top1000[
    !is.na(genes_trans_top1000) &
      genes_trans_top1000 != ""
  ]
  
  
  # ----------------------------------------------------------
  # 6.2 GENES PROTEÓMICOS
  # ----------------------------------------------------------
  
  genes_prot <- unique(
    prot_var$Gene_Symbol
  )
  
  genes_prot <- genes_prot[
    !is.na(genes_prot) &
      genes_prot != ""
  ]
  
  
  # ----------------------------------------------------------
  # 6.3 UNIVERSO COMÚN
  # ----------------------------------------------------------
  
  universo_comun <- intersect(
    universo_trans,
    universo_prot
  )
  
  universo_comun <- unique(
    universo_comun
  )
  
  
  # ----------------------------------------------------------
  # 6.4 TOP1000 EVALUABLE EN EL UNIVERSO
  # ----------------------------------------------------------
  
  genes_trans_evaluables <- intersect(
    genes_trans_top1000,
    universo_comun
  )
  
  
  # ----------------------------------------------------------
  # 6.5 GENES PROTEÓMICOS EVALUABLES
  # ----------------------------------------------------------
  
  genes_prot_evaluables <- intersect(
    genes_prot,
    universo_comun
  )
  
  
  # ----------------------------------------------------------
  # 6.6 INTERSECCIÓN OBSERVADA
  # ----------------------------------------------------------
  
  genes_comunes <- intersect(
    genes_trans_evaluables,
    genes_prot_evaluables
  )
  
  
  # ----------------------------------------------------------
  # 6.7 PARÁMETROS HIPERGEOMÉTRICOS
  # ----------------------------------------------------------
  
  N <- length(universo_comun)
  
  K <- length(genes_prot_evaluables)
  
  n <- length(genes_trans_evaluables)
  
  k <- length(genes_comunes)
  
  
  # ----------------------------------------------------------
  # 6.8 NÚMERO ESPERADO
  # ----------------------------------------------------------
  
  esperado <- ifelse(
    N > 0,
    (n * K) / N,
    NA_real_
  )
  
  
  # ----------------------------------------------------------
  # 6.9 PRUEBA HIPERGEOMÉTRICA
  # ----------------------------------------------------------
  #
  # Probabilidad de observar k o más genes compartidos.
  #
  # phyper() utiliza:
  #
  # q = k - 1
  # m = K
  # n = N-K
  # k = n
  #
  # lower.tail = FALSE
  #
  # ----------------------------------------------------------
  
  p_hipergeometrica <- phyper(
    q = k - 1,
    m = K,
    n = N - K,
    k = n,
    lower.tail = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 6.10 FOLD ENRICHMENT
  # ----------------------------------------------------------
  
  fold_enrichment <- ifelse(
    esperado > 0,
    k / esperado,
    NA_real_
  )
  
  
  # ----------------------------------------------------------
  # 6.11 FDR DE LAS TRES COMPARACIONES
  # ----------------------------------------------------------
  
  # La corrección se realizará posteriormente sobre los
  # tres p-values obtenidos.
  
  
  # ----------------------------------------------------------
  # 6.12 DIAGNÓSTICO
  # ----------------------------------------------------------
  
  cat("\n--- Tamaño de conjuntos ---\n")
  
  cat(
    "TOP1000 original:              ",
    length(genes_trans_top1000),
    "\n"
  )
  
  cat(
    "TOP1000 evaluable:             ",
    n,
    "\n"
  )
  
  cat(
    "Proteómica identificada:       ",
    length(genes_prot),
    "\n"
  )
  
  cat(
    "Proteómica evaluable:          ",
    K,
    "\n"
  )
  
  cat(
    "Universo común:                ",
    N,
    "\n"
  )
  
  cat(
    "Genes compartidos observados:  ",
    k,
    "\n"
  )
  
  cat(
    "Genes compartidos esperados:   ",
    round(esperado, 2),
    "\n"
  )
  
  cat(
    "Fold enrichment:               ",
    round(fold_enrichment, 3),
    "\n"
  )
  
  cat(
    "P hipergeométrica:             ",
    format.pval(
      p_hipergeometrica,
      digits = 4
    ),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 6.13 CREAR TABLA PROTEÓMICA
  # ----------------------------------------------------------
  
  cols_opcionales <- intersect(
    c(
      "Ensembl_ID_final",
      "Ensembl_ID",
      "Descripcion_final",
      "Descripcion",
      "Fuente_Gene"
    ),
    colnames(prot_var)
  )
  
  
  cols_fijas <- c(
    "Accesion",
    "Gene_Symbol",
    "logFC",
    "adj.P.Val",
    "P.Value",
    "Regulacion"
  )
  
  
  faltantes <- setdiff(
    cols_fijas,
    colnames(prot_var)
  )
  
  
  if (length(faltantes) > 0) {
    
    stop(
      paste(
        "Faltan columnas en la tabla proteómica:",
        paste(
          faltantes,
          collapse = ", "
        )
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # 6.14 TABLA INTEGRADA
  # ----------------------------------------------------------
  
  tabla_integrada <- prot_var %>%
    
    filter(
      Gene_Symbol %in% genes_comunes
    ) %>%
    
    select(
      all_of(
        c(
          cols_fijas[1:2],
          cols_opcionales
        )
      ),
      logFC_prot      = logFC,
      FDR_prot        = adj.P.Val,
      Pval_prot       = P.Value,
      Regulacion_prot = Regulacion
    ) %>%
    
    left_join(
      
      top1000_var %>%
        
        filter(
          Gene_Symbol %in% genes_comunes
        ) %>%
        
        select(
          Gene_Symbol,
          Score_total,
          logFC_trans     = logFC_ref,
          Direccion_trans = Direccion_consenso,
          any_of(
            c(
              "Estudio_ref",
              "N_estudios_detectado",
              "N_estudios_misma_dir"
            )
          )
        ),
      
      by = "Gene_Symbol"
    ) %>%
    
    mutate(
      Direccion_prot = case_when(logFC_prot > 0 ~ "Up", logFC_prot < 0 ~ "Down", 
                                      TRUE ~ "Indeterminado"), 
           Concordancia = case_when(Direccion_prot == "Up" & Direccion_trans == "Up" ~ "Concordante Up", 
                                    Direccion_prot == "Down" & Direccion_trans == "Down" ~ "Concordante Down", 
                                    Direccion_prot == "Up" & Direccion_trans == "Down" ~ "Discordante Prot.Up/Trans.Down", 
                                    Direccion_prot == "Down" & Direccion_trans == "Up" ~ "Discordante Prot.Down/Trans.Up", 
                                    TRUE ~ "Indeterminado"), 
      Sig_FDR = !is.na(FDR_prot) & FDR_prot < 0.05 & abs(logFC_prot) > 1, 
      Sig_nominal = !is.na(Pval_prot) & Pval_prot < 0.05 & abs(logFC_prot) > 1, 
      score_integración = Score_total + as.integer(grepl("^Concordante", Concordancia)) + as.integer(Sig_FDR))
      
  
  # ----------------------------------------------------------
  # 6.15 RESUMEN DE CONCORDANCIA
  # ----------------------------------------------------------
  
  cat("\n--- Concordancia direccional ---\n")
  
  print(
    table(
      tabla_integrada$Concordancia
    )
  )
  
  
  cat(
    "\nCandidatos concordantes:",
    sum(
      grepl(
        "Concordante",
        tabla_integrada$Concordancia
      )
    ),
    "\n"
  )
  
  
  cat(
    "Candidatos significativos FDR:",
    sum(
      tabla_integrada$Sig_FDR,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 6.16 INFORMACIÓN ESTADÍSTICA DE LA PRUEBA
  # ----------------------------------------------------------
  
  estadistica_hipergeometrica <- data.frame(
    
    Variable = nombre_variable,
    
    Universo_N = N,
    
    Proteomica_K = K,
    
    TOP1000_n = n,
    
    Interseccion_k = k,
    
    Esperado = esperado,
    
    Fold_enrichment = fold_enrichment,
    
    P_hipergeometrica = p_hipergeometrica,
    
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 6.17 RESULTADO FINAL
  # ----------------------------------------------------------
  
  return(
    list(
      
      tabla = tabla_integrada,
      
      genes_top1000 = genes_trans_evaluables,
      
      genes_proteomica = genes_prot_evaluables,
      
      genes_comunes = genes_comunes,
      
      estadistica = estadistica_hipergeometrica
    )
  )
}

# EJECUTAR LAS TRES INTEGRACIONES

resultados_integracion <- list()

for (variable in c(
  "Sexo",
  "Edad",
  "Supervivencia"
)) {
  
  resultados_integracion[[variable]] <-
    integrar_variable(
      
      top1000_var =
        top1000_std[[variable]],
      
      prot_var =
        tablas_prot[[variable]],
      
      nombre_variable =
        variable,
      
      universo_trans =
        genes_universo_trans,
      
      universo_prot =
        genes_universo_prot
    )
}

# EXTRAER TABLAS
tablas_integracion <- lapply(
  resultados_integracion,
  function(x) x$tabla)

# RESUMEN HIPERGEOMÉTRICO
resumen_hipergeometrico <- bind_rows(
  lapply(
    resultados_integracion,
    function(x) x$estadistica))

# CORRECCIÓN POR MÚLTIPLES PRUEBAS
resumen_hipergeometrico <- resumen_hipergeometrico %>%
  mutate(
    FDR_hipergeometrica =
      p.adjust(
        P_hipergeometrica,
        method = "BH"
      ),
    
    Enriquecimiento_significativo =
      FDR_hipergeometrica < 0.05
  )

print(resumen_hipergeometrico)

# CANDIDATOS MULTIVARIABLE
genes_int <- lapply(
  resultados_integracion,
  function(x) {
    
    if (is.null(x)) {
      return(character(0))
    }
    
    unique(
      x$tabla$Gene_Symbol)})

cat(
  "Sexo ∩ Edad:               ",
  length(
    intersect(
      genes_int$Sexo,
      genes_int$Edad
    )
  ),
  "\n"
)

cat(
  "Sexo ∩ Supervivencia:      ",
  length(
    intersect(
      genes_int$Sexo,
      genes_int$Supervivencia
    )
  ),
  "\n"
)

cat(
  "Edad ∩ Supervivencia:      ",
  length(
    intersect(
      genes_int$Edad,
      genes_int$Supervivencia
    )
  ),
  "\n"
)

cat(
  "Sexo ∩ Edad ∩ Supervivencia:",
  length(
    Reduce(
      intersect,
      genes_int
    )
  ),
  "\n"
)

# TABLA DE CANDIDATOS EN ≥2 VARIABLES
todos_candidatos <- bind_rows(
  
  lapply(
    names(tablas_integracion),
    function(v) {
      
      df <- tablas_integracion[[v]]
      
      if (is.null(df) ||
          nrow(df) == 0) {
        return(NULL)
      }
      
      df %>%
        
        mutate(
          Variable = v
        ) %>%
        
        select(
          Gene_Symbol,
          Variable,
          Concordancia,
          logFC_prot,
          logFC_trans,
          FDR_prot,
          Score_total,
          Sig_FDR,
          score_integración
        )
    }
  )
)


candidatos_multivariable <- todos_candidatos %>%
  
  group_by(
    Gene_Symbol
  ) %>%
  
  summarise(
    
    N_variables =
      n_distinct(Variable),
    
    Variables =
      paste(
        unique(Variable),
        collapse = " | "
      ),
    
    N_concordantes =
      sum(
        grepl(
          "Concordante",
          Concordancia
        )
      ),
    
    N_sig_FDR =
      sum(
        Sig_FDR,
        na.rm = TRUE
      ),
    
    logFC_prot_promedio =
      round(
        mean(
          logFC_prot,
          na.rm = TRUE
        ),
        3
      ),
    
    Score_integracion_max = 
      max(score_integración, 
          na.rm = TRUE),
    
    Score_trans_promedio =
      round(
        mean(
          Score_total,
          na.rm = TRUE
        ),
        2
      ),
    
    .groups = "drop"
  ) %>%
  
  filter(
    N_variables >= 2
  ) %>%
  
  arrange(
    desc(N_variables),
    desc(Score_integracion_max),
    desc(N_concordantes),
    desc(N_sig_FDR),
    
  )

# EXPORTACIÓN
write_xlsx(c(
    lapply(
      tablas_integracion,
      identity
    ),
    
    list(
      Resumen_hipergeometrico =
        resumen_hipergeometrico,
      
      Candidatos_multivariable =
        candidatos_multivariable
    )
  ),
  
  "candidatos_multiomica_integrados.xlsx")


todos_candidatos_concordantes <-bind_rows(lapply(
      names(tablas_integracion),
      function(v) {
        df <-
          tablas_integracion[[v]]
        
        if (
          is.null(df) ||
          nrow(df) == 0
        ) {
          return(NULL)
        }
        
        df %>%
          
          filter(
            grepl(
              "Concordante",
              Concordancia
            )
          ) %>%
          
          mutate(
            Variable_clinica = v
          )
      }
    )
  )


write_xlsx(todos_candidatos_concordantes,
   "candidatos_concordantes_todos.xlsx")
