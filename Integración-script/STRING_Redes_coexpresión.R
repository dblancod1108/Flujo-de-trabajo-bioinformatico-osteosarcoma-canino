# =============================================================================
#  REDES DE COEXPRESIÓN 
# =============================================================================

#Librerias
library(readr)
library(dplyr)
library(readxl)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(writexl)

# RED DE COEXPRESIÓN - CANDIDATOS CONCORDANTES SEXO
# CARGAR DATOS
ruta_excel <- "C:/Users/danie/OneDrive/Documentos/Trabajo de Grado Maestría Bioinformática/Versiones finales/Script integración/candidatos_multiomica_integrados.xlsx"

candidatos <- read_excel(ruta_excel, sheet = "Sexo") %>%
  select(Gene_Symbol, score_integración, logFC_prot,
         Concordancia, FDR_prot) %>%
  filter(!is.na(Gene_Symbol) & !is.na(score_integración)) %>%
  filter(grepl("Concordante", Concordancia)) %>%
  group_by(Gene_Symbol) %>%
  slice_max(order_by = abs(logFC_prot), n = 1, with_ties = FALSE) %>%
  ungroup()

genes_semilla <- candidatos$Gene_Symbol

cat("Total:", length(genes_semilla), "\n")
print(genes_semilla)

# Archivo de edges de STRING
# Archivo TSV descargado de STRING ("as tabular text output" — edges recíprocos A-B y B-A)
ruta_string <- "C:/Users/danie/OneDrive/Documentos/Trabajo de Grado Maestría Bioinformática/Versiones finales/Script integración/string_network_Sexo.tsv"
edges_raw <- read_tsv(ruta_string, show_col_types = FALSE)

cat("Filas:", nrow(edges_raw), "\n")
cat("Columnas:", paste(colnames(edges_raw), collapse = ", "), "\n")
print(head(edges_raw, 3))

# LIMPIAR EDGES DE STRING
# Identificar columnas automáticamente
col_n1 <- intersect(c("node1", "protein1", "#node1"),
                    colnames(edges_raw))[1]
col_n2 <- intersect(c("node2", "protein2"),
                    colnames(edges_raw))[1]
col_sc <- intersect(c("combined_score", "score"),
                    colnames(edges_raw))[1]

cat("\nColumnas identificadas — Nodo1:", col_n1,
    "| Nodo2:", col_n2, "| Score:", col_sc, "\n")

# Función para eliminar prefijo taxonómico generado por defecto en STRING (ej: "9615.GAPDH" → "GAPDH")
limpiar_id <- function(x) gsub("^[0-9]+\\.", "", x)

# Limpiar y normalizar
edges_raw2 <- edges_raw %>%
  rename(nodo1 = all_of(col_n1),
         nodo2 = all_of(col_n2),
         score = all_of(col_sc)) %>%
  mutate(nodo1 = limpiar_id(nodo1),
         nodo2 = limpiar_id(nodo2))

score_max <- max(edges_raw2$score, na.rm = TRUE)
edges_raw2 <- edges_raw2 %>%
  mutate(score_norm = ifelse(score_max > 1, score / 1000, score))

cat("Score máximo:", score_max,
    "→", ifelse(score_max > 1, "normalizado ÷1000", "ya en 0-1"), "\n")

# Filtrar: score ≥ 0.40 (medium confidence)
# Se uso 0.40 para capturar más conexiones de coexpresión
umbral <- 0.40

edges_clean <- edges_raw2 %>%
  filter(score_norm >= umbral) %>%
  filter(nodo1 != nodo2) %>%
  # Eliminar duplicados (A-B y B-A son el mismo edge no dirigido)
  mutate(par_min = pmin(nodo1, nodo2),
         par_max = pmax(nodo1, nodo2)) %>%
  distinct(par_min, par_max, .keep_all = TRUE) %>%
  select(nodo1, nodo2, score_norm)

cat("\nEdges tras filtro (score ≥", umbral, "):", nrow(edges_clean), "\n")

# Si no hay edges, sugerir bajar el umbral
if (nrow(edges_clean) == 0) {
  cat("\n⚠ Sin edges con score ≥", umbral, "\n")
  cat("Scores disponibles en el archivo:\n")
  print(summary(edges_raw2$score_norm))
  cat("→ Prueba con umbral = 0.15 o revisa que el archivo TSV\n")
  cat("  corresponde a la red de Sexo descargada de STRING.\n")
  stop("Sin edges para construir la red.")
}

cat("Rango de scores:",
    round(min(edges_clean$score_norm), 3), "—",
    round(max(edges_clean$score_norm), 3), "\n")

# CONSTRUCCIÓN DEL GRAFO
# Tabla de nodos: todos los genes semilla
# (incluyendo los que no tienen edges — quedarán aislados)
nodos_df <- data.frame(name = genes_semilla,
                       stringsAsFactors = FALSE) %>%
  left_join(candidatos %>%
              select(name = Gene_Symbol, score_integración,
                     logFC_prot, Concordancia, FDR_prot),
            by = "name") %>%
  mutate(
    Score_norm = (score_integración - min(score_integración, na.rm = TRUE)) /
      (max(score_integración, na.rm = TRUE) -
         min(score_integración, na.rm = TRUE)),
    Sig_FDR    = FDR_prot < 0.05 & abs(logFC_prot) > 1,
    # Clasificar conectividad para la visualización
    Conectado  = name %in% c(edges_clean$nodo1, edges_clean$nodo2)
  )

# Construcción del grafo con todos los nodos (incluso aislados)
red <- graph_from_data_frame(
  d        = edges_clean,
  directed = FALSE,
  vertices = nodos_df
)

V(red)$Grado <- degree(red)

cat("Nodos totales:             ", vcount(red), "\n")
cat("Nodos con conexiones:      ",
    sum(degree(red) > 0), "\n")
cat("Nodos aislados (grado = 0):",
    sum(degree(red) == 0), "\n")
cat("Edges:                     ", ecount(red), "\n")
cat("Densidad:                  ",
    round(edge_density(red), 3), "\n")

# Mostrar qué genes quedaron aislados
aislados <- V(red)$name[degree(red) == 0]
if (length(aislados) > 0) {
  cat("\nGenes sin conexión en STRING (score ≥", umbral, "):\n")
  print(aislados)
  cat("→ No tienen evidencia de coexpresión con los otros genes\n")
  cat("  de la lista a este umbral de confianza.\n")
}

# MÉTRICAS DE CENTRALIDAD
metricas <- data.frame(
  Gene        = V(red)$name,
  Grado       = degree(red),
  Betweenness = round(betweenness(red, normalized = TRUE), 4),
  Closeness   = round(closeness(red, normalized = TRUE), 4),
  Eigenvector = round(eigen_centrality(red)$vector, 4)
) %>%
  left_join(candidatos %>%
              select(Gene = Gene_Symbol, Concordancia,
                     logFC_prot, score_integración, FDR_prot),
            by = "Gene") %>%
  arrange(desc(Grado))

print(metricas %>%
        select(Gene, Grado, Betweenness, Concordancia,
               logFC_prot, score_integración))

# VISUALIZACIÓN CON ggplot2 
# Calcular coordenadas circulares para TODOS los nodos
n       <- vcount(red)
angulos <- seq(0, 2 * pi, length.out = n + 1)[1:n]

coords_nodos <- data.frame(
  name  = V(red)$name,
  grado = V(red)$Grado,
  x     = cos(angulos),
  y     = sin(angulos),
  stringsAsFactors = FALSE
) %>%
  left_join(nodos_df, by = "name")

# Coordenadas de edges para geom_segment
coords_edges <- edges_clean %>%
  left_join(coords_nodos %>% select(name, x, y),
            by = c("nodo1" = "name")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(coords_nodos %>% select(name, x, y),
            by = c("nodo2" = "name")) %>%
  rename(x2 = x, y2 = y)

# Colores según dirección de concordancia
colores_dir <- c(
  "Concordante Up"   = "#E8527A",   
  "Concordante Down" = "#2980B9"    
)

# Forma del nodo: conectado vs aislado
formas_nodo <- c(
  "Conectado" = 17,   # triángulo relleno
  "Aislado"   = 1     # círculo vacío
)

coords_nodos <- coords_nodos %>%
  mutate(Estado = ifelse(Conectado, "Conectado", "Aislado"))

p_red <- ggplot() +
  
  #Edges: líneas entre nodos conectados 
  geom_segment(
    data = coords_edges,
    aes(x     = x1, y = y1,
        xend  = x2, yend = y2,
        alpha = score_norm,
        linewidth = score_norm),
    color = "grey50"
  ) +
  scale_alpha_continuous(
    range = c(0.3, 0.9),
    name  = "Score STRING"
  ) +
  scale_linewidth_continuous(
    range = c(0.4, 2.0),
    name  = "Score STRING"
  ) +
  
  # Nodos 
  geom_point(
    data  = coords_nodos,
    aes(x     = x,
        y     = y,
        color = Concordancia,
        size  = Score_norm,
        shape = Estado),
    alpha = 0.92
  ) +
  scale_color_manual(
    values   = colores_dir,
    na.value = "grey60",
    name     = "Dirección\nconcordante"
  ) +
  scale_size_continuous(
    range = c(4, 13),
    name  = "Score\nintegración"
  ) +
  scale_shape_manual(
    values = formas_nodo,
    name   = "Conectividad"
  ) +
  
  # Etiquetas: todas las proteínas 
  geom_text(
    data = coords_nodos,
    aes(x     = x * 1.22,
        y     = y * 1.22,
        label = name,
        hjust = ifelse(x >  0.1, 0,
                       ifelse(x < -0.1, 1, 0.5)),
        vjust = ifelse(y >  0.1, 0,
                       ifelse(y < -0.1, 1, 0.5)),
        # Negrita para los conectados, normal para aislados
        fontface = ifelse(Estado == "Conectado", "bold", "plain")
    ),
    size  = 3.2,
    color = "grey15"
  ) +
  
  coord_fixed(xlim = c(-1.6, 1.6),
              ylim = c(-1.6, 1.6)) +
  
  theme_void() +
  theme(
    plot.title      = element_text(face  = "bold", size = 13,
                                   hjust = 0.5,
                                   margin = margin(b = 6)),
    plot.subtitle   = element_text(size  = 9, color = "grey35",
                                   hjust = 0.5,
                                   margin = margin(b = 4)),
    plot.caption    = element_text(size  = 7.5, color = "grey55",
                                   hjust = 1,
                                   margin = margin(t = 6)),
    legend.position = "right",
    legend.title    = element_text(size = 8.5, face = "bold"),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(20, 20, 20, 20)
  ) +
  labs(
    title    = "Red de coexpresión STRING — Candidatos concordantes Sexo",
    subtitle = paste0(
      vcount(red), " genes | ",
      ecount(red), " conexiones (score STRING ≥ ", umbral, ") | ",
      sum(degree(red) > 0), " conectados — ",
      sum(degree(red) == 0), " aislados"
    ),
    caption = paste0(
      "▲ Conectado  ○ Aislado (sin evidencia de coexpresión con otros genes de la lista)\n",
      "Rojo = Up en Hembra | Azul = Down en Hembra | ",
      "Tamaño proporcional al Score de integración multiómica\n",
      "Fuente de coexpresión: STRING v12.0 | Canis lupus familiaris (taxón 9615)"
    )
  )

# GUARDAR GRAFO EN FORMATO png
ggsave("red_coexpresion_STRING_Sexo.png",
       plot   = p_red,
       width  = 12,
       height = 10,
       dpi    = 300,
       bg     = "white")

ggsave("red_coexpresion_STRING_Sexo.pdf",
       plot   = p_red,
       width  = 12,
       height = 10,
       device = cairo_pdf)

#  RED DE COEXPRESIÓN STRING — CANDIDATOS CONCORDANTES SUPERVIVENCIA
# CARGAR DATOS
ruta_string_sup <- "C:/Users/danie/OneDrive/Documentos/Trabajo de Grado Maestría Bioinformática/Versiones finales/Script integración/string_network_Supervivencia.tsv" 

# Candidatos concordantes de Supervivencia
candidatos_supervivencia <- read_excel(ruta_excel, sheet = "Supervivencia") %>%
  select(Gene_Symbol, score_integración, logFC_prot,
         Concordancia, FDR_prot) %>%
  filter(!is.na(Gene_Symbol) & !is.na(score_integración)) %>%
  filter(grepl("Concordante", Concordancia)) %>%
  group_by(Gene_Symbol) %>%
  slice_max(order_by = abs(logFC_prot), n = 1, with_ties = FALSE) %>%
  ungroup()

genes_semilla_sup <- candidatos_supervivencia$Gene_Symbol

cat("Total:", length(genes_semilla_sup), "\n")
print(candidatos_supervivencia %>%
        select(Gene_Symbol, Concordancia,
               logFC_prot, score_integración) %>%
        arrange(desc(score_integración)))

if (length(genes_semilla_sup) < 2) {
  stop("Menos de 2 genes concordantes — no es posible construir la red.")
}

# Archivo TSV de STRING 
edges_raw_sup <- read_tsv(ruta_string_sup, show_col_types = FALSE)

cat("Filas:", nrow(edges_raw_sup), "\n")
cat("Columnas:", paste(colnames(edges_raw_sup), collapse = ", "), "\n")
print(head(edges_raw_sup, 3))

# Verificar que el archivo tiene la estructura correcta
# Debe tener columnas tipo: node1/protein1, node2/protein2, combined_score
# Si ves: #node, identifier, x_position → es el archivo equivocado
if (any(c("x_position", "y_position") %in% colnames(edges_raw_sup))) {
  stop(paste(
    "⚠ Archivo incorrecto: cargaste 'network coordinates' en vez de",
    "'tabular text output'.\n",
    "En STRING → Exports → selecciona 'as tabular text output'\n",
    "(el que dice 'lists reciprocal edges: A-B, B-A')"
  ))
}

# Identificar columnas automáticamente
col_n1_sup <- intersect(c("node1", "protein1", "#node1"),
                        colnames(edges_raw_sup))[1]
col_n2_sup <- intersect(c("node2", "protein2"),
                        colnames(edges_raw_sup))[1]
col_sc_sup <- intersect(c("combined_score", "score"),
                        colnames(edges_raw_sup))[1]

cat("\nColumnas — Nodo1:", col_n1_sup,
    "| Nodo2:", col_n2_sup,
    "| Score:", col_sc_sup, "\n")

# Limpiar prefijo taxonómico
limpiar_id <- function(x) gsub("^[0-9]+\\.", "", x)

edges_raw2_sup <- edges_raw_sup %>%
  rename(nodo1 = all_of(col_n1_sup),
         nodo2 = all_of(col_n2_sup),
         score = all_of(col_sc_sup)) %>%
  mutate(nodo1 = limpiar_id(nodo1),
         nodo2 = limpiar_id(nodo2))

score_max_sup  <- max(edges_raw2_sup$score, na.rm = TRUE)
edges_raw2_sup <- edges_raw2_sup %>%
  mutate(score_norm = ifelse(score_max_sup > 1,
                             score / 1000,
                             score))

cat("Score máximo:", score_max_sup, "→",
    ifelse(score_max_sup > 1,
           "normalizado ÷1000",
           "ya en escala 0-1"), "\n")
cat("Resumen de scores:\n")
print(summary(edges_raw2_sup$score_norm))

# Filtrar por umbral
umbral_sup <- 0.40

edges_clean_sup <- edges_raw2_sup %>%
  filter(score_norm >= umbral_sup) %>%
  filter(nodo1 != nodo2) %>%
  mutate(par_min = pmin(nodo1, nodo2),
         par_max = pmax(nodo1, nodo2)) %>%
  distinct(par_min, par_max, .keep_all = TRUE) %>%
  select(nodo1, nodo2, score_norm)

if (nrow(edges_clean_sup) == 0) {
  cat("⚠ Sin edges con score ≥", umbral_sup,
      "— bajando umbral a 0.15\n")
  umbral_sup <- 0.15
  edges_clean_sup <- edges_raw2_sup %>%
    filter(score_norm >= umbral_sup) %>%
    filter(nodo1 != nodo2) %>%
    mutate(par_min = pmin(nodo1, nodo2),
           par_max = pmax(nodo1, nodo2)) %>%
    distinct(par_min, par_max, .keep_all = TRUE) %>%
    select(nodo1, nodo2, score_norm)
}

cat("\nEdges tras filtro (score ≥", umbral_sup, "):",
    nrow(edges_clean_sup), "\n")

# CONSTRUIR EL GRAFO
nodos_df_sup <- data.frame(name = genes_semilla_sup,
                           stringsAsFactors = FALSE) %>%
  left_join(candidatos_supervivencia %>%
              select(name = Gene_Symbol,
                     score_integración, logFC_prot,
                     Concordancia, FDR_prot),
            by = "name") %>%
  mutate(
    Score_norm = (score_integración -
                    min(score_integración, na.rm = TRUE)) /
      (max(score_integración, na.rm = TRUE) -
         min(score_integración, na.rm = TRUE)),
    Sig_FDR   = FDR_prot < 0.05 & abs(logFC_prot) > 1,
    Conectado = name %in% c(edges_clean_sup$nodo1,
                            edges_clean_sup$nodo2)
  )

red_sup <- graph_from_data_frame(
  d        = if (nrow(edges_clean_sup) > 0) edges_clean_sup else
    data.frame(nodo1      = character(),
               nodo2      = character(),
               score_norm = numeric()),
  directed = FALSE,
  vertices = nodos_df_sup
)

V(red_sup)$Grado <- degree(red_sup)

cat("Nodos totales:             ", vcount(red_sup), "\n")
cat("Nodos con conexiones:      ", sum(degree(red_sup) > 0), "\n")
cat("Nodos aislados (grado = 0):", sum(degree(red_sup) == 0), "\n")
cat("Edges:                     ", ecount(red_sup), "\n")
if (ecount(red_sup) > 0) {
  cat("Densidad:                  ",
      round(edge_density(red_sup), 3), "\n")}

aislados_sup <- V(red_sup)$name[degree(red_sup) == 0]
if (length(aislados_sup) > 0) {
  cat("\nGenes sin conexión (score ≥", umbral_sup, "):\n")
  print(aislados_sup)}

# MÉTRICAS DE CENTRALIDAD 
metricas_sup <- data.frame(
  Gene        = V(red_sup)$name,
  Grado       = degree(red_sup),
  Betweenness = round(betweenness(red_sup, normalized = TRUE), 4),
  Closeness   = round(closeness(red_sup, normalized = TRUE), 4),
  Eigenvector = round(eigen_centrality(red_sup)$vector, 4)
) %>%
  left_join(candidatos_supervivencia %>%
              select(Gene = Gene_Symbol, Concordancia,
                     logFC_prot, score_integración, FDR_prot),
            by = "Gene") %>%
  arrange(desc(Grado))

print(metricas_sup %>%
        select(Gene, Grado, Concordancia,
               logFC_prot, score_integración))

# VISUALIZACIÓN CON ggplot2 
n_sup       <- vcount(red_sup)
angulos_sup <- seq(0, 2 * pi, length.out = n_sup + 1)[1:n_sup]

coords_nodos_sup <- data.frame(
  name  = V(red_sup)$name,
  grado = V(red_sup)$Grado,
  x     = cos(angulos_sup),
  y     = sin(angulos_sup),
  stringsAsFactors = FALSE
) %>%
  left_join(nodos_df_sup, by = "name") %>%
  mutate(Estado = ifelse(Conectado, "Conectado", "Aislado"))

# Coordenadas de edges (solo si hay edges)
if (ecount(red_sup) > 0) {
  coords_edges_sup <- edges_clean_sup %>%
    left_join(coords_nodos_sup %>% select(name, x, y),
              by = c("nodo1" = "name")) %>%
    rename(x1 = x, y1 = y) %>%
    left_join(coords_nodos_sup %>% select(name, x, y),
              by = c("nodo2" = "name")) %>%
    rename(x2 = x, y2 = y)
} else {
  coords_edges_sup <- data.frame(
    x1 = numeric(), y1 = numeric(),
    x2 = numeric(), y2 = numeric(),
    score_norm = numeric()
  )
}

colores_dir_sup <- c(
  "Concordante Up"   = "#D0021B",
  "Concordante Down" = "#417505")

formas_nodo_sup <- c("Conectado" = 17, "Aislado" = 1)

p_red_sup <- ggplot() +
  
  # Edges 
  {
    if (nrow(coords_edges_sup) > 0)
      geom_segment(
        data = coords_edges_sup,
        aes(x = x1, y = y1, xend = x2, yend = y2,
            alpha     = score_norm,
            linewidth = score_norm),
        color = "grey50"
      )
  } +
  {
    if (nrow(coords_edges_sup) > 0)
      list(
        scale_alpha_continuous(range = c(0.3, 0.9),
                               name  = "Score STRING"),
        scale_linewidth_continuous(range = c(0.4, 2.0),
                                   name  = "Score STRING")
      )
  } +
  
  # Nodos 
  geom_point(
    data  = coords_nodos_sup,
    aes(x     = x,
        y     = y,
        color = Concordancia,
        size  = Score_norm,
        shape = Estado),
    alpha = 0.92
  ) +
  scale_color_manual(
    values   = colores_dir_sup,
    na.value = "grey60",
    name     = "Dirección\nconcordante"
  ) +
  scale_size_continuous(
    range = c(4, 13),
    name  = "Score\nintegración"
  ) +
  scale_shape_manual(
    values = formas_nodo_sup,
    name   = "Conectividad"
  ) +
  
  # Etiquetas 
  geom_text(
    data = coords_nodos_sup,
    aes(x        = x * 1.22,
        y        = y * 1.22,
        label    = name,
        hjust    = ifelse(x >  0.1, 0,
                          ifelse(x < -0.1, 1, 0.5)),
        vjust    = ifelse(y >  0.1, 0,
                          ifelse(y < -0.1, 1, 0.5)),
        fontface = ifelse(Estado == "Conectado",
                          "bold", "plain")),
    size  = 3.2,
    color = "grey15"
  ) +
  
  coord_fixed(xlim = c(-1.6, 1.6),
              ylim = c(-1.6, 1.6)) +
  
  theme_void() +
  theme(
    plot.title      = element_text(face  = "bold", size = 13,
                                   hjust = 0.5,
                                   margin = margin(b = 6)),
    plot.subtitle   = element_text(size  = 9,  color = "grey35",
                                   hjust = 0.5,
                                   margin = margin(b = 4)),
    plot.caption    = element_text(size  = 7.5, color = "grey55",
                                   hjust = 1,
                                   margin = margin(t = 6)),
    legend.position = "right",
    legend.title    = element_text(size = 8.5, face = "bold"),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(20, 20, 20, 20)
  ) +
  labs(
    title    = "Red de coexpresión STRING — Candidatos concordantes Supervivencia",
    subtitle = paste0(
      vcount(red_sup), " genes | ",
      ecount(red_sup), " conexiones (score ≥ ", umbral_sup, ") | ",
      sum(degree(red_sup) > 0), " conectados — ",
      sum(degree(red_sup) == 0), " aislados"
    ),
    caption  = paste0(
      "▲ Conectado  ○ Aislado (sin evidencia de coexpresión)\n",
      "Rojo = Up en Supervivencia Corta | ",
      "Verde = Down en Supervivencia Corta (Up en Larga)\n",
      "Fuente: STRING v12.0 | Canis lupus familiaris (taxón 9615)"
    )
  )

# GUARDAR 
ggsave("red_coexpresion_STRING_Supervivencia.png",
       plot   = p_red_sup,
       width  = 12,
       height = 10,
       dpi    = 300,
       bg     = "white")

ggsave("red_coexpresion_STRING_Supervivencia.pdf",
       plot   = p_red_sup,
       width  = 12,
       height = 10,
       device = cairo_pdf,
       bg     = "white")

#  RED DE COEXPRESIÓN STRING — CANDIDATOS CONCORDANTES EDAD
# CARGAR DATOS
ruta_string_edad <- "C:/Users/danie/OneDrive/Documentos/Trabajo de Grado Maestría Bioinformática/Versiones finales/Script integración/string_network_Edad.tsv"  

candidatos_edad <- read_excel(ruta_excel, sheet = "Edad") %>%
  select(Gene_Symbol, score_integración, logFC_prot,
         Concordancia, FDR_prot) %>%
  filter(!is.na(Gene_Symbol) & !is.na(score_integración)) %>%
  filter(grepl("Concordante", Concordancia)) %>%
  group_by(Gene_Symbol) %>%
  slice_max(order_by = abs(logFC_prot), n = 1,
            with_ties = FALSE) %>%
  ungroup()

genes_semilla_edad <- candidatos_edad$Gene_Symbol

cat("Total:", length(genes_semilla_edad), "\n")
print(candidatos_edad %>%
        select(Gene_Symbol, Concordancia,
               logFC_prot, score_integración) %>%
        arrange(desc(score_integración)))

if (length(genes_semilla_edad) < 2) {
  stop("Menos de 2 genes concordantes — no es posible construir la red.")}

# Archivo TSV de STRING 
edges_raw_edad <- read_tsv(ruta_string_edad, show_col_types = FALSE)

cat("Filas:", nrow(edges_raw_edad), "\n")
cat("Columnas:", paste(colnames(edges_raw_edad), collapse = ", "), "\n")
print(head(edges_raw_edad, 3))

# Verificar que el archivo tiene la estructura correcta de edges
if (any(c("x_position", "y_position") %in%
        colnames(edges_raw_edad))) {
  stop(paste(
    "⚠ Archivo incorrecto: cargaste 'network coordinates'",
    "en vez de 'tabular text output'.\n",
    "En STRING → Exports → selecciona",
    "'as tabular text output'\n",
    "(el que dice 'lists reciprocal edges: A-B, B-A')"
  ))
}

# LIMPIAR EDGES DE STRING
col_n1_edad <- intersect(c("node1", "protein1", "#node1"),
                         colnames(edges_raw_edad))[1]
col_n2_edad <- intersect(c("node2", "protein2"),
                         colnames(edges_raw_edad))[1]
col_sc_edad <- intersect(c("combined_score", "score"),
                         colnames(edges_raw_edad))[1]

cat("\nColumnas — Nodo1:", col_n1_edad,
    "| Nodo2:", col_n2_edad,
    "| Score:", col_sc_edad, "\n")

limpiar_id <- function(x) gsub("^[0-9]+\\.", "", x)

edges_raw2_edad <- edges_raw_edad %>%
  rename(nodo1 = all_of(col_n1_edad),
         nodo2 = all_of(col_n2_edad),
         score = all_of(col_sc_edad)) %>%
  mutate(nodo1 = limpiar_id(nodo1),
         nodo2 = limpiar_id(nodo2))

score_max_edad  <- max(edges_raw2_edad$score, na.rm = TRUE)
edges_raw2_edad <- edges_raw2_edad %>%
  mutate(score_norm = ifelse(score_max_edad > 1,
                             score / 1000,
                             score))

cat("Score máximo:", score_max_edad, "→",
    ifelse(score_max_edad > 1,
           "normalizado ÷1000",
           "ya en escala 0-1"), "\n")
cat("Resumen de scores:\n")
print(summary(edges_raw2_edad$score_norm))

# Filtrar por umbral 
umbral_edad <- 0.40

edges_clean_edad <- edges_raw2_edad %>%
  filter(score_norm >= umbral_edad) %>%
  filter(nodo1 != nodo2) %>%
  mutate(par_min = pmin(nodo1, nodo2),
         par_max = pmax(nodo1, nodo2)) %>%
  distinct(par_min, par_max, .keep_all = TRUE) %>%
  select(nodo1, nodo2, score_norm)

# Si no hay edges con 0.40, bajar automáticamente a 0.15
if (nrow(edges_clean_edad) == 0) {
  cat("\n⚠ Sin edges con score ≥", umbral_edad,
      "— bajando umbral a 0.15\n")
  umbral_edad <- 0.15
  edges_clean_edad <- edges_raw2_edad %>%
    filter(score_norm >= umbral_edad) %>%
    filter(nodo1 != nodo2) %>%
    mutate(par_min = pmin(nodo1, nodo2),
           par_max = pmax(nodo1, nodo2)) %>%
    distinct(par_min, par_max, .keep_all = TRUE) %>%
    select(nodo1, nodo2, score_norm)
}

cat("\nEdges tras filtro (score ≥", umbral_edad, "):",
    nrow(edges_clean_edad), "\n")

if (nrow(edges_clean_edad) == 0) {
  cat("⚠ Sin edges con ningún umbral.\n")
  cat("  Todos los genes de Edad están aislados en STRING.\n")}

# CONSTRUIR EL GRAFO
nodos_df_edad <- data.frame(name = genes_semilla_edad,
                            stringsAsFactors = FALSE) %>%
  left_join(candidatos_edad %>%
              select(name = Gene_Symbol,
                     score_integración, logFC_prot,
                     Concordancia, FDR_prot),
            by = "name") %>%
  mutate(
    Score_norm = (score_integración -
                    min(score_integración, na.rm = TRUE)) /
      (max(score_integración, na.rm = TRUE) -
         min(score_integración, na.rm = TRUE)),
    Sig_FDR   = FDR_prot < 0.05 & abs(logFC_prot) > 1,
    Conectado = name %in% c(edges_clean_edad$nodo1,
                            edges_clean_edad$nodo2)
  )

red_edad <- graph_from_data_frame(
  d        = if (nrow(edges_clean_edad) > 0) edges_clean_edad else
    data.frame(nodo1      = character(),
               nodo2      = character(),
               score_norm = numeric()),
  directed = FALSE,
  vertices = nodos_df_edad
)

V(red_edad)$Grado <- degree(red_edad)

cat("Nodos totales:             ", vcount(red_edad), "\n")
cat("Nodos con conexiones:      ", sum(degree(red_edad) > 0), "\n")
cat("Nodos aislados (grado = 0):", sum(degree(red_edad) == 0), "\n")
cat("Edges:                     ", ecount(red_edad), "\n")
if (ecount(red_edad) > 0) {
  cat("Densidad:                  ",
      round(edge_density(red_edad), 3), "\n")}

aislados_edad <- V(red_edad)$name[degree(red_edad) == 0]
if (length(aislados_edad) > 0) {
  cat("\nGenes sin conexión (score ≥", umbral_edad, "):\n")
  print(aislados_edad)}

# MÉTRICAS DE CENTRALIDAD
metricas_edad <- data.frame(
  Gene        = V(red_edad)$name,
  Grado       = degree(red_edad),
  Betweenness = round(betweenness(red_edad, normalized = TRUE), 4),
  Closeness   = round(closeness(red_edad,   normalized = TRUE), 4),
  Eigenvector = round(eigen_centrality(red_edad)$vector, 4)
) %>%
  left_join(candidatos_edad %>%
              select(Gene = Gene_Symbol, Concordancia,
                     logFC_prot, score_integración, FDR_prot),
            by = "Gene") %>%
  arrange(desc(Grado))

print(metricas_edad %>%
        select(Gene, Grado, Concordancia,
               logFC_prot, score_integración))

# VISUALIZACIÓN CON ggplot2
n_edad       <- vcount(red_edad)
angulos_edad <- seq(0, 2 * pi, length.out = n_edad + 1)[1:n_edad]

coords_nodos_edad <- data.frame(
  name  = V(red_edad)$name,
  grado = V(red_edad)$Grado,
  x     = cos(angulos_edad),
  y     = sin(angulos_edad),
  stringsAsFactors = FALSE
) %>%
  left_join(nodos_df_edad, by = "name") %>%
  mutate(Estado = ifelse(Conectado, "Conectado", "Aislado"))

# Coordenadas de edges
if (ecount(red_edad) > 0) {
  coords_edges_edad <- edges_clean_edad %>%
    left_join(coords_nodos_edad %>% select(name, x, y),
              by = c("nodo1" = "name")) %>%
    rename(x1 = x, y1 = y) %>%
    left_join(coords_nodos_edad %>% select(name, x, y),
              by = c("nodo2" = "name")) %>%
    rename(x2 = x, y2 = y)
} else {
  coords_edges_edad <- data.frame(
    x1 = numeric(), y1 = numeric(),
    x2 = numeric(), y2 = numeric(),
    score_norm = numeric()
  )
}

colores_dir_edad <- c(
  "Concordante Up"   = "#50B83C",   # Geriátrico
  "Concordante Down" = "#F5A623"    # Adulto
)

formas_nodo_edad <- c("Conectado" = 17, "Aislado" = 1)

p_red_edad <- ggplot() +
  
  # Edges 
  {
    if (nrow(coords_edges_edad) > 0)
      geom_segment(
        data = coords_edges_edad,
        aes(x = x1, y = y1, xend = x2, yend = y2,
            alpha     = score_norm,
            linewidth = score_norm),
        color = "grey50"
      )
  } +
  {
    if (nrow(coords_edges_edad) > 0)
      list(
        scale_alpha_continuous(range = c(0.3, 0.9),
                               name  = "Score STRING"),
        scale_linewidth_continuous(range = c(0.4, 2.0),
                                   name  = "Score STRING")
      )
  } +
  
  # Nodos 
  geom_point(
    data  = coords_nodos_edad,
    aes(x     = x,
        y     = y,
        color = Concordancia,
        size  = Score_norm,
        shape = Estado),
    alpha = 0.92
  ) +
  scale_color_manual(
    values   = colores_dir_edad,
    na.value = "grey60",
    name     = "Dirección\nconcordante"
  ) +
  scale_size_continuous(
    range = c(4, 13),
    name  = "Score\nintegración"
  ) +
  scale_shape_manual(
    values = formas_nodo_edad,
    name   = "Conectividad"
  ) +
  
  # Etiquetas 
  geom_text(
    data = coords_nodos_edad,
    aes(x        = x * 1.22,
        y        = y * 1.22,
        label    = name,
        hjust    = ifelse(x >  0.1, 0,
                          ifelse(x < -0.1, 1, 0.5)),
        vjust    = ifelse(y >  0.1, 0,
                          ifelse(y < -0.1, 1, 0.5)),
        fontface = ifelse(Estado == "Conectado",
                          "bold", "plain")),
    size  = 3.2,
    color = "grey15"
  ) +
  
  coord_fixed(xlim = c(-1.6, 1.6),
              ylim = c(-1.6, 1.6)) +
  
  theme_void() +
  theme(
    plot.title      = element_text(face  = "bold", size = 13,
                                   hjust = 0.5,
                                   margin = margin(b = 6)),
    plot.subtitle   = element_text(size  = 9,  color = "grey35",
                                   hjust = 0.5,
                                   margin = margin(b = 4)),
    plot.caption    = element_text(size  = 7.5, color = "grey55",
                                   hjust = 1,
                                   margin = margin(t = 6)),
    legend.position = "right",
    legend.title    = element_text(size = 8.5, face = "bold"),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(20, 20, 20, 20)
  ) +
  labs(
    title    = "Red de coexpresión STRING — Candidatos concordantes Edad",
    subtitle = paste0(
      vcount(red_edad), " genes | ",
      ecount(red_edad), " conexiones (score ≥ ", umbral_edad, ") | ",
      sum(degree(red_edad) > 0), " conectados — ",
      sum(degree(red_edad) == 0), " aislados"
    ),
    caption  = paste0(
      "▲ Conectado  ○ Aislado (sin evidencia de coexpresión)\n",
      "Verde = Up en Geriátrico | Naranja = Down en Geriátrico (Up en Adulto)\n",
      "Fuente: STRING v12.0 | Canis lupus familiaris (taxón 9615)"
    )
  )

# GUARDAR
ggsave("red_coexpresion_STRING_Edad.png",
       plot   = p_red_edad,
       width  = 12,
       height = 10,
       dpi    = 300,
       bg     = "white")

ggsave("red_coexpresion_STRING_Edad.pdf",
       plot   = p_red_edad,
       width  = 12,
       height = 10,
       device = cairo_pdf)

