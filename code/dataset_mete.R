# Setup Workspace
setwd("~/uni/2025-2026/non param/progetto/clorofilla/Non-parametric-Statistics-2025-2026")

# Libraries
library(tidyverse)
# 1. Carica le librerie necessarie
library(readxl)

# 2. Carica i dataset
dataset_ordinato_distanza <- read.csv("dataset/bordo_con_distanze.csv")
dataset_mete <- read_excel("dataset/Localita_balneari_Nord_Adriatico.xlsx")


mete_top <- c("Grado", "Lignano Sabbiadoro", "Bibione", "Caorle", 
              "Lido di Jesolo", "Lido di Venezia", "Sottomarina (Chioggia)", 
              "Milano Marittima", "Cervia", "Cesenatico", "Rimini")

# 2. Creiamo il dataset filtrato
# Assumendo che il tuo dataset si chiami 'pos_mete'
dataset_mete <- dataset_mete[
  dataset_mete$Località %in% mete_top, ]

# 3. Funzione per trovare il punto più vicino
# Assumiamo che le colonne si chiamino 'Lat' e 'Lon' in entrambi i file.
# Se hanno nomi diversi (es. 'latitude', 'longitude'), cambiali qui sotto.

trova_piu_vicino <- function(lat_target, lon_target, df_riferimento) {
  # Calcolo della distanza Euclidea semplice (ottima per distanze brevi)
  distanze <- sqrt((df_riferimento$Lat - lat_target)^2 + (df_riferimento$Lon - lon_target)^2)
  
  # Restituisce l'indice della riga con la distanza minima
  return(which.min(distanze))
}

# 4. Trova gli indici corrispondenti
# Applichiamo la funzione a ogni riga di dataset_mete
indici_corrispondenti <- sapply(1:nrow(dataset_mete), function(i) {
  trova_piu_vicino(dataset_mete$Lat[i], dataset_mete$Lon[i], dataset_ordinato_distanza)
})

# 5. Crea il nuovo dataset selezionando le righe trovate
dataset_approssimato <- dataset_ordinato_distanza[indici_corrispondenti, c(1:2, 5)]



# Opzionale: Aggiungi i nomi delle località originali per chiarezza
dataset_approssimato$Localita_Originale <- dataset_mete$Località # cambia 'Localita' col nome reale della colonna


plot(dataset_ordinato_distanza$Lon, dataset_ordinato_distanza$Lat)
points(dataset_mete$Lon, dataset_mete$Lat, col = "blue")
points(dataset_approssimato$Lon, dataset_approssimato$Lat, col = "red")

library(ggplot2)

# Creiamo un dataframe di supporto per le linee
df_lines <- data.frame(
  x_start = dataset_mete$Lon,
  y_start = dataset_mete$Lat,
  x_end = dataset_approssimato$Lon,
  y_end = dataset_approssimato$Lat
)
library(ggplot2)

# 1. Prepare the line data (ensure columns match your dataset names)
df_lines <- data.frame(
  x_start = dataset_mete$Lon,
  y_start = dataset_mete$Lat,
  x_end = dataset_approssimato$Lon,
  y_end = dataset_approssimato$Lat
)

# 2. Create the enhanced plot
ggplot() +
  # Background: Border points (lighter grey for better contrast)
  geom_point(data = dataset_ordinato_distanza, aes(x = Lon, y = Lat), 
             color = "grey8", size = 1, alpha = 0.4) +
  
  # Connection lines: Showing the displacement between real and approximated points
  geom_segment(data = df_lines, aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               color = "red", linetype = "solid", linewidth = 1) +
  
  # Real Locations (Mapped to color for the legend)
  geom_point(data = dataset_mete, aes(x = Lon, y = Lat, color = "Real Locations"), 
             size = 2.5) +
  
  # Approximated Points (Mapped to color for the legend)
  geom_point(data = dataset_approssimato, aes(x = Lon, y = Lat, color = "Approximated Points"), 
             size = 2.5) +
  
  # Customizing Legend and Colors
  scale_color_manual(values = c("Real Locations" = "#0073C2FF", 
                                "Approximated Points" = "#EFC000FF")) +
  
  # Styling and English Labels
  theme_minimal() +
  coord_fixed(ratio = 1.3) +
  labs(
    title = "Coordinate Matching: Real vs. Approximated Points",
    subtitle = "Displacement from original beach locations to the nearest border points",
    x = "Longitude", 
    y = "Latitude",
    color = "Legend"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )
# 6. Salva il risultato
write.csv(dataset_approssimato, "dataset/dataset_mete_approssimate.csv", row.names = FALSE)

# Visualizza i primi risultati
head(dataset_approssimato)
