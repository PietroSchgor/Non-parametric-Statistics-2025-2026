rm(list=ls())
graphics.off()
cat("\014")

setwd("~/uni/2025-2026/non param/progetto/clorofilla/Non-parametric-Statistics-2025-2026")



# 1. Caricamento e Preparazione Dati -------------------------------------------
library(tidyverse)

df <- read.csv("dataset/dataset_MERGIATO")

df$Date <- as.Date(df$Date)

library(sf)

# Prepara i dati come prima
coords_clean <- df %>% 
  filter(as.Date(Date) == as.Date("2022-05-31")) %>% 
  dplyr::select(Lon, Lat) %>%
  na.omit()

# Trasforma in oggetto spaziale (sf)
punti_sf <- st_as_sf(coords_clean, coords = c("Lon", "Lat"), crs = 4326)

# Unisci i punti
punti_uniti <- st_union(punti_sf)

# Calcola il bordo concavo
# ratio: 1 = Convex Hull, 0.1 = Molto aderente/concavo
bordo_concavo <- st_concave_hull(punti_uniti, ratio = 0.05)

# Visualizza
plot(st_geometry(bordo_concavo), border = "red", lwd = 2)
points(st_coordinates(punti_sf), col = "grey")

# Estrae le coordinate come matrice
coordinate_bordo <- st_coordinates(bordo_concavo)

# Converte in dataframe e pulisce (st_coordinates aggiunge colonne extra L1, L2)
bordo_df <- as.data.frame(coordinate_bordo) %>%
  dplyr::select(Lon = X, Lat = Y)


summary(bordo_df)


bordo_df_clean <- bordo_df %>%
  filter(Lat >= 44.10, Lon <= 13.90)

# 1. Disegna il confine completo
plot(bordo_df_clean, col = "blue") # "l" sta per linea

point0 <- bordo_df_clean %>%
  filter(Lon <= 12.586, Lat <= 44.110)
point0 <- point0[1,]
# 2. Aggiungi il primo punto sopra il grafico esistente
points(point0, col = "red", pch = 19, cex = 1.5)

point0_index <- which(abs(bordo_df_clean$Lat - 44.104) < 0.0001 & 
                        abs(bordo_df_clean$Lon - 12.583) < 0.0001)
point0_index <- as.numeric(point0_index)
# collasso in 1-D ---------------------------------------------------------


# 2. Funzione per calcolare l'ordinamento spaziale
reorder_points <- function(df, starting_index) {
  # Trasforma in matrice per velocità
  mat <- as.matrix(df[, c("Lon", "Lat")])
  n <- nrow(mat)
  
  unvisited <- 1:n
  ordered_indices <- numeric(n)
  
  # Partiamo dal primo punto
  current_idx <- starting_index
  ordered_indices[1] <- current_idx
  unvisited <- unvisited[unvisited != current_idx]
  
  for (i in 2:n) {
    # Calcola distanza euclidea tra il punto corrente e tutti i non visitati
    distances <- sqrt(colSums((t(mat[unvisited, , drop=FALSE]) - mat[current_idx, ])^2))
    
    # Trova il più vicino
    next_idx_in_unvisited <- which.min(distances)
    current_idx <- unvisited[next_idx_in_unvisited]
    
    # Aggiorna liste
    ordered_indices[i] <- current_idx
    unvisited <- unvisited[-next_idx_in_unvisited]
  }
  
  return(df[ordered_indices, ])
}

# 3. Applica l'ordinamento
bordo_ordinato <- reorder_points(bordo_df_clean, point0_index)


# 4. Visualizza per verifica
plot(bordo_ordinato$Lon, bordo_ordinato$Lat, type = "b", pch = 16, col = "blue",
     main = "Bordo Ordinato (Nearest Neighbor)")
points(bordo_ordinato[1, , drop = FALSE], col = "red", pch = 19, cex = 1.5)


# 1. Calcolo delle distanze euclidee tra punti consecutivi
# Calcoliamo la differenza tra ogni punto e il precedente
lon_diff <- diff(bordo_ordinato$Lon)
lat_diff <- diff(bordo_ordinato$Lat)

# Calcoliamo la distanza euclidea per ogni segmento: sqrt(deltaLon^2 + deltaLat^2)
# Aggiungiamo uno 0 all'inizio perché il primo punto non ha un "precedente"
distanze_segmenti <- c(0, sqrt(lon_diff^2 + lat_diff^2))

# 2. Costruzione del dataset finale
dataset_ordinato_distanza <- bordo_ordinato %>%
  mutate(
    distanza_segmento = distanze_segmenti,        # Distanza dal punto precedente
    distanza_cumulata = cumsum(distanze_segmenti) # Somma totale dal punto di partenza
  )

# 3. Visualizzazione delle prime righe per controllo
head(dataset_ordinato_distanza)
dataset_ordinato_distanza <- dataset_ordinato_distanza %>%
  distinct(Lat, Lon, .keep_all = TRUE)
# Conta quante righe sono "copie" di valori già visti
sum(duplicated(dataset_ordinato_distanza$distanza_cumulata))

dataset_ordinato_distanza$Pos_ID <- 1:dim(dataset_ordinato_distanza)[1]

#### plot ####
library(ggplot2)
library(viridis) # Per una scala colori leggibile

# 1. Grafico del percorso spaziale
plot_mappa <- ggplot(dataset_ordinato_distanza, aes(x = Lon, y = Lat)) +
  # Disegna la linea che unisce i punti seguendo l'ordine
  geom_path(aes(color = distanza_cumulata), size = 1.2, 
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Aggiunge i punti
  geom_point(alpha = 0.5) +
  # Evidenzia il punto di partenza (point0) in rosso
  geom_point(data = dataset_ordinato_distanza[1, ], color = "red", size = 4, shape = 18) +
  # Scala cromatica per la distanza
  scale_color_viridis_c(option = "plasma", name = "Distanza Totale") +
  theme_minimal() +
  labs(title = "Analisi del Bordo: Percorso e Distanza",
       subtitle = "Il punto rosso indica l'inizio; le frecce mostrano la direzione dell'ordinamento",
       x = "Longitudine", y = "Latitudine")

# 2. Grafico dell'andamento della distanza
plot_distanza <- ggplot(dataset_ordinato_distanza, aes(x = 1:nrow(dataset_ordinato_distanza), y = distanza_cumulata)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point() +
  theme_light() +
  labs(title = "Progressione della Distanza Euclidea",
       x = "Indice del Punto (Sequenza)", 
       y = "Distanza Cumulata")

# Visualizza i grafici
print(plot_mappa)
print(plot_distanza)



# Write CSV ---------------------------------------------------------------

# Se vuoi salvare il risultato in un file CSV:
#write.csv(dataset_ordinato_distanza, "dataset/bordo_con_distanze.csv", row.names = FALSE)


# Merge data_weekly -------------------------------------------------------

data_weekly <- read.csv("dataset/data_weekly.csv")

data_pos_x_weekly <- data_weekly %>%
                    inner_join(dataset_ordinato_distanza, by = c("Lat", "Lon"))

chl_pos_x_weekly <- data_pos_x_weekly %>%
  select(2:5, 14:16)

# 
# write.csv(data_pos_x_weekly, "dataset/bordo_x_weekly.csv", row.names = FALSE)
# 
# write.csv(chl_pos_x_weekly, "dataset/chl_bordo_x_weekly.csv", row.names = FALSE)
