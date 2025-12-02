# Pipeline R: Functional Data Analysis + Functional Kriging
# Dataset: /mnt/data/dataset_MERGIATO
# Scopo: prevedere Clorofilla (Chl) sulla costa a +14 giorni.

# libraries ---------------------------------------------------------------


# 0) Requisiti: installare pacchetti se necessario
required <- c("tidyverse", "lubridate", "fda", "refund", "mgcv", "gstat", "sp", "spacetime", "sf", "fields")
new.packages <- required[!(required %in% installed.packages()[,"Package"]) ]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(lubridate)
library(fda)
library(refund)   # pffr, fpca.sc
library(mgcv)
library(gstat)
library(sp)
library(spacetime)
library(sf)
library(fields)   # for spatial plotting



# dataset loading ---------------------------------------------------------


# 1) Carica dataset (usa il path locale caricato)
path <- "dataset/dataset_MERGIATO"
# il file dovrebbe essere in formato csv; se non lo fosse, adattare
df <- read.csv(path, stringsAsFactors = FALSE)

# 2) Pulizia minima
df <- df %>%
  mutate(Date = as.Date(Date)) %>%
  arrange(Lon, Lat, Date)

sum(is.na(df))


# 3) Costruisci matrici per ogni location
# identifichiamo le posizioni uniche (lon, lat)
library(data.table)

# Converti in data.table se non lo è già
setDT(df)

locations <- unique(df[, .(Lon, Lat)])
locations$loc_id <- seq_len(nrow(locations))

library(ggplot2)
library(maps) # Necessario per i confini geografici

# Recuperiamo i dati della mappa mondiale
world_map <- map_data(map('italy', fill = TRUE, col = 1:10))

ggplot() +
  # 1. Disegna la mappa di sfondo
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white") +
  # 2. Aggiungi i tuoi punti
  geom_point(data = locations, aes(x = Lon, y = Lat),
             color = "red", size = 2, alpha = 0.7) +
  # 3. Migliora l'estetica
  theme_minimal() +
  labs(title = "Distribuzione dei Punti Unici", x = "Longitudine", y = "Latitudine") +
  coord_fixed(1.3) # Mantiene le proporzioni corrette della mappa

df <- df %>% left_join(locations[], by = c("Lon","Lat"))

# sequenza temporale completa (globalmente)
dates <- sort(unique(df$Date))
T <- as.numeric(length(dates))

# costruisci matrici: righe = locazioni, colonne = tempo (ordinato)
make_matrix <- function(varname){
  m <- matrix(NA, nrow = nrow(locations), ncol = T)
  colnames(m) <- as.character(dates)
  for(i in seq_len(nrow(locations))){
    tmp <- df %>% filter(loc_id == i) %>% arrange(Date)
    if(nrow(tmp)>0){
      idx <- match(as.character(tmp$Date), colnames(m))
      m[i, idx] <- tmp[[varname]]
    }
  }
  return(m)
}

chl_mat <- make_matrix("Chl")
temp_mat <- make_matrix("Temp")
sal_mat  <- make_matrix("Salinity")
solar_mat <- make_matrix("Solar_Flux")
heat_mat <- make_matrix("Heat_Flux")
wind_u_mat <- make_matrix("uo")
wind_v_mat <- make_matrix("vo")

time <- seq(1, 366)

par(mfrow=c(1,1))
par(mar = c(4, 4, 2, 1))
# Prima il grafico
matplot(time, t(chl_mat)[, 1:30], type = 'l', lwd = 3, main = "Chlorofilla curves",
        xlab = "day", ylab = " ")


# smoothing ---------------------------------------------------------------


# 4) Smoothing funzionale (B-splines)
# definisci base temporale (parametro: nbasis da adattare)
nbasis <- min(35, floor(T/2))
time = seq(1:T)
basis <- create.bspline.basis(rangeval = c(1, T), nbasis = nbasis)

chl_fd <- Data2fd(y=t(chl_mat), argvals=time,  basisobj=basis)
temp_fd <- Data2fd(y=t(temp_mat), argvals=time, basisobj=basis)
sal_fd  <- Data2fd(y=t(temp_mat), argvals=time, basisobj=basis)



# FPCA --------------------------------------------------------------------


# FPCA (number of components nharm)
nharm <- 6
chl_fpca <- pca.fd(chl_fd, nharm=nharm, centerfns=TRUE)
# scores: nloc x nharm
# scree plot
par(mfrow=c(1,2))
plot(chl_fpca$values[1:15], xlab='j', ylab='Eigenvalues')
plot(cumsum(chl_fpca$values)[1:15]/sum(chl_fpca$values), xlab='j', ylab='CPV', ylim=c(0.8, 1))

cumsum(chl_fpca$values)[3]/sum(chl_fpca$values)
# 0.9585892
chl_fpca$varprop[1:3]

# eigenfunctions
par(mfrow = c(1,3))
plot(chl_fpca$harmonics[1,],col=1,ylab='FPC1')
abline(h=0,lty=2)
plot(chl_fpca$harmonics[2,],col=2,ylab='FPC2')
plot(chl_fpca$harmonics[3,],col=3,ylab='FPC3')

par(mfrow = c(1,3))
plot.pca.fd(chl_fpca, nx=100, pointplot=TRUE, harm=c(1,2,3), expand=0, cycle=FALSE)
#