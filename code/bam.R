
setwd("~/uni/2025-2026/non param/progetto/clorofilla/Non-parametric-Statistics-2025-2026")

library(tidyverse)
library(mgcv)
library(lubridate)
library(parallel)

# 1. Caricamento e Preparazione Dati -------------------------------------------
df <- read_csv("dataset/dataset_MERGIATO")

data_proc <- df %>%
  mutate(
    Date = ymd(Date),
    doy = yday(Date),                   # Giorno dell'anno (1-365)
    time_numeric = as.numeric(Date),    # Tempo continuo per il trend
    # Creiamo un ID univoco per ogni coordinata (per il plot successivo)
    loc_id = paste("Lat", round(Lat, 2), "- Lon", round(Lon, 2))
  ) %>%
  drop_na(Chl, Lat, Lon, Temp, Salinity) %>%
  arrange(Date) # È fondamentale che i dati siano ordinati temporalmente

# prendo un subset per velocizzare i tempi
data_proc <- data_proc[which(data_proc$Lon <= 13), ]
data_proc <- data_proc[which(data_proc$Lat <= 45), ]

# 2. Creazione Training e Validation Set (Ultima Settimana) --------------------
# Identifichiamo la data massima
num_days2predict = 3

max_date <- max(data_proc$Date)
cutoff_date <- max_date - days(num_days2predict)

cat("Data finale nel dataset:", as.character(max_date), "\n")
cat("Data di taglio (inizio validation):", as.character(cutoff_date), "\n")

# Splitting Temporale
train_set <- data_proc %>% filter(Date <= cutoff_date)
val_set   <- data_proc %>% filter(Date > cutoff_date)

cat("Righe Training Set:", nrow(train_set), "\n")
cat("Righe Validation Set:", nrow(val_set), "\n")

# 3. Addestramento del Modello BAM ---------------------------------------------
# Rileva i core per parallelizzare
n_cores <- detectCores() - 1

cat("Inizio addestramento BAM su", n_cores, "core...\n")

model_bam <- bam(
  Chl ~ 
    # Interazione Spazio-Tempo (la "curva" che cambia per location)
    # k è ridotto leggermente per sicurezza, puoi alzarlo se hai tanta RAM
    te(Lon, Lat, time_numeric, d=c(2,1), k=c(15, 15)) +
    
    # Stagionalità annuale ciclica
    s(doy, bs = "cc", k = 10) +
    
    # Effetti ambientali (Covariate)
    s(Temp, k=8) + 
    s(Salinity, k=8) + 
    s(MLD, k=8) + 
    s(Solar_Flux, k=8),
  
  data = train_set,
  method = "fREML",
  discrete = TRUE,     # Ottimizzazione essenziale per velocità
  nthreads = n_cores,  # Calcolo parallelo
  chunk.size = 10000,
  family =Gamma(link="log") #scat(link = "identity") Gamma(link="log"), #scat
  #rho = .7
)

cat("Modello addestrato.\n")

# 4. Predizione sul Validation Set ---------------------------------------------
preds <- predict(model_bam, newdata = val_set, type = "response", se.fit = TRUE)

# Aggiungiamo le predizioni al dataframe di validazione
val_results <- val_set %>%
  mutate(
    Predicted_Chl = preds$fit,
    SE = preds$se.fit,
    # Calcolo intervallo di confidenza al 95%
    Lower_CI = Predicted_Chl - 1.96 * SE,
    Upper_CI = Predicted_Chl + 1.96 * SE
  )

# Calcolo metriche di errore generali
rmse <- sqrt(mean((val_results$Chl - val_results$Predicted_Chl)^2))
mae <- mean(abs(val_results$Chl - val_results$Predicted_Chl))

cat("Risultati sul Validation Set:\n")
cat("RMSE:", round(rmse, 4), "\n")
cat("MAE:", round(mae, 4), "\n")

# 5. Visualizzazione (Plot Reale vs Predetto) ----------------------------------

# SELEZIONE CAMPIONE:
# Siccome ci sono troppe coordinate, ne prendiamo 6 a caso per vedere le curve
set.seed(10) # Per riproducibilità
sample_locations <- sample(unique(val_results$loc_id), 6)

plot_data <- val_results %>%
  filter(loc_id %in% sample_locations)

# Creazione del grafico
ggplot(plot_data, aes(x = Date)) +
  # Intervallo di confidenza (area grigia)
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "red", alpha = 0.1) +
  
  # Linea Reale (Nera solida)
  geom_line(aes(y = Chl, color = "Reale (Osservato)"), linewidth = 1) +
  geom_point(aes(y = Chl, color = "Reale (Osservato)"), size= 2) +
  
  # Linea Predetta (Rossa tratteggiata)
  geom_line(aes(y = Predicted_Chl, color = "Predetto (Modello)"), linewidth = 1, linetype = "dashed") +
  
  # Suddivisione per coordinate geografiche
  facet_wrap(~loc_id, scales = "free_y") +
  
  # Stile
  theme_minimal() +
  scale_color_manual(values = c("Reale (Osservato)" = "black", "Predetto (Modello)" = "red")) +
  labs(
    title = "Validazione Modello Clorofilla: Reale vs Predetto",
    subtitle = paste("Ultima settimana - RMSE:", round(rmse, 3)),
    y = "Concentrazione Clorofilla",
    x = "Data",
    color = "Legenda"
  ) +
  theme(legend.position = "bottom")

#### diagnostica ####

# Imposta una grafica a 2x2 per vedere 4 grafici insieme
par(mfrow = c(2, 2))
gam.check(model_bam)
par(mfrow = c(1, 1)) # Ripristina grafica

# Controlla la concurvità globale
concurvity(model_bam, full = TRUE)

# Calcola i residui sul training set
residui <- residuals(model_bam)

# Mappa dei residui medi per zona
ggplot(train_set, aes(x = Lon, y = Lat, color = residui)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Mappa dei Residui (Rosso = Sottostima, Blu = Sovrastima)")

# Calcola l'autocorrelazione dei residui
par(mfrow=c(1,1))
acf(residuals(model_bam), main = "Autocorrelazione dei Residui")





