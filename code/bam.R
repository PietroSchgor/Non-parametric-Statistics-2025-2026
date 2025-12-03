
setwd("~/uni/2025-2026/non param/progetto/clorofilla/Non-parametric-Statistics-2025-2026")

library(tidyverse)
library(mgcv)
library(lubridate)
library(parallel)

# 1. Caricamento e Preparazione Dati -------------------------------------------
df <- read_csv("dataset/dataset_MERGIATO")



library(ggplot2)
library(dplyr)
library(viridis) # Per le scale di colori scientifiche

# 1. Calcola la media per ogni pixel (Lat/Lon) su tutto l'anno
mappa_media <- df %>%
  group_by(Lon, Lat) %>%
  summarise(Mean_Chl = mean(Chl, na.rm = TRUE))

# 2. Plotta la media
ggplot(mappa_media, aes(x = Lon, y = Lat, fill = Mean_Chl)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "Media Chl") + # 'magma' va dal nero al bianco/giallo
  coord_fixed(ratio = 1.3) +
  labs(title = "Distribuzione Media della Clorofilla (2022-2023)") +
  theme_minimal()

#### 2. Logica di filtraggio ####

summary(df$Chl)
# Vogliamo raggruppare per posizione e contare i giorni con Chl > 0.5
posizioni_target <- df %>%
  group_by(Lat, Lon) %>%
  summarise(
    # La condizione (Chl > 0.5) restituisce TRUE/FALSE.
    # sum() conta i TRUE come 1. na.rm=TRUE gestisce eventuali dati mancanti.
    giorni_sopra_soglia = sum(Chl > 0.1, na.rm = TRUE),
    .groups = 'drop' # Sgruppa il risultato per evitare warning
  ) %>%
  mutate(loc_id = row_number()) %>%
  # Manteniamo solo chi ha 200 o più giorni
  filter(giorni_sopra_soglia >= 366)
  
  


# Stampa quante posizioni sono state trovate
cat("Numero di posizioni trovate:", nrow(posizioni_target), "\n")


library(gridExtra)
library(ggplot2)

# --- Grafico 1: Assegnalo a p1 ---
p1 <- ggplot(mappa_media, aes(x = Lon, y = Lat, fill = Mean_Chl)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "Media Chl") +
  coord_fixed(ratio = 1.3, xlim = c(12.25, 14), ylim = c(44, 45.75)) +
  labs(title = "Distribuzione Media") +
  theme_minimal() +
  theme(legend.position = "bottom") # Opzionale: sposta la legenda sotto

# --- Grafico 2: Assegnalo a p2 ---
p2 <- ggplot(posizioni_target, aes(x = Lon, y = Lat, fill = 1)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "") +
  coord_fixed(ratio = 1.3, xlim = c(12.25, 14), ylim = c(44, 45.75)) +
  labs(title = "Posizioni filtrate") +
  theme_minimal() +
  theme(legend.position = "bottom")

# --- Uniscili ---
# Affiancati orizzontalmente:
grid.arrange(p1, p2, ncol = 2)

# Facciamo un inner join tra il dataset originale e le posizioni filtrate
data_proc <- df %>%
  inner_join(posizioni_target, by = c("Lat", "Lon")) %>%
  mutate(
    Date = ymd(Date),
    doy = yday(Date),                   # Giorno dell'anno (1-365)
    time_numeric = as.numeric(Date),    # Tempo continuo per il trend
    # Creiamo un ID univoco per ogni coordinata (per il plot successivo)
  ) %>%
  drop_na(Chl, Lat, Lon, Temp, Salinity) %>%
  arrange(Date) # È fondamentale che i dati siano ordinati temporalmente

summary(df$Chl)
summary(data_proc$Chl)

# prendo un subset per velocizzare i tempi

# summary(df$Lon)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12.29   12.75   13.12   13.12   13.46   13.96 
#data_proc <- data_proc[which(data_proc$Lon <= 13), ]
# summary(df$Lat)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 44.02   44.35   44.73   44.76   45.15   45.73
#data_proc <- data_proc[which(data_proc$Lat <= 44.73), ]

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
  facet_wrap(~loc_id, scales = "fixed") +
  
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






# NUOVO BAM AR -------------------------------------------------------------------

#### --- 1. Feature Engineering: Creazione Lag --- ####
# Ordiniamo per posizione e data per essere sicuri che il 'lag' prenda il giorno prima
data_lagged <- data_proc %>%
  arrange(Lon, Lat, Date) %>%
  group_by(Lon, Lat) %>%
  mutate(
    # Variabili Autoregressive (Valori del giorno precedente)
    Chl_lag1      = lag(Chl, 1),
    Temp_lag1     = lag(Temp, 1),
    Salinity_lag1 = lag(Salinity, 1),
    MLD_lag1      = lag(MLD, 1),
    Solar_lag1    = lag(Solar_Flux, 1),
    
    # Manteniamo le info geografiche/temporali per il modello
    # (Non serve laggare Lat, Lon, Doy perché li conosciamo per il futuro)
  ) %>%
  ungroup() %>%
  # Rimuoviamo la prima riga di ogni serie temporale (che ora ha NA nei lag)
  drop_na(Chl_lag1, Temp_lag1, Salinity_lag1, MLD_lag1, Solar_lag1)

# --- 2. Split Training / Validation ---
# Usiamo i dati 'laggati'
train_set <- data_lagged %>% filter(Date <= cutoff_date)
val_set   <- data_lagged %>% filter(Date > cutoff_date) # Questo ci serve solo per confronto (Ground Truth)

cat("Nuovo Training Set:", nrow(train_set), "righe\n")


#### --- 3. Addestramento Modello Autoregressivo --- ####

weights_vec <- (exp(train_set$Chl*.3) - 1)*.5
summary(weights_vec)
summary(train_set$Chl)

model_bam_ar <- bam(
  Chl ~ 
    # 1. Autoregressione sulla Clorofilla (Fondamentale)
    s(Chl_lag1, k=20) +
    
    # 2. Autoregressione sulle Covariate Ambientali
    s(Temp_lag1, k=20) + 
    #s(Salinity_lag1, k=8) + 
    #s(MLD_lag1, k=8) + 
    #s(Solar_lag1, k=15) +
    
    # 3. Componenti Spazio-Temporali (Queste NON vanno laggate, perché sappiamo che giorno è domani)
    te(Lon, Lat, doy, d=c(2,1), bs =c("tp", "cc"), k=c(20, 20)),
  
  data = train_set,
  method = "fREML",
  discrete = TRUE,
  nthreads = n_cores,
  chunk.size = 10000,
  weights = weights_vec,
  family = tw(link = "log") 
)

summary(model_bam_ar)

#### --- 4. Forecasting Ricorsivo a 3 step --- ####

# Prendiamo l'ultimo giorno del training set: questo è il nostro punto di partenza
last_observed_day <- train_set %>%
  filter(Date == max(Date))

# Lista per salvare i risultati dei 3 giorni
forecast_results <- list()

# Copiamo i dati di partenza per iniziare il loop
input_data <- last_observed_day

cat("Inizio previsione ricorsiva per 3 giorni...\n")

for (i in 1:3) {
  
  # A. Determiniamo la data che stiamo per prevedere
  target_date <- unique(last_observed_day$Date) + days(i)
  
  # B. Aggiorniamo le variabili temporali (Date, doy, time_numeric)
  # Queste le conosciamo con certezza
  input_data <- input_data %>%
    mutate(
      Date = target_date,
      doy = yday(target_date),
      time_numeric = as.numeric(target_date)
    )
  
  # C. PREDIZIONE
  # Usiamo il modello sui dati 'input_data' (che contengono i lag del giorno precedente)
  prediction <- predict(model_bam_ar, newdata = input_data, type = "response")
  
  # D. Salviamo il risultato
  # Aggiungiamo una colonna per indicare che giorno di forecast è (1, 2 o 3)
  day_result <- input_data %>%
    mutate(
      Pred_Chl = prediction,
      Forecast_Day = i
    ) %>%
    select(Date, Lon, Lat, Pred_Chl, Forecast_Day) # Teniamo solo l'essenziale
  
  forecast_results[[i]] <- day_result
  
  # E. PREPARAZIONE PER IL GIRO SUCCESSIVO (Il passo ricorsivo)
  # Per prevedere il giorno successivo, la 'Pred_Chl' di oggi diventa la 'Chl_lag1' di domani.
  # Per le variabili fisiche (Temp, Salinity...), non avendole predette, usiamo la Persistenza:
  # manteniamo i valori lag attuali (Temp_lag1 rimane quello che era).
  # Se avessi un modello meteo, aggiorneresti anche Temp_lag1.
  
  input_data <- input_data %>%
    mutate(
      Chl_lag1 = prediction # <--- Qui avviene la magia ricorsiva
      # Temp_lag1 = Temp_lag1 (Assunzione di persistenza meteo)
    )
}

# Uniamo i 3 giorni in un unico dataframe
forecast_final <- bind_rows(forecast_results)

#### --- 5. Visualizzazione: Confronto con il Validation Set Reale --- ####

# Uniamo le previsioni con i dati veri del validation set per vedere l'errore
comparison <- val_set %>%
  select(Date, Lon, Lat, Chl_Real = Chl) %>% # La Chl vera
  inner_join(forecast_final, by = c("Date", "Lon", "Lat"))

# Calcoliamo l'errore (RMSE) per giorno di previsione
metrics <- comparison %>%
  group_by(Forecast_Day) %>%
  summarise(
    RMSE = sqrt(mean((Chl_Real - Pred_Chl)^2)),
    MAE = mean(abs(Chl_Real - Pred_Chl))
  )

print(metrics)


# 1. Selezioniamo N posizioni casuali (es. 9 posizioni)
set.seed(1) # Per rendere la scelta riproducibile
n_posizioni <- 9

# Estraiamo le coordinate uniche presenti nel dataset di confronto
unique_locs <- comparison %>%
  select(Lat, Lon) %>%
  distinct()

# Se abbiamo meno di 9 posizioni, ne prendiamo quante ce ne sono
n_sample <- min(n_posizioni, nrow(unique_locs))
sampled_locs <- unique_locs %>% sample_n(n_sample)

# 2. Filtriamo i dati per queste posizioni
plot_subset <- comparison %>%
  inner_join(sampled_locs, by = c("Lat", "Lon")) %>%
  # Creiamo un'etichetta leggibile per il titolo del grafico
  mutate(Location_Label = paste0("Lat: ", round(Lat, 2), " | Lon: ", round(Lon, 2)))

# 3. Plot con Facet Wrap
ggplot(plot_subset, aes(x = Forecast_Day)) +
  
  # Linea REALE (Continua + Punti pieni)
  geom_line(aes(y = Chl_Real, color = "Reale"), linewidth = 1) +
  geom_point(aes(y = Chl_Real, color = "Reale"), size = 2) +
  
  # Linea PREDETTA (Tratteggiata + Punti vuoti)
  geom_line(aes(y = Pred_Chl, color = "Predetto"), linetype = "dashed", linewidth = 1) +
  geom_point(aes(y = Pred_Chl, color = "Predetto"), shape = 1, size = 2) +
  
  # Divisione in griglia
  facet_wrap(~Location_Label, scales = "free_y", ncol = 3) +
  
  # Stile e Etichette
  scale_color_manual(values = c("Reale" = "black", "Predetto" = "red")) +
  scale_x_continuous(breaks = 1:3) + # Forza l'asse X a mostrare solo 1, 2, 3
  labs(
    title = paste("Forecast Ricorsivo a 3 Giorni - Campione di", n_sample, "posizioni"),
    subtitle = "Confronto tra Clorofilla Reale (Nero) e Predetta (Rosso)",
    y = "Clorofilla (mg/m^3)",
    x = "Giorni nel futuro (1 = Domani, 3 = Tra 3 giorni)",
    color = "Legenda"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgrey", color = NA), # Sfondo etichette grigio
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

#### diagnostica ####
# Imposta una grafica a 2x2 per vedere 4 grafici insieme
par(mfrow = c(2, 2))
gam.check(model_bam_ar)
par(mfrow = c(1, 1)) # Ripristina grafica

# Controlla la concurvità globale
concurvity(model_bam_ar, full = TRUE)

# Calcola i residui sul training set
residui <- residuals(model_bam_ar)

# Mappa dei residui medi per zona
ggplot(train_set, aes(x = Lon, y = Lat, color = residui)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Mappa dei Residui (Rosso = Sottostima, Blu = Sovrastima)")

# Calcola l'autocorrelazione dei residui
par(mfrow=c(1,1))
acf(residuals(model_bam_ar), main = "Autocorrelazione dei Residui")
