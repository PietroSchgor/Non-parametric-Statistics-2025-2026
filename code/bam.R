
setwd("~/uni/2025-2026/non param/progetto/clorofilla/Non-parametric-Statistics-2025-2026")

library(tidyverse)
library(mgcv)
library(lubridate)
library(parallel)

# 1. Caricamento e Preparazione Dati -------------------------------------------
df <- read.csv("dataset/dataset_MERGIATO")

summary(df)
summary(df$Chl)
Chl <- log((df$Chl +1)*1.5)
summary(Chl)

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
# summary(df$Chl)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.02871  0.11508  0.25781  0.42489  0.54657 11.23414 

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
  coord_fixed(ratio = 1.3, xlim = c(12, 14), ylim = c(44, 46)) +
  labs(title = "Distribuzione Media") +
  theme_minimal() +
  theme(legend.position = "bottom") # Opzionale: sposta la legenda sotto

# --- Grafico 2: Assegnalo a p2 ---
p2 <- ggplot(posizioni_target, aes(x = Lon, y = Lat, fill = 1)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "") +
  coord_fixed(ratio = 1.3, xlim = c(12, 14), ylim = c(44, 46)) +
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
    doy = yday(Date),# Giorno dell'anno (1-365)
    dow = as.factor(wday(Date)),
    time_numeric = as.numeric(Date),    # Tempo continuo per il trend
    # Creiamo un ID univoco per ogni coordinata (per il plot successivo)
  ) %>%
  arrange(Date) # È fondamentale che i dati siano ordinati temporalmente

summary(df$Chl)
summary(data_proc$Chl)
hist(data_proc$Chl, xlim = c(0, 5))
Chl <- log((data_proc$Chl*5 + 1))
summary(Chl)
hist(Chl, xlim = c(0, 5))
#data_proc$Chl <- Chl


# prendo un subset per velocizzare i tempi

# summary(df$Lon)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12.29   12.75   13.12   13.12   13.46   13.96 
#data_proc <- data_proc[which(data_proc$Lon <= 13), ]
# summary(df$Lat)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 44.02   44.35   44.73   44.76   45.15   45.73
data_45 <- data_proc[which(data_proc$Lat >= 45), ]


par(mfrow=c(1,2))
summary(data_45$Chl)
hist(data_45$Chl, xlim = c(0, 4))
Chl <- log((data_45$Chl*5 +1))
summary(Chl)
hist(Chl, xlim = c(0, 4))
par(mfrow=c(1,1))

data_45$Chl_LOG <- Chl



# 2. Creazione Training e Validation Set (Ultima Settimana) --------------------
# Identifichiamo la data massima
num_days2predict = 3

max_date <- max(data_45$Date)
cutoff_date <- max_date - days(num_days2predict)

cat("Data finale nel dataset:", as.character(max_date), "\n")
cat("Data di taglio (inizio validation):", as.character(cutoff_date), "\n")

# Splitting Temporale
train_set <- data_45 %>% filter(Date <= cutoff_date)
val_set   <- data_45 %>% filter(Date > cutoff_date)

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
  family = scat(link = "identity") #Gamma(link="log") Gamma(link="log"), #scat
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

dev.off()
# Mappa dei residui medi per zona
ggplot(train_set, aes(x = Lon, y = Lat, color = residui)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Mappa dei Residui (Rosso = Sottostima, Blu = Sovrastima)")

# Calcola l'autocorrelazione dei residui
par(mfrow=c(1,1))
acf(residuals(model_bam), main = "Autocorrelazione dei Residui")






# NUOVO BAM AR -------------------------------------------------------------------

#### --- 1. Feature Engineering: Creazione Lag --- ####
# Ordiniamo per posizione e data per essere sicuri che il 'lag' prenda il giorno prima
data_lagged <- data_45 %>%
  arrange(Lon, Lat, Date) %>%
  group_by(Lon, Lat) %>%
  mutate(
    # Variabili Autoregressive (Valori del giorno precedente)
    Chl_lag1      = lag(Chl_LOG, 1),
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
cat("Nuovo Validation Set:", nrow(val_set), "righe\n")


#### --- 3. Addestramento Modello Autoregressivo --- ####

# weights_vec <- (exp(train_set$Chl*.3) - 1)*.5
# summary(weights_vec)
# summary(train_set$Chl)

model_bam_ar <- bam(
  Chl ~ 
    # --- MODIFICA FONDAMENTALE ---
    dow +                    # Effetto base (intercetta) per ogni giorno
    s(Chl_lag1, k=15, by = dow) +  # Curva della Chl specifica per ogni giorno
    # -----------------------------
  
  # 2. Autoregressione sulle Covariate Ambientali (commentate nel tuo codice)
  # s(Salinity_lag1, k=8) + 
  # ...
  s(Temp_lag1, k=15, bs="cc") +
  
  # 3. Componenti Spazio-Temporali
  te(Lon, Lat, doy, d=c(2,1), bs =c("tp", "cc"), k=c(10, 20)),
  
  data = train_set,
  method = "fREML",
  discrete = TRUE,
  nthreads = n_cores,
  chunk.size = 10000,
  #weights = weights_vec,
  family = scat(link = "identity") 
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

dev.off()
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





# ARIMA - predizione ricorsiva -------------------------------------------------------------------

# --- 1. Caricamento Librerie ---
library(forecast)    # Il motore principale per ARIMA
library(tseries)     # Per i test di stazionarietà
library(tidyverse)   # Per manipolazione dati
library(imputeTS)    # FONDAMENTALE: Per riempire i buchi nei dati della clorofilla

# coordinate bibione: 45.630, 13.045
pos_bibione <- posizioni_target %>% 
  filter(near(Lat, 45.630, tol = .03), near(Lon, 13.044, tol = .03))
data_bibione <- data_45 %>%
  inner_join(pos_bibione, by = c("Lat", "Lon"))

# --- 0. Pre-processing e Ordinamento ---
# È cruciale che i dati siano ordinati per tempo prima di tagliare
# Assumiamo che la colonna della data si chiami 'Date' o 'time' (adattalo al tuo df)
data_bibione <- data_bibione %>% 
  arrange(Date) 

# Definiamo le dimensioni dei set
n_valid <- 7   # Ultima settimana
n_train <- 30  # 30 giorni precedenti

# Calcoliamo gli indici di taglio
tot_rows <- nrow(data_bibione)
start_valid <- tot_rows - n_valid + 1
end_train <- start_valid - 1
start_train <- end_train - n_train + 1

# Controlliamo di avere abbastanza dati
if (start_train < 1) {
  stop("Errore: Il dataset è troppo corto per prendere 30 giorni di training + 7 di validation.")
}

# --- 1. Split Training e Validation ---
train_set <- data_bibione[start_train:end_train, ]
valid_set <- data_bibione[start_valid:tot_rows, ]

print(paste("Training set:", min(train_set$Date), "al", max(train_set$Date), 
            "(", nrow(train_set), "giorni)"))
print(paste("Validation set:", min(valid_set$Date), "al", max(valid_set$Date), 
            "(", nrow(valid_set), "giorni)"))

# --- 2. Creazione Time Series (Training) ---
# Creiamo la TS solo sui 30 giorni di training
ts_train <- ts(train_set$Chl, frequency = 7)

# --- 3. Addestramento Modello (Auto-ARIMA) ---
fit <- auto.arima(ts_train, 
                  seasonal = TRUE, 
                  stepwise = FALSE, 
                  approximation = FALSE,
                  trace = FALSE) # Metto false per pulizia, mettilo a TRUE se vuoi vedere i log

print("Modello selezionato:")
print(summary(fit))

# --- 4. Previsione (Forecasting) ---
# Prevediamo per la lunghezza del validation set (h = 7)
forecast_result <- forecast(fit, h = n_valid)

# --- 5. Confronto e Validazione ---

# Creiamo un dataframe per il confronto numerico
confronto <- data.frame(
  Data = valid_set$Date,
  Reale = valid_set$Chl,
  Previsto = as.numeric(forecast_result$mean),
  Lower_95 = as.numeric(forecast_result$lower[,2]), # Intervallo confidenza 95%
  Upper_95 = as.numeric(forecast_result$upper[,2])
)

print("Confronto puntuale (Ultimi 7 giorni):")
print(confronto)

# Calcolo metriche di errore (RMSE, MAE)
# accuracy() confronta automaticamente la previsione con i dati reali (x)
metriche <- accuracy(forecast_result, valid_set$Chl)
print("Metriche di Accuratezza (Guarda la riga 'Test set'):")
print(metriche)

# --- 6. Visualizzazione Grafica ---
# Uniamo tutto per un grafico chiaro con ggplot2

# Creiamo un df unico per il plot
plot_data <- bind_rows(
  train_set %>% select(Date, Chl) %>% mutate(Tipo = "Training"),
  valid_set %>% select(Date, Chl) %>% mutate(Tipo = "Reale (Validation)")
)

# Aggiungiamo le previsioni al grafico
ggplot() +
  # Linea dei dati storici (Training)
  geom_line(data = plot_data %>% filter(Tipo == "Training"), 
            aes(x = Date, y = Chl, color = "Training"), size = 1) +
  
  # Linea dei dati reali futuri (Validation)
  geom_line(data = plot_data %>% filter(Tipo == "Reale (Validation)"), 
            aes(x = Date, y = Chl, color = "Reale (Target)"), size = 1) +
  
  # Linea della Previsione
  geom_line(data = confronto, 
            aes(x = Data, y = Previsto, color = "Previsione Modello"), size = 1, linetype = "dashed") +
  
  # Area di confidenza (l'ombra grigia)
  geom_ribbon(data = confronto, 
              aes(x = Data, ymin = Lower_95, ymax = Upper_95), 
              fill = "blue", alpha = 0.1) +
  
  labs(title = "Predizione Clorofilla-a: Confronto Training vs Validation",
       subtitle = "Addestrato su 30 giorni, testato sugli ultimi 7",
       y = "Chl-a (mg/m3)", x = "Data",
       color = "Legenda") +
  scale_color_manual(values = c("Training" = "black", 
                                "Reale (Target)" = "green", 
                                "Previsione Modello" = "red")) +
  theme_minimal()



# conformal prediction ----------------------------------------------------


# --- 1. Preparazione dei Dati (Split in 3 parti) ---
# Diciamo che hai 100 giorni totali.
# Giorni 1-60: Training (Addestramento modello)
# Giorni 61-90: Calibration (Calcolo soglia di errore)
# Giorni 91-100: Test (Previsione futura reale)

n_total <- nrow(data_bibione)
n_test <- 10        # Ultimi 10 giorni
n_calib <- 30       # 30 giorni per calibrare l'intervallo
n_train <- 30

# Creiamo i 3 dataset

test_data  <- data_bibione[(n_total - n_test + 1):n_total, ]
calib_data <- data_bibione[(n_total - n_test - n_train - n_calib+1):(n_total - n_test - n_train), ]
train_data <- data_bibione[(n_total - n_test - n_train+1):(n_total - n_test), ]

# --- 2. Addestramento (Training) ---
# Usiamo un modello ARIMA, ma potresti usare qualsiasi cosa (XGBoost, Regressione, ecc.)
model <- auto.arima(ts(train_data$Chl, frequency = 7))

# --- 3. Calibrazione (Conformal Step) ---
# Prevediamo sul Calibration set
# Nota: Dobbiamo fare forecast iterativo o one-step-ahead per essere onesti
# Per semplicità qui facciamo un forecast su tutta la finestra calibrazione
calib_forecast_obj <- forecast(model, h = n_calib)
calib_pred <- as.numeric(calib_forecast_obj$mean)

# Calcoliamo i "Non-Conformity Scores" (Errori Assoluti)
calib_scores <- abs(calib_data$Chl - calib_pred)

# Scegliamo la confidenza desiderata (es. 90%)
alpha <- 0.10 
# Calcoliamo il quantile (1 - alpha) degli errori
# Questo Q è la "distanza di sicurezza" che copre il 90% degli errori passati
q_hat <- quantile(calib_scores, probs = (1 - alpha), names = FALSE)

# Correzione per campioni piccoli (opzionale ma raccomandata: (n+1)/n)
q_hat_adjusted <- q_hat * (1 + 1/length(calib_scores))

print(paste("Margine di errore Conformal (Q): +/-", round(q_hat_adjusted, 3)))

# --- 4. Previsione Futura (Test Set) ---
# Prevediamo i valori futuri
test_forecast_obj <- forecast(model, h = n_test)
test_pred <- as.numeric(test_forecast_obj$mean)

# --- 5. Costruzione Intervalli Conformal ---
# Costruiamo l'intervallo sommando e sottraendo Q
conformal_intervals <- data.frame(
  Data = test_data$Date,
  Reale = test_data$Chl,
  Predizione = test_pred,
  Lower_CP = test_pred - q_hat_adjusted,
  Upper_CP = test_pred + q_hat_adjusted
)

# --- 6. Visualizzazione ---
ggplot(conformal_intervals, aes(x = Data)) +
  geom_line(aes(y = Reale, color = "Reale"), size = 1) +
  geom_line(aes(y = Predizione, color = "Predizione"), linetype = "dashed") +
  # L'intervallo di confidenza Conformal
  geom_ribbon(aes(ymin = Lower_CP, ymax = Upper_CP), fill = "purple", alpha = 0.2) +
  labs(title = "Previsione con Intervalli Conformal Prediction (90%)",
       subtitle = "La fascia viola garantisce la copertura indipendentemente dalla distribuzione",
       y = "Clorofilla") +
  theme_minimal()

hist(data_bibione$Chl)



# XGBoost - predizione a un giorno!!! --------------------------------------------------------------------

library(xgboost)
library(caret) # Utile per creare dummy variables o split, qui lo usiamo per visualizzare

# --- 1. Feature Engineering (Il passo più importante) ---
# XGBoost non capisce le date. Dobbiamo trasformare il tempo in numeri.
# E dobbiamo dirgli esplicitamente qual era la clorofilla ieri (Lag).

data_ml <- data_45 %>%
  arrange(Date) %>% # Fondamentale ordinare!
  group_by(Lat, Lon) %>% # Se hai più punti, calcola i lag separatamente per punto!
  mutate(
    # A. Variabili Temporali Cicliche
    DOY = yday(Date),           # Giorno dell'anno (1-365) -> Cattura stagionalità
    Month = month(Date),        # Mese
    Year = year(Date),
    
    # B. Lag Features (Memoria storica)
    # "Quanto era la clorofilla ieri?"
    Chl_Lag1 = lag(Chl, 1),     
    # "Quanto era 3 giorni fa?"
    Chl_Lag3 = lag(Chl, 3),
    # "Quanto era una settimana fa?"
    Chl_Lag7 = lag(Chl, 7),
    
    # C. Rolling Statistics (Tendenze recenti)
    # Media mobile degli ultimi 3 giorni (richiede pacchetto zoo, o manuale)
    # Qui usiamo una media semplice dei lag per non caricare troppe librerie
    Chl_Mean_3d = (Chl_Lag1 + lag(Chl, 2) + Chl_Lag3) / 3
  ) %>%
  ungroup() %>%
  # Rimuoviamo le prime righe che ora hanno NA (a causa dei lag)
  drop_na()

# --- 2. Preparazione dei Dati (Train / Test Split) ---
# Non fare sample random! Taglia per data.
data_ml <- data_ml %>% arrange(Date)

# Usiamo gli ultimi 15 giorni come Test
cutoff_date <- max(data_ml$Date) - 15
cutoff_training_date <- cutoff_date - 30

train_set <- data_ml %>% filter(Date <= cutoff_date)
test_set  <- data_ml %>% filter(Date > cutoff_date)

# Definiamo le colonne che il modello userà per imparare (Features)
# NOTA: NON includiamo 'Date' o 'Chl' (target) qui.
features <- c("Lat", "Lon", "DOY", "Month", 
              "Chl_Lag1", "Chl_Lag3", "Chl_Lag7", "Chl_Mean_3d")

# Creiamo le matrici ottimizzate per XGBoost (xgb.DMatrix)
dtrain <- xgb.DMatrix(data = as.matrix(train_set[, features]), 
                      label = train_set$Chl)

dtest  <- xgb.DMatrix(data = as.matrix(test_set[, features]), 
                      label = test_set$Chl)

# --- 3. Addestramento del Modello ---
# Parametri base (da ottimizzare in seguito con Cross Validation)
params <- list(
  objective = "reg:squarederror", # Vogliamo minimizzare l'errore numerico
  eta = 0.1,                      # Learning rate (più basso è più preciso ma lento)
  max_depth = 6,                  # Profondità degli alberi (se troppo alto -> overfitting)
  subsample = 0.8,                # Usa l'80% dei dati per ogni albero (evita overfitting)
  colsample_bytree = 0.8          # Usa l'80% delle colonne per ogni albero
)

# Addestriamo
print("Addestramento in corso...")
model_xgb <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 500,                  # Numero di alberi
  watchlist = list(train = dtrain, test = dtest),
  print_every_n = 50,
  early_stopping_rounds = 20      # Si ferma se non migliora per 20 giri
)

# --- 4. Predizione e Valutazione ---
# Prevediamo sul test set
preds <- predict(model_xgb, dtest)

# Uniamo i risultati per vedere come è andata
results <- test_set %>%
  select(Date, Lat, Lon, Chl) %>%
  mutate(Predicted = preds,
         Error = Chl - Predicted)

# Calcoliamo RMSE
rmse <- sqrt(mean(results$Error^2))
print(paste("RMSE Finale:", round(rmse, 3)))

# --- 5. Feature Importance (Interpretazione) ---
# Quali variabili sono state più utili?
importance_matrix <- xgb.importance(feature_names = features, model = model_xgb)
xgb.plot.importance(importance_matrix)

# --- 6. Grafico Previsione ---
# Prendiamo un punto specifico per il grafico (es. il primo Lat/Lon disponibile)
one_location <- results %>% 
  inner_join(pos_bibione, by = c("Lat", "Lon"))


ggplot(one_location, aes(x = Date)) +
  geom_line(aes(y = Chl, color = "Reale"), size = 1) +
  geom_line(aes(y = Predicted, color = "XGBoost"), linetype = "dashed", size = 1) +
  labs(title = "XGBoost - predizione a 1 giorno: Reale vs Predetto", 
       subtitle = paste("RMSE:", round(rmse, 3))) +
  theme_minimal()


