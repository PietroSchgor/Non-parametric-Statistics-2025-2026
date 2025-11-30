# 1. Installazione e Caricamento Librerie --------------------------------------

setwd("~/uni/2025-2026/non param/progetto/clorofilla/Non-parametric-Statistics-2025-2026")



if(!require(tidyverse)) install.packages("tidyverse")
if(!require(mgcv)) install.packages("mgcv")      # Il motore per la regressione non parametrica (GAM)
if(!require(lubridate)) install.packages("lubridate")
if(!require(viridis)) install.packages("viridis") # Per i colori dei grafici

library(tidyverse)
library(mgcv)
library(lubridate)
library(viridis)

# 2. Caricamento e Preparazione Dati -------------------------------------------
# Leggiamo il dataset
df <- read_csv("dataset/dataset_MERGIATO")

# Preprocessing
data_clean <- df %>%
  mutate(
    Date = ymd(Date),
    # Creiamo variabili temporali numeriche per le curve
    # Day of Year (1-365) cattura la stagionalità ciclica
    doy = yday(Date), 
    # Time Numeric cattura il trend a lungo termine (es. aumento temperature negli anni)
    time_numeric = as.numeric(Date),
    # Creiamo un ID unico per ogni punto griglia (per visualizzazione dopo)
    loc_id = paste(round(Lon,3), round(Lat,3), sep="_")
  ) %>%
  # Rimuoviamo eventuali righe con NA nel target o coordinate
  drop_na(Chl, Lat, Lon)

# Visualizziamo quanti punti unici (curve) abbiamo
n_locations <- n_distinct(data_clean$loc_id)
cat("Numero di stazioni/punti griglia da modellare:", n_locations, "\n")


data_clean <- data_clean[which(data_clean$Lon <= 13), ]
data_clean <- data_clean[which(data_clean$Lat <= 45), ]

# Visualizziamo quanti punti unici (curve) abbiamo
n_locations <- n_distinct(data_clean$loc_id)
cat("Numero di stazioni/punti griglia da modellare:", n_locations, "\n")
# 3. Definizione del Modello Non Parametrico (GAM) -----------------------------
# Qui modelliamo la Clorofilla come una funzione complessa.
# Formula:
# te(Lon, Lat, time_numeric): Interazione Spazio-Tempo. 
#     Modella una curva temporale diversa per ogni coordinata geografica.
#     k=... controlla la "flessibilità" (wiggles) della curva.
# s(doy, bs="cc"): Spline ciclica per la stagionalità annuale (l'inizio anno si raccorda con la fine).
# s(Temp), s(Salinity), ...: Effetti non lineari delle variabili ambientali.
help("gam")
model_gam <- gam(
  Chl ~ 
    # 1. La "Curva" Spazio-Temporale (Lat, Lon interagiscono con il Tempo)
    te(Lon, Lat, time_numeric, d=c(2,1), k=c(20, 10)) +
    
    # 2. Stagionalità ciclica (Giorno dell'anno)
    s(doy, bs = "cc", k = 10) +
    
    # 3. Regressori fisici (Relazioni non lineari)
    s(Temp, k=5) + 
    s(Salinity, k=5) +
    s(MLD, k=5) +
    s(Solar_Flux, k=5),
  
  data = data_clean,
  method = "REML", # Restricted Maximum Likelihood per evitare overfitting
  family = gaussian(link = "identity") # O Gamma(link="log") se Chl è molto asimmetrica
)

# Sommario del modello
summary(model_gam)

# Visualizzazione degli effetti parziali (come le variabili influenzano la Chl)
plot(model_gam, pages = 1, scheme = 1, shade = TRUE, main = "Effetti Non Parametrici")



# ------------------------------------------------------------------------------
# MODELLO OTTIMIZZATO (BAM)
# ------------------------------------------------------------------------------

library(parallel) # Per rilevare i core della CPU

# ... (parte di caricamento dati identica a prima) ...
# Assicurati che data_clean sia caricato come nel codice precedente

# Rileva i core disponibili (lasciamone uno libero per il sistema)
n_cores <- detectCores() - 1

# Differenze rispetto a gam():
# 1. bam(): Big Additive Model, specifico per grandi dataset.
# 2. discrete = TRUE: Discretizza le covariate per velocizzare enormemente l'algebra lineare.
# 3. nthreads = n_cores: Usa il calcolo parallelo (multi-core).
# 4. chunk.size: Gestisce la memoria a blocchi per non saturare la RAM.

model_bam <- bam(
  Chl ~ 
    # Interazione tensoriale (Spazio-Tempo)
    te(Lon, Lat, time_numeric, d=c(2,1), k=c(20, 10)) +
    
    # Stagionalità ciclica
    s(doy, bs = "cc", k = 10) +
    
    # Regressori fisici
    s(Temp, k=5) + 
    s(Salinity, k=5) +
    s(MLD, k=5) +
    s(Solar_Flux, k=5),
  
  data = data_clean,
  method = "fREML",      # fREML è molto più veloce di REML per i grandi dataset
  discrete = TRUE,       # <--- IL VERO BOOST DI PRESTAZIONI
  nthreads = n_cores,    # <--- USA TUTTA LA CPU
  chunk.size = 10000     # Elabora 10k righe alla volta per risparmiare RAM
)

summary(model_bam)

# ... (il resto del codice per le predizioni rimane uguale, usa predict() su model_bam)
# 4. Forecasting (Predizione) --------------------------------------------------
# Per fare forecasting, dobbiamo creare un dataset "futuro".
# Esempio: Prediciamo i prossimi 7 giorni per tutte le location.

# Troviamo l'ultima data nel dataset
last_date <- max(data_clean$Date)

# Generiamo date future
future_dates <- seq(last_date + 1, last_date + 7, by = "day")

# Creiamo il dataframe per la predizione
# NOTA: Per un vero forecast, dovresti avere le previsioni meteo (Temp, Salinity) future.
# Qui, per dimostrazione, usiamo la media delle ultime osservazioni per le variabili fisiche
# (o un approccio 'persistence').
last_obs <- data_clean %>%
  group_by(Lon, Lat) %>%
  slice_max(Date, n=1) %>%
  select(-Date, -doy, -time_numeric, -Chl) # Teniamo le coordinate e le variabili fisiche costanti

future_data <- expand_grid(
  Date = future_dates,
  last_obs
) %>%
  mutate(
    doy = yday(Date),
    time_numeric = as.numeric(Date)
  )

# Eseguiamo la predizione
predictions <- predict(model_bam, newdata = future_data, type = "response", se.fit = TRUE)

future_data <- future_data %>%
  mutate(
    Chl_Pred = predictions$fit,
    Chl_SE = predictions$se.fit, # Errore standard della predizione
    Lower_CI = Chl_Pred - 1.96 * Chl_SE,
    Upper_CI = Chl_Pred + 1.96 * Chl_SE
  )



# 5. Visualizzazione delle Curve Predette --------------------------------------
# Prendiamo 4 location a caso per vedere le "curve" (Storico + Forecast)

sample_locs <- unique(data_clean$loc_id)[1:4]

plot_data <- data_clean %>%
  filter(loc_id %in% sample_locs) %>%
  select(Date, Lon, Lat, Chl) %>%
  mutate(Type = "Osservato") %>%
  bind_rows(
    future_data %>% 
      mutate(loc_id = paste(round(Lon,3), round(Lat,3), sep="_")) %>%
      filter(loc_id %in% sample_locs) %>%
      select(Date, Lon, Lat, Chl = Chl_Pred) %>%
      mutate(Type = "Predetto")
  )

ggplot(plot_data, aes(x = Date, y = Chl, color = Type)) +
  geom_line(size = 1) +
  geom_point(data = plot_data %>% filter(Type == "Osservato"), alpha=0.3, size=0.5) +
  facet_wrap(~Lat + Lon, scales = "free_y", labeller = label_both) +
  theme_minimal() +
  labs(title = "Forecasting delle curve di Clorofilla",
       subtitle = "Regressione Non Parametrica (GAM) per coordinate",
       y = "Chlorophyll-a", x = "Data") +
  scale_color_manual(values = c("black", "red"))



##### check with real data #####



# Troviamo l'ultima data nel dataset
last_date <- max(data_clean$Date)

# Generiamo date future
future_dates <- seq(last_date + 1, last_date + 7, by = "day")

# Creiamo il dataframe per la predizione
# NOTA: Per un vero forecast, dovresti avere le previsioni meteo (Temp, Salinity) future.
# Qui, per dimostrazione, usiamo la media delle ultime osservazioni per le variabili fisiche
# (o un approccio 'persistence').
last_obs <- data_clean %>%
  group_by(Lon, Lat) %>%
  slice_max(Date, n=1) %>%
  select(-Date, -doy, -time_numeric, -Chl) # Teniamo le coordinate e le variabili fisiche costanti

future_data <- expand_grid(
  Date = future_dates,
  last_obs
) %>%
  mutate(
    doy = yday(Date),
    time_numeric = as.numeric(Date)
  )

# Eseguiamo la predizione
predictions <- predict(model_bam, newdata = future_data, type = "response", se.fit = TRUE)

future_data <- future_data %>%
  mutate(
    Chl_Pred = predictions$fit,
    Chl_SE = predictions$se.fit, # Errore standard della predizione
    Lower_CI = Chl_Pred - 1.96 * Chl_SE,
    Upper_CI = Chl_Pred + 1.96 * Chl_SE
  )



# 5. Visualizzazione delle Curve Predette --------------------------------------
# Prendiamo 4 location a caso per vedere le "curve" (Storico + Forecast)

sample_locs <- unique(data_clean$loc_id)[1:4]

plot_data <- data_clean %>%
  filter(loc_id %in% sample_locs) %>%
  select(Date, Lon, Lat, Chl) %>%
  mutate(Type = "Osservato") %>%
  bind_rows(
    future_data %>% 
      mutate(loc_id = paste(round(Lon,3), round(Lat,3), sep="_")) %>%
      filter(loc_id %in% sample_locs) %>%
      select(Date, Lon, Lat, Chl = Chl_Pred) %>%
      mutate(Type = "Predetto")
  )

ggplot(plot_data, aes(x = Date, y = Chl, color = Type)) +
  geom_line(size = 1) +
  geom_point(data = plot_data %>% filter(Type == "Osservato"), alpha=0.3, size=0.5) +
  facet_wrap(~Lat + Lon, scales = "free_y", labeller = label_both) +
  theme_minimal() +
  labs(title = "Forecasting delle curve di Clorofilla",
       subtitle = "Regressione Non Parametrica (GAM) per coordinate",
       y = "Chlorophyll-a", x = "Data") +
  scale_color_manual(values = c("black", "red"))

