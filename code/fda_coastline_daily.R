
setwd("~/uni/2025-2026/non param/progetto/clorofilla/Non-parametric-Statistics-2025-2026")



# 1. Caricamento e Preparazione Dati -------------------------------------------
library(tidyverse)

# data <- read.csv("dataset/dataset_MERGIATO")
# dataset_ordinato_distanza <- read.csv("dataset/bordo_con_distanze.csv")
# data_pos_x_daily <- data %>%
#   inner_join(dataset_ordinato_distanza, by = c("Lat", "Lon"))
# 
# chl_pos_x_daily <- data_pos_x_daily %>%
#   dplyr::select(Date, Lat, Lon, Chl, distanza_cumulata)
# 
# chl_pos_x_daily$Date <- as.Date(chl_pos_x_daily$Date)
# 
# x <- sort(unique(chl_pos_x_daily$distanza_cumulata))
# 
# 
# plot(chl_pos_x_daily$Lon, chl_pos_x_daily$Lat)
# 
# 
# 
# matrice_df <- chl_pos_x_daily %>%
#   dplyr::select(Date, distanza_cumulata, Chl) %>%
#   # Ordina prima i dati per sicurezza
#   arrange(Date, distanza_cumulata) %>%
#   # Trasforma in formato largo
#   pivot_wider(names_from = distanza_cumulata, values_from = Chl)
# 
# # Per convertirla in una vera matrice numerica (togliendo la colonna Date)
# matrice_finale <- as.matrix(matrice_df[,-1])
# rownames(matrice_finale) <- matrice_df$Date
# 
# sum(is.na(matrice_finale))
# # # Salvataggio standard (separatore: virgola, decimale: punto)
# write.csv(matrice_finale, "dataset/matrix_X_fda_daily.csv", row.names = TRUE)


matrice_finale <- read.csv("dataset/matrix_X_fda_daily.csv", row.names = 1)
matrice_finale <- as.matrix(matrice_finale)

dataset_ordinato_distanza <- read.csv("dataset/bordo_con_distanze.csv")
x <- sort(unique(dataset_ordinato_distanza$distanza_cumulata))
# 1. Carica il pacchetto (se non lo hai: install.packages("fields"))
library(fields)
library(viridis)

# 2. Definisci i colori
n_giorni <- nrow(matrice_finale)
colori <- viridis(n_giorni)

# 3. Crea il plot principale (lascia un po' di spazio a destra per la legenda)
par(mar = c(5, 4, 4, 6)) # Aumenta il margine destro (quarto valore)

matplot(x, t(matrice_finale), 
        type = 'l', lty = 1, lwd = 2, col = colori,
        xlab = "Distanza Cumulata", ylab = "Chl")

# 4. Aggiungi la color bar
# Usiamo Date per definire i limiti della barra
image.plot(legend.only = TRUE, 
           zlim = range(as.numeric(rownames(matrice_finale))), 
           col = colori, 
           legend.lab = "day", 
           legend.line = 3)

# fda ---------------------------------------------------------------------
library(fda)

# Definizione del range della distanza (asse x)
range_dist <- range(x)

# Creazione della base di B-splines
# nbasis: numero di funzioni base (regola empirica: numero punti / 4, o in base alla complessità)
# Creazione della base usando i valori reali di x come nodi
basis_obj <- create.bspline.basis(
  rangeval = range(x), 
  breaks = x,           # Fissa i nodi sui tuoi dati reali
  norder = 4                  # Spline cubico (default)
)
# Nota: nbasis sarà ora uguale a: length(x_knots) + norder - 2

# Nota: x sono i tuoi punti di distanza_cumulata
# t(matrice_finale) mette le giorni come curve separate
chl_fd <- Data2fd(argvals = x, y = t(matrice_finale), basisobj = basis_obj)

# 1. Plot principale delle curve
plot(chl_fd, col = viridis(nrow(matrice_finale)), lty = 1, lwd = 2,
     xlab = "Distanza Cumulata", ylab = "Chl (Smooth)",
     main = "FDA: B-spline con Knots evidenziati")

# 2. Estrai i nodi dall'oggetto base
knots <- basis_obj$params # Nodi interni
range_knots <- basis_obj$rangeval # Nodi ai bordi (inizio e fine)
tutti_i_knots <- c(range_knots[1], knots, range_knots[2])

# 3. Aggiungi i nodi al grafico
# Usiamo linee tratteggiate grigie per non appesantire il disegno
abline(v = tutti_i_knots, col = "gray80", lty = 3)

# Opzionale: aggiungi dei piccoli segmenti alla base per renderli più visibili
rug(tutti_i_knots, col = "red", lwd = 1.5)


# Dividi lo schermo in 2 colonne
par(mfrow = c(1, 2))

# Grafico 1: Dati Reali
matplot(x, t(matrice_finale), 
        type = "l", lty = 1, lwd = 1, col = colori,
        xlab = "Distanza", ylab = "Chl", main = "Dati Originali (Raw)")

# Grafico 2: FDA Fit
plot(chl_fd, col = colori, lty = 1, lwd = 2,
     xlab = "Distanza", ylab = "Chl", main = "FDA Smooth (B-splines)")
abline(v = tutti_i_knots, col = "gray90", lty = 3)
rug(tutti_i_knots, col = "red", lwd = 1.5)
# Ripristina il layout a 1x1
par(mfrow = c(1, 1))


#### FPCA ####

# Calcoliamo le prime 4 componenti principali
n_comp <- 4
chl_pca <- pca.fd(chl_fd, nharm = n_comp)

# Estraiamo la proporzione di varianza spiegata
var_spiegata <- chl_pca$varprop
names(var_spiegata) <- paste0("PC", 1:n_comp)
print(round(var_spiegata, 3))

# Plot cumulativo
plot(cumsum(var_spiegata), type = "b", pch = 19, 
     xlab = "Numero di PC", ylab = "Varianza Cumulata Spiegata",
     main = "Scree Plot FPCA")

# Visualizziamo le prime 2 componenti principali
par(mfrow = c(2, 2))
plot(chl_pca, nx = 100, cex.main = 0.8) 
par(mfrow = c(1, 1))



# prediction model --------------------------------------------------------

# Prendiamo gli scores delle PC (es. le prime 3)
scores_all <- chl_pca$scores[, 1:4]
n_tot <- nrow(scores_all)
n_val <- 7  # Numero di giorni per la validazione

# Split: Training (passato) e Validation (futuro recente)
train_set <- scores_all[1:(n_tot - n_val), ]
val_set   <- scores_all[(n_tot - n_val + 1):n_tot, ]

#### VAR ####
library(vars)

# Fit sul training set
fit_var_train <- VAR(train_set, p = 7, type = "both") # p=2 come esempio di lag

# Previsione per le 5 settimane di validazione
pred_backtest <- predict(fit_var_train, n.ahead = n_val)

par(mfrow = c(4, 1), mar = c(4, 4, 2, 1))

for(i in 1:4) {
  nome_pc <- paste0("PC", i)
  plot(1:n_tot, scores_all[, i], type = "l", col = "black", lwd = 1,
       main = paste("Backtesting", nome_pc), xlab = "Settimana", ylab = "Score")
  
  # Aggiungiamo la previsione in rosso
  lines((n_tot - n_val + 1):n_tot, pred_backtest$fcst[[i]][, "fcst"], 
        col = "red", lwd = 2, type = "b", pch = 19)
  
  # Evidenziamo l'area di validazione
  abline(v = n_tot - n_val + 0.5, col = "blue", lty = 2)
}
#### Plot curve predette
# Crea la matrice estraendo la colonna "fcst" da ogni elemento della lista
scores_predetti <- sapply(pred_backtest$fcst, function(x) x[, "fcst"])

# Se per caso i dati sono trasposti, forziamo la struttura [n_val x n_pcs]
if(nrow(scores_predetti) != n_val) scores_predetti <- t(scores_predetti)

# Verifica la struttura
print(dim(scores_predetti)) # Dovrebbe essere [7 x 4]


# Coefficienti delle armoniche e della media
coef_armoniche <- chl_pca$harmonics$coefs[, 1:4] # Assicurati di prendere solo le prime 4
coef_media <- as.vector(chl_pca$meanfd$coefs)

# Ricostruzione (Matrice di coefficienti per le 7 nuove curve)
coef_predetti <- coef_armoniche %*% t(scores_predetti)
coef_finali <- coef_predetti + coef_media

# Creazione dell'oggetto funzionale
curve_predette_fd <- fd(coef_finali, basis_obj)

library(viridis)

# Numero di step ahead (pari alle righe di scores_predetti)
n_val <- nrow(scores_predetti)

# Definiamo i colori per essere sicuri che corrispondano tra plot e legenda
colori_predizione <- viridis(n_val)

# Plot delle curve funzionali predette
par(mfrow = c(1,1))
plot(curve_predette_fd, 
     col = colori_predizione, 
     lty = 1, 
     lwd = 2,
     main = "Profili Spaziali Predetti (Set di Validazione)",
     xlab = "Distanza / Posizione", 
     ylab = "Valore Funzionale")

# Creazione delle etichette per la legenda
etichette_legenda <- paste("Step +", 1:n_val)

# Aggiunta della legenda
# 'topright' la posiziona in alto a destra, puoi cambiarla in 'topleft' se copre i dati
legend("topright", 
       legend = etichette_legenda, 
       col = colori_predizione, 
       lwd = 2, 
       lty = 1, 
       cex = 0.8,      # Dimensione del testo
       bty = "n",      # Rimuove il bordo della legenda
       title = "Orizzonte Previsivo")


par(mfrow = c(4, 2), mar = c(4, 4, 2, 1)) # mar riduce i margini per far stare tutto

for (i in 1:n_val) {
  # Plot della curva reale (Nera)
  plot(chl_fd[(n_tot - n_val + i)], 
       col = "black", 
       lwd = 2, 
       ylim = c(0, 5), 
       main = paste("Validazione: Step", i))
  
  # Aggiunta della curva predetta (Rossa tratteggiata)
  plot(curve_predette_fd[i], 
       col = "red", 
       lty = 2, 
       lwd = 2, 
       add = TRUE)
  
  # Aggiunta della legenda specifica per ogni riquadro
  legend("topright", 
         legend = c("Reale", paste("Predetto +", i)), 
         col = c("black", "red"), 
         lty = c(1, 2), 
         lwd = 2, 
         cex = 0.7,      # Dimensione ridotta per stare nei piccoli riquadri
         bty = "n")      # Rimuove il bordo per pulizia visiva
}

# Definiamo le curve reali per il periodo di validazione/test
curve_reali <- chl_fd[(n_tot - n_val + 1):n_tot]


calcola_rmse_giornaliero <- function(curve_reali, curve_predette) {
  # Calcola la differenza tra gli oggetti funzionali
  diff_fd <- curve_reali - curve_predette
  
  # Calcola l'integrale del quadrato della differenza (L2 norm squared)
  # inprod() calcola la matrice degli integrali incrociati. 
  # A noi interessano solo gli elementi sulla diagonale (giorno 1 vs giorno 1, ecc.)
  integrali_quadrati <- diag(inprod(diff_fd, diff_fd))
  
  # Restituisce la radice quadrata (RMSE)
  return(sqrt(integrali_quadrati))
}


err_var   <- calcola_rmse_giornaliero(curve_reali, curve_predette_fd)
# 2. Plot dell'errore
plot(1:length(err_var), err_var, 
     type = "b",               # Linea con punti
     pch = 19,                 # Pallino pieno
     col = "red",              # Colore rosso per il VAR
     lwd = 2,                  # Spessore linea
     xlab = "Giorno di Previsione (Step-ahead)", 
     ylab = "RMSE Funzionale",
     main = "Andamento Errore di Previsione - Modello VAR",
     panel.first = grid(),
     ylim = c(0, 3))    

# Opzionale: Aggiungi i valori numerici sopra i punti per precisione
text(1:length(err_var), err_var, labels = round(err_var, 4), pos = 3, cex = 0.8, col = "darkred")

par(mfrow=c(1,1))

#### SARIMA ####
library(forecast)


# --- 1. Preparazione Scores ---
scores_all <- chl_pca$scores[, 1:4]
n_tot <- nrow(scores_all)
n_val <- 7  

train_set <- scores_all[1:(n_tot - n_val), ]
val_set   <- scores_all[(n_tot - n_val + 1):n_tot, ]

# --- 2. Predizione con SARIMA ---
# Creiamo una matrice vuota per ospitare le previsioni delle 4 PC
scores_predetti <- matrix(NA, nrow = n_val, ncol = 4)
colnames(scores_predetti) <- paste0("PC", 1:4)

for(i in 1:4) {
  # auto.arima trova il miglior modello (inclusa stagionalità se presente)
  # Se conosci la frequenza (es. 365 per dati giornalieri), usa: ts(train_set[,i], frequency = ...)
  fit_sarima <- auto.arima(train_set[, i], seasonal = TRUE)
  
  # Predizione n-step ahead
  pred_obj <- forecast(fit_sarima, h = n_val)
  scores_predetti[, i] <- as.numeric(pred_obj$mean)
}

# Verifica dimensioni: deve essere [7 x 4]
print(dim(scores_predetti))

# --- 3. Ricostruzione Funzionale ---
# Coefficienti delle armoniche e della media
coef_armoniche <- chl_pca$harmonics$coefs[, 1:4] 
coef_media <- as.vector(chl_pca$meanfd$coefs)

# Ricostruzione (Matrice di coefficienti per le 7 nuove curve)
# scores_predetti è [7 x 4], coef_armoniche è [nbasis x 4]
coef_predetti <- coef_armoniche %*% t(scores_predetti)
coef_finali <- coef_predetti + coef_media

# Creazione dell'oggetto funzionale
curve_predette_fd <- fd(coef_finali, basis_obj)

# --- 4. Plot Backtesting Scores ---
par(mfrow = c(4, 1), mar = c(4, 4, 2, 1))
for(i in 1:4) {
  plot(1:n_tot, scores_all[, i], type = "l", col = "black", 
       main = paste("SARIMA Backtesting PC", i), xlab = "Tempo", ylab = "Score")
  lines((n_tot - n_val + 1):n_tot, scores_predetti[, i], 
        col = "blue", lwd = 2, type = "b", pch = 19) # Blu per SARIMA
  abline(v = n_tot - n_val + 0.5, col = "red", lty = 2)
}

# --- 5. Plot Confronto Profili Funzionali ---
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
for (i in 1:n_val) {
  # Curva reale
  plot(chl_fd[(n_tot - n_val + i)], 
       col = "black", lwd = 2, ylim = c(0, 5), 
       main = paste("SARIMA: Step", i))
  
  # Curva predetta
  plot(curve_predette_fd[i], 
       col = "blue", lty = 2, lwd = 2, add = TRUE)
  
  legend("topright", 
         legend = c("Reale", "SARIMA"), 
         col = c("black", "blue"), lty = c(1, 2), 
         cex = 0.7, bty = "n")
}

err_var   <- calcola_rmse_giornaliero(curve_reali, curve_predette_fd)
# 2. Plot dell'errore
plot(1:length(err_var), err_var, 
     type = "b",               # Linea con punti
     pch = 19,                 # Pallino pieno
     col = "blue",              # Colore rosso per il VAR
     lwd = 2,                  # Spessore linea
     xlab = "Giorno di Previsione (Step-ahead)", 
     ylab = "RMSE Funzionale",
     main = "Andamento Errore di Previsione - Modello VAR",
     panel.first = grid(),
     ylim = c(0, 3))     # Aggiunge una griglia sotto i punti

# Opzionale: Aggiungi i valori numerici sopra i punti per precisione
text(1:length(err_var), err_var, labels = round(err_var, 4), pos = 3, cex = 0.8, col = "darkred")

par(mfrow=c(1,1))

#### XGBoost ####
library(xgboost)

# --- 1. Preparazione Dati ---
scores_all <- chl_pca$scores[, 1:4]
n_tot <- nrow(scores_all)
n_val <- 7
n_lag <- 7 # Usiamo gli ultimi 7 giorni per prevedere il prossimo

# Funzione per creare il dataset con lag
prepare_xgb_data <- function(serie, lags) {
  X <- embed(serie, lags + 1)
  return(list(y = X[, 1], x = X[, -1]))
}

# --- 2. Training e Predizione Ricorsiva ---
scores_predetti <- matrix(NA, nrow = n_val, ncol = 4)

for(p in 1:4) {
  serie_pc <- scores_all[, p]
  train_data <- serie_pc[1:(n_tot - n_val)]
  
  # Preparazione training set
  data_split <- prepare_xgb_data(train_data, n_lag)
  
  # Fit del modello XGBoost
  bst <- xgboost(data = data_split$x, 
                 label = data_split$y, 
                 nrounds = 50, 
                 objective = "reg:squarederror",
                 verbose = 0)
  
  # Predizione Ricorsiva
  # Partiamo dagli ultimi dati del training per prevedere il primo step di validazione
  current_input <- rev(tail(train_data, n_lag)) 
  
  for(s in 1:n_val) {
    # Predizione dello step successivo
    pred <- predict(bst, matrix(current_input, nrow = 1))
    scores_predetti[s, p] <- pred
    
    # Aggiornamento input: togliamo il più vecchio, aggiungiamo la nuova predizione
    current_input <- c(pred, current_input[-n_lag])
  }
}

# --- 3. Ricostruzione Funzionale ---
coef_armoniche <- chl_pca$harmonics$coefs[, 1:4] 
coef_media <- as.vector(chl_pca$meanfd$coefs)

coef_predetti <- coef_armoniche %*% t(scores_predetti)
coef_finali <- coef_predetti + coef_media
curve_predette_fd <- fd(coef_finali, basis_obj)

# --- 4. Plot Backtesting Scores ---
par(mfrow = c(4, 1), mar = c(4, 4, 2, 1))
for(i in 1:4) {
  plot(1:n_tot, scores_all[, i], type = "l", col = "black", 
       main = paste("XGBoost Backtesting PC", i))
  lines((n_tot - n_val + 1):n_tot, scores_predetti[, i], 
        col = "darkgreen", lwd = 2, type = "b", pch = 19)
  abline(v = n_tot - n_val + 0.5, col = "blue", lty = 2)
}

# --- 5. Plot Confronto Profili Funzionali ---
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
for (i in 1:n_val) {
  plot(chl_fd[(n_tot - n_val + i)], col = "black", lwd = 2, ylim = c(0, 5), 
       main = paste("XGBoost: Step", i))
  plot(curve_predette_fd[i], col = "darkgreen", lty = 2, lwd = 2, add = TRUE)
  
  legend("topright", legend = c("Reale", "XGBoost"), 
         col = c("black", "darkgreen"), lty = c(1, 2), cex = 0.7, bty = "n")
}



calcola_rmse_giornaliero <- function(curve_reali, curve_predette) {
  # Calcola la differenza tra gli oggetti funzionali
  diff_fd <- curve_reali - curve_predette
  
  # Calcola l'integrale del quadrato della differenza (L2 norm squared)
  # inprod() calcola la matrice degli integrali incrociati. 
  # A noi interessano solo gli elementi sulla diagonale (giorno 1 vs giorno 1, ecc.)
  integrali_quadrati <- diag(inprod(diff_fd, diff_fd))
  
  # Restituisce la radice quadrata (RMSE)
  return(sqrt(integrali_quadrati))
}


err_var   <- calcola_rmse_giornaliero(curve_reali, curve_predette_fd)
# 2. Plot dell'errore
plot(1:length(err_var), err_var, 
     type = "b",               # Linea con punti
     pch = 19,                 # Pallino pieno
     col = "darkgreen",              # Colore rosso per il VAR
     lwd = 2,                  # Spessore linea
     xlab = "Giorno di Previsione (Step-ahead)", 
     ylab = "RMSE Funzionale",
     main = "Andamento Errore di Previsione - Modello VAR",
     panel.first = grid(),
     ylim = c(0, 3))    

# Opzionale: Aggiungi i valori numerici sopra i punti per precisione
text(1:length(err_var), err_var, labels = round(err_var, 4), pos = 3, cex = 0.8, col = "darkred")

par(mfrow=c(1,1))
#### MULTIVARIATE XGBoost ####
library(xgboost)
library(viridis)

# --- 1. Preparazione Dati Multivariati ---
scores_all <- chl_pca$scores[, 1:4]
n_tot <- nrow(scores_all)
n_val <- 7
n_lag <- 7 

# Funzione per creare il dataset di training (X = tutti i lag di tutte le PC)
prepare_multivariate_train <- function(data, lags) {
  X_list <- list()
  for(l in 1:lags) {
    # Lag l della matrice intera
    X_list[[l]] <- data[(lags - l + 1):(nrow(data) - l), ]
  }
  X <- do.call(cbind, X_list)
  Y <- data[(lags + 1):nrow(data), ]
  return(list(X = X, Y = Y))
}

# --- 2. Training Unico (No Rolling) ---
train_data <- scores_all[1:(n_tot - n_val), ]
d_train <- prepare_multivariate_train(train_data, n_lag)

# Alleniamo 4 modelli (uno per ogni PC target), ma ognuno vede i lag di tutte le PC
models_list <- list()
for (p in 1:4) {
  models_list[[p]] <- xgboost(data = d_train$X, 
                              label = d_train$Y[, p], 
                              nrounds = 50, 
                              objective = "reg:squarederror", 
                              verbose = 0)
}

# --- 3. Predizione Ricorsiva (Orizzonte 7 giorni) ---
scores_predetti_xgb <- matrix(NA, nrow = n_val, ncol = 4)

# L'input iniziale sono gli ultimi n_lag giorni del training set
# Lo trasformiamo in un vettore riga che contiene [Lag1_PC1..4, Lag2_PC1..4, ...]
current_history <- train_data[(nrow(train_data) - n_lag + 1):nrow(train_data), ]

for (s in 1:n_val) {
  # Prepariamo il vettore di input invertendo l'ordine per avere i lag corretti (Lag 1, Lag 2...)
  # La struttura deve essere identica a quella usata in d_train$X
  input_row <- matrix(as.vector(t(current_history[nrow(current_history):1, ])), nrow = 1)
  
  # Prediciamo le 4 PC per lo step s
  step_preds <- numeric(4)
  for (p in 1:4) {
    step_preds[p] <- predict(models_list[[p]], input_row)
  }
  
  # Salviamo il risultato
  scores_predetti_xgb[s, ] <- step_preds
  
  # Aggiorniamo la storia: togliamo il giorno più vecchio e aggiungiamo la previsione appena fatta
  current_history <- rbind(current_history[-1, ], step_preds)
}

# --- 4. Ricostruzione Funzionale ---
coef_pred_xgb <- chl_pca$harmonics$coefs[, 1:4] %*% t(scores_predetti_xgb)
coef_finali_xgb <- coef_pred_xgb + as.vector(chl_pca$meanfd$coefs)
curve_xgb_fd <- fd(coef_finali_xgb, basis_obj)

# --- 5. Plot Risultati ---
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
for (i in 1:n_val) {
  plot(chl_fd[n_tot - n_val + i], col = "black", lwd = 2, ylim = c(0, 5),
       main = paste("Multi-XGBoost Statico: Giorno", i))
  plot(curve_xgb_fd[i], col = "darkgreen", lty = 2, lwd = 2, add = TRUE)
  
  legend("topright", legend = c("Reale", "XGBoost Statico"), 
         col = c("black", "darkgreen"), lty = c(1, 2), cex = 0.7, bty = "n")
}

err_var   <- calcola_rmse_giornaliero(curve_reali, curve_xgb_fd)
# 2. Plot dell'errore
plot(1:length(err_var), err_var, 
     type = "b",               # Linea con punti
     pch = 19,                 # Pallino pieno
     col = "darkgreen",              # Colore rosso per il VAR
     lwd = 2,                  # Spessore linea
     xlab = "Giorno di Previsione (Step-ahead)", 
     ylab = "RMSE Funzionale",
     main = "Andamento Errore di Previsione - Modello VAR",
     panel.first = grid(),
     ylim = c(0, 3))    

# Opzionale: Aggiungi i valori numerici sopra i punti per precisione
text(1:length(err_var), err_var, labels = round(err_var, 4), pos = 3, cex = 0.8, col = "darkred")

par(mfrow=c(1,1))

#### VECM ####

library(tsDyn)
library(fda)

# --- 1. Preparazione Dati ---
scores_all <- chl_pca$scores[, 1:4]
n_tot <- nrow(scores_all)
n_val <- 7
n_train <- n_tot - n_val

train_set <- scores_all[1:n_train, ]
val_set   <- scores_all[(n_train + 1):n_tot, ]

# --- 2. Fit del Modello VECM ---
# lags: numero di ritardi (equivalente a p-1 nel VAR)
# r: rango di cointegrazione (quante relazioni di lungo periodo esistono)
# Iniziamo con r=1 (ipotesi conservativa di un equilibrio comune)
fit_vecm <- VECM(train_set, lag = 7, r = 3, estim = "ML")

# --- 3. Predizione Ricorsiva ---
# La funzione predict di tsDyn genera automaticamente la previsione n-step ahead
pred_vecm_obj <- predict(fit_vecm, n.ahead = n_val)

# Estraiamo i punteggi predetti (matrice 7x4)
scores_predetti_vecm <- as.matrix(pred_vecm_obj)

# --- 4. Ricostruzione Funzionale ---
coef_armoniche <- chl_pca$harmonics$coefs[, 1:4]
coef_media <- as.vector(chl_pca$meanfd$coefs)

# coef_predetti = [nbasis x 4] * [4 x 7] = [nbasis x 7]
coef_pred_vecm <- coef_armoniche %*% t(scores_predetti_vecm)
coef_finali_vecm <- coef_pred_vecm + coef_media

curve_vecm_fd <- fd(coef_finali_vecm, basis_obj)

# --- 5. Plot di Validazione ---
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
for (i in 1:n_val) {
  idx_reale <- n_train + i
  plot(chl_fd[idx_reale], col = "black", lwd = 2, ylim = c(0, 5),
       main = paste("VECM: Giorno", i))
  
  # Curva predetta in blu scuro
  plot(curve_vecm_fd[i], col = "darkblue", lty = 2, lwd = 2, add = TRUE)
  
  legend("topright", legend = c("Reale", "VECM"), 
         col = c("black", "darkblue"), lty = c(1, 2), cex = 0.7, bty = "n")
}
err_var   <- calcola_rmse_giornaliero(curve_reali, curve_vecm_fd)
# 2. Plot dell'errore
plot(1:length(err_var), err_var, 
     type = "b",               # Linea con punti
     pch = 19,                 # Pallino pieno
     col = "darkblue",              # Colore rosso per il VAR
     lwd = 2,                  # Spessore linea
     xlab = "Giorno di Previsione (Step-ahead)", 
     ylab = "RMSE Funzionale",
     main = "Andamento Errore di Previsione - Modello VAR",
     panel.first = grid(),
     ylim = c(0, 3))    

# Opzionale: Aggiungi i valori numerici sopra i punti per precisione
text(1:length(err_var), err_var, labels = round(err_var, 4), pos = 3, cex = 0.8, col = "darkred")

par(mfrow=c(1,1))


library(urca)

# Converti in oggetto Time Series (frequenza 1 se giornalieri senza stagionalità forte)
train_ts <- ts(train_set)

# Esegui il test di Johansen
# K deve essere almeno 2. Usiamo 7 per coerenza con i tuoi lag precedenti.
johan_test <- ca.jo(train_ts, type = "trace", K = 7, spec = "transitory")

summary(johan_test)




