
setwd("~/uni/2025-2026/non param/progetto/clorofilla/Non-parametric-Statistics-2025-2026")



# 1. Caricamento e Preparazione Dati -------------------------------------------
library(tidyverse)

# data <- read.csv("dataset/dataset_MERGIATO")
# dataset_ordinato_distanza <- read.csv("dataset/pos_con_distanze.csv")
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
# # Salvataggio standard (separatore: virgola, decimale: punto)
# write.csv(matrice_finale, "dataset/matrix_X_fda_daily.csv", row.names = TRUE)
matrice_finale <- read.csv("dataset/matrix_X_fda_daily.csv")

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
           legend.lab = "Settimana", 
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

par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

for(i in 1:3) {
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

# Plot delle 7 curve predette (sfumatura viridis)
plot(curve_predette_fd, col = viridis(n_val), lty = 1, lwd = 2,
     main = "Profili Spaziali Predetti (Set di Validazione)")

# Reale (l'ultimo pezzo del dataset originale)
plot(chl_fd[(n_tot - n_val + 1)], col = "black", lwd = 2)
# Predetta
plot(curve_predette_fd[1], col = "red", lty = 2, add = TRUE)















#### XGBoost ####
library(xgboost)

# Parametri
n_val <- 7       # Ultime 5 giorni per validazione
n_lag <- 7      # Usiamo le 2 giorni precedenti per predire la successiva
n_pcs <- ncol(scores_all)

# Funzione per creare i lag
create_lags <- function(data, lags) {
  X <- embed(data, lags + 1)
  y <- X[, 1:ncol(data)]                # Target (tutte le PC al tempo t)
  features <- X[, (ncol(data)+1):ncol(X)] # Predittori (PC ai tempi t-1, t-2...)
  return(list(y = y, X = features))
}

data_lagged <- create_lags(scores_all, n_lag)

# split
n_rows <- nrow(data_lagged$X)
train_idx <- 1:(n_rows - n_val)

X_train <- data_lagged$X[train_idx, ]
y_train <- data_lagged$y[train_idx, ]

X_val   <- data_lagged$X[(n_rows - n_val + 1):n_rows, ]
y_val   <- data_lagged$y[(n_rows - n_val + 1):n_rows, ]


modelli_xgb <- list()
previsioni_val <- matrix(NA, nrow = n_val, ncol = n_pcs)

for(i in 1:n_pcs) {
  # Creazione DMatrix (formato ottimizzato per XGBoost)
  dtrain <- xgb.DMatrix(data = X_train, label = y_train[, i])
  dval   <- xgb.DMatrix(data = X_val, label = y_val[, i])
  
  # Fit del modello
  params <- list(objective = "reg:squarederror", eta = 0.1, max_depth = 4)
  
  fit <- xgb.train(params = params, 
                   data = dtrain, 
                   nrounds = 100, 
                   watchlist = list(train = dtrain, eval = dval),
                   early_stopping_rounds = 10,
                   print_every_n = 0)
  
  modelli_xgb[[i]] <- fit
  
  # Predizione sul set di validazione
  previsioni_val[, i] <- predict(fit, X_val)
}


par(mfrow = c(n_pcs, 1), mar = c(4, 4, 2, 1))

for(i in 1:n_pcs) {
  plot(y_val[, i], type = "b", pch = 19, col = "black", 
       main = paste("Validation PC", i), xlab = "Settimana di test", ylab = "Score")
  lines(previsioni_val[, i], type = "b", pch = 18, col = "orange", lwd = 2)
  legend("topright", legend = c("Reale", "XGBoost"), col = c("black", "orange"), lty = 1, bty = "n")
}


