rm(list=ls())
graphics.off()
cat("\014")

# Setup Workspace
setwd("~/uni/2025-2026/non param/progetto/clorofilla/Non-parametric-Statistics-2025-2026")

# Libraries
library(tidyverse)
library(fields)
library(viridis)
library(fda)
library(vars)
library(forecast)
library(xgboost)
library(tsDyn)
library(urca)

# ------------------------------------------------------------------------------
# 1. Data Loading and Pre-processing
# ------------------------------------------------------------------------------

data <- read.csv("dataset/dataset_clorofilla_filtrata_3anni.csv")
data$Lon <- round(data$Lon, 3)

# Arrotonda la colonna Latitudine
data$Lat <- round(data$Lat, 3)

dataset_ordinato_distanza <- read.csv("dataset/bordo_con_distanze.csv")
data_pos_x_daily <- data %>%
  inner_join(dataset_ordinato_distanza, by = c("Lat", "Lon"))

chl_pos_x_daily <- data_pos_x_daily %>%
  dplyr::select(Date, Lat, Lon, Chl, distanza_cumulata)

chl_pos_x_daily$Date <- as.Date(chl_pos_x_daily$Date)

x <- sort(unique(chl_pos_x_daily$distanza_cumulata))


plot(chl_pos_x_daily$Lon, chl_pos_x_daily$Lat)



matrice_df <- chl_pos_x_daily %>%
  dplyr::select(Date, distanza_cumulata, Chl) %>%
  # Ordina prima i dati per sicurezza
  arrange(Date, distanza_cumulata) %>%
  # Trasforma in formato largo
  pivot_wider(names_from = distanza_cumulata, values_from = Chl)

# Per convertirla in una vera matrice numerica (togliendo la colonna Date)
matrice_finale <- as.matrix(matrice_df[,-1])
rownames(matrice_finale) <- matrice_df$Date

sum(is.na(matrice_finale))
# # Salvataggio standard (separatore: virgola, decimale: punto)
write.csv(matrice_finale, "dataset/matrix_X_fda_daily.csv", row.names = TRUE)

# Load the daily Chlorophyll matrix and coordinate distances
matrice_finale <- read.csv("dataset/matrix_X_fda_daily.csv", row.names = 1)
matrice_finale <- as.matrix(matrice_finale)

dataset_ordinato_distanza <- read.csv("dataset/bordo_con_distanze.csv")
x <- sort(unique(dataset_ordinato_distanza$distanza_cumulata))

pos_mete <- read.csv("dataset/dataset_mete_approssimate.csv")
pos_mete <- pos_mete[pos_mete$Localita_Originale != "Milano Marittima", ]
# 1. Definiamo la lista delle località più rinomate
mete_top <- c("Grado", "Lignano Sabbiadoro", "Bibione", "Caorle", 
              "Lido di Jesolo", "Lido di Venezia", "Sottomarina (Chioggia)", 
              "Cervia", "Cesenatico", "Rimini")

# 2. Creiamo il dataset filtrato
# Assumendo che il tuo dataset si chiami 'pos_mete'
pos_mete <- pos_mete[pos_mete$Localita_Originale %in% mete_top, ]
# Plot dello sfondo
plot(dataset_ordinato_distanza$Lon, dataset_ordinato_distanza$Lat, col = "black", pch = 20)

# Aggiungi i punti
points(pos_mete$Lon, pos_mete$Lat, col = "red", pch = 16)

# Aggiungi i NOMI
text(pos_mete$Lon, pos_mete$Lat, 
     labels = pos_mete$Localita_Originale, 
     pos = 3,      # Posizione: 1=sotto, 2=sinistra, 3=sopra, 4=destra
     cex = 0.7,    # Dimensione del testo
     col = "red")  # Colore del testo# Visualizing Raw Data

n_giorni <- nrow(matrice_finale)
colori <- viridis(n_giorni)

# 1. Imposta i margini (lasciamo spazio a destra per la color bar)
par(mar = c(7, 4, 4, 6)) # Aumentato il primo valore per far stare i nomi ruotati

# 2. Plot principale
matplot(x, t(matrice_finale), 
        type = 'l', lty = 1, lwd = 2, col = colori,
        xlab = "", ylab = "Chl Content",
        main = "Raw Daily Chlorophyll Profiles",
        xaxt = "n")

# 3. Aggiungi linee verticali tratteggiate per ogni località
abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 2)

# 4. Aggiungi i nomi delle località sull'asse X
# las = 2 ruota il testo di 90 gradi per evitare sovrapposizioni
axis(1, at = pos_mete$distanza_cumulata, 
     labels = pos_mete$Localita_Originale, 
     las = 2, cex.axis = 0.7, col.ticks = "red")

# 5. Aggiungi la color bar per i giorni (richiede il pacchetto fields)
image.plot(legend.only = TRUE, 
           zlim = range(as.numeric(rownames(matrice_finale))), 
           col = colori, 
           legend.lab = "Day Index", 
           legend.line = 3)



# ------------------------------------------------------------------------------
# 2. Functional Data Analysis (Smoothing)
# ------------------------------------------------------------------------------

# Define B-spline Basis
basis_obj <- create.bspline.basis(
  rangeval = range(x), 
  breaks = x,           # Placing knots at real data points
  norder = 4            # Cubic splines
)

# Convert discrete data to functional objects
chl_fd <- Data2fd(argvals = x, y = t(matrice_finale), basisobj = basis_obj)

# Comparison between Raw Data and FDA Smoothing
par(mfrow = c(1, 2))
matplot(x, t(matrice_finale), type = "l", lty = 1, lwd = 1, col = colori,
        xlab = "Distance", ylab = "Chl", main = "Original Data (Raw)")
plot(chl_fd, col = colori, lty = 1, lwd = 2,
     xlab = "Distance", ylab = "Chl", main = "FDA Smooth (B-splines)")
# Highlight knots on the plot
tutti_i_knots <- c(basis_obj$rangeval[1], basis_obj$params, basis_obj$rangeval[2])
abline(v = tutti_i_knots, col = "gray90", lty = 3)
rug(tutti_i_knots, col = "red", lwd = 1.5)
par(mfrow = c(1, 1))

# ------------------------------------------------------------------------------
# 3. Functional Principal Component Analysis (FPCA)
# ------------------------------------------------------------------------------

n_comp <- 4
chl_pca <- pca.fd(chl_fd, nharm = n_comp)

# Explained Variance Analysis
var_spiegata <- chl_pca$varprop
names(var_spiegata) <- paste0("PC", 1:n_comp)
print("Explained Variance by PC:")
print(round(var_spiegata, 3))

# Scree Plot (Cumulative)
# xaxt = "n" nasconde l'asse x automatico
plot(cumsum(var_spiegata), type = "b", pch = 19, 
     xlab = "Number of PCs", ylab = "Cumulative Proportion",
     main = "FPCA Scree Plot",
     xaxt = "n",
     xlim = c(0.7, 4.3),
     ylim = c(0.7, 1.0)) 

# Aggiungi l'asse X manuale solo con interi
# side = 1 (basso), at = posizioni (1, 2, 3, 4)
axis(side = 1, at = 1:n_comp, labels = 1:n_comp)

# Aggiungi una griglia per leggere meglio (opzionale ma consigliato)
grid()

# AGGIUNTA DELLE FLAG (Etichette)
text(x = 1:n_comp,          # Posizione orizzontale (1, 2, 3, 4)
     y = cumsum(var_spiegata),           # Posizione verticale (il valore del punto)
     labels = round(cumsum(var_spiegata), 3), # Il testo da scrivere (arrotondato a 3 decimali)
     pos = 3,               # 3 significa "sopra il punto" (1=sotto, 2=sinistra, 4=destra)
     cex = 0.8,             # Grandezza del testo (0.8 è un po' più piccolo del normale)
     col = "black")          # Colore del testo

# ------------------------------------------------------------------------------
# 4. Forecasting Preparation
# ------------------------------------------------------------------------------

# Extract PC scores for modeling
scores_all <- chl_pca$scores[, 1:4]
n_tot <- nrow(scores_all)
n_val <- 7  # Validation horizon (1 week)

# Split into Training and Validation sets
train_set <- scores_all[1:(n_tot - n_val), ]
val_set   <- scores_all[(n_tot - n_val + 1):n_tot, ]

# Define ground truth curves for error calculation
curve_reali <- chl_fd[(n_tot - n_val + 1):n_tot]

# Function to calculate functional RMSE day by day
calcola_rmse_giornaliero <- function(curve_reali, curve_predette) {
  diff_fd <- curve_reali - curve_predette
  integrali_quadrati <- diag(inprod(diff_fd, diff_fd))
  return(sqrt(integrali_quadrati))
}

# ------------------------------------------------------------------------------
# 5. Model: Vector Autoregression (VAR)
# ------------------------------------------------------------------------------

fit_var_train <- VAR(train_set, p = 7, type = "both")
pred_var_obj <- predict(fit_var_train, n.ahead = n_val)

# Extract predicted scores
scores_var <- sapply(pred_var_obj$fcst, function(x) x[, "fcst"])
if(nrow(scores_var) != n_val) scores_var <- t(scores_var)

# Reconstruct functional profiles
coef_pred_var <- chl_pca$harmonics$coefs[, 1:4] %*% t(scores_var)
curve_var_fd <- fd(coef_pred_var + as.vector(chl_pca$meanfd$coefs), basis_obj)

err_var <- calcola_rmse_giornaliero(curve_reali, curve_var_fd)

# ------------------------------------------------------------------------------
# 6. Model: SARIMA (Univariate)
# ------------------------------------------------------------------------------

scores_sarima <- matrix(NA, nrow = n_val, ncol = 4)
for(i in 1:4) {
  fit_sarima <- auto.arima(train_set[, i], seasonal = TRUE)
  scores_sarima[, i] <- as.numeric(forecast(fit_sarima, h = n_val)$mean)
}

coef_pred_sarima <- chl_pca$harmonics$coefs[, 1:4] %*% t(scores_sarima)
curve_sarima_fd <- fd(coef_pred_sarima + as.vector(chl_pca$meanfd$coefs), basis_obj)

err_sarima <- calcola_rmse_giornaliero(curve_reali, curve_sarima_fd)

# ------------------------------------------------------------------------------
# 7. Model: Multivariate XGBoost (Recursive)
# ------------------------------------------------------------------------------

# Data preparation for supervised learning (lagged features)
prepare_multivariate_train <- function(data, lags) {
  X_list <- list()
  for(l in 1:lags) { X_list[[l]] <- data[(lags - l + 1):(nrow(data) - l), ] }
  return(list(X = do.call(cbind, X_list), Y = data[(lags + 1):nrow(data), ]))
}

d_train <- prepare_multivariate_train(train_set, 7)
xgb_models <- lapply(1:4, function(p) {
  xgboost(data = d_train$X, label = d_train$Y[, p], nrounds = 50, objective = "reg:squarederror", verbose = 0)
})

# Recursive prediction loop
scores_xgb <- matrix(NA, nrow = n_val, ncol = 4)
curr_hist <- train_set[(nrow(train_set) - 7 + 1):nrow(train_set), ]
for (s in 1:n_val) {
  input_row <- matrix(as.vector(t(curr_hist[nrow(curr_hist):1, ])), nrow = 1)
  step_p <- sapply(xgb_models, predict, input_row)
  scores_xgb[s, ] <- step_p
  curr_hist <- rbind(curr_hist[-1, ], step_p)
}

curve_xgb_fd <- fd((chl_pca$harmonics$coefs[, 1:4] %*% t(scores_xgb)) + as.vector(chl_pca$meanfd$coefs), basis_obj)
err_xgb <- calcola_rmse_giornaliero(curve_reali, curve_xgb_fd)

# ------------------------------------------------------------------------------
# 8. Model: Vector Error Correction Model (VECM)
# ------------------------------------------------------------------------------

# Johansen Test for Cointegration Rank
johan_test <- ca.jo(ts(train_set), type = "trace", K = 7, spec = "transitory")
summary(johan_test)

# Fitting VECM with rank r=3
fit_vecm <- VECM(train_set, lag = 7, r = 3, estim = "ML")
scores_vecm <- as.matrix(predict(fit_vecm, n.ahead = n_val))

curve_vecm_fd <- fd((chl_pca$harmonics$coefs[, 1:4] %*% t(scores_vecm)) + as.vector(chl_pca$meanfd$coefs), basis_obj)
err_vecm <- calcola_rmse_giornaliero(curve_reali, curve_vecm_fd)

# ------------------------------------------------------------------------------
# 9. Error Comparison and Visualization
# ------------------------------------------------------------------------------

# Unified Error Plot
par(mfrow = c(1, 1))
plot(1:n_val, err_var, type = "b", pch = 19, col = "red", ylim = c(0, 3), lwd = 2,
     xlab = "Forecasting Horizon (Days Ahead)", ylab = "Functional RMSE",
     main = "Model Comparison: Forecasting Error Growth", panel.first = grid())

lines(1:n_val, err_sarima, type = "b", pch = 17, col = "blue", lwd = 2)
lines(1:n_val, err_xgb, type = "b", pch = 15, col = "darkgreen", lwd = 2)
lines(1:n_val, err_vecm, type = "b", pch = 18, col = "darkblue", lwd = 2)

legend("topleft", legend = c("VAR", "SARIMA", "XGBoost", "VECM"), 
       col = c("red", "blue", "darkgreen", "darkblue"), 
       lty = 1, pch = c(19, 17, 15, 18), bty = "n", cex = 0.8)

# Final step visualization (VAR)
par(mfrow = c(4, 2))
for (i in 1:n_val) {
  plot(chl_fd[n_tot - n_val + i], col = "black", lwd = 2, ylim = c(0, 5),
       main = paste("VECM Validation: Day", i))
  plot(curve_xgb_fd[i], col = "red", lty = 2, lwd = 2, add = TRUE)
}
par(mfrow = c(1, 1))


# ------------------------------------------------------------------------------
# 10. Summary Performance Table
# ------------------------------------------------------------------------------

# Calculate Mean Functional RMSE for each model
mean_rmse <- data.frame(
  Model = c("VAR ", "SARIMA ", "XGBoost", "VECM "),
  Avg_Functional_RMSE = c(mean(err_var), mean(err_sarima), mean(err_xgb), mean(err_vecm))
)

# Sort by performance (lower is better)
mean_rmse <- mean_rmse[order(mean_rmse$Avg_Functional_RMSE), ]

# Print the results
print("--- Final Model Comparison (Average Functional RMSE) ---")
print(mean_rmse)

# Optional: Horizontal Bar Plot for visual comparison
barplot(mean_rmse$Avg_Functional_RMSE, 
        names.arg = mean_rmse$Model, 
        col = viridis(nrow(mean_rmse)),
        main = "Average Model Error Comparison",
        xlab = "Mean Functional RMSE (L2 Norm)",
        horiz = TRUE, las = 1, cex.names = 0.7)




# VAR Hp testing ----------------------------------------------------------

fit_var_train <- VAR(train_set, p = 14, type = "both")
summary(fit_var_train)
pred_var_obj <- predict(fit_var_train, n.ahead = n_val)

# Extract predicted scores
scores_var <- sapply(pred_var_obj$fcst, function(x) x[, "fcst"])
if(nrow(scores_var) != n_val) scores_var <- t(scores_var)

# Reconstruct functional profiles
coef_pred_var <- chl_pca$harmonics$coefs[, 1:4] %*% t(scores_var)
curve_var_fd <- fd(coef_pred_var + as.vector(chl_pca$meanfd$coefs), basis_obj)

err_var <- calcola_rmse_giornaliero(curve_reali, curve_var_fd)

# Test sui residui del modello vincente
var_check <- serial.test(fit_var_train, lags.pt = 12, type = "PT.asymptotic")
print(var_check)

# Estrazione dei residui
residui_var <- resid(fit_var_train)

# Plot ACF per i residui di ogni PC
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for(i in 1:4) {
  acf(residui_var[, i], main = paste("ACF Residui - PC", i), col = "red", lwd = 2)
}


par(mfrow = c(2, 2))
for(i in 1:4) {
  hist(residui_var[, i], breaks = 20, prob = TRUE, 
       main = paste("Distribuzione Residui PC", i), col = "lightgray")
  curve(dnorm(x, mean = mean(residui_var[, i]), sd = sd(residui_var[, i])), 
        add = TRUE, col = "red", lwd = 2)
}
par(mfrow = c(1, 1))





# Unified Error Plot
par(mfrow = c(1, 1))
plot(1:n_val, err_var, type = "b", pch = 19, col = "red", ylim = c(0, 3), lwd = 2,
     xlab = "Forecasting Horizon (Days Ahead)", ylab = "Functional RMSE",
     main = "Forecasting Error Growth", panel.first = grid())
text(1:length(err_var), err_var, labels = round(err_var, 4), pos = 3, cex = 0.8, col = "darkred")

# Final step visualization (VAR)
par(mfrow = c(4, 2))
for (i in 1:n_val) {
  plot(chl_fd[n_tot - n_val + i], col = "black", lwd = 2, ylim = c(0, 5),
       main = paste("VAR Validation: Day", i))
  plot(curve_var_fd[i], col = "red", lty = 2, lwd = 2, add = TRUE)
}
plot(1:n_val, err_var, type = "b", pch = 19, col = "red", ylim = c(0, 3), lwd = 2,
     xlab = "Forecasting Horizon (Days Ahead)", ylab = "Functional RMSE",
     main = "Forecasting Error Growth", panel.first = grid())
text(1:length(err_var), err_var, labels = round(err_var, 4), pos = 3, cex = 0.8, col = "darkred")
par(mfrow = c(1, 1))

