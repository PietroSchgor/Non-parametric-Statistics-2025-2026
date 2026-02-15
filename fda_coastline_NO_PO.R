rm(list=ls()); graphics.off(); cat("\014")
# Setup Workspace
#setwd("~/uni/2025-2026/non param/progetto/clorofilla/Non-parametric-Statistics-2025-2026")
setwd("~/Documents/Nonparametric/Project/Non-parametric-Statistics-2025-2026")

{# Libraries
library(tidyverse)
library(fields)
library(viridis)
library(fda)
library(vars)
library(forecast)
library(xgboost)
library(tsDyn)
library(urca)
library(ggrepel)}

#### 1. Data Loading and Pre-processing #####
# 
# data <- read.csv("dataset/dataset_clorofilla_filtrata_3anni.csv")
# 
# # Arrotonda la colonna Longitudine
# data$Lon <- round(data$Lon, 3)
# 
# # Arrotonda la colonna Latitudine
# data$Lat <- round(data$Lat, 3)
# 
# data$Date <- as.Date(data$Date)
# 
# # Verifica il risultato
# head(data)
# dataset_ordinato_distanza <- read.csv("dataset/bordo_con_distanze.csv")
# data_pos_x_daily <- data %>%
#   inner_join(dataset_ordinato_distanza, by = c("Lat", "Lon")) %>%
#   filter(Lat >= 45)
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
# matrice_finale_Po <- as.matrix(matrice_df[,-1])
# 
# # 1. Converti la colonna in vero oggetto Date (controlla il tuo formato!)
# matrice_df$Date <- as.Date(matrice_df$Date, format = "%Y-%m-%d") 
# 
# # 2. Solo ora assegna alle rownames
# rownames(matrice_finale_Po) <- as.character(matrice_df$Date)
# 
# sum(is.na(matrice_finale_Po))
# # # Salvataggio standard (separatore: virgola, decimale: punto)
# write.csv(matrice_finale_Po, "dataset/matrix_Po_fda_daily.csv", row.names = TRUE)

{# Load the daily Chlorophyll matrix and coordinate distances
matrice_finale_Po <- read.csv("dataset/matrix_Po_fda_daily.csv", row.names = 1)
matrice_finale_Po <- as.matrix(matrice_finale_Po)

dataset_ordinato_distanza <- read.csv("dataset/bordo_con_distanze.csv")
dataset_ordinato_distanza <- dataset_ordinato_distanza %>%
  filter (Lat >= 45)
x <- sort(unique(dataset_ordinato_distanza$distanza_cumulata))

pos_mete <- read.csv("dataset/dataset_mete_approssimate.csv")
# 1. Definiamo la lista delle località più rinomate
mete_top <- c("Grado", "Lignano Sabbiadoro", "Bibione", "Caorle", 
              "Lido di Jesolo", "Lido di Venezia", "Sottomarina (Chioggia)")
# 2. Creiamo il dataset filtrato
# Assumendo che il tuo dataset si chiami 'pos_mete'
pos_mete <- pos_mete[pos_mete$Localita_Originale %in% mete_top, ]
pos_mete$Localita_Originale[pos_mete$Localita_Originale == "Sottomarina (Chioggia)"] <- "Sottomarina"
# Plot dello sfondo
plot(dataset_ordinato_distanza$Lon, dataset_ordinato_distanza$Lat, col = "black", pch = 20)

# Aggiungi i punti
points(pos_mete$Lon, pos_mete$Lat, col = "red", pch = 16)

# Aggiungi i NOMI
text(pos_mete$Lon, pos_mete$Lat, 
     labels = pos_mete$Localita_Originale, 
     pos = 3,      # Posizione: 1=sotto, 2=sinistra, 3=sopra, 4=destra
     cex = 0.7,    # Dimensione del testo
     col = "red" 
     )  # Colore del testo# Visualizing Raw Data
}


# 1. Converti le righe in formato Data (assicurati che il formato sia corretto)
date_oggetti <- as.Date(rownames(matrice_finale_Po))

# 2. Definisci la lunghezza del ciclo (usi 265 come richiesto, o 365 per l'anno)
ciclo <- 365 

#### RAW CYCLIC CHLOROPHYLL PROFILES ####
{# 3. Calcola l'indice ciclico per ogni data
# Usiamo il numero di giorni trascorsi da una data di riferimento
giorni_da_inizio <- as.numeric(date_oggetti - min(date_oggetti))
indici_ciclici <- (giorni_da_inizio %% ciclo) + 1

# 4. Crea una tavolozza base di 265 colori
tavolozza_base <- viridis(ciclo)

# 5. Assegna a ogni riga il colore corrispondente al suo punto nel ciclo
colori_ciclici <- tavolozza_base[indici_ciclici]

par(mar = c(8, 4, 4, 6)) # Margine destro ridotto a 6, bastano per la legend
matplot(x, t(matrice_finale_Po), 
        type = 'l', lty = 1, lwd = 2, # lwd ridotto per chiarezza se ci sono molte linee
        col = colori_ciclici,
        xaxt = "n",
        xlab = "", ylab = "Chlorophyll [mg/m^3]",
        main = paste("Raw Chlorophyll Profiles (Cycle:", ciclo, "days)"))

# Aggiunta località e linee (come nel tuo codice)
abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 2)
axis(1, at = pos_mete$distanza_cumulata, labels = pos_mete$Localita_Originale, 
     las = 2, cex.axis = 0.7, col.ticks = "darkgray")

# Legenda: ora rappresenta la posizione all'interno del ciclo
image.plot(legend.only = TRUE,
           zlim = c(1, 365),            # Ora la legenda va da 1 a 365 (Day of the Year)
           col = tavolozza_base,           # Usa solo i 365 colori base
           legend.lab = "Day of the Year", 
           legend.line = 2.5,           # Avvicina la scritta alla barra per non farla uscire
           smallplot = c(0.88, 0.90, 0.25, 0.85)) # FONDAMENTALE: (x_min, x_max, y_min, y_max) in %
}

#### RAW DAILY CHLOROPHYLL PROFILES ####
{n_giorni <- nrow(matrice_finale_Po)
colori <- viridis(n_giorni)

par(mar = c(8, 4, 4, 6)) # Margine destro ridotto a 6, bastano per la legend
matplot(x, t(matrice_finale_Po),
        type = 'l', lty = 1, lwd = 2, col = colori,
        xaxt = "n",
        xlab = "", ylab = "Chlorophyll [mg/m^3]",
        main = "Raw Daily Chlorophyll Profiles")

#3. Aggiungi linee verticali tratteggiate per ogni località
abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 2)

# 4. Aggiungi i nomi delle località sull'asse X
# las = 2 ruota il testo di 90 gradi per evitare sovrapposizioni
axis(1, at = pos_mete$distanza_cumulata,
     labels = pos_mete$Localita_Originale,
     las = 2, cex.axis = 0.7, col.ticks = "darkgray")

# Add a color bar for days
# Sostituisci il blocco della legenda con questo:
image.plot(legend.only = TRUE, 
           zlim = c(1, n_giorni),       # Range dal giorno 1 all'ultimo giorno assoluto
           col = colori,                # Colori continui
           legend.lab = "Day Index", 
           legend.line = 2.5,           # Avvicina il testo alla barra
           smallplot = c(0.88, 0.90, 0.25, 0.85)) # FONDAMENTALE: blocca la legenda dentro lo schermo
}

#### Functional Data Analysis (Smoothing) ####

# Define B-spline Basis
basis_obj <- create.bspline.basis(
  rangeval = range(x), 
  breaks = x,           # Placing knots at real data points
  norder = 4            # Cubic splines
)

# Convert discrete data to functional objects
chl_fd <- Data2fd(argvals = x, y = t(matrice_finale_Po), basisobj = basis_obj)

# Comparison between Raw Data and FDA Smoothing

{par(mfrow = c(2, 1), mar = c(7, 4, 3, 2))

# --- PRIMO GRAFICO: Dati Grezzi ---
matplot(x, t(matrice_finale_Po), type = "l", lty = 1, lwd = 1, col = colori,
        xlab = "", ylab = "Chlorophyll [mg/m^3]", main = "Original Raw Data",
        xaxt = "n") # Nasconde l'asse numerico

# Aggiunge linee guida per le città e l'asse personalizzato
abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 3)
axis(1, at = pos_mete$distanza_cumulata, 
     labels = pos_mete$Localita_Originale, 
     las = 2, cex.axis = 0.8, col.ticks = "darkgray")

# --- SECONDO GRAFICO: Dati Smussati (FDA) ---
plot(chl_fd, col = colori, lty = 1, lwd = 2,
     xlab = "", ylab = "Chlorophyll [mg/m^3]", main = "FDA Smooth (B-splines)",
     xaxt = "n") # Nasconde l'asse numerico

# Evidenzia i knots (i nodi della spline)
tutti_i_knots <- c(basis_obj$rangeval[1], basis_obj$params, basis_obj$rangeval[2])
abline(v = tutti_i_knots, col = "gray90", lty = 3)
rug(tutti_i_knots, col = "red", lwd = 1.5)

# Aggiunge linee guida per le città e l'asse personalizzato anche qui
abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 3)
axis(1, at = pos_mete$distanza_cumulata, 
     labels = pos_mete$Localita_Originale, 
     las = 2, cex.axis = 0.8, col.ticks = "darkgray")

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)}

#### RAW VS SMOOTHED ####
{library(fields)
  
  # 1. IMPOSTAZIONE MARGINI (Lasciamo 7 linee di margine ESTERNO a destra)
  par(mfrow = c(2, 1), 
      mar = c(7, 4, 3, 2),  
      oma = c(0, 0, 0, 7))  
  
  # --- PRIMO GRAFICO: Dati Grezzi ---
  matplot(x, t(matrice_finale_Po), type = "l", lty = 1, lwd = 1, col = colori,
          xlab = "", ylab = "Chlorophyll [mg/m^3]", main = "Original Raw Data",
          xaxt = "n") 
  
  abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 3)
  axis(1, at = pos_mete$distanza_cumulata, 
       labels = pos_mete$Localita_Originale, 
       las = 2, cex.axis = 0.8, col.ticks = "darkgray")
  
  
  # --- SECONDO GRAFICO: Dati Smussati (FDA) ---
  plot(chl_fd, col = colori, lty = 1, lwd = 2,
       xlab = "", ylab = "Chlorophyll [mg/m^3]", main = "FDA Smooth (B-splines)",
       xaxt = "n") 
  
  tutti_i_knots <- c(basis_obj$rangeval[1], basis_obj$params, basis_obj$rangeval[2])
  abline(v = tutti_i_knots, col = "gray90", lty = 3)
  rug(tutti_i_knots, col = "red", lwd = 1.5)
  
  abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 3)
  axis(1, at = pos_mete$distanza_cumulata, 
       labels = pos_mete$Localita_Originale, 
       las = 2, cex.axis = 0.8, col.ticks = "darkgray")
  
  
  # 3. IL TRUCCO PER LA LEGENDA GLOBALE
  # Diciamo a R di usare l'intero schermo (fig = c(0,1,0,1)) senza cancellare i grafici sotto (new = TRUE)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  
  # Disegniamo un grafico completamente vuoto e trasparente
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
  
  # Ora attacchiamo la color bar a questo grafico globale!
  image.plot(legend.only = TRUE, 
             zlim = c(1, nrow(matrice_finale_Po)), 
             col = colori, 
             legend.lab = "Day Index", 
             legend.line = 2.5, 
             smallplot = c(0.92, 0.94, 0.15, 0.85)) # (x_min, x_max, y_min, y_max)
  
  # Ripristino i parametri per non sballare i grafici successivi
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))}

#### ANALISI OUTLIERS ####
# Installazione e caricamento del pacchetto necessario
# install.packages("roahd")
library(roahd)

# 1. Valutazione dell'oggetto funzionale su una griglia di punti
# Usiamo i tuoi punti 'x' (distanza cumulata)
grid_points <- x
matrice_valutata <- t(eval.fd(grid_points, chl_fd)) # Deve essere (n_curve x n_punti)

# 2. Creazione dell'oggetto fData (formato richiesto da roahd)
chl_fData <- fData(grid_points, matrice_valutata)

# Impostiamo il layout per vedere entrambi i grafici
par(mfrow = c(1, 1))

# --- OUTLIERGRAM ---
# Rileva gli outlier di SHAPE (forma)
# Mostra la parabola teorica della relazione tra MBD e MEI
outliergram(chl_fData, main = "Outliergram (Shape Outliers)", xlab = "MEI", ylab = "MBD")
# 1. Definizione della parabola teorica

# 1. Calcolo degli indici MEI e MBD per i tuoi dati
mei_reali <- MEI(chl_fData)
mbd_reali <- MBD(chl_fData)

# 2. Otteniamo l'identificazione ufficiale degli outlier dal pacchetto roahd
# (così siamo sicuri che i punti rossi corrispondano ai calcoli del modello)
out_info <- outliergram(chl_fData, display = FALSE)
is_outlier_shape <- 1:length(mei_reali) %in% out_info$ID_outliers

#### OUTLIERGRAM ####

# 3. Plot dell'Outliergram Completo
# Disegniamo prima la struttura teorica

{# 1. Definizione corretta di n e dei coefficienti
n <- ncol(chl_fd$coefs) 
a0 <- -2 / (n * (n - 1))
a1 <- 2 * (n + 1) / (n - 1)
a2 <- a0

# 2. Calcolo della parabola teorica
mei_seq <- seq(0, 1, length.out = 100)
mbd_teorico <- a0 + a1 * mei_seq + a2 * (n^2) * mei_seq^2

# CORREZIONE 1: Calcolo della VERA distanza della linea rossa
# Calcoliamo dove "dovrebbe" essere ogni punto sulla parabola
mbd_teorico_punti <- a0 + a1 * mei_reali + a2 * (n^2) * mei_reali^2

# La distanza è quanto il punto reale cade sotto la parabola
distanze <- mbd_teorico_punti - mbd_reali

# La soglia matematica per gli outlier di forma è Q3 + 1.5 * IQR delle distanze
Q3 <- quantile(distanze, 0.75)
IQR_dist <- IQR(distanze)
shift_soglia <- Q3 + 1.5 * IQR_dist

# CORREZIONE 2: Calcolo dei limiti dell'asse Y (ylim)
# Diciamo a R di inquadrare la finestra dal punto più basso al punto più alto
min_y <- min(mbd_teorico - shift_soglia, mbd_reali) - 0.02
max_y <- max(mbd_teorico, mbd_reali) + 0.02

# Preparazione Colori e Mesi
date_tutte <- as.Date(rownames(matrice_finale_Po))
mesi_tutti <- as.numeric(format(date_tutte, "%m"))
colori_mesi_palette <- palette.colors(12, "Paired")
nomi_mesi <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# 3. Plot della struttura (NOTA L'AGGIUNTA DI ylim!)
plot(mei_seq, mbd_teorico, type = "l", col = "blue", lwd = 2, lty = 2,
     xlab = "MEI", ylab = "MBD", main = "Outliergram: Shape Outliers per Month",
     ylim = c(min_y, max_y)) # <--- FONDAMENTALE PER NON TAGLIARE LA LINEA ROSSA

# 4. Disegniamo la linea rossa usando lo shift matematicamente corretto
lines(mei_seq, mbd_teorico - shift_soglia, col = "red", lwd = 2)

# 5. Sovrapponiamo i punti normali (Grigi)
points(mei_reali[!is_outlier_shape], mbd_reali[!is_outlier_shape], 
       pch = 20, col = alpha("gray80", 0.5))

# 6. Sovrapponiamo gli Outlier colorati per mese
mesi_out <- mesi_tutti[is_outlier_shape]

points(mei_reali[is_outlier_shape], mbd_reali[is_outlier_shape], 
       pch = 19, 
       col = colori_mesi_palette[mesi_out], 
       cex = 1.5)
grid()

# # 7. Etichette dei mesi sopra i pallini
# text(mei_reali[is_outlier_shape], mbd_reali[is_outlier_shape], 
#      labels = nomi_mesi[mesi_out], 
#      pos = 3, cex = 0.7, col = "black", font = 2)

# 8. Doppia Legenda
legend("topleft", 
       legend = c("Theoretical Limit", "Threshold", "Dati Points"),
       col = c("blue", "red", "gray80"), 
       lty = c(2, 1, NA), pch = c(NA, NA, 20), bty = "o", bg = "white", cex = 0.8)

mesi_presenti <- sort(unique(mesi_out))
legend("bottomright", 
       legend = nomi_mesi[mesi_presenti],
       col = colori_mesi_palette[mesi_presenti], 
       pch = 19, title = "Outlier Months", bty = "o", bg = "white", cex = 0.8, ncol = 2)
}

#### FUNCTIONAL BOXPLOT ####
{
# 1. Allarghiamo il margine inferiore (Bottom = 8) per farci stare le scritte verticali
par(mar = c(8, 4, 4, 2)) 

# 2. Creazione del Functional Boxplot (nascondendo l'asse X automatico)
roahd::fbplot(chl_fData, 
              main = "Functional Boxplot (Magnitude Outliers)", 
              xlab = "",                            # Togliamo la scritta "Distance"
              ylab = "Chlorophyll [mg/m^3]", 
              xaxt = "n")                           # FONDAMENTALE: nasconde i numeri di default

# 3. Aggiungiamo le linee verticali tratteggiate in corrispondenza delle località
abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 3)

# 4. Disegniamo l'asse X personalizzato con i nomi delle città
axis(1, 
     at = pos_mete$distanza_cumulata, 
     labels = pos_mete$Localita_Originale, 
     las = 2,                # Ruota il testo di 90 gradi (in verticale)
     cex.axis = 0.8,         # Rimpicciolisce un po' i nomi per evitare sovrapposizioni
     col.ticks = "darkgray") # Colore delle tacchette

# Ripristiniamo i margini standard per i grafici successivi
par(mar = c(5, 4, 4, 2) + 0.1)
}

#### TEMPORAL DISTRIBUTION OUTLIERS ####

# Estrazione indici outlier di forma (dall'outliergram)
# Nota: outliergram() non restituisce direttamente gli indici in modo semplice, 
# ma possiamo usare la funzione sottostante:
outliers_shape <- outliergram(chl_fData, display = FALSE)$ID_outliers

# Estrazione indici outlier di magnitudo (dal boxplot)
outliers_magnitude <- fbplot(chl_fData, display = FALSE)$ID_outliers

# Unione degli outlier totali
tutti_outliers <- unique(c(outliers_shape, outliers_magnitude))

cat("Indici degli outlier rilevati:", tutti_outliers, "\n")
cat("Date corrispondenti agli outlier:", rownames(matrice_finale_Po)[tutti_outliers], "\n")

library(ggplot2)
library(scales) # Necessaria per la formattazione avanzata delle date

# 1. Preparazione dei dati (English translation of categories)
date_tutte <- as.Date(rownames(matrice_finale_Po))
n_tot_giorni <- length(date_tutte)

categoria_outlier <- rep("Normal", n_tot_giorni)
categoria_outlier[outliers_shape] <- "Shape Outlier"
categoria_outlier[outliers_magnitude] <- "Magnitude Outlier"

# Identificazione giorni con entrambi i tipi di outlier
entrambi <- intersect(outliers_shape, outliers_magnitude)
if(length(entrambi) > 0) {
  categoria_outlier[entrambi] <- "Both (Shape & Magnitude)"
}

# Creazione del dataframe
df_outliers_time <- data.frame(
  Date = date_tutte,
  Type = factor(categoria_outlier, 
                levels = c("Normal", "Shape Outlier", "Magnitude Outlier", "Both (Shape & Magnitude)")),
  Value = 1
)

# 2. Plot Temporale (Timeline Plot)
ggplot(df_outliers_time, aes(x = Date, y = Type, color = Type)) +
  geom_hline(aes(yintercept = Type), color = "gray90", linewidth = 0.5) +
  
  # Punti per gli outlier (escludiamo i Normal per pulizia visiva)
  geom_point(data = subset(df_outliers_time, Type != "Normal"), 
             size = 4, alpha = 0.8) +
  
  # Segmenti verticali (rug-style)
  geom_segment(data = subset(df_outliers_time, Type != "Normal"),
               aes(x = Date, xend = Date, y = 0.5, yend = Type), 
               linetype = "dotted", alpha = 0.5) +
  
  # --- ASSE X: MESI IN INGLESE ---
  scale_x_date(
    date_breaks = "1 month", 
    date_labels = "%b %Y" # Es: "Jan 2024"
  ) +
  
  # Palette colori
  scale_color_manual(values = c("Normal" = "gray95", 
                                "Shape Outlier" = "orange", 
                                "Magnitude Outlier" = "firebrick", 
                                "Both (Shape & Magnitude)" = "purple")) +
  
  # Titoli Tradotti
  labs(title = "Temporal Distribution of Functional Outliers",
       subtitle = "Comparison between Shape and Magnitude anomalies",
       x = "Month", 
       y = "Outlier Type", 
       color = "Category") +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotazione per leggibilità
  )


# 1. Filtriamo i dati per escludere i giorni "Normal"
df_only_outliers <- subset(df_outliers_time, Type != "Normal")

# 2. Plot: Stacked Histogram of Outliers
ggplot(df_only_outliers, aes(x = Date, fill = Type)) +
  # Creiamo l'istogramma con barre impilate
  # binwidth = 30 raggruppa circa per mese
  geom_histogram(binwidth = 30, position = "stack", alpha = 0.8, color = "white") +
  
  # --- ASSE X: MESI IN INGLESE ---
  scale_x_date(
    date_breaks = "1 month", 
    date_labels = "%b %Y"
  ) +
  
  # Palette colori coerente con il Timeline Plot
  scale_fill_manual(values = c("Shape Outlier" = "orange", 
                               "Magnitude Outlier" = "firebrick", 
                               "Both (Shape & Magnitude)" = "purple")) +
  
  # Titoli Tradotti in Inglese
  labs(title = "Monthly Frequency of Functional Outliers",
       subtitle = "Stacked counts of Shape and Magnitude anomalies",
       x = "Month", 
       y = "Outlier Count", 
       fill = "Category") +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1), # Ruotiamo le etichette per leggibilità
    
    # --- CENTRATURA DEI TITOLI ---
    plot.title = element_text(hjust = 0.5, face = "bold"), # Centra e mette in grassetto il titolo
    plot.subtitle = element_text(hjust = 0.5)              # Centra il sottotitolo
  )


#### Split e FPCA Rigorosa (SOLO su Training) ####

# Trova l'indice corrispondente a una data specifica
date_disponibili <- as.Date(rownames(matrice_finale_Po))
start_val <- which(date_disponibili == as.Date("2024-05-01")) # Esempio
# Scegli il giorno di inizio della validazione (es. giorno 100)
# Se vuoi l'ultima settimana, lasceresti: start_val <- n_tot - n_val + 1
n_val     <- 7

# Indici per il training (tutto ciò che viene PRIMA della settimana scelta)
# Nota: in una serie temporale, di solito non si usa ciò che viene "dopo" per predire il "prima"
indici_train <- 1:(start_val - 1)
indici_val   <- start_val:(start_val + n_val - 1)

# Split dei coefficienti
coef_train <- chl_fd$coefs[, indici_train]
coef_val   <- chl_fd$coefs[, indici_val]

# Ricostruzione oggetti FD
train_fd <- fd(coef_train, basis_obj)
val_fd   <- fd(coef_val, basis_obj)

# Aggiorna i metadati (nomi delle date) per non fare confusione nei plot
train_fd$fdnames$reps <- chl_fd$fdnames$reps[indici_train]
val_fd$fdnames$reps   <- chl_fd$fdnames$reps[indici_val]

# FPCA calcolata SOLO sul Training
n_comp <- 4
chl_pca_train <- pca.fd(train_fd, nharm = n_comp)

# Estrarre gli score per il training
train_set <- chl_pca_train$scores[, 1:n_comp]


# Explained Variance Analysis
var_spiegata <- chl_pca_train$varprop
names(var_spiegata) <- paste0("PC", 1:n_comp)
print("Explained Variance by PC:")
print(round(var_spiegata, 3))

# Calcoliamo la varianza cumulata per comodità
cum_var <- cumsum(var_spiegata)

{
# 1. Plot Principale 
# xaxt = "n" toglie i decimali sull'asse X
# ylim allarga lo spazio verticale (dal 80% fino a poco sopra il 100%)
plot(cum_var, type = "b", pch = 19, 
     xlab = "Number of PCs", ylab = "Cumulative Proportion",
     main = "FPCA Scree Plot",
     xaxt = "n", 
     xlim = c(0.8, 4.2),
     ylim = c(0.80, 1.02)) 

# 2. Asse X personalizzato (solo numeri interi da 1 a 4)
axis(1, at = 1:n_comp)

# Aggiungiamo una griglia di sfondo per aiutare la lettura
grid()

# 3. Aggiunta delle Flag (Etichette)
text(x = 1:n_comp - 0.1,                  # Posizione asse X (1, 2, 3, 4)
     y = cum_var,                   # Posizione asse Y (l'altezza del punto)
     labels = round(cum_var, 3),    # Testo da stampare (arrotondato a 3 decimali)
     pos = 3,                       # pos = 3 mette il testo esattamente *sopra* il punto
     cex = 0.8,                     # Rimpicciolisce un po' il font
     col = "black",                  # Colore della flag
     font = 2)                      # font = 2 mette il grassetto
}

#### PLOT QUATTRO COMPONENTI PRINCIPALI ####

# Visualizziamo le prime 4 componenti principali
par(mfrow = c(2, 2))
plot(chl_pca_train, nx = 100, cex.main = 0.8) 
par(mfrow = c(1, 1))

{# Impostiamo il layout 2x2 e allarghiamo il margine inferiore
  par(mfrow = c(2, 2), mar = c(7, 4, 3, 1))
  
  # IL SEGRETO: Creiamo una finta griglia di 100 punti equidistanti (nx = 100)
  x_grid <- seq(min(x), max(x), length.out = 100)
  
  # 1. Valutiamo la curva media sulla griglia a 100 punti
  mean_val <- as.vector(eval.fd(x_grid, chl_pca_train$meanfd))
  
  # 2. Valutiamo le 4 curve delle armoniche sulla stessa griglia
  harm_val <- eval.fd(x_grid, chl_pca_train$harmonics)
  
  # 3. Fattore di scala
  # IL FIX È QUI: Moltiplichiamo per 2 per avere +/- 2 Deviazioni Standard!
  fac <- 2 * sqrt(chl_pca_train$values[1:4]) 
  
  
  # CICLO DI PLOT PERSONALIZZATO
  for (i in 1:4) {
    
    perc_var <- round(var_spiegata[i] * 100, 1)
    
    # Calcoliamo l'effetto "Più" e "Meno" usando i 100 punti puliti
    y_plus  <- mean_val + fac[i] * harm_val[, i]
    y_minus <- mean_val - fac[i] * harm_val[, i]
    
    # 1. Disegniamo la curva media (linea continua nera) su x_grid
    plot(x_grid, mean_val, type = "l", lwd = 2, col = "black",
         ylim = range(c(y_plus, y_minus, mean_val)), # Inquadra perfettamente tutto
         xaxt = "n", xlab = "", ylab = "Chlorophyll [mg/m^3]",
         main = paste0("PC", i, ": ", perc_var, "% of Variance"), 
         cex.main = 1)
    
    # 2. Aggiungiamo i pallini "+" e "-" 
    points(x_grid, y_plus, pch = "+", col = "blue", cex = 0.8)
    points(x_grid, y_minus, pch = "-", col = "red", cex = 0.8)
    
    # 3. Aggiungiamo le linee verticali di riferimento per le città 
    abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 3)
    
    # 4. Aggiungiamo l'asse X personalizzato con le città 
    axis(1, 
         at = pos_mete$distanza_cumulata, 
         labels = pos_mete$Localita_Originale, 
         las = 2, 
         cex.axis = 0.7, 
         col.ticks = "darkgray")
  }
  
  # Ripristiniamo i parametri grafici standard alla fine
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)}

# 1. Estrarre i coefficienti
coef_val   <- val_fd$coefs                # Dimensione: [nbasis x 7]
coef_media <- chl_pca_train$meanfd$coefs  # Dimensione: [nbasis x 1]

# 2. Sottrarre la media da ogni colonna del validation
# Usiamo sweep() o una semplice sottrazione (R gestisce il riciclo se le righe coincidono)
coef_val_centered <- sweep(coef_val, 1, coef_media, "-")

# 3. Creare l'oggetto fd centrato (stessa base del training)
val_centered_fd <- fd(coef_val_centered, basis_obj)

# 4. Ora la proiezione funzionerà perfettamente
# Calcoliamo i "veri" score proiettando sulle armoniche (eigenfunctions)
val_scores_reali <- inprod(val_centered_fd, chl_pca_train$harmonics)

# Verifica: deve essere una matrice 7 x 4
print(dim(val_scores_reali))


# Funzione di Supporto per Ricostruzione Curve Predette ####

# Questa funzione usa la media e le armoniche del TRAIN per ricostruire
ricostruisci_curve <- function(scores_predetti, pca_obj) {
  # scores_predetti deve essere (n_val x n_comp)
  # armoniche: (nbasis x n_comp)
  coef_reconstr <- pca_obj$harmonics$coefs[, 1:n_comp] %*% t(scores_predetti)
  # Aggiungiamo la media del training
  curve_fd <- fd(coef_reconstr + as.vector(pca_obj$meanfd$coefs), pca_obj$harmonics$basis)
  return(curve_fd)
}

# 1. Ricostruzione delle curve di validation usando gli score REALI
# Usiamo la funzione che abbiamo definito prima
curve_val_ricostruite <- ricostruisci_curve(val_scores_reali, chl_pca_train)

# 2. Plot di confronto (Validation Reconstruction)
# Impostiamo 2 righe e 4 colonne, e allarghiamo il margine in basso a 7 per le scritte
par(mfrow = c(2, 4), mar = c(7, 3, 2, 1)) 

for (i in 1:n_val) {
  # Calcoliamo il range per avere assi coerenti
  y_lim <- range(c(eval.fd(x, val_fd[i]), eval.fd(x, curve_val_ricostruite[i])))
  
  # Plot della curva reale (nera, continua) con asse X nascosto
  plot(val_fd[i], col = "black", lwd = 2, lty = 1, ylim = y_lim,
       main = paste("Day", i, "(Validation)"),
       xlab = "", ylab = "Chl", 
       xaxt = "n") # Nascondiamo la distanza numerica
  
  # Sovrapponiamo la ricostruzione FPCA (rossa, tratteggiata)
  plot(curve_val_ricostruite[i], col = "red", lwd = 2, lty = 2, add = TRUE)
  
  # Aggiungiamo le linee verticali di riferimento per le città
  abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 3)
  
  # Aggiungiamo l'asse X personalizzato con le città ruotate
  axis(1, 
       at = pos_mete$distanza_cumulata, 
       labels = pos_mete$Localita_Originale, 
       las = 2, 
       cex.axis = 0.6,          # Leggermente più piccolo perché i grafici sono stretti
       col.ticks = "darkgray")
  
  # Aggiungiamo la legenda SOLO sul primo grafico
  if(i == 1) {
    legend("topright", 
           legend = c("True", "FPCA"),
           col = c("black", "red"), 
           lty = c(1, 2), 
           bty = "o",           # Riquadro chiuso
           bg = "white",        # Sfondo bianco per non far sovrapporre le linee!
           cex = 0.7) 
  }
}

# Ripristiniamo i parametri grafici standard (layout singolo, margini normali)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)


# 5. Model: Vector Autoregression (VAR) ####

fit_var_train <- VAR(train_set, p = 7, type = "both")
pred_var_obj <- predict(fit_var_train, n.ahead = n_val)

# Extract predicted scores
scores_var <- sapply(pred_var_obj$fcst, function(x) x[, "fcst"])
if(nrow(scores_var) != n_val) scores_var <- t(scores_var)

# Reconstruct functional profiles
coef_pred_var <- chl_pca_train$harmonics$coefs[, 1:4] %*% t(scores_var)
curve_var_fd <- fd(coef_pred_var + as.vector(chl_pca_train$meanfd$coefs), basis_obj)

# Function to calculate functional RMSE day by day
calcola_rmse_giornaliero <- function(curve_reali, curve_predette) {
  diff_fd <- curve_reali - curve_predette
  integrali_quadrati <- diag(inprod(diff_fd, diff_fd))
  return(sqrt(integrali_quadrati))
}
err_var <- calcola_rmse_giornaliero(val_fd, curve_var_fd)


# 6. Model: SARIMA (Univariate) ####

scores_sarima <- matrix(NA, nrow = n_val, ncol = 4)
for(i in 1:4) {
  fit_sarima <- auto.arima(train_set[, i], seasonal = TRUE)
  scores_sarima[, i] <- as.numeric(forecast(fit_sarima, h = n_val)$mean)
}

coef_pred_sarima <- chl_pca_train$harmonics$coefs[, 1:4] %*% t(scores_sarima)
curve_sarima_fd <- fd(coef_pred_sarima + as.vector(chl_pca_train$meanfd$coefs), basis_obj)

err_sarima <- calcola_rmse_giornaliero(val_fd, curve_sarima_fd)


# 7. Model: Multivariate XGBoost (Recursive) ####

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

curve_xgb_fd <- fd((chl_pca_train$harmonics$coefs[, 1:4] %*% t(scores_xgb)) + as.vector(chl_pca_train$meanfd$coefs), basis_obj)
err_xgb <- calcola_rmse_giornaliero(val_fd, curve_xgb_fd)


# 8. Model: Vector Error Correction Model (VECM) ####

# Johansen Test for Cointegration Rank
johan_test <- ca.jo(ts(train_set), type = "trace", K = 7, spec = "transitory")
summary(johan_test)

# Fitting VECM with rank r=3
fit_vecm <- VECM(train_set, lag = 7, r = 3, estim = "ML")
scores_vecm <- as.matrix(predict(fit_vecm, n.ahead = n_val))

curve_vecm_fd <- fd((chl_pca_train$harmonics$coefs[, 1:4] %*% t(scores_vecm)) + as.vector(chl_pca_train$meanfd$coefs), basis_obj)
err_vecm <- calcola_rmse_giornaliero(val_fd, curve_vecm_fd)


# 9. Error Comparison and Visualization ####

# Unified Error Plot
par(mfrow = c(1, 1))
plot(1:n_val, err_var, type = "b", pch = 19, col = "red", ylim = c(0, 1), lwd = 2,
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
  plot(val_fd[i], col = "black", lwd = 2, ylim = c(0, 3),
       main = paste("VAR Validation: Day", i))
  plot(curve_var_fd[i], col = "red", lty = 2, lwd = 2, add = TRUE)
}
par(mfrow = c(1, 1))


# Final step visualization (XGB)
par(mfrow = c(4, 2))
for (i in 1:n_val) {
  plot(val_fd[i], col = "black", lwd = 2, ylim = c(0, 3),
       main = paste("XGB Validation: Day", i))
  plot(curve_xgb_fd[i], col = "darkgreen", lty = 2, lwd = 2, add = TRUE)
}
par(mfrow = c(1, 1))

# Final step visualization (VECM)
par(mfrow = c(4, 2))
for (i in 1:n_val) {
  plot(val_fd[i], col = "black", lwd = 2, ylim = c(0, 3),
       main = paste("VECM Validation: Day", i))
  plot(curve_vecm_fd[i], col = "darkblue", lty = 2, lwd = 2, add = TRUE)
}
par(mfrow = c(1, 1))

{
  # Final step visualization (VECM)
  
  # Layout 2 righe x 4 colonne, margine in basso a 7 per far spazio ai nomi delle città
  par(mfrow = c(2, 4), mar = c(7, 3, 2, 1)) 
  
  for (i in 1:n_val) {
    # 1. Disegna la curva reale (nera) togliendo i numeri sull'asse X (xaxt = "n")
    plot(val_fd[i], col = "black", lwd = 2, ylim = c(0, 3),
         main = paste("VECM Validation: Day", i),
         xlab = "", ylab = "Chl", xaxt = "n")
    
    # 2. Sovrappone la curva predetta dal VECM (blu scura tratteggiata)
    plot(curve_vecm_fd[i], col = "darkblue", lty = 2, lwd = 2, add = TRUE)
    
    # 3. Aggiunge le linee verticali di riferimento
    abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 3)
    
    # 4. Aggiunge i nomi delle località ruotati di 90 gradi (las = 2)
    axis(1, 
         at = pos_mete$distanza_cumulata, 
         labels = pos_mete$Localita_Originale, 
         las = 2, 
         cex.axis = 0.6,          # Rimpiccioliamo un po' il font perché i riquadri sono stretti
         col.ticks = "darkgray")
    
    # 5. Aggiungiamo la legenda SOLO sul primo grafico per non appesantire tutto
    if(i == 1) {
      legend("topright", 
             legend = c("True", "VECM"),
             col = c("black", "darkblue"), 
             lty = c(1, 2), 
             bty = "o",           # Riquadro chiuso
             bg = "white",        # Sfondo bianco opaco per coprire le linee dietro
             cex = 0.7) 
    }
  }
  
  # Ripristino del layout singolo e dei margini standard
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) 
}


# 10. Summary Performance Table ####

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
par(mfrow=c(1,1))
barplot(mean_rmse$Avg_Functional_RMSE, 
        names.arg = mean_rmse$Model, 
        col = viridis(nrow(mean_rmse)),
        main = "Average Model Error Comparison",
        xlab = "Mean Functional RMSE (L2 Norm)",
        horiz = TRUE, las = 1, cex.names = 0.7)

{
  # 1. Creiamo una palette coerente legata al VERO nome del modello
  # (Usiamo trimws() per sicurezza, nel caso ci siano spazi nei nomi come "VAR ")
  nomi_modelli <- trimws(mean_rmse$Model)
  colori_modelli <- c("VAR" = "red", 
                      "SARIMA" = "blue", 
                      "XGBoost" = "darkgreen", 
                      "VECM" = "darkblue")[nomi_modelli]
  
  # 2. Impostiamo i margini (allargando un po' il sinistro per far leggere bene i nomi)
  par(mfrow = c(1,1), mar = c(5, 6, 4, 2))
  
  # 3. Disegniamo il Barplot
  # Salviamo il risultato nella variabile 'bp' (ci darà le coordinate Y delle barre!)
  # Allunghiamo xlim del 15% ( * 1.15 ) per fare spazio alle scritte dei numeri
  bp <- barplot(mean_rmse$Avg_Functional_RMSE, 
                names.arg = mean_rmse$Model, 
                col = colori_modelli,
                main = "Average Model Error Comparison",
                xlab = "Mean Functional RMSE (L2 Norm)",
                horiz = TRUE, 
                las = 1, 
                cex.names = 0.9,
                xlim = c(0, max(mean_rmse$Avg_Functional_RMSE) * 1.15)) 
  
  # 4. Aggiungiamo le etichette di testo (RMSE arrotondato a 2 decimali)
  # pos = 4 significa "posiziona il testo a destra delle coordinate X"
  text(x = mean_rmse$Avg_Functional_RMSE, 
       y = bp, 
       labels = round(mean_rmse$Avg_Functional_RMSE, 2), 
       pos = 4, 
       cex = 0.9, 
       font = 2,           # font = 2 rende il numeretto in grassetto
       col = "black")      # Puoi usare anche colori_modelli se vuoi i numeri colorati!
}


# # VAR Hp testing ----------------------------------------------------------
# 
# fit_var_train <- VAR(train_set, p = 14, type = "both")
# summary(fit_var_train)
# pred_var_obj <- predict(fit_var_train, n.ahead = n_val)
# 
# # Extract predicted scores
# scores_var <- sapply(pred_var_obj$fcst, function(x) x[, "fcst"])
# if(nrow(scores_var) != n_val) scores_var <- t(scores_var)
# 
# # Reconstruct functional profiles
# coef_pred_var <- chl_pca_train$harmonics$coefs[, 1:4] %*% t(scores_var)
# curve_var_fd <- fd(coef_pred_var + as.vector(chl_pca_train$meanfd$coefs), basis_obj)
# 
# err_var <- calcola_rmse_giornaliero(val_fd, curve_var_fd)
# 
# # Test sui residui del modello vincente
# var_check <- serial.test(fit_var_train, lags.pt = 12, type = "PT.asymptotic")
# print(var_check)
# 
# # Estrazione dei residui
# residui_var <- resid(fit_var_train)
# 
# # Plot ACF per i residui di ogni PC
# par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
# for(i in 1:4) {
#   acf(residui_var[, i], main = paste("ACF Residui - PC", i), col = "red", lwd = 2)
# }
# 
# 
# par(mfrow = c(2, 2))
# for(i in 1:4) {
#   hist(residui_var[, i], breaks = 20, prob = TRUE, 
#        main = paste("Distribuzione Residui PC", i), col = "lightgray")
#   curve(dnorm(x, mean = mean(residui_var[, i]), sd = sd(residui_var[, i])), 
#         add = TRUE, col = "red", lwd = 2)
# }
# par(mfrow = c(1, 1))

# Diagnostica del Modello Vincitore (VECM) ####

# 1. Estrazione dei residui del VECM
residui_vecm <- resid(fit_vecm)

# 2. Portmanteau Test (Test di Ljung-Box multivariato) per i residui del VECM
# Il pacchetto vars ha una funzione specifica serial.test, ma il VECM va prima convertito in VAR
vecm_as_var <- vec2var(johan_test, r = 3)
vecm_check <- serial.test(vecm_as_var, lags.pt = 12, type = "PT.asymptotic")
print(vecm_check)

# 3. ACF Plot for the residuals of each Component (VECM)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for(i in 1:4) {
  acf(residui_vecm[, i], 
      main = paste("ACF of VECM Residuals - PC", i), # Titolo tradotto
      col = "darkblue", 
      lwd = 2)
}
{
  # 3. ACF Plot for the residuals of each Component (VECM)
  
  # mar = c(bottom, left, top, right) -> Riportiamo top a 3 e bottom a 4
  # mgp = c(2.2, 0.8, 0) -> Distanza perfetta e sicura per le etichette degli assi
  par(mfrow = c(2, 2), 
      mar = c(4, 4, 3, 1),  
      mgp = c(2.2, 0.8, 0)) 
  
  for(i in 1:4) {
    acf(residui_vecm[, i], 
        main = paste("PC", i, "- VECM Residuals ACF"), 
        xlab = "Lag",                                  
        ylab = "ACF",                                  
        col = "darkblue", 
        lwd = 2)
  }
  
  # Ripristiniamo i parametri di default
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))
}

# 4. Plot della Distribuzione dei Residui (Istogramma + Curva Normale)
par(mfrow = c(2, 2))
for(i in 1:4) {
  hist(residui_vecm[, i], breaks = 20, prob = TRUE, 
       main = paste("Distribuzione Residui VECM - PC", i), col = "lightblue", border = "white")
  # Aggiungiamo la curva normale teorica in rosso
  curve(dnorm(x, mean = mean(residui_vecm[, i]), sd = sd(residui_vecm[, i])), 
        add = TRUE, col = "red", lwd = 2)
}
par(mfrow = c(1, 1)) # Ripristino




# # Unified Error Plot
# par(mfrow = c(1, 1))
# plot(1:n_val, err_var, type = "b", pch = 19, col = "red", ylim = c(0, 3), lwd = 2,
#      xlab = "Forecasting Horizon (Days Ahead)", ylab = "Functional RMSE",
#      main = "Forecasting Error Growth", panel.first = grid())
# text(1:length(err_var), err_var, labels = round(err_var, 4), pos = 3, cex = 0.8, col = "darkred")
# 
# # Final step visualization (VAR)
# par(mfrow = c(4, 2))
# for (i in 1:n_val) {
#   plot(val_fd[i], col = "black", lwd = 2, ylim = c(0, 5),
#        main = paste("VAR Validation: Day", i))
#   plot(curve_var_fd[i], col = "red", lty = 2, lwd = 2, add = TRUE)
# }
# plot(1:n_val, err_var, type = "b", pch = 19, col = "red", ylim = c(0, 3), lwd = 2,
#      xlab = "Forecasting Horizon (Days Ahead)", ylab = "Functional RMSE",
#      main = "Forecasting Error Growth", panel.first = grid())
# text(1:length(err_var), err_var, labels = round(err_var, 4), pos = 3, cex = 0.8, col = "darkred")
# par(mfrow = c(1, 1))

# Unified Error Plot (Singolo, se ti serve a schermo intero) ####

par(mfrow = c(1, 1))
plot(1:n_val, err_vecm, type = "b", pch = 18, col = "darkblue", 
     ylim = c(0, max(err_vecm) + 0.5), # Adatta l'asse Y automaticamente all'errore
     lwd = 2,
     xlab = "Forecasting Horizon (Days Ahead)", ylab = "Functional RMSE",
     main = "VECM Forecasting Error Growth", panel.first = grid())
text(1:length(err_vecm), err_vecm, labels = round(err_vecm, 4), pos = 3, cex = 0.8, col = "navy")


# Final Step Visualization (VECM) - LA GRIGLIA DEFINITIVA PER IL REPORT ####

# Impostiamo 2 righe e 4 colonne, allargando il margine inferiore (7) per le città
par(mfrow = c(2, 4), mar = c(7, 3, 2, 1))

# Disegniamo i primi 7 giorni con le curve e le città
for (i in 1:n_val) {
  plot(val_fd[i], col = "black", lwd = 2, ylim = c(0, 5),
       main = paste("VECM Validation: Day", i),
       xlab = "", ylab = "Chl", xaxt = "n") # Nascondiamo i numeri
  
  # Curva predetta dal VECM (blu tratteggiata)
  plot(curve_vecm_fd[i], col = "darkblue", lty = 2, lwd = 2, add = TRUE)
  
  # Linee e nomi delle città
  abline(v = pos_mete$distanza_cumulata, col = "gray80", lty = 3)
  axis(1, at = pos_mete$distanza_cumulata, labels = pos_mete$Localita_Originale, 
       las = 2, cex.axis = 0.6, col.ticks = "darkgray")
  
  # Legenda solo nel primo riquadro
  if(i == 1) {
    legend("topright", legend = c("True", "VECM"),
           col = c("black", "darkblue"), lty = c(1, 2), bty = "o", bg = "white", cex = 0.7)
  }
}

# OTTAVO RIQUADRO: Il grafico di crescita dell'errore (RMSE)

# Ripristiniamo margini normali solo per quest'ultimo riquadro
par(mar = c(5, 4, 3, 1)) 
plot(1:n_val, err_vecm, type = "b", pch = 18, col = "darkblue", 
     ylim = c(0, max(err_vecm) + 0.2), lwd = 2,
     xlab = "Days Ahead", ylab = "RMSE",
     main = "Error Growth", panel.first = grid())
text(1:length(err_vecm), err_vecm, labels = round(err_vecm, 3), pos = 3, cex = 0.8, col = "navy")

# Ripristiniamo il layout singolo e i margini di default alla fine
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

