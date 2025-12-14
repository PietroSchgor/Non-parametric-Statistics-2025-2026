rm(list = ls())
setwd("~/Documents/Nonparametric/Project/Non-parametric-Statistics-2025-2026")

library(readxl)
library(maps)
library(mapdata)

Localita_balneari_Nord_Adriatico <- read_excel("dataset/Localita_balneari_Nord_Adriatico.xlsx")

xlim_val <- c(12, 14)
ylim_val <- c(45, 46)

{
map('world', xlim = xlim_val, ylim = ylim_val, 
    col = "lightgray", fill = TRUE, bg = "lightblue") # bg=lightblue per il mare
map.axes()

points(Localita_balneari_Nord_Adriatico$Lon, 
       Localita_balneari_Nord_Adriatico$Lat, 
       col = "red", pch = 19, cex = 1)

offset_y <- 0.03 # quanto è lungo il segmentino
# segments(...): Disegna una linea nera (l'asta) che parte dalle coordinate del
# lido (y0) e va leggermente verso nord (y1 = Lat + offset_y).
segments(x0 = Localita_balneari_Nord_Adriatico$Lon, 
         y0 = Localita_balneari_Nord_Adriatico$Lat, 
         x1 = Localita_balneari_Nord_Adriatico$Lon, 
         y1 = Localita_balneari_Nord_Adriatico$Lat + offset_y, 
         col = "black", lwd = 1)
# text(...): Posiziona il testo esattamente alla fine dell'asta.
text(x = Localita_balneari_Nord_Adriatico$Lon, 
     y = Localita_balneari_Nord_Adriatico$Lat + offset_y, 
     labels = Localita_balneari_Nord_Adriatico$Località, 
     pos = 3,       # 3 = sopra
     cex = 0.7,     # Testo un po' più piccolo per non sovrapporre troppo
     col = "black",
     font = 2)      # 2 = Grassetto

title("Località Balneari Nord Adriatico")
}




{
# 1. IMPOSTAZIONI ESTETICHE REALISTICHE
# Colori stile "National Geographic"
colore_terra <- "wheat"       # Colore sabbia/grano per la terra
colore_mare  <- "azure1"      # Azzurro chiaro realistico per il mare
colore_confine <- "gray60"    # Colore sottile per i confini

# Limiti della mappa (Nord Adriatico)
xlim_val <- c(12.0, 14.0)
ylim_val <- c(45.0, 46.0)

# 2. DISEGNO DELLA MAPPA AD ALTA RISOLUZIONE
# Usa 'worldHires' se hai mapdata, altrimenti usa 'world'
database_mappa <- "worldHires" # o "world" se non hai mapdata

map(database_mappa, xlim = xlim_val, ylim = ylim_val, 
    fill = TRUE, col = colore_terra, bg = colore_mare, 
    lwd = 0.5, border = colore_confine) # Bordi più sottili

# 3. AGGIUNTA DEL RETICOLO (GRIGLIA)
# Aggiunge linee bianche tratteggiate per latitudine e longitudine (molto "nautico")
abline(h = seq(45, 46, 0.5), v = seq(12, 14, 0.5), col = "white", lty = 3)
map.axes(cex.axis = 0.8, col = "gray40") # Assi grigio scuro, meno invadenti

# 4. DISEGNO DEI PUNTI (LIDI) CON STILE "SPILLO"
points(Localita_balneari_Nord_Adriatico$Lon, 
       Localita_balneari_Nord_Adriatico$Lat, 
       pch = 21,          # 21 = cerchio riempibile
       bg = "red",        # Interno rosso
       col = "black",     # Bordo nero
       cex = 1.2)         # Leggermente più grandi

# 5. BANDIERINE REALISTICHE
offset_y <- 0.04 

# Asta della bandiera (nera sottile)
segments(x0 = Localita_balneari_Nord_Adriatico$Lon, 
         y0 = Localita_balneari_Nord_Adriatico$Lat, 
         x1 = Localita_balneari_Nord_Adriatico$Lon, 
         y1 = Localita_balneari_Nord_Adriatico$Lat + offset_y, 
         col = "black", lwd = 1)

# Etichette con un "alone" o sfondo per renderle leggibili sulla costa frastagliata
# Usiamo shadowtext (un trucco base R: scrivere prima bianco grosso sotto, poi nero sopra)
text(x = Localita_balneari_Nord_Adriatico$Lon, 
     y = Localita_balneari_Nord_Adriatico$Lat + offset_y, 
     labels = Localita_balneari_Nord_Adriatico$Località, 
     pos = 3, cex = 0.7, font = 2, col = "white", offset = 0.5) # Sfondo bianco simulato
text(x = Localita_balneari_Nord_Adriatico$Lon, 
     y = Localita_balneari_Nord_Adriatico$Lat + offset_y, 
     labels = Localita_balneari_Nord_Adriatico$Località, 
     pos = 3, cex = 0.7, font = 2, col = "black") # Testo vero e proprio

# 6. AGGIUNTA SCALA E NORD (ELEMENTI "PRO")
# Scala chilometrica in basso a destra
map.scale(x = 13.5, y = 44.1, ratio = FALSE, relwidth = 0.15, metric = TRUE, cex = 0.6)

# Freccia del Nord (disegnata a mano con arrows)
arrows(x0 = 12.2, y0 = 45.8, x1 = 12.2, y1 = 45.9, length = 0.1, lwd = 2, col = "black")
text(12.2, 45.75, "N", font = 2, cex = 1.2)

title("Mappa Lidi Nord Adriatico", col.main = "black")
}
