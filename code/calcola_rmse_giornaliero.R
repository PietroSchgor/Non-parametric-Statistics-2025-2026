
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