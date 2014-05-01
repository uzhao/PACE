akima = function(x_in, x_out, covm) {
  junk = expand.grid(x_in, x_in)
  junkzi = as.vector(covm)
  valid = !is.na(junkzi)
  junkxi = as.vector(junk[,1])[valid]
  junkyi = as.vector(junk[,2])[valid]
  junkzi = junkzi[valid]
  interp(junkxi, junkyi, junkzi, x_out, x_out)$z
}
