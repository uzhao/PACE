isregular = function(x) {
  ans = 'sparse'
  xx = unlist(x)
  f = length(xx) / length(unique(xx)) / length(x)
  if (f > 0.75) {
    ans = 'missing'
  }
  if (f == 1) {
    ans = 'regular'
  }
  ans
}
