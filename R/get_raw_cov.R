get_raw_cov = function(data, mu, regular, error) {
  n = nrow(data)
  p = ncol(data)
  data = as.matrix(sweep(data, 2, mu, '-'))
  if (regular == 'regular') {
    ans = t(data) %*% data / (n - 1)
    attr(ans, 'weight') = matrix(1, nrow = nrow(ans), ncol = ncol(ans))
  }
  else {
    ans = matrix(NA, p, p)
    w = matrix(0, p, p)
    for (i in 1:p) {
      for (j in i:p) {
        # ube cause over est lambda
        # mle cause under est lambda
        junk = mle(data[,i], data[,j])
        ans[i,j] = ans[j,i] = junk[1]
        w[i,j] = w[j,i] = junk[2]
      }
    }
    attr(ans, 'weight') = w
  }
  attr(ans, 'vary') = diag(ans)
  attr(ans, 'varyw') = diag(attr(ans, 'weight'))
  if (error) {
    diag(ans) = NA
    diag(attr(ans, 'weight')) = 0
  }
  colnames(ans) = colnames(data)
  rownames(ans) = colnames(data)
  ans
}

ube = function(x, y) {
  v = na.omit(x * y)
  if (length(v) < 2) {
    return(c(NA, 0))
  }
  else {
    j = length(v) - 1
    return(c(sum(v) / j, j + 1))
  }
}

mle = function(x, y) {
  v = na.omit(x * y)
  j = length(v)
  return(c(mean(v), j))
}
