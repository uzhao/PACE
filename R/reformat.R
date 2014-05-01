reformat = function(x, y, regular, num_bins) {
  if (regular == 'sparse') {
    ans = bin_sparse(x, y, num_bins)
  }
  else {
    tt = sort(unique(unlist(x)))
    if (regular == 'regular') {
      ans = matrix(unlist(y), nrow = length(y), ncol = length(tt), byrow = TRUE)
    }
    else {
      ans = matrix(NA, nrow = length(y), ncol = length(tt), byrow = TRUE)
      for (i in 1:length(y)) {
        ans[i, match(x[[i]], tt)] = y[[i]]
      }
    }
    colnames(ans) = tt
  }
  ans
}

bin_sparse = function(x, y, num_bins) {
  xx = unlist(x)
  yy = unlist(y)
  n = length(yy)
  # index for subject
  # -1 for easy to use in C
  ii = unlist(lapply(1:length(y), function(i) rep(i, length(y[[i]])))) - 1
  minx = min(xx)
  maxx = max(xx)
  grids = seq(minx, maxx, length.out = num_bins + 1)
  d = grids[2] - grids[1]
  data = matrix(NA, ncol = num_bins, nrow = length(y))
  subject_sum = rep(0, num_bins)
  subject_count = rep(0, num_bins)
  for (i in 1:length(y)) {
    index = as.integer(cut(x[[i]], grids, include.lowest = TRUE))
    junk = tapply(y[[i]], index, mean)
    data[i, match(names(junk), 1:num_bins)] = junk
  }
  colnames(data) = (grids[1:(length(grids)-1)] + grids[2:length(grids)]) / 2
  data
}
