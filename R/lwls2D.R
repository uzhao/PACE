lwls_seq = function(bandwidth = NULL, kernel = 2, x_in, m_in, 
                    x_out, drv = 0, degree = 1,
                    cv_mode = 0, verbose = NULL) {

  if (!is.null(verbose)) {
    cat('Bandwidth for ', verbose, ' is: ', bandwidth, '.\n', sep = "")
  }

  if (missing(x_out)) {
    x_out = x_in
  }

  cv = 0
  
  ans = matrix(NA, length(x_out), ncol(m_in))
  w = matrix(0, length(x_out), ncol(m_in))
  for (i in 1:ncol(m_in)) {
    y_in = as.numeric(m_in[,i])
    w_in = as.numeric(attr(m_in, "weight")[,i])
    junk = lwls(bandwidth, kernel, 
                x_in, y_in, w_in,
                x_out, drv, degree, cv_mode = cv_mode)
    ans[,i] = junk$output
    w[,i] = junk$weight_out
    cv = cv + junk$cv_value
  }
  attr(ans, 'weight') = w
  attr(ans, 'cv') = cv
  ans
}

bandwidth_choice_seq = function(kernel, x_in, m_in, degree, cv_mode) {
  ad_hoc_lwls = function(bw) {
    attr(lwls_seq(bw, kernel, x_in, m_in, x_in, drv = 0, degree, cv_mode), 'cv')
  }
  lower = find_min_bandwidth(x_in)
  upper = find_max_bandwidth(x_in)
  ans = optimize(ad_hoc_lwls, lower = lower, upper = upper)$minimum
  if (cv_mode == 3) {
    ans = sqrt(ans * lower)
  }
  ans  
}

