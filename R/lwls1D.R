lwls = function(bandwidth = NULL, kernel = 2, x_in, y_in,
                w_in = rep(1 / length(x_in), length(x_in)),
                x_out, drv = 0, degree = 1, cv_mode = 0,
                verbose) {
  ##kernel:
  #0 epanechnikov
  #1 rectangle
  #2 gaussian
  #3 quartic
  #4 variant of Gaussian
  
  #cv_mode
  #0 no cv calculation
  #1 ordinary cross-validation
  #2 generalized cross-validation
  #3 geometric mean of GCV and minimum bandwidth
  
  if (!missing(verbose)) {
    cat('Bandwidth for ', verbose, ' is: ', bandwidth, '.\n', sep = "")
  }

  if (missing(x_out)) {
    x_out = x_in
  }
  
  bandwidth = as.double(bandwidth)
  kernel = as.integer(kernel)
  
  valid = ((!is.na(y_in)) & (w_in != 0))
  x_in = as.double(x_in[valid])
  if (any(diff(x_in) > ((max(x_in) - min(x_in)) / 4))) {
    warning('The data is too sparse! Result could be bad.')
  }
  n_in = as.integer(length(x_in))
  y_in = as.double(y_in[valid])
  w_in = as.double(w_in[valid])

  x_out = as.double(x_out)
  n_out = as.integer(length(x_out))
  mu_out = as.double(rep(0, n_out))
  output = as.double(rep(0, n_out))
  weight_out = as.double(rep(0, n_out))

  drv = as.integer(drv)
  degree = as.integer(degree)
  cv_mode = as.integer(cv_mode)
  cv_value = as.double(0)
  .C("lwls",
     bandwidth = bandwidth, kernel = kernel, 
     x_in = x_in, y_in = y_in, w_in = w_in,
     n_in = n_in, n_out = n_out, 
     x_out = x_out, 
     mu_out = mu_out, output = output, weight_out = weight_out, 
     drv = drv, degree = degree, 
     cv_mode = cv_mode, cv_value = cv_value, 
     NAOK = TRUE, PACKAGE = "PACE")
}

find_min_bandwidth = function(x_in, number_of_points = 2) {
  d = x_in[number_of_points:length(x_in)] - x_in[1:(length(x_in) - number_of_points + 1)]
  min(d)
}

find_max_bandwidth = function(x_in) {
  max(x_in) - min(x_in)
}

bandwidth_choice_1d = function(
  kernel, x_in, y_in, w_in, degree, cv_mode) {
  
  ad_hoc_lwls = function(bw) {
    lwls(bw, kernel, x_in, y_in, w_in,
         drv = 0, degree = degree, 
         cv_mode = cv_mode)$cv_value
  }
  lower = find_min_bandwidth(x_in)
  upper = find_max_bandwidth(x_in)
  ans = optimize(ad_hoc_lwls, lower = lower, upper = upper)$minimum
  if (cv_mode == 3) {
    ans = sqrt(ans * lower)
  }
  ans
}
