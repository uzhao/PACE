FPCA = function(
  # data
  x, y, 
  # mean related
  bandwidth_mean_cv_mode = c("geo_mean", "gcv", "cv"),
  # covariance related
  bandwidth_cov_cv_mode = c("geo_mean", "gcv", "cv"),
  # choosing number of components
  no_of_components = -1,
  # method of choosing
  selection = c("FVE", "BIC", "AIC"), 
  FVE_threshold = 0.95,
  # maximum components
  maxk = 20,
  # measurements error or not
  error = TRUE,
  # grid of smoothed covariance surface
  cov_smooth_grid = 50,
  # kernel
  kernel = c("Auto", "Epanechnikov", "Rectangular", "Gaussian", "Quartic", "Gausvar"),
  # number of bins
  num_bins = NULL,
  # truncation step agreement
  # use cross validation for rho or user define its value
  # using new method in FPCscore.pdf to obtain ridge
  new_ridge = TRUE,
  rho_cv = FALSE, rho = 0,
  # abandon rate
  out_percent = 0.25, userrange = NULL) {
  
  # initialization
  num_obs = sapply(y, length)
  kernel = match.arg(kernel)
  bandwidth_mean_cv_mode = match.arg(bandwidth_mean_cv_mode)
  bandwidth_mean_cv_mode = switch(bandwidth_mean_cv_mode,
                                  "cv" = 1,
                                  "gcv" = 2,
                                  "geo_mean" = 3)
  
  bandwidth_cov_cv_mode = match.arg(bandwidth_cov_cv_mode)
  bandwidth_cov_cv_mode = switch(bandwidth_cov_cv_mode,
                                 "cv" = 1,
                                 "gcv" = 2,
                                 "geo_mean" = 3)
  selection = match.arg(selection)
  sigma = NULL
  mu_cov = NULL

  # data check
  if (sum(sapply(x, function(xx) any(is.na(xx))))) {
    stop('FPCA is aborted because x contain NA(s)!')
  }
  if ((length(y) != length(x)) | (any(sapply(y, length) != sapply(x, length)))){
    stop('FPCA is aborted because y and x don\'t match!')
  }
  if (all(num_obs == 1)) {
    stop('FPCA is aborted because the data do not contain repeated measurements!')
  }
  
  # no of components check
  if (!missing(no_of_components)) {
    no_of_components = as.integer(no_of_components)
    if (no_of_components < 0 && no_of_components != -1) {
      no_of_components = -1
      warning('no_of_components must be a positive integer!\n',
              'Will choose automatically (Default method: FVE >= 0.95)!\n', sep = "")
    }
  }

  # regular check
  sparse_flag = FALSE
  regular = isregular(x)
  if (regular == 'sparse') {
    sparse_flag = TRUE
  }
  
  # binning
  num_bins = get_num_bins(num_bins, num_obs, length(x), sparse_flag)
  data = reformat(x, y, regular, num_bins)

  if (any(is.na(data))) {
    regular = 'missing'
  }
  else {
    regular = 'regular'
  }  

  if (sparse_flag) {
    cat('The design is ', regular, ' (convert from sparse).\n', sep = "")
  }
  else {
    cat('The design is ', regular, '.\n', sep = "")
  }

  # kernel check
  if (kernel == "Auto") {
    if(regular == 'regular' & ncol(data) >= 20) {
      kernel = "Epanechnikov"
    }
    else {
      kernel = "Gaussian"
    }
  }
  cat('Use ', kernel, ' kernel.\n', sep = "")
  kernel = switch(kernel,
                  "Epanechnikov" = 0,
                  "Rectangular" = 1,
                  "Gaussian" = 2,
                  "Quartic" = 3,
                  "Gausvar" = 4)
  
  ## PACE
  # Part I: Obtain smoothed mean curve.
  cat('Part I: Obtain smoothed mean curve.\n')


  x_data = as.double(colnames(data))
  y_data = as.double(colMeans(data, na.rm = TRUE))
  w_data = as.double(colSums(!is.na(data)))
  
  x_cov = seq(min(x_data), max(x_data), length.out = cov_smooth_grid + 1)
  x_cov = (x_cov[1:(length(x_cov)-1)] + x_cov[2:length(x_cov)])/2

  bandwidth_mean = bandwidth_choice_1d(
    kernel, x_data, y_data, w_data, degree = 1, cv_mode = bandwidth_mean_cv_mode)
  
  mu_data = lwls(bandwidth_mean, kernel, x_data, y_data, w_data, x_out = x_data, 
            cv_mode = 0, verbose = "mean function")$output
  mu = lwls(bandwidth_mean, kernel, x_data, y_data, w_data, x_out = x_cov, 
            cv_mode = 0)$output

  rm(y_data)

  # Part II: Obtain smoothed covariance surface.
  cat('Part II: Obtain smoothed covariance surface.\n')
  raw_cov = get_raw_cov(data, mu_data, regular, error)
  
  png('design_plot.png')
  image(x = x_data, y = x_data, z = is.na(raw_cov), col = c('black', 'white'), main = "Design Plot")
  dev.off()

  # fill holes in cov matrix
  junk1 = attr(raw_cov, 'weight')
  junk2 = attr(raw_cov, 'vary')
  junk3 = attr(raw_cov, 'varyw')
  raw_cov = akima(x_data, x_data, raw_cov)
  attr(raw_cov, 'weight') = junk1
  attr(raw_cov, 'vary') = junk2
  attr(raw_cov, 'varyw') = junk3
  rm(junk1, junk2, junk3)

  cov_surface_smoothed = get_smoothed_cov(
    kernel, x_data, raw_cov, x_cov, bandwidth_cov_cv_mode)
  # partial_surface_smoothed = get_smoothed_partial(
  #   kernel, x_data, raw_cov, x_cov, bandwidth_cov_cv_mode)

  # Part III: Choose number of principal components functions.
  cat('Part III: Choose number of principal components functions.\n')
  if (is.numeric(out_percent)) {
    if (out_percent < 0 || out_percent > 0.25) {
      warning('out_percent must between 0 to 0.25, reset out_percent to 0.25!\n')
      out_percent = 0.25
    }
  }
  else {
    warning('out_percent input is unsupported, reset out_percent to 0.25!\n')
    out_percent = 0.25
  }

  cat('Use', out_percent, 'for leave out boundary.\n')
  userrange = as.numeric(quantile(x_data, c(out_percent, 1 - out_percent)))
  
  cat('Use userrange from', userrange[1],
      'to', userrange[2], 'for leave out boundary.\n')

  full_x_data = x_data
  cut_index = x_data >= userrange[1] & x_data <= userrange[2]
  x_data = x_data[cut_index]
  full_x_cov = x_cov
  x_cov = seq(min(x_data), max(x_data), length.out = cov_smooth_grid + 1)
  x_cov = (x_cov[1:(length(x_cov)-1)] + x_cov[2:length(x_cov)])/2

  mu_data = mu_data[cut_index]
  mu = spline(full_x_cov, mu, xout = x_cov)$y

  cov_surface_smoothed = akima(full_x_cov, x_cov, cov_surface_smoothed)
  # partial_surface_smoothed = akima(full_x_cov, x_cov, cov_surface_smoothed)

  rownames(cov_surface_smoothed) = x_cov
  colnames(cov_surface_smoothed) = x_cov
  # rownames(partial_surface_smoothed) = x_cov
  # colnames(partial_surface_smoothed) = x_cov

  rm(full_x_cov)

  data = data[,cut_index]
  kept_subjects = which(!apply(data, 1, function(x) all(is.na(x))))
  data = data[kept_subjects,]

  h = x_cov[2] - x_cov[1]
  junk = eigen(cov_surface_smoothed, symmetric = TRUE)
  index = Im(junk$values) == 0 & Re(junk$values) > 0
  eigen_values = junk$values[index]
  eigen_func_cov = junk$vectors[,index]
  rm(junk)

  maxk = min(maxk, sum(index))
  if (no_of_components > maxk) {
    warning('At most ', maxk, ' of PC can be selected!\n', sep = "")
    no_of_components = maxk
  }
  rm(index)

  lambda = h * eigen_values
  FVE = cumsum(lambda / sum(lambda))
  cat("FVE calculated from", length(lambda), "possible eigenvalues:\n ", FVE[1:min(5, maxk)], "...\n")
  if (missing(no_of_components) | no_of_components == -1) {
    if (selection == "FVE") {
      no_of_components = min(which(FVE >= FVE_threshold), maxk)
      if (no_of_components < 2) {
        no_of_components = 2
        warning("Use less than 2 components, reset to 2.\n")
      }
      cat("Best number of principal components selected by FVE: ", no_of_components, ".\n",
          "It accounts for ", FVE[no_of_components] * 100, "% of total variation (threshold = ", FVE_threshold, ").\n", sep = "")
    }
    if (selection == "BIC") {
      no_of_components = 2
      warning("No implement for BIC yet.")
    }
    if (selection == "AIC") {
      no_of_components = 2
      warning("No implement for AIC yet.")
    }
  }

  eigen_values = eigen_values[1:no_of_components]
  lambda = lambda[1:no_of_components]
  eigen_func_cov = eigen_func_cov[,1:no_of_components]

  for (i in 1:ncol(eigen_func_cov)) {
    eigen_func_cov[,i] = eigen_func_cov[,i] / sqrt(trapz(x_cov, eigen_func_cov[,i] ^ 2))
    if (eigen_func_cov[1,i] > eigen_func_cov[2,i]) {
      eigen_func_cov[,i] = -eigen_func_cov[,i]
    }
  }

  eigen_func_data = matrix(as.double(NA), ncol = ncol(eigen_func_cov), nrow = length(x_data))
  for (i in 1:ncol(eigen_func_data)) {
    eigen_func_data[,i] = spline(x = x_cov, y = eigen_func_cov[,i], xout = x_data)$y
    eigen_func_data[,i] = eigen_func_data[,i] / sqrt(trapz(x_data, eigen_func_data[,i] ^ 2))
  }

  cov_surface_fitted = matrix(as.double(0), ncol = length(x_data), nrow = length(x_data))
  for (i in 1:ncol(eigen_func_data)) {
    cov_surface_fitted = cov_surface_fitted +
      lambda[i] * eigen_func_data[,i] %*% t(eigen_func_data[,i])
  }
  cor_surface_fitted = cov2cor(cov_surface_fitted)
  rownames(cov_surface_fitted) = x_data
  colnames(cov_surface_fitted) = x_data
  rownames(cor_surface_fitted) = x_data
  colnames(cor_surface_fitted) = x_data

  png('corr_plot.png')
  image(x = x_data, y = x_data, z = cor_surface_fitted, col = heat.colors(100), main = "Correlation Plot")
  dev.off()

  # Part IV: Perform principal components analysis.
  cat('Part IV: Perform principal components analysis.\n')
  if (error) {
    # temporary
    sigma = sigma_est(bandwidth_mean* sqrt(2), kernel, full_x_data, raw_cov, cov_surface_smoothed, cut_index, x_cov)
    # sigma = sigma_est(kernel, full_x_data, raw_cov, cut_index, 
    #   bandwidth_mean_cv_mode, bandwidth_mean_cv_mode)
  }
  else {
    sigma = 0
  }

  hash = apply(data, 1, function(x) do.call(paste, as.list(as.integer(!is.na(x)))))
  unique_hash = unique(hash)
  # original getOriCurves.m

  # set up hash table to identify missing situation
  # combin same missing situation subjects to improve speed
  # not meanful when subjects are all different

  # always update sigma, means no case rho = -1 in paper

  # if needed
  # if (rho != -1)
  #   ...
  # else
  # sigma_new_2 = sigma

  if (error && new_ridge) {
    # \hat{\sigma}^2_{new,1} in FPCscore.pdf Step 1
    sigma_new_1 = pcs_est(data, mu_data, hash, unique_hash, lambda, eigen_func_data,
                          cov_surface_fitted, sigma, sigma_new_switch = TRUE)$sigma_new

    # \hat{\sigma}^2_{new,2} in FPCscore.pdf Step 2
    sigma_new_2 = pcs_est(data, mu_data, hash, unique_hash, lambda, eigen_func_data,
                          cov_surface_fitted, sigma_new_1, sigma_new_switch = TRUE)$sigma_new
    # Step 3
    if (rho_cv) {
      rho = cv_rho(data, mu_data, lambda, eigen_func_data, cov_surface_fitted)
    }
    sigma_new_2 = max(sigma_new_2, rho)
  }
  else {
    sigma_new_2 = 0
  }
  # original getScores1.m
  junk = pcs_est(data, mu_data, hash, unique_hash, lambda, eigen_func_data,
                 cov_surface_fitted, sigma_new_2, xi_var_switch = TRUE)

  if (!error) {
    sigma = NA
    rho_opt = NA
    sigmanew = NA
  }
  if (sparse_flag) {
    regular = paste(regular, '(convert from sparse)')
  }

  ans = list(
    # mean related
    x_data = x_data, mu_data = mu_data,
    x_cov = x_cov, mu_cov = mu_cov,
    # covariance(correlation) related
    raw_cov = raw_cov,
    cov_surface_smoothed = cov_surface_smoothed,
    # partial_surface_smoothed = partial_surface_smoothed,
    cov_surface_fitted = cov_surface_fitted,
    cor_surface_fitted = cor_surface_fitted,
    # PCA related
    no_of_components = no_of_components, lambda = lambda,
    eigen_func_data = eigen_func_data,
    eigen_func_cov = eigen_func_cov,
    xi_est = junk$xi_est,
    y_pred = junk$y_pred,
    xi_var = junk$xi_var,
    kept_subjects = kept_subjects,
    # etc
    sigma = sigma, regular = regular,
    AIC = NA, BIC = NA, FVE = FVE,
    # parameters for restricting the \hat{\sigma}^2_{new,2}
    rho_opt = rho, sigmanew = sigma_new_2)
  class(ans) = "FPCA"
  ans
}
