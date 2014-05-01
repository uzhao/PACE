cv_rho = function(data, mu_data, lambda, eigen_func_data, cov_surface_fitted, alpha_range = c(0.01, 0.225)) {
  # clean = function(one, full) {
  #   full[full != one]
  # }
  all_time = as.numeric(colnames(data))
  time_range = range(all_time)
  gamma = sqrt(trapz(all_time, mu_data^2) + sum(lambda)) / sqrt(time_range[2] - time_range[1])
  exist_matrix = !is.na(data)
  exist_list = apply(exist_matrix, 1, which)
  if (class(exist_list) == "list") {
    more_than_one = sapply(exist_list, length) > 1
  }
  else {
    more_than_one = apply(exist_matrix, 1, sum) > 1
  }
  data = data[more_than_one,]
  exist_list = exist_list[more_than_one]
  if (class(exist_list) == "list") {
    leave_one = sapply(exist_list, sample ,1)
  }
  else {
    leave_one = apply(exist_matrix, 1, function(x) sample(which(x), 1))
  }
  # keep = mapply(clean, leave_one, exist_list)
  y_true = rep(0, nrow(data))
  data_cv = data
  # may improve in case n is too large
  for (i in 1:nrow(data_cv)) {
    y_true[i] = data_cv[i, leave_one[i]]
    data_cv[i, leave_one[i]] = NA
  }
  #
  hash = apply(data_cv, 1, function(x) do.call(paste, as.list(as.integer(!is.na(x)))))
  unique_hash = unique(hash)
  ad_hoc_cv_func = function(rho) {
    xi_est = pcs_est(data_cv, mu_data, hash, unique_hash, lambda, eigen_func_data, cov_surface_fitted, rho)$xi_est
    pred_y = mu_data[leave_one]
    for (i in 1:length(lambda)) {
      pred_y = pred_y + xi_est * eigen_func_data[leave_one, i]
    }
    sum((y_true - pred_y)^2)
  }
  ans = optimize(ad_hoc_cv_func, alpha_range * gamma)$minimum
  cat("rho for ridge regression is:", ans, "\n")
  ans
}
