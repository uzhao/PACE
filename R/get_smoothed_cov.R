get_smoothed_cov = function(kernel, x_data, raw_cov, x_cov, bandwidth_cov_cv_mode) {
  # 1st direction
  bandwidth_cov_adhoc = bandwidth_choice_seq(
    kernel, x_data, raw_cov, degree = 1, cv_mode = bandwidth_cov_cv_mode) * 0.7
  cov_surface_smoothed = lwls_seq(bandwidth_cov_adhoc, kernel, 
                                  x_data, raw_cov, x_cov, 
                                  cv_mode = 0, verbose = "cov matrix column")
  cov_surface_smoothed = t(cov_surface_smoothed)
  attr(cov_surface_smoothed, "weight") = t(attr(cov_surface_smoothed, "weight"))

  # 2nd direction
  bandwidth_cov_adhoc = bandwidth_choice_seq(
    kernel, x_data, cov_surface_smoothed, degree = 1, cv_mode = bandwidth_cov_cv_mode) * 0.7
  cov_surface_smoothed = lwls_seq(bandwidth_cov_adhoc, kernel, 
                                  x_data, cov_surface_smoothed, x_cov, 
                                  cv_mode = 0, verbose = "cov matrix row")
  cov_surface_smoothed = t(cov_surface_smoothed)  
  attr(cov_surface_smoothed, "weight") = NULL

  # make it symmetric
  (cov_surface_smoothed + t(cov_surface_smoothed)) / 2  
}

get_smoothed_partial = function(kernel, x_data, raw_cov, x_cov, bandwidth_cov_cv_mode) {
  # 1st direction
  bandwidth_cov_adhoc = bandwidth_choice_seq(
    # after fix lwls, change degree to 2
    kernel, x_data, raw_cov, degree = 1, cv_mode = bandwidth_cov_cv_mode)
  cov_surface_smoothed = lwls_seq(bandwidth_cov_adhoc, kernel, 
                                  x_data, raw_cov, x_cov, drv = 1, degree = 2, 
                                  cv_mode = 0, verbose = "cov matrix column")
  cov_surface_smoothed = t(cov_surface_smoothed)
  attr(cov_surface_smoothed, "weight") = t(attr(cov_surface_smoothed, "weight"))

  # 2nd direction
  bandwidth_cov_adhoc = bandwidth_choice_seq(
    # after fix lwls, change degree to 2
    kernel, x_data, cov_surface_smoothed, degree = 1, cv_mode = bandwidth_cov_cv_mode)
  cov_surface_smoothed = lwls_seq(bandwidth_cov_adhoc, kernel, 
                                  x_data, cov_surface_smoothed, x_cov, drv = 1, degree = 2, 
                                  cv_mode = 0, verbose = "cov matrix row")
  cov_surface_smoothed = t(cov_surface_smoothed)  
  attr(cov_surface_smoothed, "weight") = NULL

  # make it symmetric
  (cov_surface_smoothed + t(cov_surface_smoothed)) / 2  
}
