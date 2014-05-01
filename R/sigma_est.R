# temporary
sigma_est = function(bandwidth, kernel, full_x_data, raw_cov, cov_surface_smoothed, cut_index, x_cov) {
  x_diag = full_x_data[cut_index]
  var_y = attr(raw_cov, 'vary')[cut_index]
  var_yw = attr(raw_cov, 'varyw')[cut_index]
  var_y = lwls(bandwidth, kernel, x_diag, var_y, var_yw, x_cov)$output
  
  var_x = diag(cov_surface_smoothed)
  sigma = trapz(x_diag, var_y - var_x) / (max(x_diag) - min(x_diag))
  if (sigma < 0) {
    cat('Estimated sigma is negative, reset to zero now!\n')
  }
  max(0, sigma)
}

# sigma_est = function(kernel, full_x_data, raw_cov, cut_index, diag_cv_mode, q_cv_mode) {
#   x_diag = full_x_data[cut_index]
#   # y_in
#   var_y = attr(raw_cov, 'vary')[cut_index]
#   var_yw = attr(raw_cov, 'varyw')[cut_index]
#   
#   bandwidth = bandwidth_choice_1d(
#     kernel, x_diag, var_y, var_yw, degree = 1, cv_mode = diag_cv_mode)
# 
#   # variance of y
#   var_y = lwls(bandwidth, kernel, x_diag, var_y, var_yw)$output
# 
#   # variance from x
#   s = min(which(cut_index))
#   e = max(which(cut_index))
#   var_x = qlwls(kernel, full_x_data, s, e, raw_cov, q_cv_mode)
# 
#   # variance by yself
#   sigma = trapz(x_diag, var_y - var_x) / (max(x_diag) - min(x_diag))
# 
#   if (sigma < 0) {
#     cat('Estimated sigma is negative, reset to zero now!\n')
#   }
# 
#   max(0, sigma)
# }
# 
# qlwls = function(kernel, full_x_data, s, e, raw_cov, cv_mode) {
#   k = e - s + 1
#   ss = 2 * s - 1
#   qm = matrix(NA, nrow = k, ncol = ss)
#   wm = matrix(NA, nrow = k, ncol = ss)
#   for (i in s:e) {
#     for (j in -(s-1):(s-1)) {
#       qm[i - s + 1, j + s] = raw_cov[i - j, i + j]
#       wm[i - s + 1, j + s] = attr(raw_cov, 'weight')[i - j, i + j]
#     }
#   }
# 
#   attr(qm, "weight") = wm
#   junkqm = t(qm)
#   attr(junkqm, "weight") = t(wm)
#   # after fix lwls, change degree to 2
#   bandwidth = bandwidth_choice_seq(kernel, 1:ss, junkqm, 1, cv_mode)
# 
#   # after fix lwls, change degree to 2  
#   for (i in 1:nrow(qm)) {
#     qm[i, (ss + 1) / 2] = lwls(bandwidth, kernel, 1:ss, 
#                           qm[i,], w_in = attr(qm, "weight")[i,],
#                           x_out = full_x_data[(ss + 1) / 2], 
#                           drv = 0, degree = 1, cv_mode = 0)$output
#   }
#   as.double(qm[, (ss + 1) / 2])
# }
# 
