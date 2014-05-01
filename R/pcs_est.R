pcs_est = function(data, mu, hash, unique_hash, lambda,
                   eigen_func_data, cov_surface_fitted, sigma_old,
                   # switch
                   # no switch for xi_est and y_pred since those are always needed
                   sigma_new_switch = FALSE, xi_var_switch = FALSE) {
  sigma_new = 0
  xi_est = matrix(NA, ncol = ncol(eigen_func_data), nrow = nrow(data))
  y_pred = matrix(NA, ncol = ncol(data), nrow = nrow(data))
  xi_var = list()
  # hash table acceleration does not work for very different missing pattern
  for (hash_iterator in unique_hash) {
    index = hash == hash_iterator
    flag = unlist(strsplit(hash_iterator, " ")) == 1
    part_cov_surface_fitted = cov_surface_fitted[flag, flag] + diag(sigma_old, sum(flag))
    p = ncol(part_cov_surface_fitted)
    if (sum(flag) > 1) {
      junk = .Internal(La_svd('A', part_cov_surface_fitted, double(p),
                              matrix(double(p^2),p), matrix(double(p^2),p)))
      # eval this middle_matrix since it's need for both xi_est and xi_var
      middle_matrix = diag(lambda) %*% t(eigen_func_data[flag,]) %*% t(junk$vt) %*% diag(1 / junk$d) %*% t(junk$u)
      if (sum(index) > 1) {
        junk_y = as.matrix(sweep(data[index,], 2, mu, '-')[,flag])
      }
      else {
        junk_y = matrix((data[index,] - mu)[flag], nrow = 1)
      }
      xi_est[index,] = t(middle_matrix %*% t(junk_y))
    }
    else {
      junk = 1 / part_cov_surface_fitted[1,1]
      middle_matrix = diag(lambda) %*% matrix(eigen_func_data[flag,], ncol = 1) * junk
      if (sum(index) > 1) {
        junk_y = as.matrix(sweep(data[index,], 2, mu, '-')[,flag])
      }
      else {
        junk_y = matrix((data[index,] - mu)[flag], nrow = 1)
      }
      xi_est[index,] = sapply(junk_y, function(junk_yy) middle_matrix * junk_yy)
    }
    # k: number of components
    # p: grids of output without missing (=sum(flag))

    # diag(lambda):                                  k *          k
    # t(eigen_func_data[flag,]):                     k *          p
    # pinv:                                          p *          p
    # t(junk_y):                                     p * sum(index)
    # xi_est[index,]:                                sum(index) * k

    # be careful here
    
    # if (sum(flag) > 1) {
    y_pred[index, ] = xi_est[index,] %*% t(eigen_func_data[,])
    # }
    # else {
    #   y_pred[index, flag] = xi_est[index,] %*% matrix(eigen_func_data[flag,], ncol = 1)
    # }

    # xi_est[index,]:                     sum(index) * k
    # t(eigen_func_data[flag,]):          k          * p
    # y_pred[index, flag]:                sum(index) * p
    if (sigma_new_switch) {
      sigma_new = sigma_new + sum(apply(junk_y - y_pred[index, flag], 1, function(x) mean(x^2)))
    }
    if (xi_var_switch) {
      for (ii in which(index)) {
        if (sum(flag) > 1) {
          xi_var[[ii]] = diag(lambda) - middle_matrix %*% t(diag(lambda) %*% t(eigen_func_data[flag,]))
        }
        else {
          xi_var[[ii]] = diag(lambda) - middle_matrix %*% t(diag(lambda) %*% matrix(eigen_func_data[flag,], ncol = 1))
        }
      }
    }
  }
  y_pred = sweep(y_pred, 2, mu, '+')
  if (!sigma_new_switch) {
    sigma_new = NA
  }
  if (!xi_var_switch) {
    xi_var = NA
  }
  list(xi_est = xi_est,
       y_pred = y_pred,
       sigma_new = sigma_new / nrow(data),
       xi_var = xi_var)
}
