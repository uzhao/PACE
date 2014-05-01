#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A11 d0
#define A12 d1
#define A21 d1
#define A13 d2
#define A22 d2
#define A31 d2
#define A14 d3
#define A23 d3
#define A32 d3
#define A41 d3
#define A24 d4
#define A33 d4
#define A42 d4
#define A34 d5
#define A43 d5
#define A44 d6

int lwls(
     // input stuff
     double *bandwidthP, int *kernelP, double *x_in, double *y_in, double *count_in,
     // data describing stuff
     int *n_inP, int *n_outP, 
     // output stuff
     double *x_out, double *mu_out, double *output, double *weight_out,
     // derivative and cross validation stuff
     int *drvP, int *degreeP, int *cv_modeP, double *cv_value) {
  
  double bandwidth = *bandwidthP;
  int kernel = *kernelP, n_in = *n_inP, n_out = *n_outP, 
    drv = *drvP, degree = *degreeP, cv_mode = *cv_modeP;
  
  // drv: 0 mean
  //      1 1st derivative
  //      2 2nd derivative
  
  // cv_mode: 0 non-cv
  //          1 ocv
  //          2 gcv
  //          3 geometric mean
  
  // declaration for cv
  double *h = (double *) calloc(n_out, sizeof(double));
  double sumh = 0;

  // declaration for 2d weight
  double *temp_weight = NULL;
  temp_weight = (double *) calloc(n_in, sizeof(double));

  double d0 = 0, d1 = 0, d2 = 0, d3 = 0, d4 = 0, d5 = 0, d6 = 0,
  b1 = 0, b2 = 0, b3 = 0, b4 = 0, cur_x_out = 0, d = 0, d_bw = 0, w = 0, min_bw = x_in[n_in - 1];
  int in_index = 0, out_index = 0, i = 0;
  
  if (bandwidth < 0) {
    fprintf(stderr, "Bandwidth choice for mu(t) and/or its derivative must be positive!\n");
    return 1;
  }

  // loop over out point, set output
  for (out_index = 0; out_index < n_out; out_index++) {
    d0 = d1 = d2 = d3 = d4 = d5 = d6 = b1 = b2 = b3 = b4 = 0;
    cur_x_out = x_out[out_index];
    for (in_index = 0; in_index < n_in; in_index++) {
      if (count_in[in_index] == 0) {
        continue;
      }
      if (kernel == 0 || kernel == 1 || kernel == 4) {
        // epan, rect, or quar
        if (x_in[in_index] < cur_x_out - bandwidth) {
          continue;
        }
        if (x_in[in_index] > cur_x_out + bandwidth) {
          break;
        }
      }
      else {
        // gauss, or gausvar
        if (x_in[in_index] < cur_x_out - 5 * bandwidth) {
          continue;
        }
        if (x_in[in_index] > cur_x_out + 5 * bandwidth) {
          break;
        }
      }

      d = cur_x_out - x_in[in_index];
      d_bw = d / bandwidth;

      switch (kernel) {
          // here w is the real weight
        case 0:
          w = count_in[in_index] * (1 - d_bw * d_bw);
          break;
        case 1:
          w = count_in[in_index];
          break;
        case 2:
          w = count_in[in_index] / exp(d_bw * d_bw / 2);
          break;
        case 3:
          d_bw = d_bw * d_bw;
          w = count_in[in_index] / exp(d_bw / 2) * (1.25 - 0.25 * d_bw);
          break;
        case 4:
          w = count_in[in_index] * (1 - 2 * d_bw * d_bw + d_bw * d_bw * d_bw * d_bw);
          break;
      }

      temp_weight[in_index] = w;
      d0 += w;
      d1 += w * d;
      d2 += w * d * d;
      b1 += w * y_in[in_index];
      b2 += w * d * y_in[in_index];

      if (degree > 1) {
        d3 += w * d * d * d;
        d4 += w * d * d * d * d;
        b3 += w * d * d * y_in[in_index];
      }
      if (degree > 2) {
        d5 += w * d * d * d * d * d;
        d6 += w * d * d * d * d * d * d;
        b4 += w * d * d * d * y_in[in_index];
      }
    }

    
    for (i = 0; i != n_in; i++) {
      weight_out[out_index] += temp_weight[i] * temp_weight[i] / count_in[i];
      temp_weight[i] = 0;
    }
    weight_out[out_index] = 1 / weight_out[out_index];

    // calculating cv value
    if (degree == 1) {
      // if degree == 1, we only estimate mean
      output[out_index] = (b1*d2 - b2*d1)/(- d1*d1 + d0*d2);
      if (cv_mode) {
        h[out_index] = count_in[out_index] * A22 / (A11 * A22 - A12 * A12);
        mu_out[out_index] = output[out_index];
      }
    }
    if (degree == 2) {
      // if degree == 2, we only estimate mean or 1st drv
      if (drv == 0) {
        output[out_index] = 0;
      }
      else {
        output[out_index] = -(b2*d2*d2 + b1*d1*d4 - b1*d2*d3 - b2*d0*d4 + b3*d0*d3 - b3*d1*d2)/(d4*d1*d1 - 2*d1*d2*d3 + d2*d2*d2 - d0*d4*d2 + d0*d3*d3);
      }
      if (cv_mode) {
        h[out_index] = 0;
        mu_out[out_index] = 0;
      }
    }
    if (degree == 3) {
      // if degree == 2, we only estimate mean or 2nd drv
      if (drv == 0) {
        output[out_index] = 0;
      }
      else {
      output[out_index] = 2*(b3*d0*d4*d4 - b1*d2*d4*d4 - b2*d3*d3*d3 + b1*d3*d3*d4 + b3*d2*d3*d3 + b4*d1*d3*d3 + b1*d2*d2*d6 - b4*d2*d2*d3 + b3*d1*d1*d6 - b4*d1*d1*d5 - b1*d1*d3*d6 + b1*d1*d4*d5 - b1*d2*d3*d5 + b2*d0*d3*d6 - b2*d0*d4*d5 - b2*d1*d2*d6 + b2*d1*d3*d5 + b2*d2*d3*d4 - b3*d0*d2*d6 - 2*b3*d1*d3*d4 + b4*d0*d2*d5 - b4*d0*d3*d4 + b4*d1*d2*d4)/(d6*d1*d1*d4 - d1*d1*d5*d5 - 2*d6*d1*d2*d3 + 2*d1*d2*d4*d5 + 2*d1*d3*d3*d5 - 2*d1*d3*d4*d4 + d6*d2*d2*d2 - 2*d2*d2*d3*d5 - d2*d2*d4*d4 + 3*d2*d3*d3*d4 - d0*d6*d2*d4 + d0*d2*d5*d5 - d3*d3*d3*d3 + d0*d6*d3*d3 - 2*d0*d3*d4*d5 + d0*d4*d4*d4);
      }
      if (cv_mode) {
        h[out_index] = 0;
        mu_out[out_index] = 0;
      }
    }
  }
  //end of out_index loop
  
  // gcv or geomean
  if (cv_mode > 1) {
    for (i = 0; i < n_in; i++) {
      sumh += h[i];
    }
    sumh /= n_in;
    for (i = 0; i < n_in; i++) {
      h[i] = sumh;
    }
  }
  
  // assign cv value
  if (cv_mode) {
    *cv_value = 0;
    for (i = 0; i < n_in; i++) {
      (*cv_value) += (y_in[i] - mu_out[i]) * (y_in[i] - mu_out[i]) / (1 - h[i]) / (1 - h[i]);
    }
    (*cv_value) /= n_in;
    if (isnan(*cv_value)) {
      *cv_value = INFINITY;
    }
  }
  
  free(h);
  free(temp_weight);
  h = NULL;
  temp_weight = NULL;

  return 0;
}
