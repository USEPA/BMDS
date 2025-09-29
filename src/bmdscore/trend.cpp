#include "trend.h"

#include "bmds_helper.h"

int imax(int i, int j) { return (i > j) ? i : j; }

double ***w;  // used for storing cwilcox(i,j,jk) -> w[i][j][k]
int allocated_m, allocated_n;

void w_free(int m, int n) {
  int i, j;
  for (i = m; i >= 0; i--) {
    for (j = n; j >= 0; j--) {
      if (w[i][j] != 0) free((void *)w[i][j]);
    }
    free((void *)w[i]);
  }
  free((void *)w);
  w = 0;
  allocated_m = allocated_n = 0;
}

void w_init(int m, int n) {
  int i;
  if (m > n) {
    i = n;
    n = m;
    m = i;
  }

  // zero out w if needed
  if (w && (m > allocated_m || n > allocated_n)) w_free(allocated_m, allocated_n);
  // initialize w[][] if needed
  if (!w) {
    m = imax(m, WILCOX_MAX);
    n = imax(n, WILCOX_MAX);
    w = (double ***)calloc((size_t)m + 1, sizeof(double **));
  }

  if (!w) std::cout << "error in wilcox allocation" << std::endl;
  for (i = 0; i <= m; i++) {
    w[i] = (double **)calloc((size_t)n + 1, sizeof(double *));
  }
  allocated_m = m;
  allocated_n = n;
}

double cwilcox_new(int k, int m, int n) {
  int c, i, j;
  int u = m * n;
  if (k < 0 || k > u) return (0);
  c = (int)(u / 2);
  if (k > c) k = u - k;
  if (m < n) {
    i = m;
    j = n;
  } else {
    i = n;
    j = m;
  }
  if (j == 0) return (k == 0);

  // Simplify if k is small
  if (j > 0 && k < j) return cwilcox_new(k, i, k);

  if (w[i][j] == 0) {
    w[i][j] = (double *)calloc((size_t)c + 1, sizeof(double));
    if (!w[i][j]) std::cout << "wilcox allocation error" << std::endl;
    for (int l = 0; l <= c; l++) {
      w[i][j][l] = -1;
    }
  }
  if (w[i][j][k] < 0) {
    if (j == 0)
      w[i][j][k] = (k == 0);
    else
      w[i][j][k] = cwilcox_new(k - j, i - 1, j) + cwilcox_new(k, i, j - 1);
  }
  return (w[i][j][k]);
}

int cwilcox(int k, int m, int n) {
  // reserve memory for w
  w_init(m, n);

  // call cwilcox
  double ret = cwilcox_new(k, m, n);

  // free w
  if (m > WILCOX_MAX || n > WILCOX_MAX) w_free(m, n);

  return ((int)ret);
}
