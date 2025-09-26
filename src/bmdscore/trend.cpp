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

// int cwilcox(int k, int m, int n) {
//
////	Args:
////		k(int): Sum to achieve
////		m(int): Max value for each element
////		n(int): Max number of elements
////	Returns:
////		The number of ways to obtain the sum with the given inputs
//
////	std::cout<<"Starting cwilcox with k:"<<k<<", m:"<<m<<", n:"<<n<<std::endl;
//	std::vector<std::vector<std::vector<int>>> w(m+1, std::vector<std::vector<int>>(n+1,
// std::vector<int>(1, int(BMDS_MISSING))));
//
//	int c, i, j;
//	//Calculate the total number of observations
//	int u = m*n;
//	//Check if k is too large or too small
//	if (k<0 || k>u) return 0;
//	//Calculate the half of the total observations
//	c = (int) (u/2);
//	//Adjust k if it's greater than half
//	if (k>c) k=u-k;
//	//Ensure i<j
//	if (m<n) {
//	   i=m; j=n;
//	} else {
//	   i=n; j=m;
//	}
//
//	//Check if j is 0
//	if (j==0) return (k==0);
//	//If k is less than the size of the second sample, recurse
//	if (j>0 && k<j) return cwilcox(k,i,k);
//	//Initialize array if not already done
//	if (w[i-1][j-1][0] == int(BMDS_MISSING)) {
////		std::cout<<"do stuff"<<std::endl;
//	   w[i-1][j-1] = std::vector<int>(c+1,-1);
//	}
//        //If the value has not been computed yet
////	std::cout<<"i:"<<i<<", j:"<<j<<", k:"<<k<<std::endl;
//	if (w[i-1][j-1][k-1] < 0) {
//	   if (j==0){
//	      w[i-1][j-1][k-1] = int(k==0);
//	   } else {
//	      w[i-1][j-1][k-1] = cwilcox(k-j, i-1, j) + cwilcox(k, i, j-1);
//	   }
//	}
//
////	std::cout<<"Returning w="<<w[i-1][j-1][k-1]<<std::endl;
//	return w[i-1][j-1][k-1];
//}
