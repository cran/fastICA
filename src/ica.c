#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int min_JM (int *, int *);
static int max_JM (int *, int *);

static void rowcentre_JM (double *, int, int);
static void colstandard_JM (double *, int, int);
static void rowstd_JM (double *, int, int, int);
static void transpose_mat_JM (double *, int *, int *, double *);
static void mmult_JM (double *, int, int, double *, int, int, double *);
static void orthog_mat_JM (double *, int, double *);
static void gramsch_JM (double *, int, int, int);
static void svd_JM (double *, int *, int *, double *, double *, double *);

static void Symm_logcosh_JM (double *, int, double *, int, int, double, double *, double *);
static void Symm_exp_JM (double *, int, double *, int, int, double, double *, double *);
static void Def_logcosh_JM (double *, int, double *, int, int, double, double *);
static void Def_exp_JM (double *, int, double *, int, int, double, double *);
static void calc_A_JM(double*, double*, double*, int*, int*, int*, double*, double*);
static void calc_K_JM(double*, int*, int*, double*);

#include <R_ext/Lapack.h>

static void
rowcentre_JM (double *ans, int n, int p)
{
/*  mean centres nxp matrix ans */
    double tmp;
    int i, j;
    for (i = 0; i < n; i++) {
	tmp = 0;
	for (j = 0; j < p; j++) {
	    tmp = tmp + ((double) ans[p * i + j]) / p;
	}
	for (j = 0; j < p; j++) {
	    ans[p * i + j] -= (double) tmp;
	}
    }
}

static void
colstandard_JM (double *ans, int n, int p)
{
/*  transform columns of nxp matrix ans to have zero mean and unit variance */
    double tmp[2];
    double tmp1;
    int i, j;
    for (i = 0; i < p; i++) {
	tmp[0] = 0;
	tmp[1] = 0;

	for (j = 0; j < n; j++) {
	    tmp[0] += (double) ans[p * j + i];
	    tmp[1] += ((double) ans[p * j + i]) * ((double) ans[p * j + i]);
	}

	tmp[0] = tmp[0] / n;
	tmp1 = (tmp[1] - n * (tmp[0]) * (tmp[0])) / (n - 1);

	tmp[1] = sqrt (tmp1);
	for (j = 0; j < n; j++) {
	    ans[p * j + i] =
		(double) ((((double) ans[p * j + i]) - tmp[0]) / tmp[1]);
	}
    }
}


static void
svd_JM (double *mat, int *n, int *p, double *u, double *d, double *v)
{

    /*  calculates svd decomposition of nxp matrix mat */
    /*    mat is a pointer to an nxp array of doubles */
    /*    n is a pointer to an integer specifying the no. of rows of mat */
    /*    p is a pointer to an integer specifying the no. of cols of mat */
    /*    u is a pointer to a double array of dimension (n,n) */
    /*    d is a pointer to a double array of dimension min(n,p) */
    /*    v is a pointer to a double array of dimension (p,p) */


    int info, *iwork, lwork, a, b;
    size_t iwork_size, ilwork, nn = *n, pp = *p, mm;
    double *work, *mat1, *u1, *v1;
    char jobz = 'A';

    mm = min_JM(n,p);
    iwork_size = 8 * (size_t) mm;

    a = max_JM(n,p);
    b = 4 * mm * mm + 4 * mm;
    ilwork = 3 * mm * mm + max_JM(&a, &b);
    if (ilwork > INT_MAX)
	error("svd on %d x %d exceeds Fortran indexing limits");

    work = Calloc (ilwork, double);
    iwork = Calloc (iwork_size, int);
    mat1 = Calloc (nn * pp, double);
    u1 = Calloc (nn * nn, double);
    v1 = Calloc (pp * pp, double);

    transpose_mat_JM (mat, n, p, mat1);

    lwork = ilwork;
    F77_CALL (dgesdd) (&jobz, n, p, mat1, n, d, u1, n, v1, p, work,
		       &lwork, iwork, &info);

    transpose_mat_JM (u1, n, n, u);
    transpose_mat_JM (v1, p, p, v);

    Free (mat1);
    Free (u1);
    Free (v1);
    Free (work);
    Free (iwork);
}


static void
transpose_mat_JM (double *mat, int *n, int *p, double *ans)
{
/*    transpose nxp matrix mat */
    int i, j;

    for (i = 0; i < *n; i++) {
	for (j = 0; j < *p; j++) {
	    *(ans + j * (*n) + i) = *(mat + i * (*p) + j);
	}
    }
}


static int min_JM (int *a, int *b)
{
/*  find minimum of a and b */
    int ans;

    ans = *b;
    if (*a < *b) ans = *a;

    return ans;
}

static int max_JM (int *a, int *b)
{
/*  find maximum of a and b */

    int ans;

    ans = *b;
    if (*a > *b) ans = *a;

    return ans;
}


static void
mmult_JM (double *A, int n, int p, double *B, int q, int r, double *C)
{
/*    matrix multiplication using FORTRAN BLAS routine SGEMM */
/*    A is (n*p) and B is (q*r), A*B returned to C  */

    double alpha = 1.0, beta = 0.0;
    int M, K, N;
    char transA = 'N', transB = 'N';

    if (p != q) {
	error ("Error, matrices not suitable\nfor multiplication");
    } else {
	M = n;
	K = p;
	N = r;
	F77_CALL (dgemm) (&transA, &transB, &N, &M, &K, &alpha, B, &N,
			  A, &K, &beta, C, &N);
    }
}

static void
orthog_mat_JM (double *mat, int e, double *orthog)
{
	/* take Wmat, (e*e), and return orthogonalized version to orthog_W */
	double *u, *v, *d, *temp;
	int i;
	size_t ee = e;


	u = Calloc (ee * ee, double);
	d = Calloc (ee, double);
	v = Calloc (ee * ee, double);
	temp = Calloc (ee * ee, double);

	svd_JM (mat, &e, &e, u, d, v);
	for (i = 0; i < e; i++) {
		temp[i * e + i] = 1 / (d[i]);
	}

	mmult_JM (u, e, e, temp, e, e, v);
	transpose_mat_JM (u, &e, &e, temp);
	mmult_JM (v, e, e, temp, e, e, u);
	mmult_JM (u, e, e, mat, e, e, orthog);


	Free (u);
	Free (v);
	Free (d);
	Free (temp);

}

static void
Symm_logcosh_JM (double *w_init, int e, double *data, int f, int p, double alpha, double *w_final, double *Tol)
{

	/* Function that carries out Symmetric ICA using a logcosh approximation to the neg. entropy function */

	double *mat1, *mat2, *mat3, *mat4, *mat5, *mat6;
	int i, j;
	double mean;

	if (e != f) {
		error ("error in Symm_logcosh_JM, dims dont match");
	}
	else {
		size_t es = (size_t)e * (size_t)e;
		size_t ep = (size_t)e * (size_t)p;
		mat1 = Calloc (ep, double);
		mat2 = Calloc (ep, double);
		mat3 = Calloc (es, double);
		mat4 = Calloc (es, double);
		mat5 = Calloc (es, double);
		mat6 = Calloc (es, double);

		mmult_JM (w_init, e, e, data, e, p, mat1);


		for (i = 0; i < e; i++) {
			for (j = 0; j < p; j++) {
				mat1[i * p + j] = tanh (alpha * mat1[i * p + j]);
			}
		}
		transpose_mat_JM (data, &e, &p, mat2);
		for (i = 0; i < e; i++) {
			for (j = 0; j < p; j++) {
				mat2[i * p + j] = (mat2[i * p + j]) / p;
			}
		}
		mmult_JM (mat1, e, p, mat2, p, e, mat3);
		for (i = 0; i < e; i++) {
			for (j = 0; j < p; j++) {
				mat1[i * p + j] =
					(alpha * (1 - (mat1[i * p + j]) * (mat1[i * p + j])));
			}
		}

		for (i = 0; i < e; i++) {
			mean = 0;
			for (j = 0; j < p; j++) {
				mean += ((mat1[i * p + j]) / p);
			}
			mat4[i * e + i] = mean;
		}
		mmult_JM (mat4, e, e, w_init, e, e, mat5);
		for (i = 0; i < e; i++) {
			for (j = 0; j < e; j++) {
				mat4[i * e + j] = (mat3[i * e + j] - mat5[i * e + j]);
			}
		}

		transpose_mat_JM (w_init, &e, &e, mat6);
		orthog_mat_JM (mat4, e, w_final);


		mmult_JM (w_final, e, e, mat6, e, e, mat5);
		mean = 0;
		for (i = 0; i < e; i++) {
			if (fabs (1 - fabs (mat5[i * e + i])) > mean) {
				mean = (fabs (1 - fabs (mat5[i * e + i])));
			}
		}
		*Tol = mean;
		Free (mat1);
		Free (mat2);
		Free (mat3);
		Free (mat4);
		Free (mat5);
		Free (mat6);
	}
}

static void
Def_logcosh_JM (double *w_init, int e, double *data, int f, int p, double alpha, double *w_final)
{
/* Function that carries out Deflation ICA using an logcosh approximation to the neg. entropy function */

	double *mat1, *mat2, *mat3, *mat4;
	int i, j;
	double mean;

	if (e != f) {
		error ("error in Def_logcosh_JM, dims dont match");
	}
	else {
		mat1 = Calloc (p, double);
		mat2 = Calloc ((size_t)e * (size_t)p, double);
		mat3 = Calloc (e, double);
		mat4 = Calloc (e, double);

		mmult_JM (w_init, 1, e, data, e, p, mat1);


		for (i = 0; i < p; i++) {
			mat1[i] = tanh (alpha * mat1[i]);
		}
		transpose_mat_JM (data, &e, &p, mat2);
		for (i = 0; i < e; i++) {
			for (j = 0; j < p; j++) {
				mat2[i * p + j] = (mat2[i * p + j]) / p;
			}
		}

		mmult_JM (mat1, 1, p, mat2, p, e, mat3);
		for (i = 0; i < p; i++) {
			mat1[i] = (alpha * (1 - (mat1[i]) * (mat1[i])));
		}

		mean = 0;
		for (j = 0; j < p; j++) {
			mean += ((mat1[j]) / p);
		}
		for (i = 0; i < e; i++) {
			mat4[i] = (w_init[i]) * mean;
		}
		for (i = 0; i < e; i++) {
			w_final[i] = (mat3[i] - mat4[i]);
		}

		Free (mat1);
		Free (mat2);
		Free (mat3);
		Free (mat4);

	}
}

static void
Symm_exp_JM (double *w_init, int e, double *data, int f, int p, double alpha, double *w_final, double *Tol)
		{
    /* Function that carries out Symmetric ICA using a exponential approximation to the neg. entropy function */

double *mat1, *mat2, *mat3, *mat4, *mat5, *mat0, *mat6;
int i, j;
double mean;

if (e != f) {
    error ("error in Symm_exp_JM, dims dont match");
}
else {
    size_t ep = (size_t)e * (size_t)p;
    size_t ee = (size_t)e * (size_t)e;
    mat0 = Calloc (ep, double);
    mat1 = Calloc (ep, double);
    mat2 = Calloc (ep, double);
    mat3 = Calloc (ee, double);
    mat4 = Calloc (ee, double);
    mat5 = Calloc (ee, double);
    mat6 = Calloc (ee, double);
    mmult_JM (w_init, e, e, data, e, p, mat1);
    for (i = 0; i < e; i++) {
	for (j = 0; j < p; j++) {
	    mat0[i * p + j] =
		(mat1[i * p + j]) * exp (-0.5 * (mat1[i * p + j]) *
					 (mat1[i * p + j]));
	}
    }
    transpose_mat_JM (data, &e, &p, mat2);
    for (i = 0; i < e; i++) {
	for (j = 0; j < p; j++) {
	    mat2[i * p + j] = (mat2[i * p + j]) / p;
	}
    }
    mmult_JM (mat0, e, p, mat2, p, e, mat3);
    for (i = 0; i < e; i++) {
	for (j = 0; j < p; j++) {
	    mat1[i * p + j] =
		((1 - (mat1[i * p + j]) * (mat1[i * p + j])) *
		 exp (-0.5 * (mat1 [i * p + j]) * (mat1 [i * p + j])));
	}
    }

    for (i = 0; i < e; i++) {
	mean = 0;
	for (j = 0; j < p; j++) {
	    mean += ((mat1[i * p + j]) / p);
	}
	mat4[i * e + i] = mean;
    }
    mmult_JM (mat4, e, e, w_init, e, e, mat5);
    for (i = 0; i < e; i++) {
	for (j = 0; j < e; j++) {
	    mat4[i * e + j] = (mat3[i * e + j] - mat5[i * e + j]);
	}
    }

    transpose_mat_JM (w_init, &e, &e, mat6);
    orthog_mat_JM (mat4, e, w_final);

    mmult_JM (w_final, e, e, mat6, e, e, mat5);
    mean = 0;
    for (i = 0; i < e; i++) {
	if (fabs (1 - fabs (mat5[i * e + i])) > mean) {
	    mean = (fabs (1 - fabs (mat5[i * e + i])));
	}
    }
    *Tol = mean;
    Free (mat1);
    Free (mat2);
    Free (mat3);
    Free (mat4);
    Free (mat5);
    Free (mat0);
    Free (mat6);
}
}

static void
Def_exp_JM (double *w_init, int e, double *data, int f, int p, double alpha, double *w_final)
{
    /* Function that carries out Deflation ICA using an exponential approximation to the neg. entropy function */

double *mat1, *mat2, *mat3, *mat4;
int i, j;
double mean;

if (e != f) {
    error ("error in Def_exp_JM, dims dont match");
}
else {
    mat1 = Calloc (p, double);
    mat2 = Calloc ((size_t)e * (size_t)p, double);
    mat3 = Calloc (e, double);
    mat4 = Calloc (e, double);

    mmult_JM (w_init, 1, e, data, e, p, mat1);

    for (i = 0; i < p; i++) {
	mat1[i] = ((mat1[i]) * exp (-0.5 * (mat1[i]) * (mat1[i])));
    }

    transpose_mat_JM (data, &e, &p, mat2);
    for (i = 0; i < e; i++) {
	for (j = 0; j < p; j++) {
	    mat2[i * p + j] = (mat2[i * p + j]) / p;
	}
    }

    mmult_JM (mat1, 1, p, mat2, p, e, mat3);

    mmult_JM (w_init, 1, e, data, e, p, mat1);
    for (i = 0; i < p; i++) {
	mat1[i] =
	    ((1 -
	      (mat1[i]) * (mat1[i])) * exp (-.5 * (mat1[i]) * (mat1[i])));
    }
    mean = 0;
    for (j = 0; j < p; j++) {
	mean += ((mat1[j]) / p);
    }
    for (i = 0; i < e; i++) {
	mat4[i] = (w_init[i]) * mean;
    }
    for (i = 0; i < e; i++) {
	w_final[i] = (mat3[i] - mat4[i]);
    }


    Free (mat1);
    Free (mat2);
    Free (mat3);
    Free (mat4);

}
}

static void
gramsch_JM (double *ww, int n, int m, int k)
{
int ip, jp;
double tmp;
/* do Gram-Schmidt on row k of (n*m) matrix ww */
k -= 1;
if (k > n) {
    error ("Error in gramsch");
}
else {
    for (ip = 0; ip < k; ip++) {
	tmp = 0;
	for (jp = 0; jp < m; jp++) {
	    tmp += ((ww[m * ip + jp]) * (ww[m * k + jp]));
	}
	for (jp = 0; jp < m; jp++) {
	    ww[m * k + jp] = (ww[m * k + jp] - ((ww[m * ip + jp]) * tmp));
	}
    }
}
}

static void
rowstd_JM (double *ww, int n, int m, int k)
{
/* for ww (n*m), make ||ww[k, ]|| equal 1 */
double tmp = 0;
int i;
k -= 1;
if (k > n) {
    error ("Error in rowstd");
}
else {
    for (i = 0; i < m; i++) {
	tmp += ((ww[k * m + i]) * (ww[k * m + i]));
    }
    tmp = sqrt (tmp);
    for (i = 0; i < m; i++) {
	ww[k * m + i] = ((ww[k * m + i]) / tmp);
    }
}
}


static void
calc_K_JM(double *x, int *n, int *p, double *K)
{
    int i, j;
    double *xxt, *xt, *u, *d, *v, *temp1, *temp2;
    size_t nn = *n, pp = *p;

    xxt = Calloc (nn * nn, double);
    xt = Calloc (nn * pp, double);

    /* transpose x matrix */
    transpose_mat_JM (x, n, p, xt);

    /* calculate sample covariance matrix xxt */
    mmult_JM (x, *n, *p, xt, *p, *n, xxt);
    for (i = 0; i < *n; i++) {
	for (j = 0; j < *n; j++) {
	    xxt[*n * i + j] = xxt[*n * i + j] / *p;
	}
    }
    Free (xt);

    /* calculate svd decomposition of xxt */
    u = Calloc (nn * nn, double);
    d = Calloc (nn, double);
    v = Calloc (nn * nn, double);

    svd_JM (xxt, n, n, u, d, v);

    /* calculate K matrix*/
    temp1 = Calloc (nn * nn, double);
    temp2 = Calloc (nn * nn, double);

    for (i = 0; i < *n; i++) {
	temp1[*n * i + i] = 1 / sqrt (d[i]);
    }

    transpose_mat_JM (u, n, n, temp2);
    mmult_JM (temp1, *n, *n, temp2, *n, *n, K);

    Free (temp1);
    Free (temp2);
    Free(xxt);
    Free(u);
    Free(d);
    Free(v);
}

static void
calc_A_JM(double *w, double *k, double *data,
	  int *e, int *n, int *p, double *A, double *unmixed_data)
{
    /* calculate un-mixing matrix A */
    int i;
    double *um, *umt, *umumt, *uu, *dd, *vv, *temp1, *temp2, *temp3;
    size_t nn = *n, ee = *e;

    um = Calloc (ee * nn, double);
    umt = Calloc (nn * ee, double);

    mmult_JM (w, *e, *e, k, *e, *n, um);
    mmult_JM (um, *e, *n, data, *n, *p, unmixed_data);
    transpose_mat_JM (um, e, n, umt);

    umumt = Calloc (ee * ee, double);
    mmult_JM (um, *e, *n, umt, *n, *e, umumt);

    uu = Calloc (ee * ee, double);
    dd = Calloc (ee, double);
    vv = Calloc (ee * ee, double);
    svd_JM (umumt, e, e, uu, dd, vv);

    temp1 = Calloc (ee * ee, double);
    for (i = 0; i < *e; i++) {
	temp1[*e * i + i] = 1 / (dd[i]);
    }

    temp2 = Calloc (ee * ee, double);
    temp3 = Calloc (ee * ee, double);
    transpose_mat_JM (vv, e, e, temp3);
    mmult_JM (temp3, *e, *e, temp1, *e, *e, temp2);
    transpose_mat_JM (uu, e, e, vv);
    mmult_JM (temp2, *e, *e, vv, *e, *e, uu);

    mmult_JM (umt, *n, *e, uu, *e, *e, A);

    Free(um);
    Free(umt);
    Free(umumt);
    Free(uu);
    Free(dd);
    Free(vv);
    Free(temp1);
    Free(temp2);
    Free(temp3);

}

static void
icainc_JM (double *data_matrix, double *w_matrix, int *nn, int *pp, int *ee,
	double *alpha, int *rowflag, int *colflag, int *funflag, int *maxit,
	double *lim, int *defflag, int *verbose, double *data_pre, double *Kmat1,
	double *w_final, double *ansa, double *ansx2)
{

    /* main ICA function */

    int i, j, k;
    size_t n = *nn, p = *pp, e = *ee;
    double tol;
    double *temp_w1, *temp_w2;
    double *data1, *Kmat, *temp1, *w_init;

    /* make a copy of the data matrix*/
    data1 = Calloc (n * p, double);
    for (i = 0; i < n; i++) {
	for (j = 0; j < p; j++) {
	    data_pre[i * p + j] = data_matrix[i * p + j];
	}
    }

    /* row center data matrix if required*/
    if (*rowflag == 1) {
	rowcentre_JM (data_pre, n, p);
	if (*verbose == 1)
	    Rprintf ("Centering\n");
    }

    /* standardize columns of data matrix if required*/
    if (*colflag == 1) {
	colstandard_JM (data_pre, n, p);
	Rprintf("colstandard\n");
    }

    /* calculate pre-whitening matrix Kmat */
    if (*verbose == 1)	Rprintf ("Whitening\n");
    Kmat = Calloc (n * n, double);
    calc_K_JM(data_pre, nn, pp, Kmat);

    /* pre-whiten data and reduce dimension from size n to size e */

    for (i = 0; i < e; i++) {
	for (j = 0; j < n; j++) {
	    Kmat1[i * n + j] = Kmat[i * n + j];
	}
    }
    mmult_JM (Kmat1, e, n, data_pre, n, p, data1);

    /* calculate initial (orthogonal) unmixing matrix w */
    temp1 = Calloc (e * e, double);
    w_init = Calloc (e * e, double);
    for (i = 0; i < e; i++) {
	for (j = 0; j < e; j++) {
	    temp1[i * e + j] = w_matrix[i * e + j];
	}
    }
    orthog_mat_JM (temp1, e, w_init);

    /* Main ICA code */
    if (*defflag == 0) {
	if (*funflag == 1) {
	    if (*verbose == 1)
		Rprintf("Symmetric FastICA using logcosh approx. to neg-entropy function\n");

	    i = 1;
	    Symm_logcosh_JM (w_init, e, data1, e, p, *alpha, w_final, &tol);
	    if (*verbose == 1)
		Rprintf ("Iteration %d tol=%f\n", i, tol);
	    i = 2;

	    while ((tol > (*lim)) && (i < (*maxit))) {
		Symm_logcosh_JM (w_final, e, data1, e, p, *alpha, w_final, &tol);
		if (*verbose == 1)
		    Rprintf ("Iteration %d tol=%f\n", i, tol);
		i += 1;
	    }
	}

	if (*funflag == 2) {
	    if (*verbose == 1)
		Rprintf("Symmetric FastICA using exponential approx. to neg-entropy function\n");

	    i = 1;
	    Symm_exp_JM (w_init, e, data1, e, p, *alpha, w_final, &tol);
	    if (*verbose == 1) Rprintf ("Iteration %d tol=%f\n", i, tol);

	    i = 2;
	    while ((tol > (*lim)) && (i < (*maxit))) {
		Symm_exp_JM (w_final, e, data1, e, p, *alpha, w_final, &tol);
		if (*verbose == 1) Rprintf ("Iteration %d tol=%f\n", i, tol);
		i += 1;
	    }
	}
    }

    if (*defflag == 1) {
	temp_w1 = Calloc (e, double);
	temp_w2 = Calloc (e, double);

	if (*funflag == 1) {
	    if (*verbose == 1)
		Rprintf ("Deflation FastICA using logcosh approx. to neg-entropy function\n");

	    for (i = 0; i < e; i++) {
		k = 0;
		gramsch_JM (w_init, e, e, i + 1);
		rowstd_JM (w_init, e, e, i + 1);
		tol = 1;

		while ((tol > (*lim)) && (k < (*maxit))) {
		    for (j = 0; j < e; j++) {
			temp_w1[j] = w_init[i * e + j];
		    }
		    Def_logcosh_JM (temp_w1, e, data1, e, p, *alpha, temp_w2);
		    for (j = 0; j < e; j++) {
			w_init[i * e + j] = temp_w2[j];
		    }
		    gramsch_JM (w_init, e, e, i + 1);
		    rowstd_JM (w_init, e, e, i + 1);
		    tol = 0;
		    for (j = 0; j < e; j++) {
			tol += ((temp_w1[j]) * (w_init[i * e + j]));
		    }
		    tol = (fabs (fabs (tol) - 1));
		    k += 1;
		}

		if (*verbose == 1)
		    Rprintf ("Component %d needed %d iterations tol=%f\n",
			     i + 1, k, tol);

	    }
	}
	if (*funflag == 2) {

	    if (*verbose == 1)
		Rprintf ("Deflation FastICA using exponential approx. to neg-entropy function\n");

	    for (i = 0; i < e; i++) {
		k = 0;
		gramsch_JM (w_init, e, e, i + 1);
		rowstd_JM (w_init, e, e, i + 1);
		tol = 1;

		while ((tol > (*lim)) && (k < (*maxit))) {
		    for (j = 0; j < e; j++) {
			temp_w1[j] = w_init[i * e + j];
		    }
		    Def_exp_JM (temp_w1, e, data1, e, p, *alpha, temp_w2);
		    for (j = 0; j < e; j++) {
			w_init[i * e + j] = temp_w2[j];
		    }
		    gramsch_JM (w_init, e, e, i + 1);
		    rowstd_JM (w_init, e, e, i + 1);
		    tol = 0;
		    for (j = 0; j < e; j++) {
			tol += ((temp_w1[j]) * (w_init[i * e + j]));
		    }
		    tol = (fabs (fabs (tol) - 1));
		    k += 1;
		}

		if (*verbose == 1)
		    Rprintf ("Component %d needed %d iterations tol=%f\n",
			     i + 1, k, tol);

	    }
	}
	for (i = 0; i < e; i++) {
	    for (j = 0; j < e; j++) {
		w_final[i * e + j] = w_init[i * e + j];
	    }
	}
	Free (temp_w1);
	Free (temp_w2);
    }

    /* calculate mixing matrix ansa */
    calc_A_JM(w_final, Kmat1, data_pre, ee, nn, pp, ansa, ansx2);

    Free (data1);
    Free (Kmat);
    Free (temp1);
    Free (w_init);
}

#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[] = {
    {"icainc_JM", (DL_FUNC) &icainc_JM, 18},
    {NULL, NULL, 0}
};


void
R_init_fastICA(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
