/* simfock2.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer pair[999000]	/* was [2][499500] */, kdim, ndim;
} probdata_;

#define probdata_1 probdata_

/*     ================================================================= */
/*     File: simfock2.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     July20, 2005. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Simulation of Hartree-Fock problem (version 2) */
/*     ---------------------------------------------- */
/*     Description: */
/*     Given kdim and ndim satisfying ndim <= kdim, the problem is: */

/*                  Minimize 0.5 x^T A x + b^T x */

/*                  subject to */

/*                          M orthonornal, */

/*     where x is the columnwise representation of the matrix X = M M^T, */
/*     and M in R^{kdim \times ndim}. Therefore, the problem has */
/*     ndim * kdim variables (the elements of matrix M) and */
/*     ndim + ndim * ( ndim - 1 ) / 2 nonlinear equality constraints. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock2_inip__(n, x, l, u, m, lambda, rho, equatn, 
	linear, n_in__, k_in__)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
integer *n_in__, *k_in__;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal seed;
    static integer i__, j, k;
    extern doublereal drand_();

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine must set some problem data. For achieving this */
/*     objective YOU MUST MODIFY it according to your problem. See below */
/*     where your modifications must be inserted. */

/*     Parameters of the subroutine: */

/*     On Entry: */

/*     This subroutine has no input parameters. */

/*     On Return */

/*     n        integer, */
/*              number of variables, */

/*     x        double precision x(n), */
/*              initial point, */

/*     l        double precision l(n), */
/*              lower bounds on x, */

/*     u        double precision u(n), */
/*              upper bounds on x, */

/*     m        integer, */
/*              number of constraints (excluding the bounds), */

/*     lambda   double precision lambda(m), */
/*              initial estimation of the Lagrange multipliers, */

/*     rho      double precision rho(m), */
/*              initial penalty parameters. */

/*     equatn   logical equatn(m) */
/*              for each constraint j, set equatn(j) = .true. if it is an */
/*              equality constraint of the form c_j(x) = 0, and set */
/*              equatn(j) = .false. if it is an inequality constraint of */
/*              the form c_j(x) <= 0, */

/*     linear   logical linear(m) */
/*              for each constraint j, set linear(j) = .true. if it is a */
/*              linear constraint, and set linear(j) = .false. if it is a */
/*              nonlinear constraint. */
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     EXTERNAL FUNCTIONS */
/*     Set problem data */
/*     Dimension of the space (kdim) maximum 100 */
    /* Parameter adjustments */
    --linear;
    --equatn;
    --rho;
    --lambda;
    --u;
    --l;
    --x;

    /* Function Body */
    probdata_1.kdim = *k_in__;
/*     Rank (ndim) of the desired idempotent solution: (ndim must be less than or equal to kdim) */
    probdata_1.ndim = *n_in__;
/*     Set pairs */
    k = 0;
    i__1 = probdata_1.ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = probdata_1.ndim;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ++k;
	    probdata_1.pair[(k << 1) - 2] = i__;
	    probdata_1.pair[(k << 1) - 1] = j;
	}
    }
/*     Number of variables */
    *n = probdata_1.ndim * probdata_1.kdim;
/*     Initial point */
    seed = 1324.;
    i__1 = probdata_1.ndim;
    for (j = 1; j <= i__1; ++j) {
	i__2 = probdata_1.kdim;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[(j - 1) * probdata_1.kdim + i__] = drand_(&seed);
	}
    }
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e5;
	u[i__] = 1e5;
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = probdata_1.ndim + probdata_1.ndim * (probdata_1.ndim - 1) / 2;
/*     Lagrange multipliers approximation. Most users prefer to use the */
/*     null initial Lagrange multipliers estimates. However, if the */
/*     problem that you are solving is "slightly different" from a */
/*     previously solved problem of which you know the correct Lagrange */
/*     multipliers, we encourage you to set these multipliers as initial */
/*     estimates. Of course, in this case you are also encouraged to use */
/*     the solution of the previous problem as initial estimate of the */
/*     solution. Similarly, most users prefer to use rho = 10 as initial */
/*     penalty parameters. But in the case mentioned above (good */
/*     estimates of solution and Lagrange multipliers) larger values of */
/*     the penalty parameters (say, rho = 1000) may be more useful. More */
/*     warm-start procedures are being elaborated. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lambda[i__] = 0.;
    }
/*     Initial penalty parameters */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rho[i__] = 10.;
    }
/*     For each constraint i, set equatn(i) = .true. if it is an equality */
/*     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if */
/*     it is an inequality constraint of the form c_i(x) <= 0. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	equatn[i__] = TRUE_;
    }
/*     For each constraint i, set linear(i) = .true. if it is a linear */
/*     constraint, otherwise set linear(i) = .false. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	linear[i__] = FALSE_;
    }
} /* simfock2_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock2_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal y[10000000], gy[1000000];
    static integer ny;
    extern /* Subroutine */ int simfock2_gedede__();
    static doublereal tmp;

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine must compute the objective function. For achieving */
/*     this objective YOU MUST MODIFY it according to your problem. See */
/*     below where your modifications must be inserted. */

/*     Parameters of the subroutine: */

/*     On Entry: */

/*     n        integer, */
/*              number of variables, */

/*     x        double precision x(n), */
/*              current point, */

/*     On Return */

/*     f        double precision, */
/*              objective function value at x, */

/*     flag     integer, */
/*              You must set it to any number different of 0 (zero) if */
/*              some error ocurred during the evaluation of the objective */
/*              function. (For example, trying to compute the square root */
/*              of a negative number, dividing by zero or a very small */
/*              number, etc.) If everything was o.k. you must set it */
/*              equal to zero. */
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
/*     Y = M M^t */
/*     y is vector which saves Y (kdim \times kdim) in the columnwise */
/*     representation and x is the same thing for M (kdim \times ndim) */
    i__1 = probdata_1.kdim;
    for (j = 1; j <= i__1; ++j) {
	i__2 = probdata_1.kdim;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*             Y(i,j) = row i by row j of matrix M */
	    tmp = 0.;
	    i__3 = probdata_1.ndim;
	    for (k = 1; k <= i__3; ++k) {
		tmp += x[(k - 1) * probdata_1.kdim + i__] * x[(k - 1) * 
			probdata_1.kdim + j];
	    }
	    y[(j - 1) * probdata_1.kdim + i__ - 1] = tmp;
	}
    }
/*     Now the quadratic objective function in y */
    *f = 0.;
    i__1 = probdata_1.kdim - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	*f = *f - y[(i__ - 2) * probdata_1.kdim + i__ - 1] + y[(i__ - 1) * 
		probdata_1.kdim + i__ - 1] * 2. - y[i__ * probdata_1.kdim + 
		i__ - 1];
    }
    *f = *f + y[0] * 2. - y[probdata_1.kdim] - y[(probdata_1.kdim - 2) * 
	    probdata_1.kdim + probdata_1.kdim - 1] + y[(probdata_1.kdim - 1) *
	     probdata_1.kdim + probdata_1.kdim - 1] * 2.;
    *f *= 2.;
/* Computing 2nd power */
    i__1 = probdata_1.kdim;
    ny = i__1 * i__1;
    simfock2_gedede__(&ny, y, gy);
    i__1 = probdata_1.kdim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = probdata_1.kdim;
	for (j = 1; j <= i__2; ++j) {
	    *f += gy[(i__ - 1) * probdata_1.kdim + j - 1] * y[(i__ - 1) * 
		    probdata_1.kdim + j - 1];
	}
    }
} /* simfock2_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock2_evalg__(n, x, g, flag__)
integer *n;
doublereal *x, *g;
integer *flag__;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal y[1000000], gy[1000000];
    static integer ny;
    extern /* Subroutine */ int simfock2_gedede__();
    static doublereal tmp;

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine must compute the gradient vector of the objective */
/*     function. For achieving these objective YOU MUST MODIFY it in the */
/*     way specified below. However, if you decide to use numerical */
/*     derivatives (we dont encourage this option at all!) you dont need */
/*     to modify evalg. */

/*     Parameters of the subroutine: */

/*     On Entry: */

/*     n        integer, */
/*              number of variables, */

/*     x        double precision x(n), */
/*              current point, */

/*     On Return */

/*     g        double precision g(n), */
/*              gradient vector of the objective function evaluated at x, */

/*     flag     integer, */
/*              You must set it to any number different of 0 (zero) if */
/*              some error ocurred during the evaluation of any component */
/*              of the gradient vector. (For example, trying to compute */
/*              the square root of a negative number, dividing by zero or */
/*              a very small number, etc.) If everything was o.k. you */
/*              must set it equal to zero. */
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    *flag__ = 0;
/*     Y = M M^t */
/*     y is vector which saves Y (kdim \times kdim) in the columnwise */
/*     representation and x is the same thing for M (kdim \times ndim) */
    i__1 = probdata_1.kdim;
    for (j = 1; j <= i__1; ++j) {
	i__2 = probdata_1.kdim;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*             Y(i,j) = row i by row j of matrix M */
	    tmp = 0.;
	    i__3 = probdata_1.ndim;
	    for (k = 1; k <= i__3; ++k) {
		tmp += x[(k - 1) * probdata_1.kdim + i__] * x[(k - 1) * 
			probdata_1.kdim + j];
	    }
	    y[(j - 1) * probdata_1.kdim + i__ - 1] = tmp;
	}
    }
/*     Now the derivatives of the quadratic objective function */
/*     with respect to y */
/* Computing 2nd power */
    i__1 = probdata_1.kdim;
    ny = i__1 * i__1;
    simfock2_gedede__(&ny, y, gy);
    i__1 = probdata_1.kdim - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	gy[(i__ - 2) * probdata_1.kdim + i__ - 1] += -1.;
	gy[(i__ - 1) * probdata_1.kdim + i__ - 1] += 2.;
	gy[i__ * probdata_1.kdim + i__ - 1] += -1.;
    }
    gy[0] += 2.;
    gy[(probdata_1.kdim - 2) * probdata_1.kdim + probdata_1.kdim - 1] += -1.;
    gy[(probdata_1.kdim - 1) * probdata_1.kdim + probdata_1.kdim - 1] += 2.;
    gy[probdata_1.kdim] += -1.;
    i__1 = ny;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gy[i__ - 1] *= 2.;
    }
/*     Finally, the chain rule to compute the derivatives with respect */
/*     to the elements of matrix M */
    i__1 = probdata_1.ndim;
    for (j = 1; j <= i__1; ++j) {
	i__2 = probdata_1.kdim;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    g[(j - 1) * probdata_1.kdim + i__] = 0.;
	}
    }
    i__1 = probdata_1.kdim;
    for (j = 1; j <= i__1; ++j) {
	i__2 = probdata_1.kdim;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    tmp = gy[(j - 1) * probdata_1.kdim + i__ - 1];
	    i__3 = probdata_1.ndim;
	    for (k = 1; k <= i__3; ++k) {
		g[(k - 1) * probdata_1.kdim + i__] += tmp * x[(k - 1) * 
			probdata_1.kdim + j];
		g[(k - 1) * probdata_1.kdim + j] += tmp * x[(k - 1) * 
			probdata_1.kdim + i__];
	    }
	}
    }
} /* simfock2_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock2_evalh__(n, x, hlin, hcol, hval, nnzh, flag__)
integer *n;
doublereal *x;
integer *hlin, *hcol;
doublereal *hval;
integer *nnzh, *flag__;
{
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine might compute the Hessian matrix of the objective */
/*     function. For achieving this objective YOU MAY MODIFY it according */
/*     to your problem. To modify this subroutine IS NOT MANDATORY. See */
/*     below where your modifications must be inserted. */

/*     Parameters of the subroutine: */

/*     On Entry: */

/*     n        integer, */
/*              number of variables, */

/*     x        double precision x(n), */
/*              current point, */

/*     On Return */

/*     nnzh     integer, */
/*              number of perhaps-non-null elements of the computed */
/*              Hessian, */

/*     hlin     integer hlin(nnzh), */
/*              see below, */

/*     hcol     integer hcol(nnzh), */
/*              see below, */

/*     hval     double precision hval(nnzh), */
/*              the non-null value of the (hlin(k),hcol(k)) position */
/*              of the Hessian matrix of the objective function must */
/*              be saved at hval(k). Just the lower triangular part of */
/*              Hessian matrix must be computed, */

/*     flag     integer, */
/*              You must set it to any number different of 0 (zero) if */
/*              some error ocurred during the evaluation of the Hessian */
/*              matrix of the objective funtion. (For example, trying */
/*              to compute the square root of a negative number, */
/*              dividing by zero or a very small number, etc.) If */
/*              everything was o.k. you must set it equal to zero. */
    /* Parameter adjustments */
    --x;
    --hlin;
    --hcol;
    --hval;

    /* Function Body */
    *flag__ = -1;
} /* simfock2_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock2_evalc__(n, x, ind, c__, flag__)
integer *n;
doublereal *x;
integer *ind;
doublereal *c__;
integer *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, i1, i2;

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine computes the ind-th constraint. */

/*     Parameters of the subroutine: */

/*     On Entry: */

/*     n        integer, */
/*              number of variables, */

/*     x        double precision x(n), */
/*              current point, */

/*     ind      integer, */
/*              index of the constraint to be computed, */

/*     On Return */

/*     c        double precision, */
/*              ind-th constraint evaluated at x, */

/*     flag     integer */
/*              You must set it to any number different of 0 (zero) if */
/*              some error ocurred during the evaluation of the */
/*              constraint. (For example, trying to compute the square */
/*              root of a negative number, dividing by zero or a very */
/*              small number, etc.) If everything was o.k. you must set */
/*              it equal to zero. */
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    if (*ind <= probdata_1.ndim) {
/*         Euclidian norm of column ind of matrix M equal to 1 */
	*c__ = -1.;
	i__1 = probdata_1.kdim;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = x[(*ind - 1) * probdata_1.kdim + i__];
	    *c__ += d__1 * d__1;
	}
	return 0;
    }
    if (*ind <= probdata_1.ndim + probdata_1.ndim * (probdata_1.ndim - 1) / 2)
	     {
/*         Given the number of the constraint, identify the indices */
/*         i1 and i2 of both columns of matrix M */
	i__ = *ind - probdata_1.ndim;
	i1 = probdata_1.pair[(i__ << 1) - 2];
	i2 = probdata_1.pair[(i__ << 1) - 1];
/*         Orthogonality between columns i1 and i2 of matrix M */
	*c__ = 0.;
	i__1 = probdata_1.kdim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    *c__ += x[(i1 - 1) * probdata_1.kdim + i__] * x[(i2 - 1) * 
		    probdata_1.kdim + i__];
	}
	return 0;
    }
} /* simfock2_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock2_evaljac__(n, x, ind, indjac, valjac, nnzjac, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *indjac;
doublereal *valjac;
integer *nnzjac, *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, i1, i2;

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine computes the gradient of the ind-th constraint. */

/*     Parameters of the subroutine: */

/*     On Entry: */

/*     n        integer, */
/*              number of variables, */

/*     x        double precision x(n), */
/*              current point, */

/*     ind      integer, */
/*              index of the constraint whose gradient will be computed, */

/*     On Return */

/*     nnzjac   integer, */
/*              number of perhaps-non-null elements of the computed */
/*              gradient, */

/*     indjac   integer indjac(nnzjac), */
/*              see below, */

/*     valjac   double precision valjac(nnzjac), */
/*              the non-null value of the partial derivative of the */
/*              ind-th constraint with respect to the indjac(k)-th */
/*              variable must be saved at valjac(k). */

/*     flag     integer */
/*              You must set it to any number different of 0 (zero) if */
/*              some error ocurred during the evaluation of the */
/*              constraint. (For example, trying to compute the square */
/*              root of a negative number, dividing by zero or a very */
/*              small number, etc.) If everything was o.k. you must set */
/*              it equal to zero. */
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --valjac;
    --indjac;
    --x;

    /* Function Body */
    *flag__ = 0;
    if (*ind <= probdata_1.ndim) {
/*         Euclidian norm of column ind of matrix M equal to 1 */
	*nnzjac = probdata_1.kdim;
	i__1 = probdata_1.kdim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    indjac[i__] = (*ind - 1) * probdata_1.kdim + i__;
	    valjac[i__] = x[(*ind - 1) * probdata_1.kdim + i__] * 2.;
	}
	return 0;
    }
    if (*ind <= probdata_1.ndim + probdata_1.ndim * (probdata_1.ndim - 1) / 2)
	     {
/*         Given the number of the constraint, identify the indices */
/*         i1 and i2 of both columns of matrix M */
	i__ = *ind - probdata_1.ndim;
	i1 = probdata_1.pair[(i__ << 1) - 2];
	i2 = probdata_1.pair[(i__ << 1) - 1];
/*         Orthogonality between columns i1 and i2 of matrix M */
	*nnzjac = probdata_1.kdim << 1;
	i__1 = probdata_1.kdim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    indjac[i__] = (i1 - 1) * probdata_1.kdim + i__;
	    valjac[i__] = x[(i2 - 1) * probdata_1.kdim + i__];
	    indjac[probdata_1.kdim + i__] = (i2 - 1) * probdata_1.kdim + i__;
	    valjac[probdata_1.kdim + i__] = x[(i1 - 1) * probdata_1.kdim + 
		    i__];
	}
	return 0;
    }
} /* simfock2_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock2_gedede__(ny, y, gy)
integer *ny;
doublereal *y, *gy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, jj;

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
    /* Parameter adjustments */
    --gy;
    --y;

    /* Function Body */
    i__1 = *ny;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gy[i__] = y[i__];
	for (jj = 1; jj <= 10; ++jj) {
	    j = i__ + jj;
	    if (j <= *ny) {
		gy[i__] += 1. / (doublereal) (i__ + j) * y[j];
	    }
	    j = i__ - jj;
	    if (j >= 1) {
		gy[i__] += 1. / (doublereal) (i__ + j) * y[j];
	    }
	}
    }
} /* simfock2_gedede__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock2_evalhc__(n, x, ind, hclin, hccol, hcval, nnzhc, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *hclin, *hccol;
doublereal *hcval;
integer *nnzhc, *flag__;
{
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine might compute the Hessian matrix of the ind-th */
/*     constraint. For achieving this objective YOU MAY MODIFY it */
/*     according to your problem. To modify this subroutine IS NOT */
/*     MANDATORY. See below where your modifications must be inserted. */

/*     Parameters of the subroutine: */

/*     On Entry: */

/*     n        integer, */
/*              number of variables, */

/*     x        double precision x(n), */
/*              current point, */

/*     ind      integer, */
/*              index of the constraint whose Hessian will be computed, */

/*     On Return */

/*     nnzhc    integer, */
/*              number of perhaps-non-null elements of the computed */
/*              Hessian, */

/*     hclin    integer hclin(nnzhc), */
/*              see below, */

/*     hccol    integer hccol(nnzhc), */
/*              see below, */

/*     hcval    double precision hcval(nnzhc), */
/*              the non-null value of the (hclin(k),hccol(k)) position */
/*              of the Hessian matrix of the ind-th constraint must */
/*              be saved at hcval(k). Just the lower triangular part of */
/*              Hessian matrix must be computed, */

/*     flag     integer, */
/*              You must set it to any number different of 0 (zero) if */
/*              some error ocurred during the evaluation of the Hessian */
/*              matrix of the ind-th constraint. (For example, trying */
/*              to compute the square root of a negative number, */
/*              dividing by zero or a very small number, etc.) If */
/*              everything was o.k. you must set it equal to zero. */
    /* Parameter adjustments */
    --x;
    --hclin;
    --hccol;
    --hcval;

    /* Function Body */
    *flag__ = -1;
} /* simfock2_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock2_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
integer *n;
doublereal *x;
integer *m;
doublereal *lambda, *p, *hp;
logical *goth;
integer *flag__;
{
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine might compute the product of the Hessian of the */
/*     Lagrangian times vector p (just the Hessian of the objective */
/*     function in the unconstrained or bound-constrained case). */

/*     Parameters of the subroutine: */

/*     On Entry: */

/*     n        integer, */
/*              number of variables, */

/*     x        double precision x(n), */
/*              current point, */

/*     m        integer, */
/*              number of constraints, */

/*     lambda   double precision lambda(m), */
/*              vector of Lagrange multipliers, */

/*     p        double precision p(n), */
/*              vector of the matrix-vector product, */

/*     goth     logical, */
/*              can be used to indicate if the Hessian matrices were */
/*              computed at the current point. It is set to .false. */
/*              by the optimization method every time the current */
/*              point is modified. Sugestion: if its value is .false. */
/*              then compute the Hessians, save them in a common */
/*              structure and set goth to .true.. Otherwise, just use */
/*              the Hessians saved in the common block structure, */

/*     On Return */

/*     hp       double precision hp(n), */
/*              Hessian-vector product, */

/*     goth     logical, */
/*              see above, */

/*     flag     integer, */
/*              You must set it to any number different of 0 (zero) if */
/*              some error ocurred during the evaluation of the */
/*              Hessian-vector product. (For example, trying to compute */
/*              the square root of a negative number, dividing by zero */
/*              or a very small number, etc.) If everything was o.k. you */
/*              must set it equal to zero. */
    /* Parameter adjustments */
    --hp;
    --p;
    --x;
    --lambda;

    /* Function Body */
    *flag__ = -1;
} /* simfock2_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock2_endp__(n, x, l, u, m, lambda, rho, equatn, 
	linear)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
{
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine can be used to do some extra job after the solver */
/*     has found the solution,like some extra statistics, or to save the */
/*     solution in some special format or to draw some graphical */
/*     representation of the solution. If the information given by the */
/*     solver is enough for you then leave the body of this subroutine */
/*     empty. */

/*     Parameters of the subroutine: */

/*     The paraemters of this subroutine are the same parameters of */
/*     subroutine inip. But in this subroutine there are not output */
/*     parameter. All the parameters are input parameters. */
    /* Parameter adjustments */
    --u;
    --l;
    --x;
    --linear;
    --equatn;
    --rho;
    --lambda;

    /* Function Body */
} /* simfock2_endp__ */

