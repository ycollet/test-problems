/* bratu3Db.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal phi[941192]	/* was [98][98][98] */;
    integer vind[3000]	/* was [1000][3] */;
    doublereal vval[1000];
    integer cind[2823576]	/* was [3][941192] */, np, t;
    doublereal h__, theta;
} probdata_;

#define probdata_1 probdata_

/* Table of constant values */

static doublereal c_b2 = 4.5;

/*     ================================================================= */
/*     File: bratu3Db.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     November 30, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Discretized 3D Bratu-based problem */
/*     ---------------------------------- */
/*     Description: */
/*     The discretized three-dimensional Bratu-based optimization */
/*     problem is: */

/*     Minimize \sum_{(i,j,k) \in S} ( u(i,j,k) - u_*(i,j,k) )^2 */

/*     subject to */

/*     \phi_\theta(u,i,j,k) = \phi_\theta(u_*,i,j,k), i,j,k=2,...,np-1, */

/*     where u_* was choosed as */

/*     u_*(i,j,k) = 10 q(i) q(j) q(k) */
/*                     (1-q(i)) (1-q(j)) (1-q(k)) e^{q(k)^{4.5}} */

/*     with */

/*     q(\ell) = \frac{ np - \ell }{ np - 1 } for i,j,k=1,...,np */

/*     and */

/*     \phi_\theta(v,i,j,k) = - \Delta v(i,j,k) + \theta e^{v(i,j,k)}, */

/*     \Delta v(i,j,k) =  \frac{v(i \pm 1,j,k) + */
/*                              v(i,j \pm 1,k) + */
/*                              v(i,j,k \pm 1) - 6 v(i,j,k)} {h^2}, */

/*     for i,j,k=2,...,np-1. */

/*     The number of variables is n=np^3 and the number of (equality) */
/*     constraints is m=(np-2)^3. */

/*     We set \theta=-100, h=1/(np-1), |S|=7 and the 3-uples of */
/*     indices in S are randomly selected in [1,np]^3. */

/*     The initial point is randomly generated in [0,1]^n. */

/*     Several problems can be considered varying the value of np, for */
/*     example, np = 5, 6, ..., 20. */
/*     A much more friendly description of the problem can be found in: */

/*     R. Andreani, E.G. Birgin, J.M. Martinez and M.L. Schuverdt, On */
/*     Augmented Lagrangian Methods with general lower-level constraints, */
/*     Technical Report MCDO-051015 (see www.ime.usp.br/~egbirgin/), */
/*     Department of Applied Mathematics, UNICAMP, Brazil, 2005. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int bratu3db_inip__(n, x, l, u, m, lambda, rho, equatn, 
	linear, nppar, seed)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
integer *nppar;
doublereal *seed;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(), exp();

    /* Local variables */
    static doublereal uopt[1000000]	/* was [100][100][100] */, a, b, c__;
    static integer i__, j, k;
    extern integer bratu3db_xind__();
    extern doublereal drand_();
    static integer ind;

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
/*     COMON ARRAYS */
/*     LOCAL SCALARS */
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
/*     EXTERNAL FUNCTIONS */
/*     READ PROBLEM VARIABLE DATA */
    /* Parameter adjustments */
    --linear;
    --equatn;
    --rho;
    --lambda;
    --u;
    --l;
    --x;

    /* Function Body */
    probdata_1.theta = -100.;
    probdata_1.t = 7;
/*     write(*,*) 'Number of points (np): ' */
/*     read(*,*) np */
/*     write(*,*) 'Seed for the initial point random generation: ' */
/*     read(*,*) seed */
/*     read(*,*) np,theta,t,seed */
/*     seed = 123456.0d0 * seed */
/*     DEFINE SOME PROBLEM PARAMETERS */
/*     Discretization step */
    probdata_1.np = *nppar;
    probdata_1.h__ = 1. / (probdata_1.np - 1);
/*     "Known solution" used to generate the problem data */
    i__1 = probdata_1.np;
    for (k = 1; k <= i__1; ++k) {
	i__2 = probdata_1.np;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = probdata_1.np;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		a = (real) (probdata_1.np - k) / (probdata_1.np - 1);
		b = (real) (probdata_1.np - j) / (probdata_1.np - 1);
		c__ = (real) (probdata_1.np - i__) / (probdata_1.np - 1);
		uopt[i__ + (j + k * 100) * 100 - 10101] = a * 10. * b * c__ * 
			(1. - a) * (1. - b) * (1. - c__) * exp(pow_dd(&a, &
			c_b2));
	    }
	}
    }
/*     Right hand side of the system */
    i__1 = probdata_1.np - 1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = probdata_1.np - 1;
	for (j = 2; j <= i__2; ++j) {
	    i__3 = probdata_1.np - 1;
	    for (i__ = 2; i__ <= i__3; ++i__) {
/* Computing 2nd power */
		d__1 = probdata_1.h__;
		probdata_1.phi[i__ + (j + k * 98) * 98 - 19406] = (uopt[i__ + 
			(j + (k + 1) * 100) * 100 - 10101] + uopt[i__ + (j + (
			k - 1) * 100) * 100 - 10101] + uopt[i__ + (j + 1 + k *
			 100) * 100 - 10101] + uopt[i__ + (j - 1 + k * 100) * 
			100 - 10101] + uopt[i__ + 1 + (j + k * 100) * 100 - 
			10101] + uopt[i__ - 1 + (j + k * 100) * 100 - 10101] 
			- uopt[i__ + (j + k * 100) * 100 - 10101] * 6.) / (
			d__1 * d__1) - probdata_1.theta * exp(uopt[i__ + (j + 
			k * 100) * 100 - 10101]);
	    }
	}
    }
/*     Points of the grid at which the solution is known */
    i__1 = probdata_1.t;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    probdata_1.vind[i__ + j * 1000 - 1001] = (integer) (drand_(seed) *
		     probdata_1.np) + 1;
	}
	probdata_1.vval[i__ - 1] = uopt[probdata_1.vind[i__ - 1] + (
		probdata_1.vind[i__ + 999] + probdata_1.vind[i__ + 1999] * 
		100) * 100 - 10101];
    }
/*     Indices i, j and k related to ind-th constraint */
    ind = 0;
    i__1 = probdata_1.np - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = probdata_1.np - 1;
	for (j = 2; j <= i__2; ++j) {
	    i__3 = probdata_1.np - 1;
	    for (k = 2; k <= i__3; ++k) {
		++ind;
		probdata_1.cind[ind * 3 - 3] = i__;
		probdata_1.cind[ind * 3 - 2] = j;
		probdata_1.cind[ind * 3 - 1] = k;
	    }
	}
    }
/*     Number of variables (elements of a lower-triangular ndim x ndim matrix) */
/* Computing 3rd power */
    i__1 = probdata_1.np;
    *n = i__1 * (i__1 * i__1);
/*     Lower and upper bounds */
    i__1 = probdata_1.np;
    for (k = 1; k <= i__1; ++k) {
	i__2 = probdata_1.np;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = probdata_1.np;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ind = bratu3db_xind__(&i__, &j, &k);
		l[ind] = -1e20;
		u[ind] = 1e20;
	    }
	}
    }
/*     Initial guess */
    i__1 = probdata_1.np;
    for (k = 1; k <= i__1; ++k) {
	i__2 = probdata_1.np;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = probdata_1.np;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ind = bratu3db_xind__(&i__, &j, &k);
		x[ind] = drand_(seed);
	    }
	}
    }
/*     Number of constraints (equalities plus inequalities) */
/* Computing 3rd power */
    i__1 = probdata_1.np - 2;
    *m = i__1 * (i__1 * i__1);
/*     Lagrange multipliers approximation. */
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
} /* bratu3db_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int bratu3db_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    extern integer bratu3db_xind__();
    static integer ind;

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
/*     COMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     EXTERNAL FUNCTIONS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    *f = 0.;
    i__1 = probdata_1.t;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ind = bratu3db_xind__(&probdata_1.vind[i__ - 1], &probdata_1.vind[i__ 
		+ 999], &probdata_1.vind[i__ + 1999]);
/* Computing 2nd power */
	d__1 = x[ind] - probdata_1.vval[i__ - 1];
	*f += d__1 * d__1;
    }
} /* bratu3db_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int bratu3db_evalg__(n, x, g, flag__)
integer *n;
doublereal *x, *g;
integer *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern integer bratu3db_xind__();
    static integer ind;

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
/*     COMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     EXTERNAL FUNCTIONS */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    *flag__ = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = 0.;
    }
    i__1 = probdata_1.t;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ind = bratu3db_xind__(&probdata_1.vind[i__ - 1], &probdata_1.vind[i__ 
		+ 999], &probdata_1.vind[i__ + 1999]);
	g[ind] += (x[ind] - probdata_1.vval[i__ - 1]) * 2.;
    }
} /* bratu3db_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int bratu3db_evalh__(n, x, hlin, hcol, hval, hnnz, flag__)
integer *n;
doublereal *x;
integer *hlin, *hcol;
doublereal *hval;
integer *hnnz, *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern integer bratu3db_xind__();
    static integer ind;

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

/*     hnnz     integer, */
/*              number of perhaps-non-null elements of the computed */
/*              Hessian, */

/*     hlin     integer hlin(hnnz), */
/*              see below, */

/*     hcol     integer hcol(hnnz), */
/*              see below, */

/*     hval     double precision hval(hnnz), */
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
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     EXTERNAL FUNCTIONS */
    /* Parameter adjustments */
    --x;
    --hlin;
    --hcol;
    --hval;

    /* Function Body */
    *flag__ = 0;
    *hnnz = probdata_1.t;
    i__1 = probdata_1.t;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ind = bratu3db_xind__(&probdata_1.vind[i__ - 1], &probdata_1.vind[i__ 
		+ 999], &probdata_1.vind[i__ + 1999]);
	hlin[i__] = ind;
	hcol[i__] = ind;
	hval[i__] = 2.;
    }
} /* bratu3db_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int bratu3db_evalc__(n, x, ind, c__, flag__)
integer *n;
doublereal *x;
integer *ind;
doublereal *c__;
integer *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static integer i__, j, k;
    extern integer bratu3db_xind__();
    static integer i0, i1, i2, i3, i4, i5, i6;

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
/*     COMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     EXTERNAL FUNCTIONS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    i__ = probdata_1.cind[*ind * 3 - 3];
    j = probdata_1.cind[*ind * 3 - 2];
    k = probdata_1.cind[*ind * 3 - 1];
    i0 = bratu3db_xind__(&i__, &j, &k);
    i__1 = k + 1;
    i1 = bratu3db_xind__(&i__, &j, &i__1);
    i__1 = k - 1;
    i2 = bratu3db_xind__(&i__, &j, &i__1);
    i__1 = j + 1;
    i3 = bratu3db_xind__(&i__, &i__1, &k);
    i__1 = j - 1;
    i4 = bratu3db_xind__(&i__, &i__1, &k);
    i__1 = i__ + 1;
    i5 = bratu3db_xind__(&i__1, &j, &k);
    i__1 = i__ - 1;
    i6 = bratu3db_xind__(&i__1, &j, &k);
/* Computing 2nd power */
    d__1 = probdata_1.h__;
    *c__ = (x[i1] + x[i2] + x[i3] + x[i4] + x[i5] + x[i6] - x[i0] * 6.) / (
	    d__1 * d__1) - probdata_1.theta * exp(x[i0]) - probdata_1.phi[i__ 
	    + (j + k * 98) * 98 - 19406];
} /* bratu3db_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int bratu3db_evaljac__(n, x, ind, jcvar, jcval, jcnnz, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *jcvar;
doublereal *jcval;
integer *jcnnz, *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static integer i__, j, k;
    extern integer bratu3db_xind__();
    static integer i0, i1, i2, i3, i4, i5, i6;

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

/*     jcnnz    integer, */
/*              number of perhaps-non-null elements of the computed */
/*              gradient, */

/*     jcvar    integer jcvar(jcnnz), */
/*              see below, */

/*     jcval    double precision jcval(jcnnz), */
/*              the non-null value of the partial derivative of the */
/*              ind-th constraint with respect to the jcvar(k)-th */
/*              variable must be saved at jcval(k). */

/*     flag     integer */
/*              You must set it to any number different of 0 (zero) if */
/*              some error ocurred during the evaluation of the */
/*              constraint. (For example, trying to compute the square */
/*              root of a negative number, dividing by zero or a very */
/*              small number, etc.) If everything was o.k. you must set */
/*              it equal to zero. */
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     EXTERNAL FUNCTIONS */
    /* Parameter adjustments */
    --jcval;
    --jcvar;
    --x;

    /* Function Body */
    *flag__ = 0;
    i__ = probdata_1.cind[*ind * 3 - 3];
    j = probdata_1.cind[*ind * 3 - 2];
    k = probdata_1.cind[*ind * 3 - 1];
    i0 = bratu3db_xind__(&i__, &j, &k);
    i__1 = k + 1;
    i1 = bratu3db_xind__(&i__, &j, &i__1);
    i__1 = k - 1;
    i2 = bratu3db_xind__(&i__, &j, &i__1);
    i__1 = j + 1;
    i3 = bratu3db_xind__(&i__, &i__1, &k);
    i__1 = j - 1;
    i4 = bratu3db_xind__(&i__, &i__1, &k);
    i__1 = i__ + 1;
    i5 = bratu3db_xind__(&i__1, &j, &k);
    i__1 = i__ - 1;
    i6 = bratu3db_xind__(&i__1, &j, &k);
    *jcnnz = 7;
    jcvar[1] = i1;
/* Computing 2nd power */
    d__1 = probdata_1.h__;
    jcval[1] = 1. / (d__1 * d__1);
    jcvar[2] = i2;
/* Computing 2nd power */
    d__1 = probdata_1.h__;
    jcval[2] = 1. / (d__1 * d__1);
    jcvar[3] = i3;
/* Computing 2nd power */
    d__1 = probdata_1.h__;
    jcval[3] = 1. / (d__1 * d__1);
    jcvar[4] = i4;
/* Computing 2nd power */
    d__1 = probdata_1.h__;
    jcval[4] = 1. / (d__1 * d__1);
    jcvar[5] = i5;
/* Computing 2nd power */
    d__1 = probdata_1.h__;
    jcval[5] = 1. / (d__1 * d__1);
    jcvar[6] = i6;
/* Computing 2nd power */
    d__1 = probdata_1.h__;
    jcval[6] = 1. / (d__1 * d__1);
    jcvar[7] = i0;
/* Computing 2nd power */
    d__1 = probdata_1.h__;
    jcval[7] = -6. / (d__1 * d__1) - probdata_1.theta * exp(x[i0]);
} /* bratu3db_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int bratu3db_evalhc__(n, x, ind, hclin, hccol, hcval, hcnnz, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *hclin, *hccol;
doublereal *hcval;
integer *hcnnz, *flag__;
{
    /* Builtin functions */
    double exp();

    /* Local variables */
    static integer i__, j, k;
    extern integer bratu3db_xind__();
    static integer i0;

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

/*     hcnnz    integer, */
/*              number of perhaps-non-null elements of the computed */
/*              Hessian, */

/*     hclin    integer hclin(hcnnz), */
/*              see below, */

/*     hccol    integer hccol(hcnnz), */
/*              see below, */

/*     hcval    double precision hcval(hcnnz), */
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
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     EXTERNAL FUNCTIONS */
    /* Parameter adjustments */
    --x;
    --hclin;
    --hccol;
    --hcval;

    /* Function Body */
    *flag__ = 0;
    i__ = probdata_1.cind[*ind * 3 - 3];
    j = probdata_1.cind[*ind * 3 - 2];
    k = probdata_1.cind[*ind * 3 - 1];
    i0 = bratu3db_xind__(&i__, &j, &k);
    *hcnnz = 1;
    hclin[1] = i0;
    hccol[1] = i0;
    hcval[1] = -probdata_1.theta * exp(x[i0]);
} /* bratu3db_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int bratu3db_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
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
} /* bratu3db_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int bratu3db_endp__(n, x, l, u, m, lambda, rho, equatn, 
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
} /* bratu3db_endp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
integer bratu3db_xind__(i__, j, k)
integer *i__, *j, *k;
{
    /* System generated locals */
    integer ret_val, i__1;

/*     SCALAR ARGUMENTS */
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMON ARRAYS */
/*     COMMON BLOCKS */
/* Computing 2nd power */
    i__1 = probdata_1.np;
    ret_val = i__1 * i__1 * (*i__ - 1) + probdata_1.np * (*j - 1) + *k;
    return ret_val;
} /* bratu3db_xind__ */

