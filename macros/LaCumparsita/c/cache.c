/* cache.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal ll[500000], ff[500000], mu[500000], llambda[500000], w[500000],
	     alphas, alphab, betas, betab, bwin, bwout, t;
} probdata_;

#define probdata_1 probdata_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;

/*     ================================================================= */
/*     File: cache.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     August 8th, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Cache problem */
/*     ------------- */
/*     Mathematical programming model for fixing content location cache */
/*     expiration dates at aggregation nodes in a content network, in */
/*     order to maximize the information discovered by the nodes of the */
/*     networks, while taking into account bandwidth constraints at the */
/*     aggregation nodes. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int cache_inip__(n, x, l, u, m, lambda, rho, equatn, linear)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsle(), do_lio(), e_rsle(), f_clos();
    double exp();

    /* Local variables */
    static integer i__, k;
    static doublereal tmp;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 10, 0, 0, 0 };
    static cilist io___3 = { 0, 10, 0, 0, 0 };
    static cilist io___5 = { 0, 10, 0, 0, 0 };
    static cilist io___6 = { 0, 10, 0, 0, 0 };
    static cilist io___7 = { 0, 10, 0, 0, 0 };
    static cilist io___8 = { 0, 10, 0, 0, 0 };


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
/*     COMMON BLOCKS */
/*     LOCAL SCALARS */
/*     Set problem data */
    /* Parameter adjustments */
    --linear;
    --equatn;
    --rho;
    --lambda;
    --u;
    --l;
    --x;

    /* Function Body */
    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = 13;
    o__1.ofnm = "cache-128.dat";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsle(&io___1);
    do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___3);
    do_lio(&c__5, &c__1, (char *)&probdata_1.alphas, (ftnlen)sizeof(
	    doublereal));
    do_lio(&c__5, &c__1, (char *)&probdata_1.alphab, (ftnlen)sizeof(
	    doublereal));
    do_lio(&c__5, &c__1, (char *)&probdata_1.betas, (ftnlen)sizeof(doublereal)
	    );
    do_lio(&c__5, &c__1, (char *)&probdata_1.betab, (ftnlen)sizeof(doublereal)
	    );
    do_lio(&c__5, &c__1, (char *)&probdata_1.bwin, (ftnlen)sizeof(doublereal))
	    ;
    do_lio(&c__5, &c__1, (char *)&probdata_1.bwout, (ftnlen)sizeof(doublereal)
	    );
    e_rsle();
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_rsle(&io___5);
	do_lio(&c__5, &c__1, (char *)&probdata_1.ff[i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsle();
    }
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_rsle(&io___6);
	do_lio(&c__5, &c__1, (char *)&probdata_1.llambda[i__ - 1], (ftnlen)
		sizeof(doublereal));
	e_rsle();
    }
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_rsle(&io___7);
	do_lio(&c__5, &c__1, (char *)&probdata_1.mu[i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsle();
    }
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_rsle(&io___8);
	do_lio(&c__5, &c__1, (char *)&probdata_1.ll[i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
    probdata_1.t = 0.;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tmp = probdata_1.ll[i__ - 1] * probdata_1.llambda[i__ - 1] / 
		probdata_1.mu[i__ - 1];
	probdata_1.t += probdata_1.ff[i__ - 1] * tmp;
	probdata_1.w[i__ - 1] = tmp / probdata_1.mu[i__ - 1];
    }
/*     Number of variables */
    *n = k * 3;
/*     Initial point */
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[(k << 1) + i__] = .1;
	x[i__] = (1. - exp(-probdata_1.ff[i__ - 1] * x[(k << 1) + i__])) / x[(
		k << 1) + i__];
	x[k + i__] = (1. - exp(-probdata_1.mu[i__ - 1] * x[(k << 1) + i__])) /
		 x[(k << 1) + i__];
    }
/*     Lower and upper bounds */
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = 0.;
	u[i__] = probdata_1.ff[i__ - 1];
	l[k + i__] = 0.;
	u[k + i__] = probdata_1.mu[i__ - 1];
	l[(k << 1) + i__] = 0.;
	u[(k << 1) + i__] = 1e20;
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = (k << 1) + 2;
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
    i__1 = k << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	equatn[i__] = TRUE_;
    }
    equatn[(k << 1) + 1] = FALSE_;
    equatn[(k << 1) + 2] = FALSE_;
/*     For each constraint i, set linear(i) = .true. if it is a linear */
/*     constraint, otherwise set linear(i) = .false. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	linear[i__] = FALSE_;
    }
} /* cache_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int cache_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;

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
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    k = *n / 3;
    *f = 0.;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*f -= probdata_1.w[i__ - 1] * (probdata_1.mu[i__ - 1] * x[i__] + 
		probdata_1.ff[i__ - 1] * x[k + i__] - x[i__] * x[k + i__]);
    }
    *f /= probdata_1.t;
} /* cache_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int cache_evalg__(n, x, g, flag__)
integer *n;
doublereal *x, *g;
integer *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;

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
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    *flag__ = 0;
    k = *n / 3;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = -probdata_1.w[i__ - 1] * (probdata_1.mu[i__ - 1] - x[k + i__]
		) / probdata_1.t;
	g[k + i__] = -probdata_1.w[i__ - 1] * (probdata_1.ff[i__ - 1] - x[i__]
		) / probdata_1.t;
	g[(k << 1) + i__] = 0.;
    }
} /* cache_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int cache_evalh__(n, x, hlin, hcol, hval, nnzh, flag__)
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
} /* cache_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int cache_evalc__(n, x, ind, c__, flag__)
integer *n;
doublereal *x;
integer *ind;
doublereal *c__;
integer *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static integer i__, k;

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine must compute the ind-th constraint of your */
/*     problem. For achieving this objective YOU MUST MOFIFY it */
/*     according to your problem. See below the places where your */
/*     modifications must be inserted. */

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
    k = *n / 3;
    if (*ind <= k) {
	i__ = *ind;
	*c__ = 1. - exp(-probdata_1.ff[i__ - 1] * x[(k << 1) + i__]) - x[i__] 
		* x[(k << 1) + i__];
    } else if (*ind <= k << 1) {
	i__ = *ind - k;
	*c__ = 1. - exp(-probdata_1.mu[i__ - 1] * x[(k << 1) + i__]) - x[k + 
		i__] * x[(k << 1) + i__];
    } else if (*ind == (k << 1) + 1) {
	*c__ = 0.;
	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    *c__ = *c__ + probdata_1.betas * probdata_1.ll[i__ - 1] * 
		    probdata_1.ff[i__ - 1] + probdata_1.alphab * 
		    probdata_1.ll[i__ - 1] * probdata_1.llambda[i__ - 1] / 
		    probdata_1.mu[i__ - 1] * probdata_1.ff[i__ - 1] / (x[(k <<
		     1) + i__] * probdata_1.ff[i__ - 1] + 1.);
	}
	*c__ -= probdata_1.bwin;
    } else if (*ind == (k << 1) + 2) {
	*c__ = 0.;
	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    *c__ = *c__ + probdata_1.alphas * probdata_1.ll[i__ - 1] * 
		    probdata_1.llambda[i__ - 1] / probdata_1.mu[i__ - 1] * 
		    probdata_1.ff[i__ - 1] + probdata_1.betab * probdata_1.ll[
		    i__ - 1] * probdata_1.ff[i__ - 1] / (x[(k << 1) + i__] * 
		    probdata_1.ff[i__ - 1] + 1.);
	}
	*c__ -= probdata_1.bwout;
    }
} /* cache_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int cache_evaljac__(n, x, ind, indjac, valjac, nnzjac, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *indjac;
doublereal *valjac;
integer *nnzjac, *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static integer i__, k;

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine must compute the gradient of the ind-th constraint. */
/*     For achieving these objective YOU MUST MODIFY it in the way */
/*     specified below. */

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
    k = *n / 3;
    if (*ind <= k) {
	i__ = *ind;
	*nnzjac = 2;
	indjac[1] = (k << 1) + i__;
	valjac[1] = probdata_1.ff[i__ - 1] * exp(-probdata_1.ff[i__ - 1] * x[(
		k << 1) + i__]) - x[i__];
	indjac[2] = i__;
	valjac[2] = -x[(k << 1) + i__];
    } else if (*ind <= k << 1) {
	i__ = *ind - k;
	*nnzjac = 2;
	indjac[1] = (k << 1) + i__;
	valjac[1] = probdata_1.mu[i__ - 1] * exp(-probdata_1.mu[i__ - 1] * x[(
		k << 1) + i__]) - x[k + i__];
	indjac[2] = k + i__;
	valjac[2] = -x[(k << 1) + i__];
    } else if (*ind == (k << 1) + 1) {
	*nnzjac = k;
	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    indjac[i__] = (k << 1) + i__;
/* Computing 2nd power */
	    d__1 = probdata_1.ff[i__ - 1];
/* Computing 2nd power */
	    d__2 = x[(k << 1) + i__] * probdata_1.ff[i__ - 1] + 1.;
	    valjac[i__] = -(d__1 * d__1) * probdata_1.alphab * probdata_1.ll[
		    i__ - 1] * probdata_1.llambda[i__ - 1] / probdata_1.mu[
		    i__ - 1] / (d__2 * d__2);
	}
    } else if (*ind == (k << 1) + 2) {
	*nnzjac = k;
	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    indjac[i__] = (k << 1) + i__;
/* Computing 2nd power */
	    d__1 = probdata_1.ff[i__ - 1];
/* Computing 2nd power */
	    d__2 = x[(k << 1) + i__] * probdata_1.ff[i__ - 1] + 1.;
	    valjac[i__] = -probdata_1.betab * probdata_1.ll[i__ - 1] * (d__1 *
		     d__1) / (d__2 * d__2);
	}
    }
} /* cache_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int cache_evalhc__(n, x, ind, hclin, hccol, hcval, nnzhc, 
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
} /* cache_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int cache_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
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
} /* cache_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int cache_endp__(n, x, l, u, m, lambda, rho, equatn, linear)
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
} /* cache_endp__ */

