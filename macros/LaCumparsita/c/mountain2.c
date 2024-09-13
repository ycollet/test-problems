/* mountain2.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal pi[100000], pf[100000], dmax2;
    integer np;
} probdata_;

#define probdata_1 probdata_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;

/*     ================================================================= */
/*     File: mountain2.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     April 10, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Mountain Pass problem */
/*     --------------------- */
/*     Given a surface S(x,y) and initial and final points pi, pf in R^2 */
/*     the problem consists on finding a path pi, p1, p2, ..., pN, pf */
/*     from pi to pf such that max_{1 <= k <= N} S(pk) is as small as */
/*     possible. Moreover, the distance between consecutive points in the */
/*     path must be less than or equal to a prescribed tolerance. */

/*     The nonlinear programming formulation of the problem follows: */

/*     min z */

/*     subject to */

/*     d(pi     ,p1)^2 <= dmax^2 */
/*     d(p_{k-1},pk)^2 <= dmax^2, k = 2, ..., N */
/*     d(pN     ,pf)^2 <= dmax^2 */

/*     S(pk) <= z, k = 1, 2, ..., N */

/*     where dmax is the prescribed maximum distance and d(.,.) is the */
/*     Euclidian distance. */

/*     The problem has n = 2 N + 1 variables and m = 2 N + 1 inequality */
/*     constraints, where N is the number of intermediate points in the */
/*     path. */

/*     In the implementation, pi = (-10,-10), pf = (10,10). N is a user */
/*     defined parameter. dmax is also defined by the user and must */
/*     satisfy dmax >= d(pi,pf) / (N + 1). The initial values for p1, p2, */
/*     ..., pN correspond to a perturbation of the equally spaced points */
/*     in the segment [pi,pf]. The intial value of z is z = max { S(pi), */
/*     S(p1), S(p2), ..., S(pN), S(pf) }. The considered surfaces are */

/*     Mountain 1: S(x,y) = sin( x * y ) + sin( x + y ) */

/*     Mountain 2: S(x,y) = sin(x) + cos(y) */
/*     References: */

/*     [1] J. Horák, Constrained mountain pass algorithm for the */
/*     numerical solution of semilinear elliptic problems, Numerische */
/*     Mathematik 98, pp. 251-276, 2004. */

/*     [2] J. J. Moré and T. S. Munson, Computing mountain passes and */
/*     transition states, Mathematical Programming 100, pp. 151-182, 2004. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int mountain2_inip__(n, x, l, u, m, lambda, rho, equatn, 
	linear, np_in__, dmax2_in__, seed_in__)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
integer *np_in__;
doublereal *dmax2_in__, *seed_in__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal seed, dist;
    static integer b, i__, j;
    extern doublereal drand_();
    static doublereal tmp;
    extern doublereal mountain2_fmoun__();

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
/*     COMMON BLOCKS */
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
    probdata_1.pi[0] = -10.;
    probdata_1.pi[1] = -10.;
    probdata_1.pf[0] = 10.;
    probdata_1.pf[1] = 10.;
/* Computing 2nd power */
    d__1 = probdata_1.pi[0] - probdata_1.pf[0];
/* Computing 2nd power */
    d__2 = probdata_1.pi[1] - probdata_1.pf[1];
    dist = sqrt(d__1 * d__1 + d__2 * d__2);
/*     The distance between the extremes of the path is 'dist' */
/*     Number of points in the path excluding extremes (maximum 100,000) */
    probdata_1.np = *np_in__;
    dist /= probdata_1.np + 1;
/*     The minimum distance between consecutive points in the path is 'dist' */
    probdata_1.dmax2 = *dmax2_in__;
    if (probdata_1.dmax2 < 0.) {
	probdata_1.dmax2 = dist * 2.;
    }
/* Computing 2nd power */
    d__1 = probdata_1.dmax2;
    probdata_1.dmax2 = d__1 * d__1;
/*     Seed for the initial point random generation */
    seed = *seed_in__;
/*     Number of variables */
    *n = (probdata_1.np << 1) + 1;
/*     Initial point */
/* Computing MAX */
    d__1 = mountain2_fmoun__(probdata_1.pi, &probdata_1.pi[1]), d__2 = 
	    mountain2_fmoun__(probdata_1.pf, &probdata_1.pf[1]);
    x[*n] = max(d__1,d__2);
    i__1 = probdata_1.np;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b = i__ - 1 << 1;
	for (j = 1; j <= 2; ++j) {
	    tmp = probdata_1.pi[j - 1] + i__ * (probdata_1.pf[j - 1] - 
		    probdata_1.pi[j - 1]) / (probdata_1.np + 1);
	    x[b + j] = tmp * (drand_(&seed) * .1 + 1.);
	}
/* Computing MAX */
	d__1 = x[*n], d__2 = mountain2_fmoun__(&x[b + 1], &x[b + 2]);
	x[*n] = max(d__1,d__2);
    }
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e20;
	u[i__] = 1e20;
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = probdata_1.np + 1 + probdata_1.np;
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
	equatn[i__] = FALSE_;
    }
/*     For each constraint i, set linear(i) = .true. if it is a linear */
/*     constraint, otherwise set linear(i) = .false. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	linear[i__] = FALSE_;
    }
/* L100: */
} /* mountain2_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int mountain2_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
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
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    *f = x[*n];
} /* mountain2_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int mountain2_evalg__(n, x, g, flag__)
integer *n;
doublereal *x, *g;
integer *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

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
/*     COMMON BLOCKS */
/*     LOCAL SCALARS */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    *flag__ = 0;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = 0.;
    }
    g[*n] = 1.;
} /* mountain2_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int mountain2_evalh__(n, x, hlin, hcol, hval, nnzh, flag__)
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
    *flag__ = 0;
    *nnzh = 0;
} /* mountain2_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int mountain2_evalc__(n, x, ind, c__, flag__)
integer *n;
doublereal *x;
integer *ind;
doublereal *c__;
integer *flag__;
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static integer b, b1, b2;
    extern doublereal mountain2_fmoun__();

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
/*     COMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    if (*ind == 1) {
	b1 = *ind - 1 << 1;
/* Computing 2nd power */
	d__1 = x[b1 + 1] - probdata_1.pi[0];
/* Computing 2nd power */
	d__2 = x[b1 + 2] - probdata_1.pi[1];
	*c__ = d__1 * d__1 + d__2 * d__2 - probdata_1.dmax2;
    } else if (*ind <= probdata_1.np) {
	b1 = *ind - 1 << 1;
	b2 = *ind - 2 << 1;
/* Computing 2nd power */
	d__1 = x[b1 + 1] - x[b2 + 1];
/* Computing 2nd power */
	d__2 = x[b1 + 2] - x[b2 + 2];
	*c__ = d__1 * d__1 + d__2 * d__2 - probdata_1.dmax2;
    } else if (*ind == probdata_1.np + 1) {
	b2 = *ind - 2 << 1;
/* Computing 2nd power */
	d__1 = probdata_1.pf[0] - x[b2 + 1];
/* Computing 2nd power */
	d__2 = probdata_1.pf[1] - x[b2 + 2];
	*c__ = d__1 * d__1 + d__2 * d__2 - probdata_1.dmax2;
    } else if (*ind <= probdata_1.np + 1 + probdata_1.np) {
	b = *ind - probdata_1.np - 2 << 1;
	*c__ = mountain2_fmoun__(&x[b + 1], &x[b + 2]) - x[*n];
    }
} /* mountain2_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int mountain2_evaljac__(n, x, ind, indjac, valjac, nnzjac, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *indjac;
doublereal *valjac;
integer *nnzjac, *flag__;
{
    static integer b, j, b1, b2;
    extern /* Subroutine */ int mountain2_gmoun__();

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine must compute the gradient of the constraint ind. */
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
/*     COMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --valjac;
    --indjac;
    --x;

    /* Function Body */
    *flag__ = 0;
    if (*ind == 1) {
	*nnzjac = 2;
	b1 = *ind - 1 << 1;
	for (j = 1; j <= 2; ++j) {
	    indjac[j] = b1 + j;
	    valjac[j] = (x[b1 + j] - probdata_1.pi[j - 1]) * 2.;
	}
    } else if (*ind <= probdata_1.np) {
	*nnzjac = 4;
	b1 = *ind - 1 << 1;
	b2 = *ind - 2 << 1;
	for (j = 1; j <= 2; ++j) {
	    indjac[j] = b1 + j;
	    valjac[j] = (x[b1 + j] - x[b2 + j]) * 2.;
	    indjac[j + 2] = b2 + j;
	    valjac[j + 2] = (x[b1 + j] - x[b2 + j]) * -2.;
	}
    } else if (*ind == probdata_1.np + 1) {
	*nnzjac = 2;
	b2 = *ind - 2 << 1;
	for (j = 1; j <= 2; ++j) {
	    indjac[j] = b2 + j;
	    valjac[j] = (probdata_1.pf[j - 1] - x[b2 + j]) * -2.;
	}
    } else if (*ind <= probdata_1.np + 1 + probdata_1.np) {
	*nnzjac = 3;
	b = *ind - probdata_1.np - 2 << 1;
	mountain2_gmoun__(&x[b + 1], &x[b + 2], &valjac[1], &valjac[2]);
	indjac[1] = b + 1;
	indjac[2] = b + 2;
	indjac[3] = *n;
	valjac[3] = -1.;
    }
} /* mountain2_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int mountain2_evalhc__(n, x, ind, hclin, hccol, hcval, nnzhc,
	 flag__)
integer *n;
doublereal *x;
integer *ind, *hclin, *hccol;
doublereal *hcval;
integer *nnzhc, *flag__;
{
    static integer b, j, b1, b2;
    extern /* Subroutine */ int mountain2_hmoun__();

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
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;
    --hclin;
    --hccol;
    --hcval;

    /* Function Body */
    *flag__ = 0;
    if (*ind == 1) {
	*nnzhc = 2;
	b1 = *ind - 1 << 1;
	for (j = 1; j <= 2; ++j) {
	    hclin[j] = b1 + j;
	    hccol[j] = b1 + j;
	    hcval[j] = 2.;
	}
    } else if (*ind <= probdata_1.np) {
	*nnzhc = 6;
	b1 = *ind - 1 << 1;
	b2 = *ind - 2 << 1;
	for (j = 1; j <= 2; ++j) {
	    hclin[j] = b1 + j;
	    hccol[j] = b1 + j;
	    hcval[j] = 2.;
	    hclin[j + 2] = b2 + j;
	    hccol[j + 2] = b2 + j;
	    hcval[j + 2] = 2.;
	    hclin[j + 4] = b1 + j;
	    hccol[j + 4] = b2 + j;
	    hcval[j + 4] = -2.;
	}
    } else if (*ind == probdata_1.np + 1) {
	*nnzhc = 2;
	b2 = *ind - 2 << 1;
	for (j = 1; j <= 2; ++j) {
	    hclin[j] = b2 + j;
	    hccol[j] = b2 + j;
	    hcval[j] = 2.;
	}
    } else if (*ind <= probdata_1.np + 1 + probdata_1.np) {
	*nnzhc = 3;
	b = *ind - probdata_1.np - 2 << 1;
	mountain2_hmoun__(&x[b + 1], &x[b + 2], &hcval[1], &hcval[2], &hcval[
		3]);
	hclin[1] = b + 1;
	hccol[1] = b + 1;
	hclin[2] = b + 2;
	hccol[2] = b + 2;
	hclin[3] = b + 2;
	hccol[3] = b + 1;
    }
} /* mountain2_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int mountain2_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
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
} /* mountain2_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int mountain2_endp__(n, x, l, u, m, lambda, rho, equatn, 
	linear)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsle(), e_wsle(), do_lio();

    /* Local variables */
    static doublereal fmax;
    static integer imax, b;
    static doublereal f, g[2], h__[4]	/* was [2][2] */;
    static integer i__;
    extern doublereal mountain2_fmoun__();
    extern /* Subroutine */ int mountain2_gmoun__(), mountain2_hmoun__();

    /* Fortran I/O blocks */
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, 0, 0 };
    static cilist io___28 = { 0, 6, 0, 0, 0 };
    static cilist io___30 = { 0, 6, 0, 0, 0 };
    static cilist io___31 = { 0, 6, 0, 0, 0 };
    static cilist io___32 = { 0, 6, 0, 0, 0 };
    static cilist io___33 = { 0, 6, 0, 0, 0 };
    static cilist io___34 = { 0, 6, 0, 0, 0 };
    static cilist io___35 = { 0, 6, 0, 0, 0 };


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
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMON ARRAYS */
/*     LOCAL SCALARS */
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --u;
    --l;
    --x;
    --linear;
    --equatn;
    --rho;
    --lambda;

    /* Function Body */
    s_wsle(&io___19);
    e_wsle();
    fmax = -1e99;
    i__1 = probdata_1.np;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b = i__ - 1 << 1;
	f = mountain2_fmoun__(&x[b + 1], &x[b + 2]);
	if (f > fmax) {
	    fmax = f;
	    imax = i__;
	}
    }
    b = imax - 1 << 1;
    s_wsle(&io___25);
    do_lio(&c__9, &c__1, "Maximizer on optimal path: ", (ftnlen)27);
    do_lio(&c__5, &c__1, (char *)&x[b + 1], (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&x[b + 2], (ftnlen)sizeof(doublereal));
    e_wsle();
    s_wsle(&io___26);
    do_lio(&c__9, &c__1, "Maximum: ", (ftnlen)9);
    do_lio(&c__5, &c__1, (char *)&fmax, (ftnlen)sizeof(doublereal));
    e_wsle();
    mountain2_gmoun__(&x[b + 1], &x[b + 2], g, &g[1]);
    s_wsle(&io___28);
    do_lio(&c__9, &c__1, "Gradient of the energy: ", (ftnlen)24);
    do_lio(&c__5, &c__1, (char *)&g[0], (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&g[1], (ftnlen)sizeof(doublereal));
    e_wsle();
    mountain2_hmoun__(&x[b + 1], &x[b + 2], h__, &h__[3], &h__[2]);
    s_wsle(&io___30);
    do_lio(&c__9, &c__1, "Hessian of the energy: ", (ftnlen)23);
    e_wsle();
    s_wsle(&io___31);
    do_lio(&c__9, &c__1, "Determinant: ", (ftnlen)13);
/* Computing 2nd power */
    d__2 = h__[2];
    d__1 = h__[0] * h__[3] - d__2 * d__2;
    do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsle();
    if (imax != 1) {
	b = imax - 2 << 1;
	s_wsle(&io___32);
	do_lio(&c__9, &c__1, "Previous point in the path: ", (ftnlen)28);
	do_lio(&c__5, &c__1, (char *)&x[b + 1], (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&x[b + 2], (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___33);
	do_lio(&c__9, &c__1, "Energy at previous point in the path: ", (
		ftnlen)38);
	d__1 = mountain2_fmoun__(&x[b + 1], &x[b + 2]);
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
    }
    if (imax != probdata_1.np) {
	b = imax << 1;
	s_wsle(&io___34);
	do_lio(&c__9, &c__1, "Following point in the path: ", (ftnlen)29);
	do_lio(&c__5, &c__1, (char *)&x[b + 1], (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&x[b + 2], (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___35);
	do_lio(&c__9, &c__1, "Energy at following point in the path: ", (
		ftnlen)39);
	d__1 = mountain2_fmoun__(&x[b + 1], &x[b + 2]);
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
    }
} /* mountain2_endp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
doublereal mountain2_fmoun__(x, y)
doublereal *x, *y;
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sin(), cos();

/*     SCALAR ARGUMENTS */
    ret_val = sin(*x) + cos(*y);
    return ret_val;
} /* mountain2_fmoun__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int mountain2_gmoun__(x, y, dfdx, dfdy)
doublereal *x, *y, *dfdx, *dfdy;
{
    /* Builtin functions */
    double cos(), sin();

/*     SCALAR ARGUMENTS */
    *dfdx = cos(*x);
    *dfdy = -sin(*y);
} /* mountain2_gmoun__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int mountain2_hmoun__(x, y, dfdxx, dfdyy, dfdxy)
doublereal *x, *y, *dfdxx, *dfdyy, *dfdxy;
{
    /* Builtin functions */
    double sin(), cos();

/*     SCALAR ARGUMENTS */
    *dfdxx = -sin(*x);
    *dfdyy = -cos(*y);
    *dfdxy = 0.;
} /* mountain2_hmoun__ */

