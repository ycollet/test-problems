/* packccmn.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer nite, ndim;
    doublereal iterad, objrad;
} packdata_;

#define packdata_1 packdata_

/*     ================================================================= */
/*     File: packccmn.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     April 7, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Packing fixed-dimension circular items within a fixed-dimension */
/*     --------------------------------------------------------------- */
/*     circular object maximizing the number of items */
/*     ---------------------------------------------- */

/*     We wish to place k circles of radii ri (i=1, ..., k) into a */
/*     circular object with radius R in such a way that the circles are */
/*     not overlapped. Therefore, given k, the radii ri (i=1, ..., k) */
/*     and R, the goal is to solve the problem: */

/*     Minimize \sum_{i not equal j} max(0, (ri+rj)^2 - d(pi,pj)^2 )^2 */

/*     subject to d(0,pi)^2 <= (R - ri)^2, i = 1, ..., k, */

/*     where d(.,.) is the Euclidian distance. If the objective function */
/*     value at the global minimizer of this problem is zero then the */
/*     answer of the decision problem: "Given k circular items of radii */
/*     ri (i=1, ..., k) and a circular object of radius R, whether is it */
/*     possible to locate all the items within the object or not." is YES, */
/*     otherwise, the answer is NO. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packccmn_inip__(n, x, l, u, m, lambda, rho, equatn, 
	linear)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(), sin();

    /* Local variables */
    static doublereal seed;
    static integer i__;
    static doublereal r__, t;
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
/*     COMMON SCALARS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     SET PROBLEM DATA */
/*     Dimension of the space */
    /* Parameter adjustments */
    --linear;
    --equatn;
    --rho;
    --lambda;
    --u;
    --l;
    --x;

    /* Function Body */
    packdata_1.ndim = 2;
/*     Number of identical circular items to be packed */
    packdata_1.nite = 80;
/*     Radius of the circular items to be packed */
    packdata_1.iterad = 1.;
/*     Radius of the circular object within which the items will */
/*     be packed */
    packdata_1.objrad = 10.;
/*     Number of variables */
    *n = packdata_1.nite << 1;
/*     Initial point */
    seed = 12337.;
    i__1 = packdata_1.nite;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__ = (packdata_1.objrad - packdata_1.iterad) * drand_(&seed);
	t = drand_(&seed) * 6.2831799999999998;
	x[(i__ << 1) - 1] = r__ * cos(t);
	x[i__ * 2] = r__ * sin(t);
    }
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e20;
	u[i__] = 1e20;
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = packdata_1.nite;
/*     Lagrange multipliers approximation */
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
} /* packccmn_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packccmn_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal dist;
    static integer i__, j, k;
    static doublereal fparc;
    static integer ind1, ind2;

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
/*     COMMON SCALARS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
/*     COMPUTE DENSE OVERLAPPING */
    *f = 0.;
    i__1 = packdata_1.nite;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = packdata_1.nite;
	for (j = i__ + 1; j <= i__2; ++j) {
	    dist = 0.;
	    i__3 = packdata_1.ndim;
	    for (k = 1; k <= i__3; ++k) {
		ind1 = (i__ - 1) * packdata_1.ndim + k;
		ind2 = (j - 1) * packdata_1.ndim + k;
/* Computing 2nd power */
		d__1 = x[ind1] - x[ind2];
		dist += d__1 * d__1;
	    }
/* Computing MAX */
/* Computing 2nd power */
	    d__3 = packdata_1.iterad * 2.;
	    d__1 = 0., d__2 = d__3 * d__3 - dist;
	    fparc = max(d__1,d__2);
/* Computing 2nd power */
	    d__1 = fparc;
	    *f += d__1 * d__1;
	}
    }
} /* packccmn_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packccmn_evalg__(n, x, g, flag__)
integer *n;
doublereal *x, *g;
integer *flag__;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal dist;
    static integer i__, j, k;
    static doublereal fparc, tmp;
    static integer ind1, ind2;

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
/*     COMMON SCALARS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    *flag__ = 0;
/*     COMPUTE DENSE OVERLAPPING */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = 0.;
    }
    i__1 = packdata_1.nite;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = packdata_1.nite;
	for (j = i__ + 1; j <= i__2; ++j) {
	    dist = 0.;
	    i__3 = packdata_1.ndim;
	    for (k = 1; k <= i__3; ++k) {
		ind1 = (i__ - 1) * packdata_1.ndim + k;
		ind2 = (j - 1) * packdata_1.ndim + k;
/* Computing 2nd power */
		d__1 = x[ind1] - x[ind2];
		dist += d__1 * d__1;
	    }
/* Computing MAX */
/* Computing 2nd power */
	    d__3 = packdata_1.iterad * 2.;
	    d__1 = 0., d__2 = d__3 * d__3 - dist;
	    fparc = max(d__1,d__2);
	    if (fparc != 0.) {
		i__3 = packdata_1.ndim;
		for (k = 1; k <= i__3; ++k) {
		    ind1 = (i__ - 1) * packdata_1.ndim + k;
		    ind2 = (j - 1) * packdata_1.ndim + k;
		    tmp = fparc * 4. * (x[ind1] - x[ind2]);
		    g[ind1] -= tmp;
		    g[ind2] += tmp;
		}
	    }
	}
    }
} /* packccmn_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packccmn_evalh__(n, x, hlin, hcol, hval, nnzh, flag__)
integer *n;
doublereal *x;
integer *hlin, *hcol;
doublereal *hval;
integer *nnzh, *flag__;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal dist, diagb[40000]	/* was [2][2][10000] */;
    static integer i__, j, k, l;
    static doublereal tmp;
    static integer ind1, ind2, ind3, ind4;

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
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     LOCAL SCALARS */
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;
    --hlin;
    --hcol;
    --hval;

    /* Function Body */
    *flag__ = 0;
/*     INITALIZE DIAGONAL BLOCKS */
    i__1 = packdata_1.nite;
    for (k = 1; k <= i__1; ++k) {
	i__2 = packdata_1.ndim;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = packdata_1.ndim;
	    for (i__ = j; i__ <= i__3; ++i__) {
		diagb[i__ + (j + (k << 1) << 1) - 7] = 0.;
	    }
	}
    }
/*     COMPUTE DENSE OVERLAPPING SECOND DERIVATIVES */
    *nnzh = 0;
    i__1 = packdata_1.nite;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = packdata_1.nite;
	for (j = i__ + 1; j <= i__2; ++j) {
	    dist = 0.;
	    i__3 = packdata_1.ndim;
	    for (k = 1; k <= i__3; ++k) {
		ind1 = (i__ - 1) * packdata_1.ndim + k;
		ind2 = (j - 1) * packdata_1.ndim + k;
/* Computing 2nd power */
		d__1 = x[ind1] - x[ind2];
		dist += d__1 * d__1;
	    }
/* Computing 2nd power */
	    d__1 = packdata_1.iterad * 2.;
	    if (dist <= d__1 * d__1) {
		i__3 = packdata_1.ndim;
		for (k = 1; k <= i__3; ++k) {
		    ind1 = (i__ - 1) * packdata_1.ndim + k;
		    ind2 = (j - 1) * packdata_1.ndim + k;
/* Computing 2nd power */
		    d__1 = x[ind1] - x[ind2];
/* Computing 2nd power */
		    d__2 = packdata_1.iterad * 2.;
		    tmp = d__1 * d__1 * 8. - (d__2 * d__2 - dist) * 4.;
		    if (tmp != 0.) {
/*                         H(ind1,ind1) = H(ind1,ind1) + tmp */
			diagb[k + (k + (i__ << 1) << 1) - 7] += tmp;
/*                         H(ind2,ind2) = H(ind2,ind2) + tmp */
			diagb[k + (k + (j << 1) << 1) - 7] += tmp;
/*                         H(ind2,ind1) = H(ind2,ind1) - tmp */
			++(*nnzh);
			hlin[*nnzh] = ind2;
			hcol[*nnzh] = ind1;
			hval[*nnzh] = -tmp;
		    }
		    i__4 = k - 1;
		    for (l = 1; l <= i__4; ++l) {
			ind3 = (i__ - 1) * packdata_1.ndim + l;
			ind4 = (j - 1) * packdata_1.ndim + l;
			tmp = (x[ind3] - x[ind4]) * 8. * (x[ind1] - x[ind2]);
			if (tmp != 0.) {
/*                             H(ind1,ind3) = H(ind1,ind3) + tmp */
			    diagb[k + (l + (i__ << 1) << 1) - 7] += tmp;
/*                             H(ind2,ind4) = H(ind2,ind4) + tmp */
			    diagb[k + (l + (j << 1) << 1) - 7] += tmp;
/*                             H(ind2,ind3) = H(ind2,ind3) - tmp */
			    ++(*nnzh);
			    hlin[*nnzh] = ind2;
			    hcol[*nnzh] = ind3;
			    hval[*nnzh] = -tmp;
			}
		    }
		    i__4 = packdata_1.ndim;
		    for (l = k + 1; l <= i__4; ++l) {
			ind3 = (i__ - 1) * packdata_1.ndim + l;
			ind4 = (j - 1) * packdata_1.ndim + l;
			tmp = (x[ind3] - x[ind4]) * 8. * (x[ind1] - x[ind2]);
			if (tmp != 0.) {
/*                             H(ind2,ind3) = H(ind2,ind3) - tmp */
			    ++(*nnzh);
			    hlin[*nnzh] = ind2;
			    hcol[*nnzh] = ind3;
			    hval[*nnzh] = -tmp;
			}
		    }
		}
	    }
	}
    }
/*     ADD DIAGONAL BLOCKS TO THE SPARSE STRUCTURE */
    i__1 = packdata_1.nite;
    for (k = 1; k <= i__1; ++k) {
	i__2 = packdata_1.ndim;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = packdata_1.ndim;
	    for (i__ = j; i__ <= i__3; ++i__) {
		if (diagb[i__ + (j + (k << 1) << 1) - 7] != 0.) {
		    ++(*nnzh);
		    hlin[*nnzh] = (k - 1) * packdata_1.ndim + i__;
		    hcol[*nnzh] = (k - 1) * packdata_1.ndim + j;
		    hval[*nnzh] = diagb[i__ + (j + (k << 1) << 1) - 7];
		}
	    }
	}
    }
} /* packccmn_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packccmn_evalc__(n, x, ind, c__, flag__)
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
    static integer i__;

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
/*     COMMON SCALARS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    if (0 < *ind && *ind <= packdata_1.nite) {
/* Computing 2nd power */
	d__1 = packdata_1.objrad - packdata_1.iterad;
	*c__ = -(d__1 * d__1);
	i__1 = packdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = x[packdata_1.ndim * (*ind - 1) + i__];
	    *c__ += d__1 * d__1;
	}
	return 0;
    }
} /* packccmn_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packccmn_evaljac__(n, x, ind, indjac, valjac, nnzjac, 
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
    static integer i__;

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
/*     COMMON SCALARS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --valjac;
    --indjac;
    --x;

    /* Function Body */
    *flag__ = 0;
    if (0 < *ind && *ind <= packdata_1.nite) {
	*nnzjac = packdata_1.ndim;
	i__1 = packdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    indjac[i__] = packdata_1.ndim * (*ind - 1) + i__;
	    valjac[i__] = x[packdata_1.ndim * (*ind - 1) + i__] * 2.;
	}
	return 0;
    }
} /* packccmn_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packccmn_evalhc__(n, x, ind, hclin, hccol, hcval, nnzhc, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *hclin, *hccol;
doublereal *hcval;
integer *nnzhc, *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;

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
/*     COMMON SCALARS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;
    --hclin;
    --hccol;
    --hcval;

    /* Function Body */
    *flag__ = 0;
    if (0 < *ind && *ind <= packdata_1.nite) {
	*nnzhc = packdata_1.ndim;
	i__1 = packdata_1.ndim;
	for (k = 1; k <= i__1; ++k) {
	    hclin[k] = packdata_1.ndim * (*ind - 1) + k;
	    hccol[k] = packdata_1.ndim * (*ind - 1) + k;
	    hcval[k] = 2.;
	}
	return 0;
    }
} /* packccmn_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packccmn_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
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
} /* packccmn_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packccmn_endp__(n, x, l, u, m, lambda, rho, equatn, 
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
} /* packccmn_endp__ */

