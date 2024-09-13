/* piecefit.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal uobs[2000], usol[2000], h__;
    integer p;
} ovo_;

#define ovo_1 ovo_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;

/*     ================================================================= */
/*     File: piecefit.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     July 10th, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Piecefit problem */
/*     ---------------- */
/*     We defined usol, a piecewise linear function defined on [0,4] */
/*     with kinks in 1 and 3 and discotinuity in 2. We discretized [0,4] */
/*     with n points, including extremes. We defined uobs, a 10% */
/*     perturbation of usol. The objective function is the sum of the p */
/*     smaller squares of the discretized second derivatives of x. p, */
/*     the OVO parameter, is given by the user. The restriction is that */
/*     the root mean squared deviation of the solution found with */
/*     respect to uobs is smaller than 0.2. In other words, using LOVO, */
/*     this code aims to find a piecewise harmonic one-dimensional */
/*     function. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int piecefit_inip__(n, x, l, u, m, lambda, rho, equatn, 
	linear, n_in__, p_in__)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
integer *n_in__, *p_in__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal seed;
    static integer i__;
    static doublereal t;
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
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     Number of variables (points in [0,4]) */
    /* Parameter adjustments */
    --linear;
    --equatn;
    --rho;
    --lambda;
    --u;
    --l;
    --x;

    /* Function Body */
    *n = *n_in__;
/*     p (integer close to n-2), perhaps n-3 */
    ovo_1.p = *p_in__;
/*     Generation of uobs */
    seed = 17172937.;
    ovo_1.h__ = 4. / (doublereal) (*n - 1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = (doublereal) (i__ - 1) * ovo_1.h__;
	if (t <= 1.) {
	    ovo_1.uobs[i__ - 1] = t + 1.;
	} else if (t <= 2.) {
	    ovo_1.uobs[i__ - 1] = 1. - (t - 2.);
	} else if (t <= 3.) {
	    ovo_1.uobs[i__ - 1] = 2.;
	} else {
	    ovo_1.uobs[i__ - 1] = 2. - (t - 3.);
	}
	ovo_1.usol[i__ - 1] = ovo_1.uobs[i__ - 1];
	ovo_1.uobs[i__ - 1] += ovo_1.uobs[i__ - 1] * .1 * (drand_(&seed) * 2. 
		- 1.);
    }
/*     write(*,*) 'uobs: ',(uobs(i),i=1,n) */
/*     Initial point */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = ovo_1.uobs[i__ - 1];
    }
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e20;
	u[i__] = 1e20;
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = 1;
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
	equatn[i__] = FALSE_;
    }
/*     For each constraint i, set linear(i) = .true. if it is a linear */
/*     constraint, otherwise set linear(i) = .false. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	linear[i__] = FALSE_;
    }
} /* piecefit_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int piecefit_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer imax;
    static doublereal zmax;
    static integer i__, j;
    static doublereal resid[2000];

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
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    *f = 0.;
    i__1 = *n - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = ovo_1.h__;
	resid[i__ - 1] = (x[i__ + 2] - x[i__ + 1] * 2. + x[i__]) / (d__1 * 
		d__1);
    }
    i__1 = ovo_1.p;
    for (j = 1; j <= i__1; ++j) {
	zmax = 1e99;
	i__2 = *n - 2;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if ((d__1 = resid[i__ - 1], abs(d__1)) < zmax) {
		zmax = (d__1 = resid[i__ - 1], abs(d__1));
		imax = i__;
	    }
	}
/* Computing 2nd power */
	d__1 = resid[imax - 1];
	*f += d__1 * d__1;
	resid[imax - 1] = 1e99;
    }
} /* piecefit_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int piecefit_evalg__(n, x, g, flag__)
integer *n;
doublereal *x, *g;
integer *flag__;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer imax;
    static doublereal zmax;
    static integer i__, j;
    static doublereal resid[2000];

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
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    *flag__ = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = 0.;
    }
    i__1 = *n - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = ovo_1.h__;
	resid[i__ - 1] = (x[i__ + 2] - x[i__ + 1] * 2. + x[i__]) / (d__1 * 
		d__1);
    }
    i__1 = ovo_1.p;
    for (j = 1; j <= i__1; ++j) {
	zmax = 1e99;
	i__2 = *n - 2;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if ((d__1 = resid[i__ - 1], abs(d__1)) < zmax) {
		zmax = (d__1 = resid[i__ - 1], abs(d__1));
		imax = i__;
	    }
	}
/* Computing 2nd power */
	d__1 = ovo_1.h__;
	g[imax + 2] += resid[imax - 1] * 2. / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = ovo_1.h__;
	g[imax + 1] -= resid[imax - 1] * 4. / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = ovo_1.h__;
	g[imax] += resid[imax - 1] * 2. / (d__1 * d__1);
	resid[imax - 1] = 1e99;
    }
} /* piecefit_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int piecefit_evalh__(n, x, hlin, hcol, hval, nnzh, flag__)
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
} /* piecefit_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int piecefit_evalc__(n, x, ind, c__, flag__)
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
    static doublereal z__;

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
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    z__ = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = x[i__] - ovo_1.uobs[i__ - 1];
	z__ += d__1 * d__1;
    }
    z__ /= (doublereal) (*n);
    *c__ = z__ - .01;
} /* piecefit_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int piecefit_evaljac__(n, x, ind, indjac, valjac, nnzjac, 
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
    *nnzjac = *n;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	indjac[i__] = i__;
	valjac[i__] = (x[i__] - ovo_1.uobs[i__ - 1]) * 2. / (doublereal) (*n);
    }
} /* piecefit_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int piecefit_evalhc__(n, x, ind, hclin, hccol, hcval, nnzhc, 
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
} /* piecefit_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int piecefit_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
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
} /* piecefit_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int piecefit_endp__(n, x, l, u, m, lambda, rho, equatn, 
	linear)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static doublereal zmax;
    static integer i__;
    static doublereal z__;

    /* Fortran I/O blocks */
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };


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
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
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
    zmax = 0.;
    s_wsle(&io___18);
    do_lio(&c__9, &c__1, "i, t, sol-found, true-sol, second-derivative", (
	    ftnlen)44);
    e_wsle();
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == 1 || i__ == *n) {
	    z__ = 0.;
	} else {
/* Computing 2nd power */
	    d__1 = ovo_1.h__;
	    z__ = (x[i__ + 1] - x[i__] * 2. + x[i__ - 1]) / (d__1 * d__1);
	}
	s_wsle(&io___21);
	do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	d__1 = (i__ - 1) * ovo_1.h__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&ovo_1.usol[i__ - 1], (ftnlen)sizeof(
		doublereal));
	do_lio(&c__5, &c__1, (char *)&z__, (ftnlen)sizeof(doublereal));
	e_wsle();
	zmax = max(zmax,z__);
    }
    s_wsle(&io___22);
    do_lio(&c__9, &c__1, " Maximum second derivative: ", (ftnlen)28);
    do_lio(&c__5, &c__1, (char *)&zmax, (ftnlen)sizeof(doublereal));
    e_wsle();
} /* piecefit_endp__ */

