/* pedazos4.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal t[30], y[30];
} probdata_;

#define probdata_1 probdata_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;

/*     ================================================================= */
/*     File: pedazos4.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     July 20th, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Pedazos4 problem: Fit a continuous piecewise polynomial */
/*     ------------------------------------------------------- */
/*     We generate a table (t_i, y_i) with 30 equally spaced points */
/*     between 1 and 30. The first 10 points are on a line, the second */
/*     10 on a parabola and the third 10 on a cubic. The overall */
/*     function is continuous at 10 and 20. We fit a line to the first */
/*     10, a quadratic to the second 10 and a cubic to the third ten, */
/*     with the conditions that we preserve continuity at 10 and 20. So, */
/*     we have 9 unknowns and 2 equality constraints. The objective */
/*     function is quadratic and the constraints are linear. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int pedazos4_inip__(n, x, l, u, m, lambda, rho, equatn, 
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

    /* Local variables */
    static doublereal seed;
    static integer i__;
    extern doublereal drand_();

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine sets some problem data. */

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
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     Number of variables */
    /* Parameter adjustments */
    --linear;
    --equatn;
    --rho;
    --lambda;
    --u;
    --l;
    --x;

    /* Function Body */
    *n = 9;
/*     Adjust a line, a parabola and a cubic */
    for (i__ = 1; i__ <= 10; ++i__) {
	probdata_1.t[i__ - 1] = (doublereal) i__;
	probdata_1.y[i__ - 1] = probdata_1.t[i__ - 1];
    }
    for (i__ = 11; i__ <= 20; ++i__) {
	probdata_1.t[i__ - 1] = (doublereal) i__;
/* Computing 2nd power */
	d__1 = probdata_1.t[i__ - 1] - 15.;
	probdata_1.y[i__ - 1] = d__1 * d__1 * .4;
    }
    for (i__ = 21; i__ <= 30; ++i__) {
	probdata_1.t[i__ - 1] = (doublereal) i__;
/* Computing 3rd power */
	d__1 = probdata_1.t[i__ - 1] - 25.;
	probdata_1.y[i__ - 1] = d__1 * (d__1 * d__1) * -10. / 125.;
    }
/*     Perturbation */
    seed = 17172937.;
    for (i__ = 1; i__ <= 30; ++i__) {
	probdata_1.y[i__ - 1] += probdata_1.y[i__ - 1] * .2 * (drand_(&seed) *
		 2. - 1.);
    }
/*     Initial point */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = 0.;
    }
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e20;
	u[i__] = 1e20;
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = 2;
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
	linear[i__] = TRUE_;
    }
} /* pedazos4_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int pedazos4_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__;
    static doublereal r__;

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine computes the objective function. */

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
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    *f = 0.;
    for (i__ = 1; i__ <= 30; ++i__) {
	if (probdata_1.t[i__ - 1] <= 10.) {
/* Computing 2nd power */
	    d__1 = x[1] + x[2] * probdata_1.t[i__ - 1] - probdata_1.y[i__ - 1]
		    ;
	    r__ = d__1 * d__1;
	} else if (probdata_1.t[i__ - 1] <= 20.) {
/* Computing 2nd power */
	    d__2 = probdata_1.t[i__ - 1];
/* Computing 2nd power */
	    d__1 = x[3] + x[4] * probdata_1.t[i__ - 1] + x[5] * (d__2 * d__2) 
		    - probdata_1.y[i__ - 1];
	    r__ = d__1 * d__1;
	} else {
/* Computing 2nd power */
	    d__2 = probdata_1.t[i__ - 1];
/* Computing 3rd power */
	    d__3 = probdata_1.t[i__ - 1];
/* Computing 2nd power */
	    d__1 = x[6] + x[7] * probdata_1.t[i__ - 1] + x[8] * (d__2 * d__2) 
		    + x[9] * (d__3 * (d__3 * d__3)) - probdata_1.y[i__ - 1];
	    r__ = d__1 * d__1;
	}
	*f += r__;
    }
} /* pedazos4_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int pedazos4_evalg__(n, x, g, flag__)
integer *n;
doublereal *x, *g;
integer *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal r__;

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine computes the gradient vector of the objective */
/*     function. */

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
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
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
    for (i__ = 1; i__ <= 30; ++i__) {
	if (probdata_1.t[i__ - 1] <= 10.) {
	    r__ = x[1] + x[2] * probdata_1.t[i__ - 1] - probdata_1.y[i__ - 1];
	    g[1] += r__ * 2.;
	    g[2] += r__ * 2. * probdata_1.t[i__ - 1];
	} else if (probdata_1.t[i__ - 1] <= 20.) {
/* Computing 2nd power */
	    d__1 = probdata_1.t[i__ - 1];
	    r__ = x[3] + x[4] * probdata_1.t[i__ - 1] + x[5] * (d__1 * d__1) 
		    - probdata_1.y[i__ - 1];
	    g[3] += r__ * 2.;
	    g[4] += r__ * 2. * probdata_1.t[i__ - 1];
/* Computing 2nd power */
	    d__1 = probdata_1.t[i__ - 1];
	    g[5] += r__ * 2. * (d__1 * d__1);
	} else {
/* Computing 2nd power */
	    d__1 = probdata_1.t[i__ - 1];
/* Computing 3rd power */
	    d__2 = probdata_1.t[i__ - 1];
	    r__ = x[6] + x[7] * probdata_1.t[i__ - 1] + x[8] * (d__1 * d__1) 
		    + x[9] * (d__2 * (d__2 * d__2)) - probdata_1.y[i__ - 1];
	    g[6] += r__ * 2.;
	    g[7] += r__ * 2. * probdata_1.t[i__ - 1];
/* Computing 2nd power */
	    d__1 = probdata_1.t[i__ - 1];
	    g[8] += r__ * 2. * (d__1 * d__1);
/* Computing 3rd power */
	    d__1 = probdata_1.t[i__ - 1];
	    g[9] += r__ * 2. * (d__1 * (d__1 * d__1));
	}
    }
} /* pedazos4_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int pedazos4_evalh__(n, x, hlin, hcol, hval, nnzh, flag__)
integer *n;
doublereal *x;
integer *hlin, *hcol;
doublereal *hval;
integer *nnzh, *flag__;
{
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine computes the Hessian matrix of the objective */
/*     function. */

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
    *nnzh = 19;
    hlin[1] = 1;
    hcol[1] = 1;
    hval[1] = 20.;
    hlin[2] = 2;
    hcol[2] = 1;
    hval[2] = 110.;
    hlin[3] = 2;
    hcol[3] = 2;
    hval[3] = 770.;
    hlin[4] = 3;
    hcol[4] = 3;
    hval[4] = 20.;
    hlin[5] = 4;
    hcol[5] = 3;
    hval[5] = 310.;
    hlin[6] = 4;
    hcol[6] = 4;
    hval[6] = 4970.;
    hlin[7] = 5;
    hcol[7] = 3;
    hval[7] = 4970.;
    hlin[8] = 5;
    hcol[8] = 4;
    hval[8] = 82150.;
    hlin[9] = 5;
    hcol[9] = 5;
    hval[9] = 1394666.;
    hlin[10] = 6;
    hcol[10] = 6;
    hval[10] = 20.;
    hlin[11] = 7;
    hcol[11] = 6;
    hval[11] = 510.;
    hlin[12] = 7;
    hcol[12] = 7;
    hval[12] = 13170.;
    hlin[13] = 8;
    hcol[13] = 6;
    hval[13] = 13170.;
    hlin[14] = 8;
    hcol[14] = 7;
    hval[14] = 344250.;
    hlin[15] = 8;
    hcol[15] = 8;
    hval[15] = 9102666.;
    hlin[16] = 9;
    hcol[16] = 6;
    hval[16] = 344250.;
    hlin[17] = 9;
    hcol[17] = 7;
    hval[17] = 9102666.;
    hlin[18] = 9;
    hcol[18] = 8;
    hval[18] = 243308250.;
    hlin[19] = 9;
    hcol[19] = 9;
    hval[19] = 6568950810.;
} /* pedazos4_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int pedazos4_evalc__(n, x, ind, c__, flag__)
integer *n;
doublereal *x;
integer *ind;
doublereal *c__;
integer *flag__;
{
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
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    if (*ind == 1) {
	*c__ = x[1] + x[2] * 10. - (x[3] + x[4] * 10. + x[5] * 100.);
    } else if (*ind == 2) {
	*c__ = x[3] + x[4] * 20. + x[5] * 400. - (x[6] + x[7] * 20. + x[8] * 
		400. + x[9] * 8e3);
    }
} /* pedazos4_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int pedazos4_evaljac__(n, x, ind, indjac, valjac, nnzjac, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *indjac;
doublereal *valjac;
integer *nnzjac, *flag__;
{
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
    /* Parameter adjustments */
    --valjac;
    --indjac;
    --x;

    /* Function Body */
    *flag__ = 0;
    if (*ind == 1) {
	*nnzjac = 5;
	indjac[1] = 1;
	valjac[1] = 1.;
	indjac[2] = 2;
	valjac[2] = 10.;
	indjac[3] = 3;
	valjac[3] = -1.;
	indjac[4] = 4;
	valjac[4] = -10.;
	indjac[5] = 5;
	valjac[5] = -100.;
    } else if (*ind == 2) {
	*nnzjac = 7;
	indjac[1] = 3;
	valjac[1] = 1.;
	indjac[2] = 4;
	valjac[2] = 20.;
	indjac[3] = 5;
	valjac[3] = 400.;
	indjac[4] = 6;
	valjac[4] = -1.;
	indjac[5] = 7;
	valjac[5] = -20.;
	indjac[6] = 8;
	valjac[6] = -400.;
	indjac[7] = 9;
	valjac[7] = -8e3;
    }
} /* pedazos4_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int pedazos4_evalhc__(n, x, ind, hclin, hccol, hcval, nnzhc, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *hclin, *hccol;
doublereal *hcval;
integer *nnzhc, *flag__;
{
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine computes the Hessian matrix of the ind-th */
/*     constraint. */

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
    *flag__ = 0;
    *nnzhc = 0;
} /* pedazos4_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int pedazos4_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
integer *n;
doublereal *x;
integer *m;
doublereal *lambda, *p, *hp;
logical *goth;
integer *flag__;
{
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine computes the product of the Hessian of the */
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
} /* pedazos4_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int pedazos4_endp__(n, x, l, u, m, lambda, rho, equatn, 
	linear)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static integer i__;
    static doublereal z__, tt;

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };


/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine can be used to do some extra job after the solver */
/*     has found the solution, like some extra statistics, or to save the */
/*     solution in some special format or to draw some graphical */
/*     representation of the solution. If the information given by the */
/*     solver is enough for you then leave the body of this subroutine */
/*     empty. */

/*     Parameters of the subroutine: */

/*     The paraemters of this subroutine are the same parameters of */
/*     subroutine inip. But in this subroutine there are not output */
/*     parameter. All the parameters are input parameters. */
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
    s_wsle(&io___7);
    do_lio(&c__9, &c__1, "t, observed y, computed y", (ftnlen)25);
    e_wsle();
    for (i__ = 1; i__ <= 30; ++i__) {
	if (probdata_1.t[i__ - 1] <= 10.) {
	    z__ = x[1] + x[2] * probdata_1.t[i__ - 1];
	} else if (probdata_1.t[i__ - 1] <= 20.) {
/* Computing 2nd power */
	    d__1 = probdata_1.t[i__ - 1];
	    z__ = x[3] + x[4] * probdata_1.t[i__ - 1] + x[5] * (d__1 * d__1);
	} else {
/* Computing 2nd power */
	    d__1 = probdata_1.t[i__ - 1];
/* Computing 3rd power */
	    d__2 = probdata_1.t[i__ - 1];
	    z__ = x[6] + x[7] * probdata_1.t[i__ - 1] + x[8] * (d__1 * d__1) 
		    + x[9] * (d__2 * (d__2 * d__2));
	}
	s_wsle(&io___10);
	do_lio(&c__5, &c__1, (char *)&probdata_1.t[i__ - 1], (ftnlen)sizeof(
		doublereal));
	do_lio(&c__5, &c__1, (char *)&probdata_1.y[i__ - 1], (ftnlen)sizeof(
		doublereal));
	do_lio(&c__5, &c__1, (char *)&z__, (ftnlen)sizeof(doublereal));
	e_wsle();
    }
    tt = 1.;
L10:
    if (tt <= 30.01) {
	if (tt <= 10.) {
	    z__ = x[1] + x[2] * tt;
	} else if (tt <= 20.) {
/* Computing 2nd power */
	    d__1 = tt;
	    z__ = x[3] + x[4] * tt + x[5] * (d__1 * d__1);
	} else {
/* Computing 2nd power */
	    d__1 = tt;
/* Computing 3rd power */
	    d__2 = tt;
	    z__ = x[6] + x[7] * tt + x[8] * (d__1 * d__1) + x[9] * (d__2 * (
		    d__2 * d__2));
	}
	s_wsle(&io___12);
	do_lio(&c__5, &c__1, (char *)&tt, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&z__, (ftnlen)sizeof(doublereal));
	e_wsle();
	tt += .01;
	goto L10;
    }
} /* pedazos4_endp__ */

