/* condor.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal uint[1000];
} objective_;

#define objective_1 objective_

struct {
    doublereal fon[1000], h__;
} constraints_;

#define constraints_1 constraints_

struct {
    doublereal usol[1000];
} solution_;

#define solution_1 solution_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;

/*     ================================================================= */
/*     File: condor.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     July 7th, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Condor problem */
/*     -------------- */
/*     We generate (using subroutine evalu) a function called usol as a */
/*     Fourier polynomial with up to q = 20 terms. A grid in [0, 2 pi] */
/*     with up to 1000 points is established. We generate function fi(t) */
/*     in such a way that usol is a solution of the discretization of */
/*     the differential equation u''(t) + [u'(t)]^2 + u(t) = fi(t). We */
/*     select up to 1000 observation points equally spaced in [-2, 18] */
/*     and observe the solution usol at these observation points, */
/*     obtaining the observed solution uobs. Finally, uint, the */
/*     interpolation of uobs in the grid points is computed. The */
/*     optimization problem is to obtain the solution of the */
/*     discretized differential equation which is closest, in the */
/*     least-squares sense, to uint. Instances with q = 20, 500 */
/*     observation points and 1000 grid points are challenging. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int condor_inip__(n, x, l, u, m, lambda, rho, equatn, linear,
	 n_in__, p_in__, q_in__)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
integer *n_in__, *p_in__, *q_in__;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double cos(), sin(), atan();

    /* Local variables */
    static doublereal uobs[1000];
    extern /* Subroutine */ int condor_evalu__();
    static doublereal a[20], b[20];
    static integer i__, j, p, q;
    static doublereal t[1000], z__, a0, u1, u2, pen, aux;

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
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
/*     ****************************************************************** */
/*     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR PROBLEM */
/*     DATA: */
/*     ****************************************************************** */
/*     Number of points in the grid (nmax=1000) */
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
/*     Number of observed point (pmax<=nmax) */
    p = *p_in__;
/*     Number of terms of the Fourier polynomial (qmax=20) */
    q = *q_in__;
/*     Fourier polynomial coefficients */
    i__1 = q;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__ - 1] = cos((doublereal) i__) / (doublereal) i__;
/* Computing 2nd power */
	d__1 = sin((doublereal) i__);
	b[i__ - 1] = d__1 * d__1 / (doublereal) i__;
    }
    a0 = 2.;
    z__ = 0.;
    i__1 = p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t[i__ - 1] = 20. / (doublereal) p * i__ - 2.;
    }
    i__1 = p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	condor_evalu__(&a0, a, b, &q, &t[i__ - 1], &uobs[i__ - 1], &u1, &u2, &
		z__);
    }
/*     write(*,*) 't = ',(t(i),i=1,p) */
/*     write(*,*) 'uobs = ',(uobs(i),i=1,p) */
/*     Independent term of the differential equation (fon) */
    constraints_1.h__ = atan(1.) * 8. / (doublereal) (*n - 1);
    i__1 = *n - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ = i__ * constraints_1.h__;
	condor_evalu__(&a0, a, b, &q, &z__, &aux, &u1, &u2, &
		constraints_1.fon[i__ - 1]);
    }
/*     Analytic solution at each grid point (usol) */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ = (i__ - 1) * constraints_1.h__;
	condor_evalu__(&a0, a, b, &q, &z__, &solution_1.usol[i__ - 1], &u1, &
		u2, &aux);
    }
/*     Interpolation */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ = (doublereal) (i__ - 1) * constraints_1.h__;
	if (z__ <= t[1]) {
	    pen = (uobs[1] - uobs[0]) / (t[1] - t[0]);
	    objective_1.uint[i__ - 1] = uobs[0] + pen * (z__ - t[0]);
	} else if (z__ >= t[p - 2]) {
	    pen = (uobs[p - 1] - uobs[p - 2]) / (t[p - 1] - t[p - 2]);
	    objective_1.uint[i__ - 1] = uobs[p - 2] + pen * (z__ - t[p - 2]);
	} else {
	    i__2 = p - 2;
	    for (j = 2; j <= i__2; ++j) {
		if (z__ >= t[j - 1] && z__ <= t[j]) {
		    pen = (uobs[j] - uobs[j - 1]) / (t[j] - t[j - 1]);
		    objective_1.uint[i__ - 1] = uobs[j - 1] + pen * (z__ - t[
			    j - 1]);
		}
	    }
	}
    }
/*     write(*,*) 'Interpolated points: ' */
/*     do i = 1,n */
/*         write(*,*) i,(i-1)*h,uint(i) */
/*     end do */
/*     Initial point */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = objective_1.uint[i__ - 1];
    }
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e20;
	u[i__] = 1e20;
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = *n - 2;
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
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE INIP. */
/*     ****************************************************************** */
} /* condor_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int condor_evalu__(a0, a, b, m, x, u, u1, u2, fi)
doublereal *a0, *a, *b;
integer *m;
doublereal *x, *u, *u1, *u2, *fi;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double cos(), sin();

    /* Local variables */
    static integer k;

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    *u = *a0;
    *u1 = 0.;
    *u2 = 0.;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	*u = *u + a[k] * cos(k * *x) + b[k] * sin(k * *x);
	*u1 = *u1 - a[k] * k * sin(k * *x) + b[k] * k * cos(k * *x);
	*u2 = *u2 - a[k] * k * k * cos(k * *x) - b[k] * k * k * sin(k * *x);
    }
/* Computing 2nd power */
    d__1 = *u1;
    *fi = *u2 + d__1 * d__1 + *u;
} /* condor_evalu__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int condor_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;

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
/*     COMMON ARRAYS */
/*     LOCAL SCALAR */
/*     COMMON BLOCKS */
/*     Objective function */
/*     ****************************************************************** */
/*     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR OBJECTIVE */
/*     FUNCTION: */
/*     ****************************************************************** */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    *f = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = x[i__] - objective_1.uint[i__ - 1];
	*f += d__1 * d__1;
    }
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALF. */
/*     ****************************************************************** */
} /* condor_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int condor_evalg__(n, x, g, flag__)
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
/*     COMMON ARRAYS */
/*     LOCAL SCALAR */
/*     COMMON BLOCKS */
/*     Gradient vector */
/*     ****************************************************************** */
/*     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENT */
/*     VECTOR OF YOUR OBJECTIVE FUNCTION: */
/*     ****************************************************************** */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    *flag__ = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = (x[i__] - objective_1.uint[i__ - 1]) * 2.;
    }
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALG. */
/*     ****************************************************************** */
} /* condor_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int condor_evalh__(n, x, hlin, hcol, hval, nnzh, flag__)
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
/*     ****************************************************************** */
/*     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE */
/*     HESSIAN MATRIX OF YOUR OBJECTIVE FUNCTION: */
/*     ****************************************************************** */
    /* Parameter adjustments */
    --x;
    --hlin;
    --hcol;
    --hval;

    /* Function Body */
    *flag__ = -1;
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALH. */
/*     ****************************************************************** */
} /* condor_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int condor_evalc__(n, x, ind, c__, flag__)
integer *n;
doublereal *x;
integer *ind;
doublereal *c__;
integer *flag__;
{
    /* System generated locals */
    doublereal d__1, d__2;

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
/*     COMMON ARRAYS */
/*     COMMON BLOCKS */
/*     Constraints */
/*     ****************************************************************** */
/*     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR */
/*     CONSTRAINTS: */
/*     ****************************************************************** */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
/* Computing 2nd power */
    d__1 = constraints_1.h__;
/* Computing 2nd power */
    d__2 = (x[*ind + 2] - x[*ind]) / (constraints_1.h__ * 2.);
    *c__ = (x[*ind + 2] - x[*ind + 1] * 2. + x[*ind]) / (d__1 * d__1) + d__2 *
	     d__2 + x[*ind + 1] - constraints_1.fon[*ind - 1];
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALC. */
/*     ****************************************************************** */
} /* condor_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int condor_evaljac__(n, x, ind, indjac, valjac, nnzjac, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *indjac;
doublereal *valjac;
integer *nnzjac, *flag__;
{
    /* System generated locals */
    doublereal d__1, d__2;

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
/*     COMMON ARRAYS */
/*     COMMON BLOCKS */
/*     Sparse gradient vector of the ind-th constraint */
/*     ****************************************************************** */
/*     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENTS */
/*     OF YOUR CONSTRAINTS: */
/*     ****************************************************************** */
    /* Parameter adjustments */
    --valjac;
    --indjac;
    --x;

    /* Function Body */
    *flag__ = 0;
    *nnzjac = 3;
    indjac[1] = *ind;
/* Computing 2nd power */
    d__1 = constraints_1.h__;
/* Computing 2nd power */
    d__2 = constraints_1.h__ * 2.;
    valjac[1] = 1. / (d__1 * d__1) - (x[*ind + 2] - x[*ind]) * 2. / (d__2 * 
	    d__2);
    indjac[2] = *ind + 1;
/* Computing 2nd power */
    d__1 = constraints_1.h__;
    valjac[2] = -2. / (d__1 * d__1) + 1.;
    indjac[3] = *ind + 2;
/* Computing 2nd power */
    d__1 = constraints_1.h__;
/* Computing 2nd power */
    d__2 = constraints_1.h__ * 2.;
    valjac[3] = 1. / (d__1 * d__1) + (x[*ind + 2] - x[*ind]) * 2. / (d__2 * 
	    d__2);
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALJAC. */
/*     ****************************************************************** */
} /* condor_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int condor_evalhc__(n, x, ind, hclin, hccol, hcval, nnzhc, 
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
/*     ****************************************************************** */
/*     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE */
/*     HESSIANS OF YOUR CONSTRAINTS: */
/*     ****************************************************************** */
    /* Parameter adjustments */
    --x;
    --hclin;
    --hccol;
    --hcval;

    /* Function Body */
    *flag__ = -1;
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALHC. */
/*     ****************************************************************** */
} /* condor_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int condor_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
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
/*     ****************************************************************** */
/*     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO CUMPUTE */
/*     THE PRODUCT OF THE HESSIANS OF THE LAGRANGIAN TIME VECTOR P: */
/*     ****************************************************************** */
    /* Parameter adjustments */
    --hp;
    --p;
    --x;
    --lambda;

    /* Function Body */
    *flag__ = -1;
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALHC. */
/*     ****************************************************************** */
} /* condor_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int condor_endp__(n, x, l, u, m, lambda, rho, equatn, linear)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsle(), e_wsle(), do_lio();

    /* Local variables */
    static doublereal zmax;
    static integer i__;

    /* Fortran I/O blocks */
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };


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
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     ****************************************************************** */
/*     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO COMPLEMENT THE */
/*     INFORMATION RELATED TO THE SOLUTION GIVEN BY THE SOLVER */
/*     ****************************************************************** */
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
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = zmax, d__3 = (d__1 = x[i__] - solution_1.usol[i__ - 1], abs(
		d__1));
	zmax = max(d__2,d__3);
    }
    s_wsle(&io___20);
    e_wsle();
    s_wsle(&io___21);
    do_lio(&c__9, &c__1, "Comparison between analytic and computed solution: "
	    , (ftnlen)51);
    e_wsle();
    s_wsle(&io___22);
    do_lio(&c__9, &c__1, "Analytic Computed Difference", (ftnlen)28);
    e_wsle();
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsle(&io___23);
	do_lio(&c__5, &c__1, (char *)&solution_1.usol[i__ - 1], (ftnlen)
		sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	d__1 = solution_1.usol[i__ - 1] - x[i__];
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
    }
    s_wsle(&io___24);
    do_lio(&c__9, &c__1, "Maximum error: ", (ftnlen)15);
    do_lio(&c__5, &c__1, (char *)&zmax, (ftnlen)sizeof(doublereal));
    e_wsle();
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE ENDP */
/*     ****************************************************************** */
} /* condor_endp__ */

