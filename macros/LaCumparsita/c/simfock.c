/* simfock.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer kdim, ndim;
} probdata_;

#define probdata_1 probdata_

/*     ================================================================= */
/*     File: simfock.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     May 12, 2005. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Simulation of Hartree-Fock problem */
/*     ---------------------------------- */
/*     Given kdim, ndim in 1, 2, 3, ..., ndim <= kdim, n = kdim * kdim, */
/*     x in R^n, we define mat(x) in R^{kdim x kdim} as the matrix whose */
/*     columnwise representation is x. The problem is: */

/*                  Minimize 0.5 x^T A x + b^T x */

/*                  subject to */

/*                          mat(x) mat(x) = mat(x) */

/*                          trace[mat(x)] = ndim */

/*                          mat(x) = mat(x)^T */

/*     Therefore, the problem has kdim^2 variables, kdim^2 nonlinear */
/*     equality constraints and kdim * ( kdim - 1 ) / 2 + 1 linear */
/*     equality constraints. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock_inip__(n, x, l, u, m, lambda, rho, equatn, 
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
    static integer i__, j;

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
/*     Number of variables */
/* Computing 2nd power */
    i__1 = probdata_1.kdim;
    *n = i__1 * i__1;
/*     Initial point */
    i__1 = probdata_1.kdim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = probdata_1.kdim;
	for (j = 1; j <= i__2; ++j) {
	    x[(j - 1) * probdata_1.kdim + i__] = 1.;
	}
    }
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e5;
	u[i__] = 1e5;
    }
/*     Number of constraints (equalities plus inequalities) */
/* Computing 2nd power */
    i__1 = probdata_1.kdim;
    *m = (i__1 * i__1 << 1) + 1;
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
} /* simfock_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal g[10000000];
    static integer i__, j;
    extern /* Subroutine */ int simfock_gedede__();

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
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    *f = 0.;
    i__1 = probdata_1.kdim - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	*f = *f - x[(i__ - 2) * probdata_1.kdim + i__] + x[(i__ - 1) * 
		probdata_1.kdim + i__] * 2. - x[i__ * probdata_1.kdim + i__];
    }
    *f = *f + x[1] * 2. - x[probdata_1.kdim + 1] - x[(probdata_1.kdim - 2) * 
	    probdata_1.kdim + probdata_1.kdim] + x[(probdata_1.kdim - 1) * 
	    probdata_1.kdim + probdata_1.kdim] * 2.;
    *f *= 2.;
    simfock_gedede__(n, &x[1], g);
    i__1 = probdata_1.kdim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = probdata_1.kdim;
	for (j = 1; j <= i__2; ++j) {
	    *f += g[(i__ - 1) * probdata_1.kdim + j - 1] * x[(i__ - 1) * 
		    probdata_1.kdim + j];
	}
    }
} /* simfock_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock_evalg__(n, x, g, flag__)
integer *n;
doublereal *x, *g;
integer *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int simfock_gedede__();

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
    simfock_gedede__(n, &x[1], &g[1]);
    i__1 = probdata_1.kdim - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	g[(i__ - 2) * probdata_1.kdim + i__] += -1.;
	g[(i__ - 1) * probdata_1.kdim + i__] += 2.;
	g[i__ * probdata_1.kdim + i__] += -1.;
    }
    g[1] += 2.;
    g[(probdata_1.kdim - 2) * probdata_1.kdim + probdata_1.kdim] += -1.;
    g[(probdata_1.kdim - 1) * probdata_1.kdim + probdata_1.kdim] += 2.;
    g[probdata_1.kdim + 1] += -1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] *= 2.;
    }
} /* simfock_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock_evalh__(n, x, hlin, hcol, hval, nnzh, flag__)
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
} /* simfock_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock_evalc__(n, x, ind, c__, flag__)
integer *n;
doublereal *x;
integer *ind;
doublereal *c__;
integer *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, iminus;

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
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
/* Computing 2nd power */
    i__1 = probdata_1.kdim;
    if (*ind <= i__1 * i__1) {
	j = *ind % probdata_1.kdim;
	if (j == 0) {
	    j = probdata_1.kdim;
	}
	i__ = (*ind - j) / probdata_1.kdim + 1;
	*c__ = -x[(j - 1) * probdata_1.kdim + i__];
	i__1 = probdata_1.kdim;
	for (k = 1; k <= i__1; ++k) {
	    *c__ += x[(k - 1) * probdata_1.kdim + i__] * x[(j - 1) * 
		    probdata_1.kdim + k];
	}
	return 0;
    }
/* Computing 2nd power */
    i__1 = probdata_1.kdim;
    if (*ind <= i__1 * i__1 << 1) {
/* Computing 2nd power */
	i__1 = probdata_1.kdim;
	iminus = *ind - i__1 * i__1;
	j = iminus % probdata_1.kdim;
	if (j == 0) {
	    j = probdata_1.kdim;
	}
	i__ = (iminus - j) / probdata_1.kdim + 1;
	*c__ = x[(j - 1) * probdata_1.kdim + i__] - x[(i__ - 1) * 
		probdata_1.kdim + j];
	return 0;
    }
    *c__ = -((doublereal) probdata_1.ndim);
    i__1 = probdata_1.kdim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*c__ += x[(i__ - 1) * probdata_1.kdim + i__];
    }
} /* simfock_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock_evaljac__(n, x, ind, indjac, valjac, nnzjac, 
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
    static integer i__, j, k, iminus;

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
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --valjac;
    --indjac;
    --x;

    /* Function Body */
    *flag__ = 0;
/* Computing 2nd power */
    i__1 = probdata_1.kdim;
    if (*ind <= i__1 * i__1) {
	*nnzjac = (probdata_1.kdim << 1) + 1;
	j = *ind % probdata_1.kdim;
	if (j == 0) {
	    j = probdata_1.kdim;
	}
	i__ = (*ind - j) / probdata_1.kdim + 1;
	i__1 = probdata_1.kdim;
	for (k = 1; k <= i__1; ++k) {
	    indjac[k] = (j - 1) * probdata_1.kdim + k;
	    valjac[k] = x[(k - 1) * probdata_1.kdim + i__];
	    indjac[probdata_1.kdim + k] = (k - 1) * probdata_1.kdim + i__;
	    valjac[probdata_1.kdim + k] = x[(j - 1) * probdata_1.kdim + k];
	}
	indjac[(probdata_1.kdim << 1) + 1] = (j - 1) * probdata_1.kdim + i__;
	valjac[(probdata_1.kdim << 1) + 1] = -1.;
	return 0;
    }
/* Computing 2nd power */
    i__1 = probdata_1.kdim;
    if (*ind <= i__1 * i__1 << 1) {
	*nnzjac = 2;
/* Computing 2nd power */
	i__1 = probdata_1.kdim;
	iminus = *ind - i__1 * i__1;
	j = iminus % probdata_1.kdim;
	if (j == 0) {
	    j = probdata_1.kdim;
	}
	i__ = (iminus - j) / probdata_1.kdim + 1;
	indjac[1] = (j - 1) * probdata_1.kdim + i__;
	valjac[1] = 1.;
	indjac[2] = (i__ - 1) * probdata_1.kdim + j;
	valjac[2] = -1.;
	return 0;
    }
    *nnzjac = probdata_1.kdim;
    i__1 = probdata_1.kdim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	indjac[i__] = (i__ - 1) * probdata_1.kdim + i__;
	valjac[i__] = 1.;
    }
} /* simfock_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock_gedede__(n, x, g)
integer *n;
doublereal *x, *g;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, jj;

/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = x[i__];
	for (jj = 1; jj <= 10; ++jj) {
	    j = i__ + jj;
	    if (j <= *n) {
		g[i__] += 1. / (doublereal) (i__ + j) * x[j];
	    }
	    j = i__ - jj;
	    if (j >= 1) {
		g[i__] += 1. / (doublereal) (i__ + j) * x[j];
	    }
	}
    }
} /* simfock_gedede__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock_evalhc__(n, x, ind, hclin, hccol, hcval, nnzhc, 
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
} /* simfock_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
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
} /* simfock_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int simfock_endp__(n, x, l, u, m, lambda, rho, equatn, 
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
} /* simfock_endp__ */

