/* kissing2.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer pair[999000]	/* was [2][499500] */, ndim, npun;
} probdata_;

#define probdata_1 probdata_

/*     ================================================================= */
/*     File: kissing2.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     October 30th, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Kissing Number problem */
/*     ---------------------- */
/*     Given npun and ndim, try to find npun points p_1, p_2, ..., */
/*     p_npun in the unitary sphere in R^ndim such that */

/*          dist(pi,pj) >= 1, for all i < j. */

/*     The regular formulation of this problem would be */

/*     || pi ||^2 = 1, for all i, and dist(pi,pj)^2 >= 1, for all i < j. */

/*     In this formulation, using that */

/*     dist(pi,pj)^2 = || pi ||^2 + || pj ||^2 - 2 < pi,pj > */

/*     equations are */

/*     || pi ||^2 = 1, for all i, and < pi,pj > <= 0.5, for all i < j. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int kissing2_inip__(n, x, l, u, m, lambda, rho, equatn, 
	linear, nd_in__, np_in__, seed_in__)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
integer *nd_in__, *np_in__;
doublereal *seed_in__;
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
/*     Dimension of the space (ndim) and number of points (npun) */
    /* Parameter adjustments */
    --linear;
    --equatn;
    --rho;
    --lambda;
    --u;
    --l;
    --x;

    /* Function Body */
    probdata_1.ndim = *nd_in__;
    probdata_1.npun = *np_in__;
/*     Seed for the initial point random generation */
    seed = *seed_in__;
/*     Set pairs */
    k = 0;
    i__1 = probdata_1.npun;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = probdata_1.npun;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ++k;
	    probdata_1.pair[(k << 1) - 2] = i__;
	    probdata_1.pair[(k << 1) - 1] = j;
	}
    }
/*     Number of variables */
    *n = probdata_1.ndim * probdata_1.npun;
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e20;
	u[i__] = 1e20;
    }
/*     Initial point */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = drand_(&seed) * 2. - 1.;
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = probdata_1.npun + probdata_1.npun * (probdata_1.npun - 1) / 2;
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
	rho[i__] = 1.;
    }
/*     For each constraint i, set equatn(i) = .true. if it is an equality */
/*     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if */
/*     it is an inequality constraint of the form c_i(x) <= 0. */
    i__1 = probdata_1.npun;
    for (i__ = 1; i__ <= i__1; ++i__) {
	equatn[i__] = TRUE_;
    }
    i__1 = probdata_1.npun + 1;
    i__2 = probdata_1.npun + probdata_1.npun * (probdata_1.npun - 1) / 2;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	equatn[i__] = FALSE_;
    }
/*     For each constraint i, set linear(i) = .true. if it is a linear */
/*     constraint, otherwise set linear(i) = .false. */
    i__2 = *m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	linear[i__] = FALSE_;
    }
} /* kissing2_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int kissing2_evalf__(n, x, f, flag__)
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
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    *f = 0.;
} /* kissing2_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int kissing2_evalg__(n, x, g, flag__)
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
/*     LOCAL SCALARS */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    *flag__ = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = 0.;
    }
} /* kissing2_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int kissing2_evalh__(n, x, hlin, hcol, hval, hnnz, flag__)
integer *n;
doublereal *x;
integer *hlin, *hcol;
doublereal *hval;
integer *hnnz, *flag__;
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
    /* Parameter adjustments */
    --x;
    --hlin;
    --hcol;
    --hval;

    /* Function Body */
    *flag__ = 0;
    *hnnz = 0;
} /* kissing2_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int kissing2_evalc__(n, x, ind, c__, flag__)
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
    static integer i__, i1, i2, ind1, ind2;

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
    if (1 <= *ind && *ind <= probdata_1.npun) {
	*c__ = -1.;
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = x[probdata_1.ndim * (*ind - 1) + i__];
	    *c__ += d__1 * d__1;
	}
    } else if (probdata_1.npun + 1 <= *ind && *ind <= probdata_1.npun + 
	    probdata_1.npun * (probdata_1.npun - 1) / 2) {
	i1 = probdata_1.pair[(*ind - probdata_1.npun << 1) - 2];
	i2 = probdata_1.pair[(*ind - probdata_1.npun << 1) - 1];
	*c__ = -.5;
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ind1 = probdata_1.ndim * (i1 - 1) + i__;
	    ind2 = probdata_1.ndim * (i2 - 1) + i__;
	    *c__ += x[ind1] * x[ind2];
	}
    }
} /* kissing2_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int kissing2_evaljac__(n, x, ind, jcvar, jcval, jcnnz, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *jcvar;
doublereal *jcval;
integer *jcnnz, *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, i1, i2, ind1, ind2;

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

/*     jcnnz   integer, */
/*              number of perhaps-non-null elements of the computed */
/*              gradient, */

/*     jcvar   integer jcvar(jcnnz), */
/*              see below, */

/*     jcval   double precision jcval(jcnnz), */
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
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --jcval;
    --jcvar;
    --x;

    /* Function Body */
    *flag__ = 0;
    if (1 <= *ind && *ind <= probdata_1.npun) {
	*jcnnz = probdata_1.ndim;
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jcvar[i__] = probdata_1.ndim * (*ind - 1) + i__;
	    jcval[i__] = x[probdata_1.ndim * (*ind - 1) + i__] * 2.;
	}
    } else if (probdata_1.npun + 1 <= *ind && *ind <= probdata_1.npun + 
	    probdata_1.npun * (probdata_1.npun - 1) / 2) {
	i1 = probdata_1.pair[(*ind - probdata_1.npun << 1) - 2];
	i2 = probdata_1.pair[(*ind - probdata_1.npun << 1) - 1];
	*jcnnz = probdata_1.ndim << 1;
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ind1 = probdata_1.ndim * (i1 - 1) + i__;
	    ind2 = probdata_1.ndim * (i2 - 1) + i__;
	    jcvar[i__] = ind1;
	    jcval[i__] = x[ind2];
	    jcvar[probdata_1.ndim + i__] = ind2;
	    jcval[probdata_1.ndim + i__] = x[ind1];
	}
    }
} /* kissing2_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int kissing2_evalhc__(n, x, ind, hclin, hccol, hcval, hcnnz, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *hclin, *hccol;
doublereal *hcval;
integer *hcnnz, *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, i1, i2, ind1, ind2;

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
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;
    --hclin;
    --hccol;
    --hcval;

    /* Function Body */
    *flag__ = 0;
    if (1 <= *ind && *ind <= probdata_1.npun) {
	*hcnnz = probdata_1.ndim;
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    hclin[i__] = probdata_1.ndim * (*ind - 1) + i__;
	    hccol[i__] = probdata_1.ndim * (*ind - 1) + i__;
	    hcval[i__] = 2.;
	}
    } else if (probdata_1.npun + 1 <= *ind && *ind <= probdata_1.npun + 
	    probdata_1.npun * (probdata_1.npun - 1) / 2) {
	i1 = probdata_1.pair[(*ind - probdata_1.npun << 1) - 2];
	i2 = probdata_1.pair[(*ind - probdata_1.npun << 1) - 1];
	*hcnnz = probdata_1.ndim;
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ind1 = probdata_1.ndim * (i1 - 1) + i__;
	    ind2 = probdata_1.ndim * (i2 - 1) + i__;
	    hclin[i__] = ind2;
	    hccol[i__] = ind1;
	    hcval[i__] = 1.;
	}
    }
} /* kissing2_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int kissing2_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
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
} /* kissing2_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int kissing2_endp__(n, x, l, u, m, lambda, rho, equatn, 
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
} /* kissing2_endp__ */

