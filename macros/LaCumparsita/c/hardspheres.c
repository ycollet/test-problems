/* hardspheres.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer pair[999000]	/* was [2][499500] */, ndim, np;
} probdata_;

#define probdata_1 probdata_

/*     ================================================================= */
/*     File: hardspheres.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     November 23, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Hard spheres problem */
/*     -------------------- */
/*     This problem is associated to the family of Hard-Spheres */
/*     problem. It belongs to a family of sphere packing */
/*     problems, a class of challenging problems dating from the */
/*     beginning of the seventeenth century. In the tradition of famous */
/*     problems in mathematics, the statements of these problems are */
/*     elusively simple, and have withstood the attacks of many worthy */
/*     mathematicians (e.g. Newton, Hilbert, Gregory), while most of its */
/*     instances remain open problems. Furthermore, it is related to */
/*     practical problems in chemistry, biology and physics. */
/*     The Hard-Spheres Problem is to maximize the minimum pairwise */
/*     distance between np points on a sphere in R^ndim. This problem */
/*     may be reduced to a nonlinear optimization problem that turns */
/*     out, as might be expected from the mentioned history, to be a */
/*     particularly hard, nonconvex problem, with a potentially large */
/*     number of (nonoptimal) points satisfying KKT conditions. We have */
/*     thus a class of problems indexed by the parameters ndim and np, */
/*     that provides a suitable set of test problems for evaluating */
/*     Nonlinear Programming codes. */
/*     After some algebric manipulations, we can formulate this problem */
/*     as */

/*               Minimize z */

/*               subject to */

/*               z \geq <x_i,x_j> for all different pair of indices i,j */

/*               ||x_i||^2 = 1    for all i = 1,...,NP */

/*     The goal is to find an objective value less than 0.5 (This means */
/*     that the NP points stored belong to the sphere and every distance */
/*     between two of them is greater than 1.0). */

/*     Obs: the starting point is aleatorally chosen although each */
/*     variable belongs to [-1,1]. */

/*     References: */

/*     [1] "Validation of an Augmented Lagrangian algorithm with a */
/*         Gauss-Newton Hessian approximation using a set of */
/*         Hard-Spheres problems", N. Krejic, J. M. Martinez, M. Mello */
/*         and E. A. Pilotta, Tech. Report RP 29/98, IMECC-UNICAMP, */
/*         Campinas, 1998. */

/*     [2] "Inexact-Restoration Algorithm for Constrained Optimization", */
/*         J. M. Martinez and E. A. Pilotta, Tech. Report, */
/*         IMECC-UNICAMP, Campinas, 1998. */

/*     [3] "Sphere Packings, Lattices and Groups", J. H. Conway and */
/*         N. J. C. Sloane, Springer-Verlag, NY, 1988. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hardspheres_inip__(n, x, l, u, m, lambda, rho, equatn, 
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
/*     Dimension of the space (ndim) and number of points (np) */
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
    probdata_1.np = *np_in__;
/*     Seed for the initial point random generation */
    seed = *seed_in__;
/*     Set pairs */
    k = 0;
    i__1 = probdata_1.np;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = probdata_1.np;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ++k;
	    probdata_1.pair[(k << 1) - 2] = i__;
	    probdata_1.pair[(k << 1) - 1] = j;
	}
    }
/*     Number of variables */
    *n = probdata_1.ndim * probdata_1.np + 1;
/*     Initial point */
/*     seed = 120927.0d0 */
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = drand_(&seed) * 2. - 1.;
    }
    x[*n] = drand_(&seed);
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e20;
	u[i__] = 1e20;
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = probdata_1.np + probdata_1.np * (probdata_1.np - 1) / 2;
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
    i__1 = probdata_1.np;
    for (i__ = 1; i__ <= i__1; ++i__) {
	equatn[i__] = TRUE_;
    }
    i__1 = probdata_1.np + probdata_1.np * (probdata_1.np - 1) / 2;
    for (i__ = probdata_1.np + 1; i__ <= i__1; ++i__) {
	equatn[i__] = FALSE_;
    }
/*     For each constraint i, set linear(i) = .true. if it is a linear */
/*     constraint, otherwise set linear(i) = .false. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	linear[i__] = FALSE_;
    }
} /* hardspheres_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hardspheres_evalf__(n, x, f, flag__)
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
    *f = x[*n];
} /* hardspheres_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hardspheres_evalg__(n, x, g, flag__)
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
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = 0.;
    }
    g[*n] = 1.;
} /* hardspheres_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hardspheres_evalh__(n, x, hlin, hcol, hval, hnnz, flag__)
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
} /* hardspheres_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hardspheres_evalc__(n, x, ind, c__, flag__)
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
    if (0 < *ind && *ind <= probdata_1.np) {
	*c__ = -1.;
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = x[probdata_1.ndim * (*ind - 1) + i__];
	    *c__ += d__1 * d__1;
	}
	return 0;
    }
    if (probdata_1.np < *ind && *ind <= probdata_1.np + probdata_1.np * (
	    probdata_1.np - 1) / 2) {
	i1 = probdata_1.pair[(*ind - probdata_1.np << 1) - 2];
	i2 = probdata_1.pair[(*ind - probdata_1.np << 1) - 1];
	*c__ = -x[*n];
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    *c__ += x[probdata_1.ndim * (i1 - 1) + i__] * x[probdata_1.ndim * 
		    (i2 - 1) + i__];
	}
	return 0;
    }
} /* hardspheres_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hardspheres_evaljac__(n, x, ind, jcvar, jcval, jcnnz, 
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
    static integer i__, i1, i2;

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
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --jcval;
    --jcvar;
    --x;

    /* Function Body */
    *flag__ = 0;
    if (0 < *ind && *ind <= probdata_1.np) {
	*jcnnz = probdata_1.ndim;
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jcvar[i__] = probdata_1.ndim * (*ind - 1) + i__;
	    jcval[i__] = x[probdata_1.ndim * (*ind - 1) + i__] * 2.;
	}
	return 0;
    }
    if (probdata_1.np < *ind && *ind <= probdata_1.np + probdata_1.np * (
	    probdata_1.np - 1) / 2) {
	i1 = probdata_1.pair[(*ind - probdata_1.np << 1) - 2];
	i2 = probdata_1.pair[(*ind - probdata_1.np << 1) - 1];
	*jcnnz = (probdata_1.ndim << 1) + 1;
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jcvar[i__] = probdata_1.ndim * (i1 - 1) + i__;
	    jcval[i__] = x[probdata_1.ndim * (i2 - 1) + i__];
	    jcvar[probdata_1.ndim + i__] = probdata_1.ndim * (i2 - 1) + i__;
	    jcval[probdata_1.ndim + i__] = x[probdata_1.ndim * (i1 - 1) + i__]
		    ;
	}
	jcvar[*jcnnz] = *n;
	jcval[*jcnnz] = -1.;
	return 0;
    }
} /* hardspheres_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hardspheres_evalhc__(n, x, ind, hclin, hccol, hcval, 
	hcnnz, flag__)
integer *n;
doublereal *x;
integer *ind, *hclin, *hccol;
doublereal *hcval;
integer *hcnnz, *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, i1, i2, tmp;

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
    if (1 <= *ind && *ind <= probdata_1.np) {
	*hcnnz = probdata_1.ndim;
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    hclin[i__] = probdata_1.ndim * (*ind - 1) + i__;
	    hccol[i__] = probdata_1.ndim * (*ind - 1) + i__;
	    hcval[i__] = 2.;
	}
	return 0;
    }
    if (*ind >= probdata_1.np + 1 && *ind <= probdata_1.np + probdata_1.np * (
	    probdata_1.np - 1) / 2) {
	i1 = probdata_1.pair[(*ind - probdata_1.np << 1) - 2];
	i2 = probdata_1.pair[(*ind - probdata_1.np << 1) - 1];
	if (i1 < i2) {
	    tmp = i1;
	    i1 = i2;
	    i2 = tmp;
	}
	*hcnnz = probdata_1.ndim;
	i__1 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    hclin[i__] = probdata_1.ndim * (i1 - 1) + i__;
	    hccol[i__] = probdata_1.ndim * (i2 - 1) + i__;
	    hcval[i__] = 1.;
	}
	return 0;
    }
} /* hardspheres_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hardspheres_evalhlp__(n, x, m, lambda, p, hp, goth, 
	flag__)
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
} /* hardspheres_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hardspheres_endp__(n, x, l, u, m, lambda, rho, equatn, 
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
} /* hardspheres_endp__ */

