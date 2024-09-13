/* packcrmn-feas.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer pair[999000]	/* was [2][499500] */, nite, ndim;
    doublereal iterad, objdim[3];
} packdata_;

#define packdata_1 packdata_

/*     ================================================================= */
/*     File: packcrmn-feas.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     November 15, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Packing fixed-dimension circular items within a fixed-dimension */
/*     --------------------------------------------------------------- */
/*     rectangular object maximizing the number of items */
/*     ------------------------------------------------- */

/*     FEASIBILITY PROBLEM VERSION */
/*     --------------------------- */

/*     We wish to fit k circles of radii ri (i=1, ..., k) into a */
/*     rectangular object with dimensions L and W in such a way that the */
/*     circles do not overlapped. Therefore, given k, the radii ri */
/*     (i=1, ..., k), L and W, the goal is to solve the feasibility */
/*     problem: */

/*     d( pi, pj )^2 >= ( ri + rj )^2, for all i < j, */

/*     ri <= [pi]_1 <= L - ri and ri <= [pi]_2 <= W - ri, for all i, */

/*     where d(.,.) is the Euclidian distance. */

/*     In fact, the problem is coded for any arbitrary dimension */
/*     (not only 2). */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packcrmn_feas_inip__(n, x, l, u, m, lambda, rho, equatn, 
	linear, nd_in__, nite_in__, iterad_in__, objdim_in__, seed_in__)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
integer *nd_in__, *nite_in__;
doublereal *iterad_in__, *objdim_in__, *seed_in__;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal seed;
    static integer i__, j, k;
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
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     SET PROBLEM DATA */
/*     Dimension of the space (ndim) and number of items (nite) */
    /* Parameter adjustments */
    --objdim_in__;
    --linear;
    --equatn;
    --rho;
    --lambda;
    --u;
    --l;
    --x;

    /* Function Body */
    packdata_1.ndim = *nd_in__;
    packdata_1.nite = *nite_in__;
/*     Radius of the circular items */
    packdata_1.iterad = *iterad_in__;
/*     Rectangular object dimensions */
    i__1 = packdata_1.ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	packdata_1.objdim[i__ - 1] = objdim_in__[i__];
    }
/*     Seed for the initial point random generation */
    seed = *seed_in__;
/*     Set pairs */
    k = 0;
    i__1 = packdata_1.nite;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = packdata_1.nite;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ++k;
	    packdata_1.pair[(k << 1) - 2] = i__;
	    packdata_1.pair[(k << 1) - 1] = j;
	}
    }
/*     Number of variables */
    *n = packdata_1.ndim * packdata_1.nite;
/*     Lower and upper bounds */
    i__1 = packdata_1.nite;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = packdata_1.ndim;
	for (j = 1; j <= i__2; ++j) {
	    ind = (i__ - 1) * packdata_1.ndim + j;
	    l[ind] = packdata_1.iterad;
	    u[ind] = packdata_1.objdim[j - 1] - packdata_1.iterad;
	}
    }
/*     Initial point */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = l[i__] + (u[i__] - l[i__]) * drand_(&seed);
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = packdata_1.nite * (packdata_1.nite - 1) / 2;
/*     Lagrange multipliers approximation */
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
} /* packcrmn_feas_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packcrmn_feas_evalf__(n, x, f, flag__)
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
} /* packcrmn_feas_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packcrmn_feas_evalg__(n, x, g, flag__)
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
} /* packcrmn_feas_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packcrmn_feas_evalh__(n, x, hlin, hcol, hval, hnnz, 
	flag__)
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
} /* packcrmn_feas_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packcrmn_feas_evalc__(n, x, ind, c__, flag__)
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
    i1 = packdata_1.pair[(*ind << 1) - 2];
    i2 = packdata_1.pair[(*ind << 1) - 1];
/* Computing 2nd power */
    d__1 = packdata_1.iterad;
    *c__ = d__1 * d__1 * 4.;
    i__1 = packdata_1.ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ind1 = packdata_1.ndim * (i1 - 1) + i__;
	ind2 = packdata_1.ndim * (i2 - 1) + i__;
/* Computing 2nd power */
	d__1 = x[ind1] - x[ind2];
	*c__ -= d__1 * d__1;
    }
} /* packcrmn_feas_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packcrmn_feas_evaljac__(n, x, ind, jcvar, jcval, jcnnz, 
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
    static doublereal tmp;
    static integer ind1, ind2;

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
    i1 = packdata_1.pair[(*ind << 1) - 2];
    i2 = packdata_1.pair[(*ind << 1) - 1];
    *jcnnz = packdata_1.ndim << 1;
    i__1 = packdata_1.ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ind1 = packdata_1.ndim * (i1 - 1) + i__;
	ind2 = packdata_1.ndim * (i2 - 1) + i__;
	tmp = (x[ind1] - x[ind2]) * 2.;
	jcvar[i__] = ind1;
	jcval[i__] = -tmp;
	jcvar[packdata_1.ndim + i__] = ind2;
	jcval[packdata_1.ndim + i__] = tmp;
    }
} /* packcrmn_feas_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packcrmn_feas_evalhc__(n, x, ind, hclin, hccol, hcval, 
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
    i1 = packdata_1.pair[(*ind << 1) - 2];
    i2 = packdata_1.pair[(*ind << 1) - 1];
    *hcnnz = packdata_1.ndim * 3;
    i__1 = packdata_1.ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ind1 = packdata_1.ndim * (i1 - 1) + i__;
	ind2 = packdata_1.ndim * (i2 - 1) + i__;
	hclin[i__] = ind1;
	hccol[i__] = ind1;
	hcval[i__] = -2.;
	hclin[packdata_1.ndim + i__] = ind2;
	hccol[packdata_1.ndim + i__] = ind2;
	hcval[packdata_1.ndim + i__] = -2.;
	hclin[(packdata_1.ndim << 1) + i__] = ind2;
	hccol[(packdata_1.ndim << 1) + i__] = ind1;
	hcval[(packdata_1.ndim << 1) + i__] = 2.;
    }
} /* packcrmn_feas_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packcrmn_feas_evalhlp__(n, x, m, lambda, p, hp, goth, 
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
} /* packcrmn_feas_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int packcrmn_feas_endp__(n, x, l, u, m, lambda, rho, equatn, 
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
} /* packcrmn_feas_endp__ */

