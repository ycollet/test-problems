/* ellipsoid.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal p[200000]	/* was [10][20000] */;
    integer e[100]	/* was [10][10] */, ndim, npun;
} probdata_;

#define probdata_1 probdata_

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b16 = 1.;
static doublereal c_b22 = 0.;

/*     ================================================================= */
/*     File: ellipsoid.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     November 28, 2006 (second derivatives added). */
/*     Previous: October 18, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Smallest enclosing ellipsoid problem */
/*     ------------------------------------ */
/*     Description: */
/*     Given npun fixed points pi in R^ndim, minimize the volume of the */
/*     ellipsoid, centered at the origin, that fits them. The model can */
/*     be written as: */

/*                  Minimize - \sum_{i=1,ndim} log ( l_{ii} ) */

/*                  subject to */

/*                          pi^T L L^T pi <= 1, i = 1, ..., npun, */

/*                          l_{ii} >= 0, i = 1, ..., ndim, */

/*     where L is ndim x ndim lower-triangular matrix. */

/*     Interesting cases: ndim = 3, 4, 10; npun = 10,000, 20,000. Points */
/*     pi generated using the Cauchy distribution. */

/*     The fixed npun data points are generated using Cauchy distribution. */

/*     References: */

/*     M.J. Todd and E.A. Yildirim, ``On Khachiyan's algorithm for the */
/*     computation of minimum volume enclosing ellipsoids'' (9/05). */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int ellipsoid_inip__(n, x, l, u, m, lambda, rho, equatn, 
	linear, nd_in__, np_in__)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
integer *nd_in__, *np_in__;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal seed;
    static integer i__, j, k;
    extern doublereal drand_();
    extern /* Subroutine */ int ellipsoid_cauchyd__();

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
/*     Set problem data (directly related to the number of variables */
/*     and constraints) */
/*     Dimension of the space (ndim) maximum 10 */
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
/*     Number of points (npun) maximum 20,000 */
    probdata_1.npun = *np_in__;
/*     Set the mapping between the elements of the lower-triagular matrix */
/*     L and its column-wise representation at vector x. */
    k = 0;
    i__1 = probdata_1.ndim;
    for (j = 1; j <= i__1; ++j) {
	i__2 = probdata_1.ndim;
	for (i__ = j; i__ <= i__2; ++i__) {
	    ++k;
	    probdata_1.e[i__ + j * 10 - 11] = k;
	}
    }
/*     Generate the npun points in R^ndim using the Cauchy distribution */
    seed = 123456.;
    i__1 = probdata_1.npun;
    for (j = 1; j <= i__1; ++j) {
	i__2 = probdata_1.ndim;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ellipsoid_cauchyd__(&seed, &probdata_1.p[i__ + j * 10 - 11]);
	}
    }
/*     Or, read points from a file */
/*     do j = 1,npun */
/*         read(*,*) (p(i,j),i=1,ndim) */
/*     end do */
/*     Number of variables (elements of a lower-triangular ndim x ndim matrix) */
    *n = probdata_1.ndim * (probdata_1.ndim + 1) / 2;
/*     Lower and upper bounds */
    i__1 = probdata_1.ndim;
    for (j = 1; j <= i__1; ++j) {
	l[probdata_1.e[j + j * 10 - 11]] = 1e-16;
	u[probdata_1.e[j + j * 10 - 11]] = 1e20;
	i__2 = probdata_1.ndim;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    l[probdata_1.e[i__ + j * 10 - 11]] = -1e20;
	    u[probdata_1.e[i__ + j * 10 - 11]] = 1e20;
	}
    }
/*     Initial point */
    seed = 123456.;
    i__1 = probdata_1.ndim;
    for (j = 1; j <= i__1; ++j) {
	i__2 = probdata_1.ndim;
	for (i__ = j; i__ <= i__2; ++i__) {
	    x[probdata_1.e[i__ + j * 10 - 11]] = drand_(&seed);
	}
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = probdata_1.npun;
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
} /* ellipsoid_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int ellipsoid_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log();

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
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    *f = 0.;
    i__1 = probdata_1.ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*f -= log(x[probdata_1.e[i__ + i__ * 10 - 11]]);
    }
} /* ellipsoid_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int ellipsoid_evalg__(n, x, g, flag__)
integer *n;
doublereal *x, *g;
integer *flag__;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;

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
    i__1 = probdata_1.ndim;
    for (j = 1; j <= i__1; ++j) {
	g[probdata_1.e[j + j * 10 - 11]] = -1. / x[probdata_1.e[j + j * 10 - 
		11]];
	i__2 = probdata_1.ndim;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    g[probdata_1.e[i__ + j * 10 - 11]] = 0.;
	}
    }
} /* ellipsoid_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int ellipsoid_evalh__(n, x, hlin, hcol, hval, hnnz, flag__)
integer *n;
doublereal *x;
integer *hlin, *hcol;
doublereal *hval;
integer *hnnz, *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;

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
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;
    --hlin;
    --hcol;
    --hval;

    /* Function Body */
    *flag__ = 0;
    *hnnz = probdata_1.ndim;
    i__1 = probdata_1.ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	hlin[i__] = probdata_1.e[i__ + i__ * 10 - 11];
	hcol[i__] = probdata_1.e[i__ + i__ * 10 - 11];
/* Computing 2nd power */
	d__1 = x[probdata_1.e[i__ + i__ * 10 - 11]];
	hval[i__] = 1. / (d__1 * d__1);
    }
} /* ellipsoid_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int ellipsoid_evalc__(n, x, ind, c__, flag__)
integer *n;
doublereal *x;
integer *ind;
doublereal *c__;
integer *flag__;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal y[10];

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
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
/*     Compute y = L^T p_{ind} */
    i__1 = probdata_1.ndim;
    for (j = 1; j <= i__1; ++j) {
	y[j - 1] = 0.;
	i__2 = probdata_1.ndim;
	for (i__ = j; i__ <= i__2; ++i__) {
	    y[j - 1] += x[probdata_1.e[i__ + j * 10 - 11]] * probdata_1.p[i__ 
		    + *ind * 10 - 11];
	}
    }
/*     Compute c = y^t y - 1.0d0 */
    *c__ = -1.;
    i__1 = probdata_1.ndim;
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	d__1 = y[j - 1];
	*c__ += d__1 * d__1;
    }
} /* ellipsoid_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int ellipsoid_evaljac__(n, x, ind, jcvar, jcval, jcnnz, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *jcvar;
doublereal *jcval;
integer *jcnnz, *flag__;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal y[10];

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
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --jcval;
    --jcvar;
    --x;

    /* Function Body */
    *flag__ = 0;
    i__1 = probdata_1.ndim;
    for (j = 1; j <= i__1; ++j) {
	y[j - 1] = 0.;
	i__2 = probdata_1.ndim;
	for (i__ = j; i__ <= i__2; ++i__) {
	    y[j - 1] += x[probdata_1.e[i__ + j * 10 - 11]] * probdata_1.p[i__ 
		    + *ind * 10 - 11];
	}
    }
    *jcnnz = probdata_1.ndim * (probdata_1.ndim + 1) / 2;
    i__1 = probdata_1.ndim;
    for (j = 1; j <= i__1; ++j) {
	i__2 = probdata_1.ndim;
	for (i__ = j; i__ <= i__2; ++i__) {
	    jcvar[probdata_1.e[i__ + j * 10 - 11]] = probdata_1.e[i__ + j * 
		    10 - 11];
	    jcval[probdata_1.e[i__ + j * 10 - 11]] = y[j - 1] * 2. * 
		    probdata_1.p[i__ + *ind * 10 - 11];
	}
    }
} /* ellipsoid_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int ellipsoid_evalhc__(n, x, ind, hclin, hccol, hcval, hcnnz,
	 flag__)
integer *n;
doublereal *x;
integer *ind, *hclin, *hccol;
doublereal *hcval;
integer *hcnnz, *flag__;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;

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
    *hcnnz = 0;
    i__1 = probdata_1.ndim;
    for (j = 1; j <= i__1; ++j) {
	i__2 = probdata_1.ndim;
	for (i__ = j; i__ <= i__2; ++i__) {
	    i__3 = probdata_1.ndim;
	    for (k = i__; k <= i__3; ++k) {
		++(*hcnnz);
		hclin[*hcnnz] = probdata_1.e[k + j * 10 - 11];
		hccol[*hcnnz] = probdata_1.e[i__ + j * 10 - 11];
		hcval[*hcnnz] = probdata_1.p[k + *ind * 10 - 11] * 2. * 
			probdata_1.p[i__ + *ind * 10 - 11];
	    }
	}
    }
} /* ellipsoid_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int ellipsoid_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
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
} /* ellipsoid_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int ellipsoid_endp__(n, x, l, u, m, lambda, rho, equatn, 
	linear)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
{
    extern /* Subroutine */ int ellipsoid_drawsol__();

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
/*     ellipsoid_drawsol__(n, &x[1]); */
} /* ellipsoid_endp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int ellipsoid_cauchyd__(seed, x)
doublereal *seed, *x;
{
    /* Builtin functions */
    double tan();

    /* Local variables */
    extern doublereal drand_();

/*     SCALAR ARGUMENTS */
/*     PARAMETERS */
/*     EXTERNAL FUNCTIONS */
    *x = tan(drand_(seed) * 3.1415926535898);
} /* ellipsoid_cauchyd__ */

/*     ***************************************************************** */
/*     ***************************************************************** */
/* Subroutine */ int ellipsoid_drawsol__(n, x)
integer *n;
doublereal *x;
{
    /* Format strings */
    static char fmt_10[] = "(\002beginfig(\002,i2,\002);\002/,\002u = \002,f\
20.10,\002 cm;\002)";
    static char fmt_20[] = "(\002draw fullcircle\002,/,5x,\002xscaled  \002,\
f20.10,\002u\002,/,5x,\002yscaled  \002,f20.10,\002u\002,/,5x,\002shifted \
(\002,f20.10,\002u,\002,f20.10,\002u) \002,/,5x,\002rotated  \002,f20.10,\
\002;\002)";
    static char fmt_30[] = "(\002pickup pencircle scaled 4pt;\002)";
    static char fmt_40[] = "(\002draw (\002,f20.10,\002u,\002,f20.10,\002u)\
;\002)";
    static char fmt_100[] = "(\002endfig;\002/,\002end;\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_wsfe(), do_fio(), e_wsfe();
    double sqrt(), acos(), d_sign();
    integer f_clos();

    /* Local variables */
    static doublereal xmin, ymin, xmax, ymax;
    static integer i__, j;
    static doublereal m[4]	/* was [2][2] */, alpha, v[2], scale;
    static integer nprob;
    static doublereal aa, bb, cc, dd, lambda[2];

    /* Fortran I/O blocks */
    static cilist io___25 = { 0, 10, 0, fmt_10, 0 };
    static cilist io___34 = { 0, 10, 0, fmt_20, 0 };
    static cilist io___35 = { 0, 10, 0, fmt_30, 0 };
    static cilist io___37 = { 0, 10, 0, fmt_40, 0 };
    static cilist io___38 = { 0, 10, 0, fmt_100, 0 };


/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     This subroutine generate a metapost file with the graphical */
/*     representation of the problem (just for ndim=2). */
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
/*     JUST FOR NDIM=2! */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (probdata_1.ndim != 2) {
	return 0;
    }
/*     PROBLEM ID */
    nprob = 1;
/*     SCALING */
    xmin = 1e99;
    xmax = -1e99;
    i__1 = probdata_1.npun;
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	d__1 = xmin, d__2 = probdata_1.p[j * 10 - 10];
	xmin = min(d__1,d__2);
/* Computing MAX */
	d__1 = xmax, d__2 = probdata_1.p[j * 10 - 10];
	xmax = max(d__1,d__2);
    }
    ymin = 1e99;
    ymax = -1e99;
    i__1 = probdata_1.npun;
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	d__1 = ymin, d__2 = probdata_1.p[j * 10 - 9];
	ymin = min(d__1,d__2);
/* Computing MAX */
	d__1 = ymax, d__2 = probdata_1.p[j * 10 - 9];
	ymax = max(d__1,d__2);
    }
/* Computing MIN */
    d__1 = 10. / (xmax - xmin), d__2 = 10. / (ymax - ymin);
    scale = min(d__1,d__2);
/*     DRAW THE PROBLEM */
    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = 6;
    o__1.ofnm = "sol.mp";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_wsfe(&io___25);
    do_fio(&c__1, (char *)&nprob, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&c_b16, (ftnlen)sizeof(doublereal));
    e_wsfe();
/*     ELLIPSE */
/*     COMPUTE M = L L^T */
/* Computing 2nd power */
    d__1 = x[probdata_1.e[0]];
    m[0] = d__1 * d__1;
    m[1] = x[probdata_1.e[0]] * x[probdata_1.e[1]];
/* Computing 2nd power */
    d__1 = x[probdata_1.e[1]];
/* Computing 2nd power */
    d__2 = x[probdata_1.e[11]];
    m[3] = d__1 * d__1 + d__2 * d__2;
/*     COMPUTE M EIGENVALUES */
    aa = 1.;
    bb = -(m[0] + m[3]);
/* Computing 2nd power */
    d__1 = m[1];
    cc = m[0] * m[3] - d__1 * d__1;
/* Computing 2nd power */
    d__1 = bb;
    dd = sqrt(d__1 * d__1 - aa * 4. * cc);
    lambda[0] = (-bb + dd) / (aa * 2.);
    lambda[1] = (-bb - dd) / (aa * 2.);
/*     COMPUTE ANGLE (USING THE EIGENVECTOR RELATED TO LAMBDA1) */
    if (abs(m[1]) > 1e-8) {
	v[0] = 1.;
	v[1] = -(m[0] - lambda[0]) / m[1];
    } else {
	v[1] = 1.;
	v[0] = -m[1] / (m[0] - lambda[0]);
    }
/* Computing 2nd power */
    d__1 = v[0];
/* Computing 2nd power */
    d__2 = v[1];
    alpha = v[0] / sqrt(d__1 * d__1 + d__2 * d__2);
    alpha = acos(alpha) / 3.1415926535898 * 180.;
    alpha = d_sign(&alpha, &v[1]);
/*     COMPUTE ELLIPSE SEMI-AXIS */
    lambda[0] = sqrt(1. / lambda[0]);
    lambda[1] = sqrt(1. / lambda[1]);
    s_wsfe(&io___34);
    d__1 = lambda[0] * 2. * scale;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    d__2 = lambda[1] * 2. * scale;
    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&c_b22, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&c_b22, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&alpha, (ftnlen)sizeof(doublereal));
    e_wsfe();
/*     DOTS WITHOUT LABELS */
    s_wsfe(&io___35);
    e_wsfe();
    i__1 = probdata_1.npun;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsfe(&io___37);
	d__1 = scale * probdata_1.p[i__ * 10 - 10];
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	d__2 = scale * probdata_1.p[i__ * 10 - 9];
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/*     END */
    s_wsfe(&io___38);
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*     NON-EXECUTABLE STATEMENTS */
} /* ellipsoid_drawsol__ */

