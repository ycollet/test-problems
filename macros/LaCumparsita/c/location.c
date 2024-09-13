/* location.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal ccent[2000], radii[1000];
    integer nc;
} circl_;

#define circl_1 circl_

struct {
    integer nvs[1000], np, totnvs;
} poly1_;

#define poly1_1 poly1_

struct {
    doublereal vert[26000];
} poly2_;

#define poly2_1 poly2_

struct {
    doublereal edges[39000];
} poly3_;

#define poly3_1 poly3_

struct {
    doublereal pcent[2000];
} poly4_;

#define poly4_1 poly4_

struct {
    doublereal pol[13000];
} poly5_;

#define poly5_1 poly5_

struct {
    doublereal xdisp, xscal, ydisp, yscal;
} ellip_;

#define ellip_1 ellip_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__9 = 9;

/*     ================================================================= */
/*     File: location.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     September 14th, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     Location problem */
/*     ---------------- */
/*     This is a variant of the family of location problems introduced */
/*     in [1]. In the original problem, given a set of np disjoint */
/*     polygons P1, P2, ..., Pnp in R^2 one wishes to find the point z1 */
/*     in P1 that minimizes the sum of the distances to the other */
/*     polygons. In this variant [2], we have, in addition to the np */
/*     polygons, nc circles. Moreover, there is an ellipse which has a */
/*     non empty intersection with P1 and such that z1 must be inside */
/*     the ellipse and zi, i = 2, ..., np+nc must be outside. The */
/*     detailed formulation can be found in [2]. */

/*     [1] E. G. Birgin, J. M. Martínez and M. Raydan, "Algorithm 813: */
/*     SPG - software for convex-constrained optimization", ACM */
/*     Transactions on Mathematical Software 27, pp. 340-349, 2001. */

/*     [2] R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt, */
/*     "On Augmented Lagrangian methods with general lower-level */
/*     constraints", submitted, 2005. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int location_inip__(n, x, l, u, m, lambda, rho, equatn, 
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
    integer s_rsle(), do_lio(), e_rsle(), s_wsle(), e_wsle();
    /* Subroutine */ int s_stop();

    /* Local variables */
    static doublereal rmin, rmax, step;
    static integer pnum, i__;
    static doublereal xstep, ystep;
    static integer nx, ny;
    static doublereal secmar;
    static integer inform__;
    static doublereal procit;
    static integer nvmapp, nvmipp;
    static doublereal propol;
    extern /* Subroutine */ int location_genpro__();

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 5, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };


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
/*     EXTERNAL SUBROUTINES */
/*     COMMON BLOCKS */
/*     DEFINE SOME PROBLEM PARAMETERS */
/*     GRID STEP (THIS IS JUST FOR SCALING THE PROBLEM) */
    /* Parameter adjustments */
    --linear;
    --equatn;
    --rho;
    --lambda;
    --u;
    --l;
    --x;

    /* Function Body */
    step = 1.;
/*     SECURITY MARGIN (BETWEEN 0 AND 1) TO AVOID CITIES OVERLAPPING */
    secmar = .1;
/*     DIFFERENT STEPS MAY BE USED IN THE GRID */
    xstep = step;
    ystep = step;
/*     MINIMUM AND MAXIMUM RADII FOR THE CIRCLE- AND POLYGON-CITIES */
/*     GENERATION */
    rmin = step / 2. * .5;
    rmax = (1. - secmar) * (step / 2.);
/*     READ PROBLEM VARIABLE DATA */
/*     THIS PROBLEM DATA DEFINES DIFFERENT INSTANCIES OF THE PROBLEM AND */
/*     IS RELATED TO THE NUMBER OF VARIABLES AND CONSTRAINTS WHICH */
/*     APPROXIMATELY ARE: */

/*     - THE NUMBER OF POLYGON-CITIES IS APPROXIMATELY: */
/*       NPC = NX * NY * PROCIT * PROPOL */

/*     - THE NUMBER OF CIRCLE-CITIES IS APPROXIMATELY: */
/*       NCC = MX 8 NY * PROCIT * ( 1 - PROPOL ) */

/*     - THE TOTAL NUMBER IF CITIES IS APPROXIMATELY: */
/*       NC = NPC + NCC */

/*     - NUMBER OF VARIABLES: */
/*       2 * ( NC + 1 ) */

/*     - NUMBER OF LINEAR CONSTRAINTS: */
/*       NPC * (NVMAPP - NVMIPP ) / 2 + 4 */

/*     - NUMBER OF NONLINEAR CONSTRAINTS: */
/*       NCC + NC + 1 */

/*     WHERE: */

/*     - PNUM IS A NUMBER THAT IDENTIFIES THE PROBLEM INSTANCE, */
/*     - NX IS THE NUMBER OF POINTS IN THE GRID ABSCISSA, */
/*     - NY IS THE NUMBER OF POINTS IN THE GRID ORDINATE, */
/*     - PROCIT IS THE PROBABILITY OF HAVING A CITY AT A GRID POINT, */
/*     - PROPOL IS THE PROBABILITY OF A CITY TO BE A POLYGON (THE CITIES */
/*       WHICH ARE NOT POLYGONS ARE CIRCLES), */
/*     - NVMIPP IS THE MINIMUM NUNBER OF VERTICES OF A POLYGON, */
/*     - NVMAPP IS THE MAXIMUM NUMBER OF VERTICES OF A POLYGON. */
    s_rsle(&io___7);
    do_lio(&c__3, &c__1, (char *)&pnum, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nx, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&ny, (ftnlen)sizeof(integer));
    do_lio(&c__5, &c__1, (char *)&procit, (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&propol, (ftnlen)sizeof(doublereal));
    do_lio(&c__3, &c__1, (char *)&nvmipp, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nvmapp, (ftnlen)sizeof(integer));
    e_rsle();
    if (procit < 0. || procit > 1.) {
	s_wsle(&io___15);
	do_lio(&c__9, &c__1, "PROCIT MUST BE BETWEEN ZERO AND ONE", (ftnlen)
		35);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    if (propol < 0. || propol > 1.) {
	s_wsle(&io___16);
	do_lio(&c__9, &c__1, "PROPOL MUST BE BETWEEN ZERO AND ONE", (ftnlen)
		35);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    if (nvmipp < 3) {
	s_wsle(&io___17);
	do_lio(&c__9, &c__1, "NVMIPP MUST BE GREATER THAN OR EQUAL TO 3", (
		ftnlen)41);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    if (nvmapp < nvmipp) {
	s_wsle(&io___18);
	do_lio(&c__9, &c__1, "NVMAPP MUST BE GREATER THAN OR EQUAL TO", (
		ftnlen)39);
	e_wsle();
	s_wsle(&io___19);
	do_lio(&c__9, &c__1, "NVMIPP. INCREASE NVMAPP OR REDUCE NVMIPP.", (
		ftnlen)41);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
/*     GENERATE THE PROBLEM */
    location_genpro__(&nx, &ny, &xstep, &ystep, &procit, &propol, &nvmipp, &
	    nvmapp, &rmin, &rmax, &inform__);
    if (inform__ != 0) {
	s_wsle(&io___21);
	do_lio(&c__9, &c__1, "THERE WAS AN ERROR IN THE PROBLEM GENERATION", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___22);
	do_lio(&c__9, &c__1, "(PROBABLY DUE TO MEMORY SPACE AVAILABILITY).", (
		ftnlen)44);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
/*     Number of variables */
    *n = poly1_1.np + circl_1.nc << 1;
/*     Initial point */
    i__1 = poly1_1.np;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[(i__ << 1) - 1] = poly4_1.pcent[(i__ << 1) - 2];
	x[i__ * 2] = poly4_1.pcent[(i__ << 1) - 1];
    }
    i__1 = circl_1.nc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[(poly1_1.np + i__ << 1) - 1] = circl_1.ccent[(i__ << 1) - 2];
	x[(poly1_1.np + i__) * 2] = circl_1.ccent[(i__ << 1) - 1];
    }
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e20;
	u[i__] = 1e20;
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = poly1_1.np + circl_1.nc + poly1_1.totnvs + circl_1.nc;
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
/*     Ellipse in and out */
    i__1 = poly1_1.np + circl_1.nc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	linear[i__] = FALSE_;
    }
/*     Polygons */
    i__1 = poly1_1.np + circl_1.nc + poly1_1.totnvs;
    for (i__ = poly1_1.np + circl_1.nc + 1; i__ <= i__1; ++i__) {
	linear[i__] = TRUE_;
    }
/*     Circles */
    i__1 = poly1_1.np + circl_1.nc + poly1_1.totnvs + circl_1.nc;
    for (i__ = poly1_1.np + circl_1.nc + poly1_1.totnvs + 1; i__ <= i__1; 
	    ++i__) {
	linear[i__] = FALSE_;
    }
} /* location_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int location_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static doublereal dist, diff1, diff2;
    static integer i__, ndist;

    /* Fortran I/O blocks */
    static cilist io___29 = { 0, 6, 0, 0, 0 };


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
/*     LOCAL SCALARS */
/*     INTRINSIC FUNCTIONS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    *f = 0.;
    ndist = *n / 2 - 1;
    i__1 = ndist;
    for (i__ = 1; i__ <= i__1; ++i__) {
	diff1 = x[1] - x[(i__ << 1) + 1];
	diff2 = x[2] - x[(i__ << 1) + 2];
/* Computing 2nd power */
	d__1 = diff1;
/* Computing 2nd power */
	d__2 = diff2;
	dist = sqrt(d__1 * d__1 + d__2 * d__2);
	if (dist <= 1e-4) {
	    s_wsle(&io___29);
	    do_lio(&c__9, &c__1, "ERROR IN PROBLEM DEFINITION (DIST TOO SMAL\
L)", (ftnlen)44);
	    e_wsle();
	    *flag__ = -1;
	}
	*f += dist;
    }
    *f /= ndist;
} /* location_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int location_evalg__(n, x, g, flag__)
integer *n;
doublereal *x, *g;
integer *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal dist, diff1, diff2;
    static integer i__, ndist;

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
/*     INTRINSIC FUNCTIONS */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    *flag__ = 0;
    g[1] = 0.;
    g[2] = 0.;
    ndist = *n / 2 - 1;
    i__1 = ndist;
    for (i__ = 1; i__ <= i__1; ++i__) {
	diff1 = x[1] - x[(i__ << 1) + 1];
	diff2 = x[2] - x[(i__ << 1) + 2];
/* Computing 2nd power */
	d__1 = diff1;
/* Computing 2nd power */
	d__2 = diff2;
	dist = sqrt(d__1 * d__1 + d__2 * d__2);
	g[(i__ << 1) + 1] = -(diff1 / dist) / ndist;
	g[(i__ << 1) + 2] = -(diff2 / dist) / ndist;
	g[1] -= g[(i__ << 1) + 1];
	g[2] -= g[(i__ << 1) + 2];
    }
} /* location_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int location_evalh__(n, x, hlin, hcol, hval, hnnz, flag__)
integer *n;
doublereal *x;
integer *hlin, *hcol;
doublereal *hval;
integer *hnnz, *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal dist, diff1, diff2;
    static integer i__, ndist;
    static doublereal g11, g21, g22;

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
/*     LOCAL SCALARS */
/*     INTRINSIC FUNCTIONS */
    /* Parameter adjustments */
    --x;
    --hlin;
    --hcol;
    --hval;

    /* Function Body */
    *flag__ = 0;
    ndist = *n / 2 - 1;
    *hnnz = 0;
    i__1 = ndist;
    for (i__ = 1; i__ <= i__1; ++i__) {
	diff1 = x[1] - x[(i__ << 1) + 1];
	diff2 = x[2] - x[(i__ << 1) + 2];
/* Computing 2nd power */
	d__1 = diff1;
/* Computing 2nd power */
	d__2 = diff2;
	dist = sqrt(d__1 * d__1 + d__2 * d__2);
/* Computing 2nd power */
	d__1 = dist;
	g11 = -(diff1 * diff1 / (d__1 * d__1) - 1.) / dist / ndist;
/* Computing 2nd power */
	d__1 = dist;
	g22 = -(diff2 * diff2 / (d__1 * d__1) - 1.) / dist / ndist;
/* Computing 2nd power */
	d__1 = dist;
	g21 = -(diff1 * diff2 / (d__1 * d__1)) / dist / ndist;
	hlin[*hnnz + 1] = 1;
	hcol[*hnnz + 1] = 1;
	hval[*hnnz + 1] = g11;
	hlin[*hnnz + 2] = (i__ << 1) + 1;
	hcol[*hnnz + 2] = (i__ << 1) + 1;
	hval[*hnnz + 2] = g11;
	hlin[*hnnz + 3] = (i__ << 1) + 1;
	hcol[*hnnz + 3] = 1;
	hval[*hnnz + 3] = -g11;
	hlin[*hnnz + 4] = 2;
	hcol[*hnnz + 4] = 2;
	hval[*hnnz + 4] = g22;
	hlin[*hnnz + 5] = (i__ << 1) + 2;
	hcol[*hnnz + 5] = (i__ << 1) + 2;
	hval[*hnnz + 5] = g22;
	hlin[*hnnz + 6] = (i__ << 1) + 2;
	hcol[*hnnz + 6] = 2;
	hval[*hnnz + 6] = -g22;
	hlin[*hnnz + 7] = 2;
	hcol[*hnnz + 7] = 1;
	hval[*hnnz + 7] = g21;
	hlin[*hnnz + 8] = (i__ << 1) + 2;
	hcol[*hnnz + 8] = (i__ << 1) + 1;
	hval[*hnnz + 8] = g21;
	hlin[*hnnz + 9] = (i__ << 1) + 1;
	hcol[*hnnz + 9] = 2;
	hval[*hnnz + 9] = -g21;
	hlin[*hnnz + 10] = (i__ << 1) + 2;
	hcol[*hnnz + 10] = 1;
	hval[*hnnz + 10] = -g21;
	*hnnz += 10;
    }
} /* location_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int location_evalc__(n, x, ind, c__, flag__)
integer *n;
doublereal *x;
integer *ind;
doublereal *c__;
integer *flag__;
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j;
    static doublereal p[2];

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
/*     Ellipse_in */
    if (*ind == 1) {
	p[0] = (x[(*ind << 1) - 1] - ellip_1.xdisp) / ellip_1.xscal;
	p[1] = (x[*ind * 2] - ellip_1.ydisp) / ellip_1.yscal;
/* Computing 2nd power */
	d__1 = p[0];
/* Computing 2nd power */
	d__2 = p[1];
	*c__ = d__1 * d__1 + d__2 * d__2 - 1.;
/*     Ellipse_out */
    } else if (*ind <= poly1_1.np + circl_1.nc) {
	p[0] = (x[(*ind << 1) - 1] - ellip_1.xdisp) / ellip_1.xscal;
	p[1] = (x[*ind * 2] - ellip_1.ydisp) / ellip_1.yscal;
/* Computing 2nd power */
	d__1 = p[0];
/* Computing 2nd power */
	d__2 = p[1];
	*c__ = 1. - d__1 * d__1 - d__2 * d__2;
/*     Polygons */
    } else if (*ind <= poly1_1.np + circl_1.nc + poly1_1.totnvs) {
	j = *ind - (poly1_1.np + circl_1.nc);
	i__ = (integer) poly5_1.pol[j - 1];
	*c__ = poly3_1.edges[j * 3 - 3] * x[(i__ << 1) - 1] + poly3_1.edges[j 
		* 3 - 2] * x[i__ * 2] + poly3_1.edges[j * 3 - 1];
/*     Circles */
    } else {
/* if ( ind .le. no + nc + totnvs + nc ) then */
	j = *ind - (poly1_1.np + circl_1.nc + poly1_1.totnvs);
	i__ = j + poly1_1.np;
/* Computing 2nd power */
	d__1 = x[(i__ << 1) - 1] - circl_1.ccent[(j << 1) - 2];
/* Computing 2nd power */
	d__2 = x[i__ * 2] - circl_1.ccent[(j << 1) - 1];
/* Computing 2nd power */
	d__3 = circl_1.radii[j - 1];
	*c__ = d__1 * d__1 + d__2 * d__2 - d__3 * d__3;
    }
} /* location_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int location_evaljac__(n, x, ind, jcvar, jcval, jcnnz, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *jcvar;
doublereal *jcval;
integer *jcnnz, *flag__;
{
    static integer i__, j;
    static doublereal p[2];

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
/*     LOCAL ARRAYS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --jcval;
    --jcvar;
    --x;

    /* Function Body */
    *flag__ = 0;
/*     Ellipse_in */
    if (*ind == 1) {
	p[0] = (x[(*ind << 1) - 1] - ellip_1.xdisp) / ellip_1.xscal;
	p[1] = (x[*ind * 2] - ellip_1.ydisp) / ellip_1.yscal;
	*jcnnz = 2;
	jcvar[1] = (*ind << 1) - 1;
	jcval[1] = p[0] * 2. / ellip_1.xscal;
	jcvar[2] = *ind << 1;
	jcval[2] = p[1] * 2. / ellip_1.yscal;
/*     Ellipse_out */
    } else if (*ind <= poly1_1.np + circl_1.nc) {
	p[0] = (x[(*ind << 1) - 1] - ellip_1.xdisp) / ellip_1.xscal;
	p[1] = (x[*ind * 2] - ellip_1.ydisp) / ellip_1.yscal;
	*jcnnz = 2;
	jcvar[1] = (*ind << 1) - 1;
	jcval[1] = p[0] * -2. / ellip_1.xscal;
	jcvar[2] = *ind << 1;
	jcval[2] = p[1] * -2. / ellip_1.yscal;
/*     Polygons */
    } else if (*ind <= poly1_1.np + circl_1.nc + poly1_1.totnvs) {
	j = *ind - (poly1_1.np + circl_1.nc);
	i__ = (integer) poly5_1.pol[j - 1];
	*jcnnz = 2;
	jcvar[1] = (i__ << 1) - 1;
	jcval[1] = poly3_1.edges[j * 3 - 3];
	jcvar[2] = i__ << 1;
	jcval[2] = poly3_1.edges[j * 3 - 2];
/*     Circles */
    } else {
/* if ( ind .le. no + nc + totnvs + nc ) then */
	j = *ind - (poly1_1.np + circl_1.nc + poly1_1.totnvs);
	i__ = j + poly1_1.np;
	*jcnnz = 2;
	jcvar[1] = (i__ << 1) - 1;
	jcval[1] = (x[(i__ << 1) - 1] - circl_1.ccent[(j << 1) - 2]) * 2.;
	jcvar[2] = i__ << 1;
	jcval[2] = (x[i__ * 2] - circl_1.ccent[(j << 1) - 1]) * 2.;
    }
} /* location_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int location_evalhc__(n, x, ind, hclin, hccol, hcval, hcnnz, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *hclin, *hccol;
doublereal *hcval;
integer *hcnnz, *flag__;
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static integer i__, j;

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
/*     Ellipse_in */
    if (*ind == 1) {
	*hcnnz = 2;
	hclin[1] = (*ind << 1) - 1;
	hccol[1] = (*ind << 1) - 1;
/* Computing 2nd power */
	d__1 = ellip_1.xscal;
	hcval[1] = 2. / (d__1 * d__1);
	hclin[2] = *ind << 1;
	hccol[2] = *ind << 1;
/* Computing 2nd power */
	d__1 = ellip_1.yscal;
	hcval[2] = 2. / (d__1 * d__1);
/*     Ellipse_out */
    } else if (*ind <= poly1_1.np + circl_1.nc) {
	*hcnnz = 2;
	hclin[1] = (*ind << 1) - 1;
	hccol[1] = (*ind << 1) - 1;
/* Computing 2nd power */
	d__1 = ellip_1.xscal;
	hcval[1] = -2. / (d__1 * d__1);
	hclin[2] = *ind << 1;
	hccol[2] = *ind << 1;
/* Computing 2nd power */
	d__1 = ellip_1.yscal;
	hcval[2] = -2. / (d__1 * d__1);
/*     Polygons */
    } else if (*ind <= poly1_1.np + circl_1.nc + poly1_1.totnvs) {
	*hcnnz = 0;
/*     Circles */
    } else {
/* if ( ind .le. no + nc + totnvs + nc ) then */
	j = *ind - (poly1_1.np + circl_1.nc + poly1_1.totnvs);
	i__ = j + poly1_1.np;
	*hcnnz = 2;
	hclin[1] = (i__ << 1) - 1;
	hccol[1] = (i__ << 1) - 1;
	hcval[1] = 2.;
	hclin[2] = i__ << 1;
	hccol[2] = i__ << 1;
	hcval[2] = 2.;
    }
} /* location_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int location_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
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
} /* location_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int location_endp__(n, x, l, u, m, lambda, rho, equatn, 
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
} /* location_endp__ */

/*     ***************************************************************** */
/*     ***************************************************************** */
/* Subroutine */ int location_genpro__(nx, ny, xstep, ystep, procit, propol, 
	nvmipp, nvmapp, rmin, rmax, inform__)
integer *nx, *ny;
doublereal *xstep, *ystep, *procit, *propol;
integer *nvmipp, *nvmapp;
doublereal *rmin, *rmax;
integer *inform__;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal seed, c__;
    static integer i__, j, k;
    extern doublereal drand_();
    static doublereal cx, cy, lx, ly, ux, vx, vy, uy;
    static logical inscre, insrec, outcre;
    extern /* Subroutine */ int location_gencir__(), location_genpol__();

/*     This subroutine generates a location problem (see the Figure). */
/*     First, a regular grid with nx horizontal points and ny vertical */
/*     points in the positive orthant is considered. The points of the */
/*     grid start at the origin with an horizontal distance of xstep and */
/*     a vertical distance of ystept. This grid will be the space at */
/*     which the cities (represented by polygons) will be distributed. */
/*     Before building the cities, an area (rectangle) of preservation */
/*     where almost nothing can be done, is defined. This area of */
/*     preservation will receive, after the construction of the cities, */
/*     an hydraulic plant of energy generation (to supply the energy to */
/*     the cities). Then, in the rest of the space, the cities are built. */
/*     At each point of the grid (out of the central region) a city */
/*     (represented by a polygon) will be built with probability procit. */
/*     The definition of the polygon uses variables nvmipp, nvmapp, rmin */
/*     and rmax in a way described in the genpol (generate polygon) */
/*     subroutine. To transmit the energy from the plant to the cities, a */
/*     tower inside each city and a tower inside the central region must */
/*     be built. The objective of the problem is to determine the */
/*     location of this towers in order to minimize the sum of the */
/*     distances from each city tower to the central one. */

/*     On Entry: */

/*     nx    integer, */
/*           number of horizontal points in the grid, */

/*     ny    integer, */
/*           number of vertical points in the grid, */

/*     xstep double precision, */
/*           horizontal distance between points of the grid, */

/*     ystep double precision, */
/*           vertical distance between points of the grid, */

/*     procit  double precision, */
/*           probability of defining a city at point of the grid */
/*           (0 <= procit <= 1), */

/*     propol double precision, */
/*           probability of the city been defined by a polygon (0 <= */
/*           propol <= 1), if the city is not defined by a polygon then */
/*           it is defined by a circle, */

/*     nvmipp integer, */
/*           parameter for the polygon generation (described in genpol */
/*           subroutine), */

/*     nvmapp integer, */
/*           parameter for the polygon generation (described in genpol */
/*           subroutine), */

/*     rmin  double precision, */
/*           parameter for the polygon generation (described in genpol */
/*           subroutine), */

/*     rmax  double precision, */
/*           parameter for the polygon generation (described in genpol */
/*           subroutine). */

/*     On output: */

/*     As described in the genpol subroutine, the output is saved in the */
/*     polyg common block. */
/*     PARAMETERS */
/*     SCALAR ARGUMENTS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     EXTERNAL FUNCTIONS */
/*     EXTERNAL SUBROUTINES */
/*     COMMON BLOCKS */
    *inform__ = 0;
/*     SEED FOR THE RANDOM GENERATION */
    seed = 760013.;
/*     DEFINE CENTRAL REGION */
    ellip_1.xdisp = (*nx - 1) * .2 * *xstep;
    ellip_1.ydisp = (*ny - 1) * .5 * *ystep;
    ellip_1.xscal = (*nx - 1) * .15 * *xstep;
    ellip_1.yscal = (*ny - 1) * .5 * *ystep;
/*     DEFINE CENTRAL POINT POLYGON */
    poly1_1.nvs[0] = 4;
    lx = ellip_1.xdisp + ellip_1.xscal * .5;
    ux = ellip_1.xdisp + ellip_1.xscal * 1.5;
    ly = ellip_1.ydisp - ellip_1.yscal * .75;
    uy = ellip_1.ydisp + ellip_1.yscal * .75;
    poly2_1.vert[0] = lx;
    poly2_1.vert[1] = ly;
    poly2_1.vert[2] = lx;
    poly2_1.vert[3] = uy;
    poly2_1.vert[4] = ux;
    poly2_1.vert[5] = uy;
    poly2_1.vert[6] = ux;
    poly2_1.vert[7] = ly;
    poly3_1.edges[0] = -1.;
    poly3_1.edges[1] = 0.;
    poly3_1.edges[2] = lx;
    poly3_1.edges[3] = 0.;
    poly3_1.edges[4] = 1.;
    poly3_1.edges[5] = -uy;
    poly3_1.edges[6] = 1.;
    poly3_1.edges[7] = 0.;
    poly3_1.edges[8] = -ux;
    poly3_1.edges[9] = 0.;
    poly3_1.edges[10] = -1.;
    poly3_1.edges[11] = ly;
    poly1_1.np = 1;
    poly1_1.totnvs = 4;
    poly4_1.pcent[0] = (lx + ux) / 2.;
    poly4_1.pcent[1] = (ly + uy) / 2.;
    for (i__ = 1; i__ <= 4; ++i__) {
	poly5_1.pol[i__ - 1] = 1.;
    }
/*     THESE THREE LINES BELOW SHOULD BE USED INSTEAD OF THE 27 LINES */
/*     ABOVE IF THE CENTRAL RECTANGLE WOULD NOT TO BE CONSIDERED. */
/*     nvs(1)    =       0 */
/*     np        =       1 */
/*     totnvs    =       0 */
    circl_1.nc = 0;
/*     DEFINE CITY-POLYGONS AND CITY-CIRCLES CENTERED AT THE GRID */
/*     POINTS. THE CITIES MUST BE OUTSIDE THE CENTRAL REGION AND NOT */
/*     COMPLETELY INSIDE THE CASSINIAN CURVE */
    i__1 = *nx - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	i__2 = *ny - 1;
	for (j = 0; j <= i__2; ++j) {
	    cx = i__ * *xstep;
	    cy = j * *ystep;
	    if (drand_(&seed) <= *propol) {
/* GENERATE A CITY-POLYGON */
		location_genpol__(&cx, &cy, nvmipp, nvmapp, rmin, rmax, &seed,
			 inform__);
		if (*inform__ != 0) {
		    return 0;
		}
/* TEST WHETHER THE POLYGON IS IN THE INTERIOR, IT IS */
/* OUTSIDE OR IT CUTS CENTRAL REGION */
/* 1) CITIES IN THE INTERIOR ARE FORBIDDEN TO AVOID */
/* INFEASIBLE PROBLEMS */
/* 2) CITIES THAT TOUCH THE CENTRAL REGION ARE DESIRED */
/* WITH THE AIM OF MAKING THE CENTRAL-REGION CONSTRAINT */
/* ACTIVE AT THE SOLUTION AS MANY TIMES AS POSSIBLE */
/* 3) CITIES THAT DO NOT TOUCH THE CENTRAL REGION ARE */
/* ARE INTRODUCED INTO THE PROBLEM WITH PROBABILITY PROCI */
/* IT IS NOT EASY TO DETERMINE THE PREVIOUS RELATIONS (IN */
/* OUTSIDE, HALF AND A HALF) BETWEEN THE POLYGON AND THE */
/* NON-C0NVEX CENTRAL REGION. SO WE WILL CONSIDER THAT: */
/* A POLYGON IS INSIDE THE CENTRAL REGION IF ALL ITS VERT */
/* ARE INSIDE OR IN THE BORDER OF THE CENTRAL REGION. */
/* A POLYGON IS OUTSIDE THE CENTRAL REGION IF ALL ITS VER */
/* ARE OUTSIDE OR IN THE BORDER OF THE CENTRAL REGION. */
/* A POLYGON CUTS THE CENTRAL REGION IF IT HAS AT LEAST A */
/* VERTEX INSIDE OR IN THE BORDER OF THE CENTRAL REGION A */
/* AT LEAST A VERTEX OUTSIDE OR IN THE BORDER OF THE CENT */
/* REGION. SO, TO BE HALF AND HALF IS EQUIVALENT TO NOT T */
/* INSIDE, NOR OUTISDE. */
		inscre = TRUE_;
		outcre = TRUE_;
		i__3 = poly1_1.nvs[poly1_1.np];
		for (k = 1; k <= i__3; ++k) {
		    vx = (poly2_1.vert[(poly1_1.totnvs + k << 1) - 2] - 
			    ellip_1.xdisp) / ellip_1.xscal;
		    vy = (poly2_1.vert[(poly1_1.totnvs + k << 1) - 1] - 
			    ellip_1.ydisp) / ellip_1.yscal;
/* Computing 2nd power */
		    d__1 = vx;
/* Computing 2nd power */
		    d__2 = vy;
		    c__ = d__1 * d__1 + d__2 * d__2 - 1.;
		    if (c__ > 0.) {
			inscre = FALSE_;
		    }
		    if (c__ < 0.) {
			outcre = FALSE_;
		    }
		}
		if (cx >= lx - *rmax && cx <= ux + *rmax && cy >= ly - *rmax 
			&& cy <= uy + *rmax) {
		    insrec = TRUE_;
		} else {
		    insrec = FALSE_;
		}
/* ADD THE CITY-POLYGON TO THE PROBLEM */
		if (! insrec && (! inscre && ! outcre || outcre && drand_(&
			seed) <= *procit)) {
		    ++poly1_1.np;
		    poly1_1.totnvs += poly1_1.nvs[poly1_1.np - 1];
		}
	    } else {
/* GENERATE A CITY-CIRCLE */
		location_gencir__(&cx, &cy, rmin, rmax, &seed, inform__);
		if (*inform__ != 0) {
		    return 0;
		}
/* VERIFY WHERE IT IS */
		inscre = TRUE_;
		outcre = TRUE_;
		for (k = 1; k <= 4; ++k) {
		    if (k == 1) {
			vx = (cx - circl_1.radii[circl_1.nc] - ellip_1.xdisp) 
				/ ellip_1.xscal;
			vy = (cy - ellip_1.ydisp) / ellip_1.yscal;
		    } else if (k == 2) {
			vx = (cx + circl_1.radii[circl_1.nc] - ellip_1.xdisp) 
				/ ellip_1.xscal;
			vy = (cy - ellip_1.ydisp) / ellip_1.yscal;
		    } else if (k == 3) {
			vx = (cx - ellip_1.xdisp) / ellip_1.xscal;
			vy = (cy - circl_1.radii[circl_1.nc] - ellip_1.ydisp) 
				/ ellip_1.yscal;
		    } else if (k == 4) {
			vx = (cx - ellip_1.xdisp) / ellip_1.xscal;
			vy = (cy + circl_1.radii[circl_1.nc] - ellip_1.ydisp) 
				/ ellip_1.yscal;
		    }
/* Computing 2nd power */
		    d__1 = vx;
/* Computing 2nd power */
		    d__2 = vy;
		    c__ = d__1 * d__1 + d__2 * d__2 - 1.;
		    if (c__ > 0.) {
			inscre = FALSE_;
		    }
		    if (c__ < 0.) {
			outcre = FALSE_;
		    }
		}
		if (cx >= lx - *rmax && cx <= ux + *rmax && cy >= ly - *rmax 
			&& cy <= uy + *rmax) {
		    insrec = TRUE_;
		} else {
		    insrec = FALSE_;
		}
/* ADD THE CITY-CIRCLE TO THE PROBLEM */
		if (! insrec && (! inscre && ! outcre || outcre && drand_(&
			seed) <= *procit)) {
		    ++circl_1.nc;
		}
	    }
	}
    }
} /* location_genpro__ */

/*     ***************************************************************** */
/*     ***************************************************************** */
/* Subroutine */ int location_gencir__(cx, cy, rmin, rmax, seed, inform__)
doublereal *cx, *cy, *rmin, *rmax, *seed;
integer *inform__;
{
    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    extern doublereal drand_();

    /* Fortran I/O blocks */
    static cilist io___67 = { 0, 6, 0, 0, 0 };
    static cilist io___68 = { 0, 6, 0, 0, 0 };
    static cilist io___69 = { 0, 6, 0, 0, 0 };
    static cilist io___70 = { 0, 6, 0, 0, 0 };
    static cilist io___71 = { 0, 6, 0, 0, 0 };


/*     PARAMETERS */
/*     SCALAR ARGUMENTS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     EXTERNAL FUNCTIONS */
/*     COMMON BLOCKS */
    *inform__ = 0;
/*     VERIFY SPACE AVAILABILITY FOR CIRCLES */
    if (circl_1.nc == 1000) {
	s_wsle(&io___67);
	do_lio(&c__9, &c__1, "THE MAXIMUM NUMBER OF CIRCLES WAS ACHIEVED.", (
		ftnlen)43);
	e_wsle();
	s_wsle(&io___68);
	do_lio(&c__9, &c__1, "INCREASE NCMAX OR REDUCE THE NUMBER OF GRID ", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___69);
	do_lio(&c__9, &c__1, "POINTS (NX*NY) OR THE PROBABILITY OF HAVING ", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___70);
	do_lio(&c__9, &c__1, "A CITY-CIRCLE AT A GRID POINT", (ftnlen)29);
	e_wsle();
	s_wsle(&io___71);
	do_lio(&c__9, &c__1, "(PROCIT*(1-PROPOL)).", (ftnlen)20);
	e_wsle();
	*inform__ = -1;
	return 0;
    }
/*     SAVE THE CENTER */
    circl_1.ccent[circl_1.nc * 2] = *cx;
    circl_1.ccent[(circl_1.nc << 1) + 1] = *cy;
/*     GENERATE THE RADIUS */
    circl_1.radii[circl_1.nc] = *rmin + (*rmax - *rmin) * drand_(seed);
} /* location_gencir__ */

/*     ***************************************************************** */
/*     ***************************************************************** */
/* Subroutine */ int location_genpol__(cx, cy, nvmipp, nvmapp, rmin, rmax, 
	seed, inform__)
doublereal *cx, *cy;
integer *nvmipp, *nvmapp;
doublereal *rmin, *rmax, *seed;
integer *inform__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();
    double cos(), sin();

    /* Local variables */
    extern /* Subroutine */ int location_constr__();
    static doublereal angl[13000], dist;
    static integer i__;
    static doublereal r__;
    extern doublereal drand_();
    static doublereal lseed;
    extern /* Subroutine */ int location_class__();
    static doublereal mindist;

    /* Fortran I/O blocks */
    static cilist io___72 = { 0, 6, 0, 0, 0 };
    static cilist io___73 = { 0, 6, 0, 0, 0 };
    static cilist io___74 = { 0, 6, 0, 0, 0 };
    static cilist io___75 = { 0, 6, 0, 0, 0 };
    static cilist io___76 = { 0, 6, 0, 0, 0 };
    static cilist io___78 = { 0, 6, 0, 0, 0 };
    static cilist io___79 = { 0, 6, 0, 0, 0 };
    static cilist io___80 = { 0, 6, 0, 0, 0 };
    static cilist io___81 = { 0, 6, 0, 0, 0 };
    static cilist io___82 = { 0, 6, 0, 0, 0 };
    static cilist io___83 = { 0, 6, 0, 0, 0 };
    static cilist io___84 = { 0, 6, 0, 0, 0 };
    static cilist io___85 = { 0, 6, 0, 0, 0 };


/*     This subroutine generates a polygon in R^2 with its */
/*     vertices in a sphere centered at point (cx,cy). The */
/*     number of vertices is randomly generated */
/*     satisfying nvmipp <= number of vertices  <= */
/*     nvmapp. The ratio of the sphere is also randomly */
/*     generated satisfying rmin <= ratio < rmax. The */
/*     generated polygon is stored in the common block "polyg". */

/*     On Entry: */

/*     cx    double precision, */
/*           first coordinate of the center of the sphere, */

/*     cy    double precision, */
/*           second coordinate of the center of the sphere, */

/*     nvmipp integer, */
/*           minimum number of vertices, */

/*     nvmapp integer, */
/*           maximum number of vertices, */

/*     rmin double precision, */
/*           minimum ratio of the sphere, */

/*     rmax double precision, */
/*           maximum ratio of the sphere, */

/*     seed double precision, */
/*           seed for the random generation. */

/*     On Output: */

/*     inform integer */
/*            termination parameter */

/*            0 = everything was ok */

/*          < 0 = there was not enough memory space. */

/*     The generated polygon is stored in the polyg common block */
/*     described below. */

/*     Common block polyg: */

/*     common /polyg/nvs,vert,edges,np,totnvs */

/*     This structure represents, at any time, np polygons. */
/*     Position i of array nvs indicates the number of vertices */
/*     of polygon i. Arrays vert and edges store the */
/*     vertices and edges of the polygons. */

/*     For example, if nvs(1) = 3 it indicates that the first */
/*     polygon has 3 vertices (edges). Then, if the vertices */
/*     are (x1,y1), (x2,y2) and (x3,y3), we have that vert(1) */
/*     = x1, vert(2) = y1, vert(3) = x2, vert(4) = y2, vert(5) */
/*     = x3, and vert(6) = y3. And, if the edges (written as */
/*     ax + by + c = 0) are a1 x + b1 y + c1 = 0, a2 x + b2 y */
/*     + c2 = 0, and a3 x + b3 y + c3 = 0 then edges(1) = a1, */
/*     edges(2) = b1, edges(3) = c1, edges(4) = a2, edges(5) = */
/*     b2, edges(6) = c2, edges(7) = a3, edges(8) = b3 and */
/*     edges(9) = c3. */

/*     totnvs indicates the total number of vertices */
/*     of the set of polygons. This information is used when */
/*     a new polygon is created to know the first free position */
/*     of arrays vert and edges at which the vertices and edges */
/*     of the new polygon will be saved. */

/*     Two additional details: */

/*     1) For each polygon, the vertices are ordered clockwise */
/*     and edge i corresponds to the edge between vertices */
/*     i and i+1 (0 if i=n). */

/*     2) For each edge of the form ax + bx + c = 0, constants */
/*     a, b and c are chosen in such a way that */
/*     (|a| = 1 or |b| = 1) and (a cx + b cy + c <= 0). */
/*     PARAMETERS */
/*     SCALAR ARGUMENTS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     LOCAL SCALARS */
/*     LOCAL ARRAYS */
/*     EXTERNAL FUNCTIONS */
/*     EXTERNAL SUBROUTINES */
/*     INTRINSIC FUNCTIONS */
/*     COMMON BLOCKS */
    *inform__ = 0;
/*     VERIFY SPACE AVAILABILITY FOR POLYGONS */
    if (poly1_1.np == 1000) {
	s_wsle(&io___72);
	do_lio(&c__9, &c__1, "THE MAXIMUM NUMBER OF POLYGONS WAS ACHIEVED.", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___73);
	do_lio(&c__9, &c__1, "INCREASE NPMAX OR REDUCE THE NUMBER OF GRID ", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___74);
	do_lio(&c__9, &c__1, "POINTS (NX*NY) OR THE PROBABILITY OF HAVING ", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___75);
	do_lio(&c__9, &c__1, "A CITY-POLYGON AT A GRID POINT", (ftnlen)30);
	e_wsle();
	s_wsle(&io___76);
	do_lio(&c__9, &c__1, "(PROCIT*PROPOL).", (ftnlen)16);
	e_wsle();
	*inform__ = -1;
	return 0;
    }
/*     GENERATE THE NUMBER OF VERTICES */
    lseed = *seed + 157318.;
    poly1_1.nvs[poly1_1.np] = *nvmipp + (integer) ((*nvmapp - *nvmipp + 1) * 
	    drand_(&lseed));
/*     VERIFY SPACE AVAILABILITY FOR VERTICES AND SIDES */
    if (poly1_1.totnvs + poly1_1.nvs[poly1_1.np] > 13000) {
	s_wsle(&io___78);
	do_lio(&c__9, &c__1, "THE MAXIMUM NUMBER OF POLYGONS VERTICES WAS", (
		ftnlen)43);
	e_wsle();
	s_wsle(&io___79);
	do_lio(&c__9, &c__1, "ACHIEVED. INCREASE NVSMAX OR REDUCE THE", (
		ftnlen)39);
	e_wsle();
	s_wsle(&io___80);
	do_lio(&c__9, &c__1, "AVERAGE NUMBER OF VERTICES (NVMAPP-NVMIPP)/2", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___81);
	do_lio(&c__9, &c__1, "OR THE NUMBER OF POLYGONS. TO REDUCE THE", (
		ftnlen)40);
	e_wsle();
	s_wsle(&io___82);
	do_lio(&c__9, &c__1, "NUMBER OF POLYGONS, REDUCE THE NUMBER OF", (
		ftnlen)40);
	e_wsle();
	s_wsle(&io___83);
	do_lio(&c__9, &c__1, "GRID POINTS (NX*NY) OR THE PROBABILITY OF", (
		ftnlen)41);
	e_wsle();
	s_wsle(&io___84);
	do_lio(&c__9, &c__1, "HAVING A CITY-POLYGONS AT A GRID POINT", (
		ftnlen)38);
	e_wsle();
	s_wsle(&io___85);
	do_lio(&c__9, &c__1, "(PROCIT*PROPOL).", (ftnlen)16);
	e_wsle();
	*inform__ = -1;
	return 0;
    }
/*     SAVE THE "CENTER" TO BE USED AS INITIAL POINT */
    poly4_1.pcent[poly1_1.np * 2] = *cx;
    poly4_1.pcent[(poly1_1.np << 1) + 1] = *cy;
/*     IDENTIFY THE CONSTRAINTS WITH THE POLYGON */
    i__1 = poly1_1.totnvs + poly1_1.nvs[poly1_1.np];
    for (i__ = poly1_1.totnvs + 1; i__ <= i__1; ++i__) {
	poly5_1.pol[i__ - 1] = (doublereal) (poly1_1.np + 1);
    }
/*     GENERATE THE RADIUS OF THE SPHERE */
    r__ = *rmin + (*rmax - *rmin) * drand_(seed);
/*     GENERATE ALL ANGLES SATISFYING 0 <= ANGLE_I < 2*PI */
L10:
    i__1 = poly1_1.nvs[poly1_1.np];
    for (i__ = 1; i__ <= i__1; ++i__) {
	angl[i__ - 1] = drand_(&lseed) * 6.2831853071796004;
    }
/*     CLASSIFY THE ANGLES IN DECREASING ORDER */
    location_class__(&poly1_1.nvs[poly1_1.np], angl);
/*     CONSTRUCT THE VERTICES */
    i__1 = poly1_1.nvs[poly1_1.np];
    for (i__ = 1; i__ <= i__1; ++i__) {
	poly2_1.vert[(poly1_1.totnvs + i__ << 1) - 2] = *cx + r__ * cos(angl[
		i__ - 1]);
	poly2_1.vert[(poly1_1.totnvs + i__ << 1) - 1] = *cy + r__ * sin(angl[
		i__ - 1]);
    }
/*     FOR NUMERICAL STABILITY IN THE CONSTRAINT CALCULATION, */
/*     AVOID TOO SIMILAR ANGLES (EQUIVALENT TO TOO NEAR POINTS) */
    mindist = 1e99;
    i__1 = poly1_1.totnvs + poly1_1.nvs[poly1_1.np] - 1;
    for (i__ = poly1_1.totnvs + 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = poly2_1.vert[(i__ << 1) - 2] - poly2_1.vert[i__ * 2];
/* Computing 2nd power */
	d__2 = poly2_1.vert[(i__ << 1) - 1] - poly2_1.vert[(i__ << 1) + 1];
	dist = d__1 * d__1 + d__2 * d__2;
	if (dist < mindist) {
	    mindist = dist;
	}
    }
    i__ = poly1_1.totnvs + poly1_1.nvs[poly1_1.np];
/* Computing 2nd power */
    d__1 = poly2_1.vert[(i__ << 1) - 2] - poly2_1.vert[poly1_1.totnvs * 2];
/* Computing 2nd power */
    d__2 = poly2_1.vert[(i__ << 1) - 1] - poly2_1.vert[(poly1_1.totnvs << 1) 
	    + 1];
    dist = d__1 * d__1 + d__2 * d__2;
    if (dist < mindist) {
	mindist = dist;
    }
    if (mindist < 1e-8) {
/*         write (*,fmt=*) 'DISCARDING ANGLES.' */
	goto L10;
    }
/*     CONSTRUCT THE EDGES */
    i__1 = poly1_1.totnvs + poly1_1.nvs[poly1_1.np] - 2;
    for (i__ = poly1_1.totnvs + 1; i__ <= i__1; ++i__) {
	location_constr__(&poly2_1.vert[(i__ << 1) - 2], &poly2_1.vert[(i__ <<
		 1) - 1], &poly2_1.vert[i__ * 2], &poly2_1.vert[(i__ << 1) + 
		1], &poly2_1.vert[(i__ << 1) + 2], &poly2_1.vert[(i__ << 1) + 
		3], &poly3_1.edges[i__ * 3 - 3], &poly3_1.edges[i__ * 3 - 2], 
		&poly3_1.edges[i__ * 3 - 1], inform__);
	if (*inform__ != 0) {
	    return 0;
	}
    }
    i__ = poly1_1.totnvs + poly1_1.nvs[poly1_1.np] - 1;
    location_constr__(&poly2_1.vert[(i__ << 1) - 2], &poly2_1.vert[(i__ << 1) 
	    - 1], &poly2_1.vert[i__ * 2], &poly2_1.vert[(i__ << 1) + 1], &
	    poly2_1.vert[poly1_1.totnvs * 2], &poly2_1.vert[(poly1_1.totnvs <<
	     1) + 1], &poly3_1.edges[i__ * 3 - 3], &poly3_1.edges[i__ * 3 - 2]
	    , &poly3_1.edges[i__ * 3 - 1], inform__);
    if (*inform__ != 0) {
	return 0;
    }
    i__ = poly1_1.totnvs + poly1_1.nvs[poly1_1.np];
    location_constr__(&poly2_1.vert[(i__ << 1) - 2], &poly2_1.vert[(i__ << 1) 
	    - 1], &poly2_1.vert[poly1_1.totnvs * 2], &poly2_1.vert[(
	    poly1_1.totnvs << 1) + 1], &poly2_1.vert[(poly1_1.totnvs << 1) + 
	    2], &poly2_1.vert[(poly1_1.totnvs << 1) + 3], &poly3_1.edges[i__ *
	     3 - 3], &poly3_1.edges[i__ * 3 - 2], &poly3_1.edges[i__ * 3 - 1],
	     inform__);
    if (*inform__ != 0) {
	return 0;
    }
} /* location_genpol__ */

/*     ***************************************************************** */
/*     ***************************************************************** */
/* Subroutine */ int location_constr__(x1, y1, x2, y2, x3, y3, a, b, c__, 
	inform__)
doublereal *x1, *y1, *x2, *y2, *x3, *y3, *a, *b, *c__;
integer *inform__;
{
    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Fortran I/O blocks */
    static cilist io___91 = { 0, 6, 0, 0, 0 };


/*     This subroutine computes the real constants a, b and c of */
/*     the straight line ax + by + c = 0 in R^2 defined by the */
/*     points (x1,y1) and (x2,y2); such that the point (x3,y3) */
/*     satisfies the constraint ax + by + c <= 0. */

/*     On Entry: */

/*     x1    double precision, */
/*           first coordinate of point (x1,y1), */

/*     y1    double precision, */
/*           second coordinate of point (x1,y1), */

/*     x2    double precision, */
/*           first coordinate of point (x2,y2), */

/*     y2    double precision, */
/*           second coordinate of point (x2,y2), */

/*     x3    double precision, */
/*           first coordinate of point (x3,y3), */

/*     y3    double precision, */
/*           second coordinate of point (x3,y3). */

/*     On Return */

/*     a,b,c double precision */
/*           the desired constants. */
/*     SCALAR ARGUMENTS */
    if (*x1 == *x2 && *y1 == *y2) {
	s_wsle(&io___91);
	do_lio(&c__9, &c__1, "ERROR IN FUNCTION CONSTRAINT: X1=X2 AND Y1=Y2", 
		(ftnlen)45);
	e_wsle();
	*inform__ = -1;
	return 0;
    }
    if (*y1 != *y2) {
	*a = 1.;
	*b = -(*x2 - *x1) / (*y2 - *y1);
	*c__ = -(*x1 + *b * *y1);
    } else {
	*a = 0.;
	*b = 1.;
	*c__ = -(*y1);
    }
    if (*a * *x3 + *b * *y3 + *c__ > 0.) {
	*a = -(*a);
	*b = -(*b);
	*c__ = -(*c__);
    }
} /* location_constr__ */

/*     ***************************************************************** */
/*     ***************************************************************** */
/* Subroutine */ int location_class__(n, x)
integer *n;
doublereal *x;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal xmax;
    static integer i__, j;
    static doublereal aux;
    static integer pos;

/*     This subroutine classifies the elements of a vector in */
/*     decreasing order, i.e., on output: x(1) >= x(2) >= */
/*     ... >= x(n). */

/*     On Entry: */

/*     n     integer, */
/*           number of elements of the vector to be classified, */

/*     x     double precision x(n), */
/*           vector to be classified. */

/*     On Return */

/*     x     double precision x(n), */
/*           classified vector. */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xmax = x[i__];
	pos = i__;
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    if (x[j] > xmax) {
		xmax = x[j];
		pos = j;
	    }
	}
	if (pos != i__) {
	    aux = x[i__];
	    x[i__] = x[pos];
	    x[pos] = aux;
	}
    }
} /* location_class__ */

