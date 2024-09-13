/* genpack-cc-mina.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer nite, ndim;
    doublereal iterad;
} packdata_;

#define packdata_1 packdata_

struct {
    integer start[100040004]	/* was [10002][10002] */, next[10000], br[
	    20000]	/* was [2][10000] */, nbr, nreg;
} partdata_;

#define partdata_1 partdata_

struct {
    doublereal diagb[40000]	/* was [2][2][10000] */;
} diblocks_;

#define diblocks_1 diblocks_

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b33 = 0.;

/*     ================================================================= */
/*     File: genpack-cc-mina.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     July 2nd, 2006. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     PACKING CIRCLES INTO CIRCLES */
/*     ---------------------------- */
/*     Given a fixed number of fixed-size identical circles (called */
/*     items), the objective is to find the smallest circle (called */
/*     object) that fits the  items without overlapping. The model */
/*     follows: */

/*     Mininimize R */
/*     subject to d(0,ci)^2 <= (R - r)^2, for all i */
/*                R >= r */
/*                sum_{i < j} max{ 0, (2r)^2 - d(ci,cj)^2 }^2 = 0 */

/*     (The variables are the centers ci in R^2 of the circular items */
/*     and the radius R of the circular object that fits the items.) */

/*     [1] E. G. Birgin and F. N. C. Sobral, "A tool based on */
/*     nonlinear programming models for packing 2D and 3D circular */
/*     items", Computers and Operations Research, 2006 */
/*     (DOI: 10.1016/j.cor.2006.11.002). */

/*     [2] J. M. Martínez and D. P. Ronconi, "Optimizing the */
/*     Packing of Cylinders into a Rectangular Container: A Nonlinear */
/*     Approach", European Journal of Operational Research 160, pp. */
/*     19-33, 2005. */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int genpack_cc_mina_inip__(n, x, l, u, m, lambda, rho, 
	equatn, linear, iterad_in__, nite_in__, seed_in__)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
doublereal *iterad_in__;
integer *nite_in__;
doublereal *seed_in__;
{
    /* Format strings */
    static char fmt_7000[] = "(/,\002 Welcome to GENPACK !!!\002,//,\002 Thi\
s program uses ALGENCAN and the strategy\002,\002 described in:\002,//,\002 \
E. G. Birgin and F. N. C. Sobral, A tool based on\002,\002 nonlinear program\
ming\002,/,\002 models for packing 2D and\002,\002 3D circular items, Comput\
ers and Operations\002,/,\002 Research, 2006 (DOI: 10.1016/j.cor.2006.11.002\
).\002,/,/,\002 (See also: E. G. Birgin, J. M. Mart\355nez and\002,\002 D. P\
. Ronconi, Optimizing the\002,/,\002 Packing of\002,\002 Cylinders into a Ci\
rcular Container: A Nonlinear\002,\002 Approach,\002,/,\002 European Journal\
 of Operational\002,\002 Research 160, pp. 19-33, 2005.)\002,/,/,\002 to sol\
ve the problem of packing a fixed number of\002,\002 (fixed-size) identica\
l\002,/,\002 circles into a circular region, minimizing the area\002,\002 of\
 the region.\002,/)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(), e_wsfe();
    double sqrt();

    /* Local variables */
    extern integer ceil_();
    static doublereal seed, dist, objradub;
    static integer i__, j;
    extern doublereal drand_();
    static doublereal tmp;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, fmt_7000, 0 };


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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     FUNCTIONS */
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
    s_wsfe(&io___1);
    e_wsfe();
/*     Read problem data */
/*     Enter the circular items radius (iterad) */
/*     and the number of items (nite) */
    packdata_1.iterad = *iterad_in__;
    packdata_1.nite = *nite_in__;
/*     Seed for the initial point random generation */
    seed = *seed_in__;
/*     Initialize regions */
/* Computing 2nd power */
    d__1 = packdata_1.iterad;
    objradub = sqrt(d__1 * d__1 * 2. * packdata_1.nite);
    d__1 = objradub / packdata_1.iterad;
    partdata_1.nreg = ceil_(&d__1);
    i__1 = partdata_1.nreg + 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	i__2 = partdata_1.nreg + 1;
	for (j = 0; j <= i__2; ++j) {
	    partdata_1.start[i__ + j * 10002] = 0;
	}
    }
/*     Number of variables */
    *n = packdata_1.ndim * packdata_1.nite + 1;
/*     Initial point */
/* Computing 2nd power */
    d__1 = packdata_1.iterad;
    x[*n] = sqrt(d__1 * d__1 * packdata_1.nite * (drand_(&seed) + 1.));
    tmp = x[*n] - packdata_1.iterad;
    i__ = 1;
L10:
    if (i__ <= packdata_1.nite) {
	dist = 0.;
	i__1 = packdata_1.ndim;
	for (j = 1; j <= i__1; ++j) {
	    x[packdata_1.ndim * i__ + 1 - j] = -tmp + tmp * 2. * drand_(&seed)
		    ;
/* Computing 2nd power */
	    d__1 = x[packdata_1.ndim * i__ + 1 - j];
	    dist += d__1 * d__1;
	}
/* Computing 2nd power */
	d__1 = tmp;
	if (dist <= d__1 * d__1) {
	    ++i__;
	}
	goto L10;
    }
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e20;
	u[i__] = 1e20;
    }
    l[*n] = packdata_1.iterad;
/*     Number of constraints (equalities plus inequalities) */
    *m = packdata_1.nite + 1;
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
    i__1 = *m - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	equatn[i__] = FALSE_;
    }
    equatn[*m] = TRUE_;
/*     For each constraint i, set linear(i) = .true. if it is a linear */
/*     constraint, otherwise set linear(i) = .false. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	linear[i__] = FALSE_;
    }
    return 0;
} /* genpack_cc_mina_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int genpack_cc_mina_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    *f = x[*n];
} /* genpack_cc_mina_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int genpack_cc_mina_evalg__(n, x, g, flag__)
integer *n;
doublereal *x, *g;
integer *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
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
} /* genpack_cc_mina_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int genpack_cc_mina_evalh__(n, x, hlin, hcol, hval, nnzh, 
	flag__)
integer *n;
doublereal *x;
integer *hlin, *hcol;
doublereal *hval;
integer *nnzh, *flag__;
{
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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
    /* Parameter adjustments */
    --x;
    --hlin;
    --hcol;
    --hval;

    /* Function Body */
    *flag__ = 0;
    *nnzh = 0;
} /* genpack_cc_mina_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int genpack_cc_mina_evalc__(n, x, ind, c__, flag__)
integer *n;
doublereal *x;
integer *ind;
doublereal *c__;
integer *flag__;
{
    extern doublereal belong_(), overlap_();

/*     This subroutine must compute the ind-th constraint of your problem. */
/*     For achieving this objective YOU MUST MOFIFY it according to your */
/*     problem. See below the places where your modifications must be */
/*     inserted. */

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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     COMMON BLOCKS */
/*     EXTERNAL FUNCTIONS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    if (1 <= *ind && *ind <= packdata_1.nite) {
	*c__ = belong_(n, &x[1], ind);
    } else if (*ind == packdata_1.nite + 1) {
	*c__ = overlap_(n, &x[1]);
    }
} /* genpack_cc_mina_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int genpack_cc_mina_evaljac__(n, x, ind, indjac, valjac, 
	nnzjac, flag__)
integer *n;
doublereal *x;
integer *ind, *indjac;
doublereal *valjac;
integer *nnzjac, *flag__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int goverlap_();
    static integer i__;
    extern /* Subroutine */ int gbelong_();

/*     This subroutine must compute the gradient of the constraint i. For */
/*     achieving these objective YOU MUST MODIFY it in the way specified */
/*     below. */

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
/*              the non-null value of the partial derivative of the ind-th */
/*              constraint with respect to the indjac(k)-th variable must */
/*              be saved at valjac(k). */

/*     flag     integer */
/*              You must set it to any number different of 0 (zero) if */
/*              some error ocurred during the evaluation of the */
/*              constraint. (For example, trying to compute the square */
/*              root of a negative number, dividing by zero or a very */
/*              small number, etc.) If everything was o.k. you must set */
/*              it equal to zero. */
/*     COMMON SCALARS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --valjac;
    --indjac;
    --x;

    /* Function Body */
    *flag__ = 0;
    if (1 <= *ind && *ind <= packdata_1.nite) {
	gbelong_(n, &x[1], ind, &indjac[1], &valjac[1], nnzjac);
    } else if (*ind == packdata_1.nite + 1) {
	goverlap_(n, &x[1], &valjac[1]);
	*nnzjac = *n - 1;
	i__1 = *nnzjac;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    indjac[i__] = i__;
	}
    }
} /* genpack_cc_mina_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int genpack_cc_mina_evalhc__(n, x, ind, hclin, hccol, hcval, 
	nnzhc, flag__)
integer *n;
doublereal *x;
integer *ind, *hclin, *hccol;
doublereal *hcval;
integer *nnzhc, *flag__;
{
    extern /* Subroutine */ int hoverlap_(), hbelong_();

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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;
    --hclin;
    --hccol;
    --hcval;

    /* Function Body */
    *flag__ = 0;
    if (1 <= *ind && *ind <= packdata_1.nite) {
	hbelong_(n, &x[1], ind, &hclin[1], &hccol[1], &hcval[1], nnzhc);
    } else if (*ind == packdata_1.nite + 1) {
	hoverlap_(n, &x[1], &hclin[1], &hccol[1], &hcval[1], nnzhc);
    }
} /* genpack_cc_mina_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int genpack_cc_mina_evalhlp__(n, x, m, lambda, p, hp, goth, 
	flag__)
integer *n;
doublereal *x;
integer *m;
doublereal *lambda, *p, *hp;
logical *goth;
integer *flag__;
{
/*     This subroutine might compute the product of the Hessian of the */
/*     Lagrangian times vector p. For achieving this objective YOU MAY */
/*     MODIFY it according to your problem. To modify this subroutine IS */
/*     NOT MANDATORY. See below where your modifications must be */
/*     inserted. */

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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     ****************************************************************** */
/*     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO COMPUTE */
/*     THE PRODUCT OF THE HESSIAN OF THE LAGRANGIAN TIMES A VECTOR */
/*     ****************************************************************** */
    /* Parameter adjustments */
    --hp;
    --p;
    --x;
    --lambda;

    /* Function Body */
    *flag__ = -1;
} /* genpack_cc_mina_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int genpack_cc_mina_endp__(n, x, l, u, m, lambda, rho, 
	equatn, linear)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
{
    /* Format strings */
    static char fmt_1[] = "(/,\002 Radius                        = \002,f12.\
8,/,\002 Maximum overlapping violation = \002,3x,1p,d9.1,/,\002 Maximum allo\
cation violation  = \002,3x,1p,d9.1)";

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static doublereal valoc, vover;
    extern doublereal maxaloc_(), maxover_();
    extern /* Subroutine */ int genpack_cc_mina_drawsol__(), 
	    genpack_cc_mina_savesol__();

    /* Fortran I/O blocks */
    static cilist io___12 = { 0, 6, 0, fmt_1, 0 };


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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     FUNCTIONS */
    /* Parameter adjustments */
    --u;
    --l;
    --x;
    --linear;
    --equatn;
    --rho;
    --lambda;

    /* Function Body */
    vover = maxover_(n, &x[1]);
    valoc = maxaloc_(n, &x[1]);
    s_wsfe(&io___12);
    do_fio(&c__1, (char *)&x[*n], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&vover, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&valoc, (ftnlen)sizeof(doublereal));
    e_wsfe();
/*     genpack_cc_mina_drawsol__(&packdata_1.nite, &packdata_1.iterad, n, &x[1],  */
/* 	    "solution.mp", (ftnlen)11); */
/*     genpack_cc_mina_savesol__(&packdata_1.nite, &packdata_1.iterad, & */
/* 	    packdata_1.ndim, n, &x[1], &vover, &valoc, "solution.dat", ( */
/* 	    ftnlen)12); */
} /* genpack_cc_mina_endp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
doublereal maxover_(n, x)
integer *n;
doublereal *x;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    extern /* Subroutine */ int classify_();
    static integer i__, j, k, l;
    extern doublereal mnrover_();

/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    classify_(n, &x[1]);
    ret_val = 0.;
    i__1 = partdata_1.nbr;
    for (k = 1; k <= i__1; ++k) {
	i__ = partdata_1.br[(k << 1) - 2];
	j = partdata_1.br[(k << 1) - 1];
	l = partdata_1.start[i__ + j * 10002];
L10:
	if (l != 0) {
/* Computing MAX */
	    d__1 = ret_val, d__2 = mnrover_(n, &x[1], &l, &partdata_1.next[l 
		    - 1]);
	    ret_val = max(d__1,d__2);
/* Computing MAX */
	    d__1 = ret_val, d__2 = mnrover_(n, &x[1], &l, &partdata_1.start[
		    i__ - 1 + (j + 1) * 10002]);
	    ret_val = max(d__1,d__2);
/* Computing MAX */
	    d__1 = ret_val, d__2 = mnrover_(n, &x[1], &l, &partdata_1.start[
		    i__ + (j + 1) * 10002]);
	    ret_val = max(d__1,d__2);
/* Computing MAX */
	    d__1 = ret_val, d__2 = mnrover_(n, &x[1], &l, &partdata_1.start[
		    i__ + 1 + (j + 1) * 10002]);
	    ret_val = max(d__1,d__2);
/* Computing MAX */
	    d__1 = ret_val, d__2 = mnrover_(n, &x[1], &l, &partdata_1.start[
		    i__ + 1 + j * 10002]);
	    ret_val = max(d__1,d__2);
	    l = partdata_1.next[l - 1];
	    goto L10;
	}
    }
    return ret_val;
} /* maxover_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
doublereal mnrover_(n, x, i__, lstart)
integer *n;
doublereal *x;
integer *i__, *lstart;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal dist;
    static integer j, k, ind1, ind2;

/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = 0.;
    j = *lstart;
L10:
    if (j != 0) {
	dist = 0.;
	i__1 = packdata_1.ndim;
	for (k = 1; k <= i__1; ++k) {
	    ind1 = (*i__ - 1) * packdata_1.ndim + k;
	    ind2 = (j - 1) * packdata_1.ndim + k;
/* Computing 2nd power */
	    d__1 = x[ind1] - x[ind2];
	    dist += d__1 * d__1;
	}
	dist = sqrt(dist);
/* Computing MAX */
	d__1 = ret_val, d__2 = packdata_1.iterad * 2. - dist;
	ret_val = max(d__1,d__2);
	j = partdata_1.next[j - 1];
	goto L10;
    }
    return ret_val;
} /* mnrover_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
doublereal maxaloc_(n, x)
integer *n;
doublereal *x;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal dist;
    static integer i__, k;

/*     COMMON SCALARS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = 0.;
    i__1 = packdata_1.nite;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dist = 0.;
	i__2 = packdata_1.ndim;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    d__1 = x[(i__ - 1) * packdata_1.ndim + k];
	    dist += d__1 * d__1;
	}
	dist = sqrt(dist);
/* Computing MAX */
	d__1 = ret_val, d__2 = dist - (x[*n] - packdata_1.iterad);
	ret_val = max(d__1,d__2);
    }
    return ret_val;
} /* maxaloc_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int genpack_cc_mina_drawsol__(nite, iterad, n, x, solfile, 
	solfile_len)
integer *nite;
doublereal *iterad;
integer *n;
doublereal *x;
char *solfile;
ftnlen solfile_len;
{
    /* Format strings */
    static char fmt_10[] = "(\002beginfig(1);\002/,\002u = \002,f20.10,\002 \
cm;\002)";
    static char fmt_20[] = "(\002draw fullcircle\002,/,\002     xscaled  \
\002,f20.10,\002u\002,/,\002     yscaled  \002,f20.10,\002u\002,/,\002     s\
hifted (\002,f20.10,\002u,\002,f20.10,\002u);\002)";
    static char fmt_30[] = "(\002endfig;\002,/,\002end;\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_wsfe(), do_fio(), e_wsfe(), f_clos();

    /* Local variables */
    static integer i__;
    static doublereal scale;

    /* Fortran I/O blocks */
    static cilist io___26 = { 0, 10, 0, fmt_10, 0 };
    static cilist io___27 = { 0, 10, 0, fmt_20, 0 };
    static cilist io___29 = { 0, 10, 0, fmt_20, 0 };
    static cilist io___30 = { 0, 10, 0, fmt_30, 0 };


/*     This subroutine generate a metapost file with the */
/*     graphical representation of the problem. */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = 11;
    o__1.ofnm = solfile;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/*     SCALING */
    scale = 10. / (x[*n] * 2.);
    s_wsfe(&io___26);
    do_fio(&c__1, (char *)&scale, (ftnlen)sizeof(doublereal));
    e_wsfe();
/*     CIRCULAR OBJECT */
    s_wsfe(&io___27);
    d__1 = x[*n] * 2.;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    d__2 = x[*n] * 2.;
    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&c_b33, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&c_b33, (ftnlen)sizeof(doublereal));
    e_wsfe();
/*     CIRCULAR ITEMS */
    i__1 = *nite;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsfe(&io___29);
	d__1 = *iterad * 2.;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	d__2 = *iterad * 2.;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&x[(i__ << 1) - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&x[i__ * 2], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    s_wsfe(&io___30);
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*     NON-EXECUTABLE STATEMENTS */
} /* genpack_cc_mina_drawsol__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int genpack_cc_mina_savesol__(nite, iterad, ndim, n, x, 
	vover, valoc, solfile, solfile_len)
integer *nite;
doublereal *iterad;
integer *ndim, *n;
doublereal *x, *vover, *valoc;
char *solfile;
ftnlen solfile_len;
{
    /* Format strings */
    static char fmt_10[] = "(/,\002Circular object radius        = \002,f16.\
10,/,\002Circular items radius         = \002,f16.10,/,\002Number of packed \
items        = \002,8x,i8,/,\002Maximum overlapping violation = \002,7x,1p,d\
9.1,/,\002Maximum allocation violation  = \002,7x,1p,d9.1,//,\002Items locat\
ion:\002,/)";
    static char fmt_20[] = "(i8,3(1x,1p,d16.8))";

    /* System generated locals */
    integer i__1, i__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_wsfe(), do_fio(), e_wsfe(), f_clos();

    /* Local variables */
    static integer i__, j;

    /* Fortran I/O blocks */
    static cilist io___31 = { 0, 10, 0, fmt_10, 0 };
    static cilist io___33 = { 0, 10, 0, fmt_20, 0 };


/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = 12;
    o__1.ofnm = solfile;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/*     PROBLEM DATA */
    s_wsfe(&io___31);
    do_fio(&c__1, (char *)&x[*n], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*iterad), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*nite), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*vover), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*valoc), (ftnlen)sizeof(doublereal));
    e_wsfe();
/*     ITEMS LOCATION */
    i__1 = *nite;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsfe(&io___33);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__2 = *ndim;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&x[(i__ - 1) * *ndim + j], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsfe();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*     NON-EXECUTABLE STATEMENTS */
} /* genpack_cc_mina_savesol__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
doublereal overlap_(n, x)
integer *n;
doublereal *x;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int classify_();
    static integer i__, j, k, l;
    extern doublereal nrdist_();

/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     CLASSIFY ITEMS BY REGIONS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    classify_(n, &x[1]);
/*     COMPUTE SPARSE OVERLAPPING */
    ret_val = 0.;
    i__1 = partdata_1.nbr;
    for (k = 1; k <= i__1; ++k) {
	i__ = partdata_1.br[(k << 1) - 2];
	j = partdata_1.br[(k << 1) - 1];
	l = partdata_1.start[i__ + j * 10002];
L10:
	if (l != 0) {
	    ret_val = ret_val + nrdist_(n, &x[1], &l, &partdata_1.next[l - 1])
		     + nrdist_(n, &x[1], &l, &partdata_1.start[i__ - 1 + (j + 
		    1) * 10002]) + nrdist_(n, &x[1], &l, &partdata_1.start[
		    i__ + (j + 1) * 10002]) + nrdist_(n, &x[1], &l, &
		    partdata_1.start[i__ + 1 + (j + 1) * 10002]) + nrdist_(n, 
		    &x[1], &l, &partdata_1.start[i__ + 1 + j * 10002]);
	    l = partdata_1.next[l - 1];
	    goto L10;
	}
    }
    return ret_val;
} /* overlap_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int goverlap_(n, x, g)
integer *n;
doublereal *x, *g;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int classify_();
    static integer i__, j, k, l;
    extern /* Subroutine */ int gnrdist_();

/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     CLASSIFY ITEMS BY REGIONS */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    classify_(n, &x[1]);
/*     COMPUTE SPARSE OVERLAPPING */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = 0.;
    }
    i__1 = partdata_1.nbr;
    for (k = 1; k <= i__1; ++k) {
	i__ = partdata_1.br[(k << 1) - 2];
	j = partdata_1.br[(k << 1) - 1];
	l = partdata_1.start[i__ + j * 10002];
L10:
	if (l != 0) {
	    gnrdist_(n, &x[1], &g[1], &l, &partdata_1.next[l - 1]);
	    gnrdist_(n, &x[1], &g[1], &l, &partdata_1.start[i__ - 1 + (j + 1) 
		    * 10002]);
	    gnrdist_(n, &x[1], &g[1], &l, &partdata_1.start[i__ + (j + 1) * 
		    10002]);
	    gnrdist_(n, &x[1], &g[1], &l, &partdata_1.start[i__ + 1 + (j + 1) 
		    * 10002]);
	    gnrdist_(n, &x[1], &g[1], &l, &partdata_1.start[i__ + 1 + j * 
		    10002]);
	    l = partdata_1.next[l - 1];
	    goto L10;
	}
    }
} /* goverlap_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hoverlap_(n, x, hlin, hcol, hval, nnzh)
integer *n;
doublereal *x;
integer *hlin, *hcol;
doublereal *hval;
integer *nnzh;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int classify_();
    static integer i__, j, k, l;
    extern /* Subroutine */ int hnrdist_();

/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     CLASSIFY ITEMS BY REGIONS */
    /* Parameter adjustments */
    --x;
    --hlin;
    --hcol;
    --hval;

    /* Function Body */
    classify_(n, &x[1]);
/*     INITALIZE DIAGONAL BLOCKS */
    i__1 = packdata_1.nite;
    for (k = 1; k <= i__1; ++k) {
	i__2 = packdata_1.ndim;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = packdata_1.ndim;
	    for (i__ = j; i__ <= i__3; ++i__) {
		diblocks_1.diagb[i__ + (j + (k << 1) << 1) - 7] = 0.;
	    }
	}
    }
/*     COMPUTE SPARSE OVERLAPPING SECOND DERIVATIVES */
    *nnzh = 0;
    i__1 = partdata_1.nbr;
    for (k = 1; k <= i__1; ++k) {
	i__ = partdata_1.br[(k << 1) - 2];
	j = partdata_1.br[(k << 1) - 1];
	l = partdata_1.start[i__ + j * 10002];
L10:
	if (l != 0) {
	    hnrdist_(n, &x[1], &hlin[1], &hcol[1], &hval[1], nnzh, &l, &
		    partdata_1.next[l - 1]);
	    hnrdist_(n, &x[1], &hlin[1], &hcol[1], &hval[1], nnzh, &l, &
		    partdata_1.start[i__ - 1 + (j + 1) * 10002]);
	    hnrdist_(n, &x[1], &hlin[1], &hcol[1], &hval[1], nnzh, &l, &
		    partdata_1.start[i__ + (j + 1) * 10002]);
	    hnrdist_(n, &x[1], &hlin[1], &hcol[1], &hval[1], nnzh, &l, &
		    partdata_1.start[i__ + 1 + (j + 1) * 10002]);
	    hnrdist_(n, &x[1], &hlin[1], &hcol[1], &hval[1], nnzh, &l, &
		    partdata_1.start[i__ + 1 + j * 10002]);
	    l = partdata_1.next[l - 1];
	    goto L10;
	}
    }
/*     ADD DIAGONAL BLOCKS TO THE SPARSE STRUCTURE */
    i__1 = packdata_1.nite;
    for (k = 1; k <= i__1; ++k) {
	i__2 = packdata_1.ndim;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = packdata_1.ndim;
	    for (i__ = j; i__ <= i__3; ++i__) {
		if (diblocks_1.diagb[i__ + (j + (k << 1) << 1) - 7] != 0.) {
		    ++(*nnzh);
		    hlin[*nnzh] = (k - 1) * packdata_1.ndim + i__;
		    hcol[*nnzh] = (k - 1) * packdata_1.ndim + j;
		    hval[*nnzh] = diblocks_1.diagb[i__ + (j + (k << 1) << 1) 
			    - 7];
		}
	    }
	}
    }
} /* hoverlap_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
doublereal nrdist_(n, x, i__, lstart)
integer *n;
doublereal *x;
integer *i__, *lstart;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static doublereal dist;
    static integer j, k, ind1, ind2;

/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = 0.;
    j = *lstart;
L10:
    if (j != 0) {
	dist = 0.;
	i__1 = packdata_1.ndim;
	for (k = 1; k <= i__1; ++k) {
	    ind1 = (*i__ - 1) * packdata_1.ndim + k;
	    ind2 = (j - 1) * packdata_1.ndim + k;
/* Computing 2nd power */
	    d__1 = x[ind1] - x[ind2];
	    dist += d__1 * d__1;
	}
/* Computing MAX */
/* Computing 2nd power */
	d__3 = packdata_1.iterad * 2.;
	d__1 = 0., d__2 = d__3 * d__3 - dist;
	dist = max(d__1,d__2);
/* Computing 2nd power */
	d__1 = dist;
	ret_val += d__1 * d__1;
	j = partdata_1.next[j - 1];
	goto L10;
    }
    return ret_val;
} /* nrdist_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int gnrdist_(n, x, g, i__, lstart)
integer *n;
doublereal *x, *g;
integer *i__, *lstart;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal dist;
    static integer j, k;
    static doublereal tmp;
    static integer ind1, ind2;

/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    j = *lstart;
L10:
    if (j != 0) {
	dist = 0.;
	i__1 = packdata_1.ndim;
	for (k = 1; k <= i__1; ++k) {
	    ind1 = (*i__ - 1) * packdata_1.ndim + k;
	    ind2 = (j - 1) * packdata_1.ndim + k;
/* Computing 2nd power */
	    d__1 = x[ind1] - x[ind2];
	    dist += d__1 * d__1;
	}
/* Computing MAX */
/* Computing 2nd power */
	d__3 = packdata_1.iterad * 2.;
	d__1 = 0., d__2 = d__3 * d__3 - dist;
	dist = max(d__1,d__2);
	if (dist != 0.) {
	    i__1 = packdata_1.ndim;
	    for (k = 1; k <= i__1; ++k) {
		ind1 = (*i__ - 1) * packdata_1.ndim + k;
		ind2 = (j - 1) * packdata_1.ndim + k;
		tmp = dist * 4. * (x[ind1] - x[ind2]);
		g[ind1] -= tmp;
		g[ind2] += tmp;
	    }
	}
	j = partdata_1.next[j - 1];
	goto L10;
    }
} /* gnrdist_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hnrdist_(n, x, hlin, hcol, hval, nnzh, i__, lstart)
integer *n;
doublereal *x;
integer *hlin, *hcol;
doublereal *hval;
integer *nnzh, *i__, *lstart;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal dist;
    static integer j, k, l, col, lin;
    static doublereal tmp;
    static integer ind1, ind2, ind3, ind4;

/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;
    --hlin;
    --hcol;
    --hval;

    /* Function Body */
    j = *lstart;
L10:
    if (j != 0) {
	if (j > *i__) {
	    col = *i__;
	    lin = j;
	} else {
	    col = j;
	    lin = *i__;
	}
	dist = 0.;
	i__1 = packdata_1.ndim;
	for (k = 1; k <= i__1; ++k) {
	    ind1 = (col - 1) * packdata_1.ndim + k;
	    ind2 = (lin - 1) * packdata_1.ndim + k;
/* Computing 2nd power */
	    d__1 = x[ind1] - x[ind2];
	    dist += d__1 * d__1;
	}
/* Computing 2nd power */
	d__1 = packdata_1.iterad * 2.;
	if (dist <= d__1 * d__1) {
	    i__1 = packdata_1.ndim;
	    for (k = 1; k <= i__1; ++k) {
		ind1 = (col - 1) * packdata_1.ndim + k;
		ind2 = (lin - 1) * packdata_1.ndim + k;
/* Computing 2nd power */
		d__1 = x[ind1] - x[ind2];
/* Computing 2nd power */
		d__2 = packdata_1.iterad * 2.;
		tmp = d__1 * d__1 * 8. - (d__2 * d__2 - dist) * 4.;
		if (tmp != 0.) {
/*                     H(ind1,ind1) = H(ind1,ind1) + tmp */
		    diblocks_1.diagb[k + (k + (col << 1) << 1) - 7] += tmp;
/*                     H(ind2,ind2) = H(ind2,ind2) + tmp */
		    diblocks_1.diagb[k + (k + (lin << 1) << 1) - 7] += tmp;
/*                     H(ind2,ind1) = H(ind2,ind1) - tmp */
		    ++(*nnzh);
		    hlin[*nnzh] = ind2;
		    hcol[*nnzh] = ind1;
		    hval[*nnzh] = -tmp;
		}
		i__2 = k - 1;
		for (l = 1; l <= i__2; ++l) {
		    ind3 = (col - 1) * packdata_1.ndim + l;
		    ind4 = (lin - 1) * packdata_1.ndim + l;
		    tmp = (x[ind3] - x[ind4]) * 8. * (x[ind1] - x[ind2]);
		    if (tmp != 0.) {
/*                         H(ind1,ind3) = H(ind1,ind3) + tmp */
			diblocks_1.diagb[k + (l + (col << 1) << 1) - 7] += 
				tmp;
/*                         H(ind2,ind4) = H(ind2,ind4) + tmp */
			diblocks_1.diagb[k + (l + (lin << 1) << 1) - 7] += 
				tmp;
/*                         H(ind2,ind3) = H(ind2,ind3) - tmp */
			++(*nnzh);
			hlin[*nnzh] = ind2;
			hcol[*nnzh] = ind3;
			hval[*nnzh] = -tmp;
		    }
		}
		i__2 = packdata_1.ndim;
		for (l = k + 1; l <= i__2; ++l) {
		    ind3 = (col - 1) * packdata_1.ndim + l;
		    ind4 = (lin - 1) * packdata_1.ndim + l;
		    tmp = (x[ind3] - x[ind4]) * 8. * (x[ind1] - x[ind2]);
		    if (tmp != 0.) {
/*                         H(ind2,ind3) = H(ind2,ind3) - tmp */
			++(*nnzh);
			hlin[*nnzh] = ind2;
			hcol[*nnzh] = ind3;
			hval[*nnzh] = -tmp;
		    }
		}
	    }
	}
	j = partdata_1.next[j - 1];
	goto L10;
    }
} /* hnrdist_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int classify_(n, x)
integer *n;
doublereal *x;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k;
    extern /* Subroutine */ int region_();

/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
/*     CLEAN-UP THE START STRUCTURE */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    i__1 = partdata_1.nbr;
    for (k = 1; k <= i__1; ++k) {
	i__ = partdata_1.br[(k << 1) - 2];
	j = partdata_1.br[(k << 1) - 1];
	partdata_1.start[i__ + j * 10002] = 0;
    }
/*     FILL-IN THE START STRUCTURE AGAIN */
    partdata_1.nbr = 0;
    i__1 = packdata_1.nite;
    for (k = 1; k <= i__1; ++k) {
	region_(&partdata_1.nreg, &x[(k - 1) * packdata_1.ndim + 1], &x[(k - 
		1) * packdata_1.ndim + 2], &i__, &j);
	if (partdata_1.start[i__ + j * 10002] == 0) {
	    ++partdata_1.nbr;
	    partdata_1.br[(partdata_1.nbr << 1) - 2] = i__;
	    partdata_1.br[(partdata_1.nbr << 1) - 1] = j;
	}
	partdata_1.next[k - 1] = partdata_1.start[i__ + j * 10002];
	partdata_1.start[i__ + j * 10002] = k;
    }
} /* classify_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int region_(nreg, x, y, i__, j)
integer *nreg;
doublereal *x, *y;
integer *i__, *j;
{
    /* System generated locals */
    integer i__1;

/*     COMMON SCALARS */
/*     SCALAR ARGUMENTS */
/*     COMMON BLOCKS */
    *i__ = (integer) ((*x + *nreg * packdata_1.iterad) / (packdata_1.iterad * 
	    2.)) + 1;
/* Computing MIN */
    i__1 = max(1,*i__);
    *i__ = min(i__1,*nreg);
    *j = (integer) ((*y + *nreg * packdata_1.iterad) / (packdata_1.iterad * 
	    2.)) + 1;
/* Computing MIN */
    i__1 = max(1,*j);
    *j = min(i__1,*nreg);
} /* region_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
integer ceil_(x)
doublereal *x;
{
    /* System generated locals */
    integer ret_val;

/*     SCALAR ARGUMENTS */
    if (*x == (doublereal) ((integer) (*x))) {
	ret_val = (integer) (*x);
    } else {
	ret_val = (integer) (*x + 1);
    }
    return ret_val;
} /* ceil_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
doublereal belong_(n, x, ind)
integer *n;
doublereal *x;
integer *ind;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__;

/*     COMMON SCALARS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;

    /* Function Body */
/* Computing 2nd power */
    d__1 = x[*n] - packdata_1.iterad;
    ret_val = -(d__1 * d__1);
    i__1 = packdata_1.ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = x[packdata_1.ndim * (*ind - 1) + i__];
	ret_val += d__1 * d__1;
    }
    return ret_val;
} /* belong_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int gbelong_(n, x, ind, gind, gval, nnzg)
integer *n;
doublereal *x;
integer *ind, *gind;
doublereal *gval;
integer *nnzg;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/*     COMMON SCALARS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --gval;
    --gind;
    --x;

    /* Function Body */
    *nnzg = packdata_1.ndim + 1;
    i__1 = packdata_1.ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gind[i__] = packdata_1.ndim * (*ind - 1) + i__;
	gval[i__] = x[packdata_1.ndim * (*ind - 1) + i__] * 2.;
    }
    gind[*nnzg] = *n;
    gval[*nnzg] = (x[*n] - packdata_1.iterad) * -2.;
} /* gbelong_ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int hbelong_(n, x, ind, hlin, hcol, hval, nnzh)
integer *n;
doublereal *x;
integer *ind, *hlin, *hcol;
doublereal *hval;
integer *nnzh;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;

/*     COMMON SCALARS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
/*     COMMON BLOCKS */
    /* Parameter adjustments */
    --x;
    --hlin;
    --hcol;
    --hval;

    /* Function Body */
    *nnzh = packdata_1.ndim + 1;
    i__1 = packdata_1.ndim;
    for (k = 1; k <= i__1; ++k) {
	hlin[k] = packdata_1.ndim * (*ind - 1) + k;
	hcol[k] = packdata_1.ndim * (*ind - 1) + k;
	hval[k] = 2.;
    }
    hlin[*nnzh] = *n;
    hcol[*nnzh] = *n;
    hval[*nnzh] = -2.;
} /* hbelong_ */

