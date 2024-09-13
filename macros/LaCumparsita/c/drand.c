/* drand.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b2 = 1.;

doublereal drand_(ix)
doublereal *ix;
{
    /* Initialized data */

    static doublereal a = 16807.;
    static doublereal b15 = 32768.;
    static doublereal b16 = 65536.;
    static doublereal p = 2147483647.;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double d_mod();

    /* Local variables */
    static doublereal xalo, k, leftlo, fhi, xhi;

    xhi = *ix / b16;
    xhi -= d_mod(&xhi, &c_b2);
    xalo = (*ix - xhi * b16) * a;
    leftlo = xalo / b16;
    leftlo -= d_mod(&leftlo, &c_b2);
    fhi = xhi * a + leftlo;
    k = fhi / b15;
    k -= d_mod(&k, &c_b2);
    *ix = xalo - leftlo * b16 - p + (fhi - k * b15) * b16 + k;
    if (*ix < 0.) {
	*ix += p;
    }
    ret_val = *ix * 4.656612875e-10;
    return ret_val;
} /* drand_ */

