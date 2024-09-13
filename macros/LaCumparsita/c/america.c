/* america.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal put[264]	/* was [132][2] */, area[17];
    integer nbord[17], nplin[2244]	/* was [17][132] */, nstates, npun;
} probdata_;

#define probdata_1 probdata_

/* Table of constant values */

static logical c_false = FALSE_;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;

/*     ================================================================= */
/*     File: america.f */
/*     ================================================================= */
/*     ================================================================= */
/*     Module: Subroutines that define the problem */
/*     ================================================================= */
/*     Last update of any of the component of this module: */
/*     May 20th, 2005. */
/*     Users are encouraged to download periodically updated versions of */
/*     this code at the COLLECTION home page: */

/*     www.ime.usp.br/~egbirgin/collection/ */

/*     and periodically updated versions of the TANGO Project solvers at */
/*     the TANGO home page: */

/*     www.ime.usp.br/~egbirgin/tango/ */
/*     ================================================================= */
/*     America problem */
/*     --------------- */
/*     The objective is to draw a map of America in which the countries */
/*     appear with areas that are proportional to their real values. The */
/*     unknowns are $132$ points in $\R^2$, which define the boundaries */
/*     of $17$ countries. Each point is assigned to the boundary of one */
/*     or more countries. The computed area of each country is calculated */
/*     as a function of its boundary points using Green's formula. The */
/*     constraints of the problem are: */

/*                   Computed area = True area */

/*     for each country. As initial approximation we took the coordinates */
/*     of the~$132$ points in the New York Times map of America which, of */
/*     course, do not fit with true areas. The objective function is */
/*     $\frac{1}{2}\sum_{j=1}^{132} \|P_j - Q_j\|_2^2$, where */
/*     $P_1,\ldots,P_{132}$ are the unknowns and $Q_1,\ldots,Q_{132}$ are */
/*     the initial approximation. */

/*     Considered "countries": */

/*      1 Cuba */
/*      2 Canada */
/*      3 USA */
/*      4 Mexico */
/*      5 America Central */
/*      6 Colombia */
/*      7 Venezuela */
/*      8 Guianas */
/*      9 Brasil */
/*     10 Ecuador */
/*     11 Peru */
/*     12 Bolivia */
/*     13 Paraguay */
/*     14 Chile */
/*     15 Argentina */
/*     16 Alaska */
/*     17 Uruguay */
/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int america_inip__(n, x, l, u, m, lambda, rho, equatn, 
	linear)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern doublereal america_comparea__();
    static integer i__;
    static doublereal falsearea, totarea;
    extern /* Subroutine */ int america_drawsol__();

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
/*     EXTERNAL FUNCTIONS */
/*     ****************************************************************** */
/*     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR PROBLEM */
/*     DATA: */
/*     ****************************************************************** */
/*     Set problem data */
/*     Number of states and number of points */
    /* Parameter adjustments */
    --linear;
    --equatn;
    --rho;
    --lambda;
    --u;
    --l;
    --x;

    /* Function Body */
    probdata_1.nstates = 17;
    probdata_1.npun = 132;
/*     Number of boundary points for each state and indices of the */
/*     boundary points */
    probdata_1.nbord[0] = 4;
    probdata_1.nplin[0] = 49;
    probdata_1.nplin[17] = 50;
    probdata_1.nplin[34] = 51;
    probdata_1.nplin[51] = 52;
    probdata_1.nbord[1] = 21;
    probdata_1.nplin[1] = 83;
    probdata_1.nplin[18] = 35;
    probdata_1.nplin[35] = 64;
    probdata_1.nplin[52] = 63;
    probdata_1.nplin[69] = 84;
    probdata_1.nplin[86] = 36;
    probdata_1.nplin[103] = 99;
    probdata_1.nplin[120] = 100;
    probdata_1.nplin[137] = 37;
    probdata_1.nplin[154] = 38;
    probdata_1.nplin[171] = 39;
    probdata_1.nplin[188] = 101;
    probdata_1.nplin[205] = 102;
    probdata_1.nplin[222] = 40;
    probdata_1.nplin[239] = 41;
    probdata_1.nplin[256] = 42;
    probdata_1.nplin[273] = 127;
    probdata_1.nplin[290] = 81;
    probdata_1.nplin[307] = 82;
    probdata_1.nplin[324] = 43;
    probdata_1.nplin[341] = 46;
    probdata_1.nbord[2] = 22;
    probdata_1.nplin[2] = 33;
    probdata_1.nplin[19] = 108;
    probdata_1.nplin[36] = 109;
    probdata_1.nplin[53] = 110;
    probdata_1.nplin[70] = 111;
    probdata_1.nplin[87] = 112;
    probdata_1.nplin[104] = 34;
    probdata_1.nplin[121] = 94;
    probdata_1.nplin[138] = 62;
    probdata_1.nplin[155] = 98;
    probdata_1.nplin[172] = 61;
    probdata_1.nplin[189] = 60;
    probdata_1.nplin[206] = 96;
    probdata_1.nplin[223] = 97;
    probdata_1.nplin[240] = 36;
    probdata_1.nplin[257] = 84;
    probdata_1.nplin[274] = 63;
    probdata_1.nplin[291] = 64;
    probdata_1.nplin[308] = 35;
    probdata_1.nplin[325] = 120;
    probdata_1.nplin[342] = 56;
    probdata_1.nplin[359] = 121;
    probdata_1.nbord[3] = 26;
    probdata_1.nplin[3] = 86;
    probdata_1.nplin[20] = 113;
    probdata_1.nplin[37] = 114;
    probdata_1.nplin[54] = 31;
    probdata_1.nplin[71] = 32;
    probdata_1.nplin[88] = 87;
    probdata_1.nplin[105] = 115;
    probdata_1.nplin[122] = 88;
    probdata_1.nplin[139] = 103;
    probdata_1.nplin[156] = 89;
    probdata_1.nplin[173] = 85;
    probdata_1.nplin[190] = 34;
    probdata_1.nplin[207] = 112;
    probdata_1.nplin[224] = 111;
    probdata_1.nplin[241] = 110;
    probdata_1.nplin[258] = 109;
    probdata_1.nplin[275] = 108;
    probdata_1.nplin[292] = 33;
    probdata_1.nplin[309] = 118;
    probdata_1.nplin[326] = 104;
    probdata_1.nplin[343] = 116;
    probdata_1.nplin[360] = 119;
    probdata_1.nplin[377] = 105;
    probdata_1.nplin[394] = 106;
    probdata_1.nplin[411] = 117;
    probdata_1.nplin[428] = 107;
    probdata_1.nbord[4] = 8;
    probdata_1.nplin[4] = 30;
    probdata_1.nplin[21] = 29;
    probdata_1.nplin[38] = 91;
    probdata_1.nplin[55] = 90;
    probdata_1.nplin[72] = 32;
    probdata_1.nplin[89] = 31;
    probdata_1.nplin[106] = 92;
    probdata_1.nplin[123] = 93;
    probdata_1.nbord[5] = 9;
    probdata_1.nplin[5] = 28;
    probdata_1.nplin[22] = 25;
    probdata_1.nplin[39] = 23;
    probdata_1.nplin[56] = 22;
    probdata_1.nplin[73] = 55;
    probdata_1.nplin[90] = 95;
    probdata_1.nplin[107] = 24;
    probdata_1.nplin[124] = 29;
    probdata_1.nplin[141] = 30;
    probdata_1.nbord[6] = 9;
    probdata_1.nplin[6] = 22;
    probdata_1.nplin[23] = 21;
    probdata_1.nplin[40] = 20;
    probdata_1.nplin[57] = 19;
    probdata_1.nplin[74] = 122;
    probdata_1.nplin[91] = 123;
    probdata_1.nplin[108] = 24;
    probdata_1.nplin[125] = 95;
    probdata_1.nplin[142] = 55;
    probdata_1.nbord[7] = 5;
    probdata_1.nplin[7] = 20;
    probdata_1.nplin[24] = 126;
    probdata_1.nplin[41] = 18;
    probdata_1.nplin[58] = 17;
    probdata_1.nplin[75] = 19;
    probdata_1.nbord[8] = 22;
    probdata_1.nplin[8] = 53;
    probdata_1.nplin[25] = 16;
    probdata_1.nplin[42] = 15;
    probdata_1.nplin[59] = 10;
    probdata_1.nplin[76] = 8;
    probdata_1.nplin[93] = 7;
    probdata_1.nplin[110] = 6;
    probdata_1.nplin[127] = 65;
    probdata_1.nplin[144] = 5;
    probdata_1.nplin[161] = 68;
    probdata_1.nplin[178] = 76;
    probdata_1.nplin[195] = 48;
    probdata_1.nplin[212] = 129;
    probdata_1.nplin[229] = 47;
    probdata_1.nplin[246] = 125;
    probdata_1.nplin[263] = 17;
    probdata_1.nplin[280] = 18;
    probdata_1.nplin[297] = 126;
    probdata_1.nplin[314] = 20;
    probdata_1.nplin[331] = 21;
    probdata_1.nplin[348] = 22;
    probdata_1.nplin[365] = 23;
    probdata_1.nbord[9] = 4;
    probdata_1.nplin[9] = 27;
    probdata_1.nplin[26] = 26;
    probdata_1.nplin[43] = 25;
    probdata_1.nplin[60] = 28;
    probdata_1.nbord[10] = 10;
    probdata_1.nplin[10] = 128;
    probdata_1.nplin[27] = 59;
    probdata_1.nplin[44] = 12;
    probdata_1.nplin[61] = 13;
    probdata_1.nplin[78] = 16;
    probdata_1.nplin[95] = 53;
    probdata_1.nplin[112] = 23;
    probdata_1.nplin[129] = 25;
    probdata_1.nplin[146] = 26;
    probdata_1.nplin[163] = 27;
    probdata_1.nbord[11] = 7;
    probdata_1.nplin[11] = 14;
    probdata_1.nplin[28] = 9;
    probdata_1.nplin[45] = 10;
    probdata_1.nplin[62] = 15;
    probdata_1.nplin[79] = 16;
    probdata_1.nplin[96] = 13;
    probdata_1.nplin[113] = 12;
    probdata_1.nbord[12] = 6;
    probdata_1.nplin[12] = 54;
    probdata_1.nplin[29] = 7;
    probdata_1.nplin[46] = 8;
    probdata_1.nplin[63] = 10;
    probdata_1.nplin[80] = 9;
    probdata_1.nplin[97] = 70;
    probdata_1.nbord[13] = 16;
    probdata_1.nplin[13] = 1;
    probdata_1.nplin[30] = 11;
    probdata_1.nplin[47] = 132;
    probdata_1.nplin[64] = 131;
    probdata_1.nplin[81] = 14;
    probdata_1.nplin[98] = 12;
    probdata_1.nplin[115] = 59;
    probdata_1.nplin[132] = 130;
    probdata_1.nplin[149] = 57;
    probdata_1.nplin[166] = 58;
    probdata_1.nplin[183] = 2;
    probdata_1.nplin[200] = 73;
    probdata_1.nplin[217] = 77;
    probdata_1.nplin[234] = 71;
    probdata_1.nplin[251] = 75;
    probdata_1.nplin[268] = 74;
    probdata_1.nbord[14] = 20;
    probdata_1.nplin[14] = 1;
    probdata_1.nplin[31] = 2;
    probdata_1.nplin[48] = 73;
    probdata_1.nplin[65] = 71;
    probdata_1.nplin[82] = 72;
    probdata_1.nplin[99] = 75;
    probdata_1.nplin[116] = 74;
    probdata_1.nplin[133] = 69;
    probdata_1.nplin[150] = 67;
    probdata_1.nplin[167] = 66;
    probdata_1.nplin[184] = 4;
    probdata_1.nplin[201] = 6;
    probdata_1.nplin[218] = 7;
    probdata_1.nplin[235] = 54;
    probdata_1.nplin[252] = 70;
    probdata_1.nplin[269] = 9;
    probdata_1.nplin[286] = 14;
    probdata_1.nplin[303] = 131;
    probdata_1.nplin[320] = 132;
    probdata_1.nplin[337] = 11;
    probdata_1.nbord[15] = 8;
    probdata_1.nplin[15] = 43;
    probdata_1.nplin[32] = 80;
    probdata_1.nplin[49] = 44;
    probdata_1.nplin[66] = 78;
    probdata_1.nplin[83] = 124;
    probdata_1.nplin[100] = 45;
    probdata_1.nplin[117] = 79;
    probdata_1.nplin[134] = 46;
    probdata_1.nbord[16] = 5;
    probdata_1.nplin[16] = 4;
    probdata_1.nplin[33] = 3;
    probdata_1.nplin[50] = 5;
    probdata_1.nplin[67] = 65;
    probdata_1.nplin[84] = 6;
/*     Real area of each state */
    probdata_1.area[0] = 115.;
    probdata_1.area[1] = 9922.;
    probdata_1.area[2] = 7885.;
    probdata_1.area[3] = 1973.;
    probdata_1.area[4] = 543.;
    probdata_1.area[5] = 1139.;
    probdata_1.area[6] = 912.;
    probdata_1.area[7] = 470.;
    probdata_1.area[8] = 8512.;
    probdata_1.area[9] = 461.;
    probdata_1.area[10] = 1285.;
    probdata_1.area[11] = 1099.;
    probdata_1.area[12] = 407.;
    probdata_1.area[13] = 752.;
    probdata_1.area[14] = 2767.;
    probdata_1.area[15] = 1478.;
    probdata_1.area[16] = 187.;
/*     Putative coordinates of each boundary point (initial approximation) */
    probdata_1.put[0] = 11.7;
    probdata_1.put[132] = 5.7;
    probdata_1.put[1] = 12.1;
    probdata_1.put[133] = 5.6;
    probdata_1.put[2] = 13.5;
    probdata_1.put[134] = 8.3;
    probdata_1.put[3] = 13.1;
    probdata_1.put[135] = 8.4;
    probdata_1.put[4] = 13.7;
    probdata_1.put[136] = 8.6;
    probdata_1.put[5] = 13.2;
    probdata_1.put[137] = 9.;
    probdata_1.put[6] = 13.6;
    probdata_1.put[138] = 9.4;
    probdata_1.put[7] = 13.4;
    probdata_1.put[139] = 10.;
    probdata_1.put[8] = 12.5;
    probdata_1.put[140] = 10.;
    probdata_1.put[9] = 13.1;
    probdata_1.put[141] = 10.4;
    probdata_1.put[10] = 11.7;
    probdata_1.put[142] = 6.4;
    probdata_1.put[11] = 11.6;
    probdata_1.put[143] = 10.6;
    probdata_1.put[12] = 11.6;
    probdata_1.put[144] = 11.2;
    probdata_1.put[13] = 12.;
    probdata_1.put[145] = 9.9;
    probdata_1.put[14] = 12.7;
    probdata_1.put[146] = 11.2;
    probdata_1.put[15] = 11.8;
    probdata_1.put[147] = 11.6;
    probdata_1.put[16] = 13.8;
    probdata_1.put[148] = 13.5;
    probdata_1.put[17] = 13.5;
    probdata_1.put[149] = 13.2;
    probdata_1.put[18] = 12.9;
    probdata_1.put[150] = 14.1;
    probdata_1.put[19] = 12.8;
    probdata_1.put[151] = 13.7;
    probdata_1.put[20] = 12.4;
    probdata_1.put[152] = 13.3;
    probdata_1.put[21] = 11.9;
    probdata_1.put[153] = 13.3;
    probdata_1.put[22] = 11.7;
    probdata_1.put[154] = 12.7;
    probdata_1.put[23] = 11.6;
    probdata_1.put[155] = 14.6;
    probdata_1.put[24] = 11.1;
    probdata_1.put[156] = 12.8;
    probdata_1.put[25] = 10.8;
    probdata_1.put[157] = 12.5;
    probdata_1.put[26] = 10.5;
    probdata_1.put[158] = 12.5;
    probdata_1.put[27] = 10.7;
    probdata_1.put[159] = 13.2;
    probdata_1.put[28] = 10.9;
    probdata_1.put[160] = 14.2;
    probdata_1.put[29] = 10.8;
    probdata_1.put[161] = 14.1;
    probdata_1.put[30] = 9.2;
    probdata_1.put[162] = 15.;
    probdata_1.put[31] = 9.4;
    probdata_1.put[163] = 15.3;
    probdata_1.put[32] = 6.6;
    probdata_1.put[164] = 17.6;
    probdata_1.put[33] = 8.4;
    probdata_1.put[165] = 16.6;
    probdata_1.put[34] = 6.4;
    probdata_1.put[166] = 20.;
    probdata_1.put[35] = 12.5;
    probdata_1.put[167] = 19.5;
    probdata_1.put[36] = 13.8;
    probdata_1.put[168] = 20.6;
    probdata_1.put[37] = 11.8;
    probdata_1.put[169] = 22.4;
    probdata_1.put[38] = 11.2;
    probdata_1.put[170] = 20.5;
    probdata_1.put[39] = 10.;
    probdata_1.put[171] = 22.2;
    probdata_1.put[40] = 11.5;
    probdata_1.put[172] = 23.;
    probdata_1.put[41] = 11.;
    probdata_1.put[173] = 24.6;
    probdata_1.put[42] = 5.8;
    probdata_1.put[174] = 23.8;
    probdata_1.put[43] = 3.5;
    probdata_1.put[175] = 23.8;
    probdata_1.put[44] = 3.;
    probdata_1.put[176] = 21.5;
    probdata_1.put[45] = 5.5;
    probdata_1.put[177] = 22.;
    probdata_1.put[46] = 15.7;
    probdata_1.put[178] = 12.2;
    probdata_1.put[47] = 15.3;
    probdata_1.put[179] = 10.8;
    probdata_1.put[48] = 10.4;
    probdata_1.put[180] = 16.2;
    probdata_1.put[49] = 11.;
    probdata_1.put[181] = 15.8;
    probdata_1.put[50] = 11.3;
    probdata_1.put[182] = 16.;
    probdata_1.put[51] = 10.4;
    probdata_1.put[183] = 16.3;
    probdata_1.put[52] = 11.6;
    probdata_1.put[184] = 12.;
    probdata_1.put[53] = 13.;
    probdata_1.put[185] = 9.3;
    probdata_1.put[54] = 12.;
    probdata_1.put[186] = 13.8;
    probdata_1.put[55] = 6.;
    probdata_1.put[187] = 18.7;
    probdata_1.put[56] = 11.5;
    probdata_1.put[188] = 5.7;
    probdata_1.put[57] = 11.8;
    probdata_1.put[189] = 5.5;
    probdata_1.put[58] = 11.4;
    probdata_1.put[190] = 10.6;
    probdata_1.put[59] = 10.7;
    probdata_1.put[191] = 17.4;
    probdata_1.put[60] = 10.6;
    probdata_1.put[192] = 16.6;
    probdata_1.put[61] = 10.4;
    probdata_1.put[193] = 17.2;
    probdata_1.put[62] = 11.;
    probdata_1.put[194] = 19.1;
    probdata_1.put[63] = 10.;
    probdata_1.put[195] = 20.;
    probdata_1.put[64] = 13.5;
    probdata_1.put[196] = 8.8;
    probdata_1.put[65] = 13.2;
    probdata_1.put[197] = 7.7;
    probdata_1.put[66] = 12.6;
    probdata_1.put[198] = 7.7;
    probdata_1.put[67] = 14.2;
    probdata_1.put[199] = 9.7;
    probdata_1.put[68] = 12.3;
    probdata_1.put[200] = 6.7;
    probdata_1.put[69] = 13.1;
    probdata_1.put[201] = 9.6;
    probdata_1.put[70] = 12.2;
    probdata_1.put[202] = 5.2;
    probdata_1.put[71] = 12.5;
    probdata_1.put[203] = 5.3;
    probdata_1.put[72] = 12.2;
    probdata_1.put[204] = 5.4;
    probdata_1.put[73] = 12.1;
    probdata_1.put[205] = 5.6;
    probdata_1.put[74] = 12.2;
    probdata_1.put[206] = 5.4;
    probdata_1.put[75] = 15.;
    probdata_1.put[207] = 10.;
    probdata_1.put[76] = 11.9;
    probdata_1.put[208] = 5.3;
    probdata_1.put[77] = 2.7;
    probdata_1.put[209] = 23.;
    probdata_1.put[78] = 4.;
    probdata_1.put[210] = 22.;
    probdata_1.put[79] = 4.5;
    probdata_1.put[211] = 24.;
    probdata_1.put[80] = 9.;
    probdata_1.put[212] = 23.5;
    probdata_1.put[81] = 7.3;
    probdata_1.put[213] = 23.7;
    probdata_1.put[82] = 6.;
    probdata_1.put[214] = 21.;
    probdata_1.put[83] = 12.;
    probdata_1.put[215] = 19.5;
    probdata_1.put[84] = 8.5;
    probdata_1.put[216] = 16.;
    probdata_1.put[85] = 7.7;
    probdata_1.put[217] = 15.8;
    probdata_1.put[86] = 9.3;
    probdata_1.put[218] = 15.5;
    probdata_1.put[87] = 9.9;
    probdata_1.put[219] = 16.;
    probdata_1.put[88] = 9.;
    probdata_1.put[220] = 15.6;
    probdata_1.put[89] = 10.2;
    probdata_1.put[221] = 15.3;
    probdata_1.put[90] = 10.2;
    probdata_1.put[222] = 14.5;
    probdata_1.put[91] = 9.7;
    probdata_1.put[223] = 14.7;
    probdata_1.put[92] = 10.;
    probdata_1.put[224] = 14.4;
    probdata_1.put[93] = 9.3;
    probdata_1.put[225] = 17.2;
    probdata_1.put[94] = 11.5;
    probdata_1.put[226] = 14.;
    probdata_1.put[95] = 11.3;
    probdata_1.put[227] = 17.8;
    probdata_1.put[96] = 11.4;
    probdata_1.put[228] = 18.4;
    probdata_1.put[97] = 10.5;
    probdata_1.put[229] = 16.6;
    probdata_1.put[98] = 13.1;
    probdata_1.put[230] = 19.5;
    probdata_1.put[99] = 12.4;
    probdata_1.put[231] = 19.9;
    probdata_1.put[100] = 11.2;
    probdata_1.put[232] = 20.2;
    probdata_1.put[101] = 10.;
    probdata_1.put[233] = 21.5;
    probdata_1.put[102] = 9.5;
    probdata_1.put[234] = 16.;
    probdata_1.put[103] = 7.1;
    probdata_1.put[235] = 16.3;
    probdata_1.put[104] = 6.8;
    probdata_1.put[236] = 17.5;
    probdata_1.put[105] = 6.9;
    probdata_1.put[237] = 17.5;
    probdata_1.put[106] = 7.8;
    probdata_1.put[238] = 15.9;
    probdata_1.put[107] = 6.8;
    probdata_1.put[239] = 17.6;
    probdata_1.put[108] = 7.4;
    probdata_1.put[240] = 17.4;
    probdata_1.put[109] = 7.7;
    probdata_1.put[241] = 17.5;
    probdata_1.put[110] = 8.1;
    probdata_1.put[242] = 17.1;
    probdata_1.put[111] = 8.3;
    probdata_1.put[243] = 17.3;
    probdata_1.put[112] = 8.5;
    probdata_1.put[244] = 15.3;
    probdata_1.put[113] = 8.9;
    probdata_1.put[245] = 15.3;
    probdata_1.put[114] = 9.7;
    probdata_1.put[246] = 15.6;
    probdata_1.put[115] = 7.2;
    probdata_1.put[247] = 16.3;
    probdata_1.put[116] = 7.4;
    probdata_1.put[248] = 16.6;
    probdata_1.put[117] = 6.8;
    probdata_1.put[249] = 16.8;
    probdata_1.put[118] = 7.;
    probdata_1.put[250] = 16.8;
    probdata_1.put[119] = 6.1;
    probdata_1.put[251] = 19.5;
    probdata_1.put[120] = 6.2;
    probdata_1.put[252] = 18.2;
    probdata_1.put[121] = 12.5;
    probdata_1.put[253] = 14.5;
    probdata_1.put[122] = 12.1;
    probdata_1.put[254] = 14.5;
    probdata_1.put[123] = 2.75;
    probdata_1.put[255] = 22.15;
    probdata_1.put[124] = 14.6;
    probdata_1.put[256] = 12.8;
    probdata_1.put[125] = 13.;
    probdata_1.put[257] = 13.2;
    probdata_1.put[126] = 10.2;
    probdata_1.put[258] = 23.3;
    probdata_1.put[127] = 10.7;
    probdata_1.put[259] = 11.25;
    probdata_1.put[128] = 15.3;
    probdata_1.put[260] = 11.3;
    probdata_1.put[129] = 11.4;
    probdata_1.put[261] = 7.5;
    probdata_1.put[130] = 11.65;
    probdata_1.put[262] = 9.;
    probdata_1.put[131] = 11.6;
    probdata_1.put[263] = 7.25;
/*     Number of variables */
    *n = probdata_1.npun << 1;
/*     Initial point */
    i__1 = probdata_1.npun;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[(i__ << 1) - 1] = probdata_1.put[i__ - 1];
	x[i__ * 2] = probdata_1.put[i__ + 131];
    }
/*     america_drawsol__(n, &x[1], &c_false, "america-inisol.tex", (ftnlen)18); */
/*     Scale areas */
    totarea = 0.;
    i__1 = probdata_1.nstates;
    for (i__ = 1; i__ <= i__1; ++i__) {
	totarea += probdata_1.area[i__ - 1];
    }
    falsearea = 0.;
    i__1 = probdata_1.nstates;
    for (i__ = 1; i__ <= i__1; ++i__) {
	falsearea += america_comparea__(n, &x[1], &i__);
    }
    i__1 = probdata_1.nstates;
    for (i__ = 1; i__ <= i__1; ++i__) {
	probdata_1.area[i__ - 1] = probdata_1.area[i__ - 1] * falsearea / 
		totarea;
    }
/*     Lower and upper bounds */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l[i__] = -1e3;
	u[i__] = 1e3;
    }
/*     Number of constraints (equalities plus inequalities) */
    *m = probdata_1.nstates;
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
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE INIP. */
/*     ****************************************************************** */
} /* america_inip__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int america_evalf__(n, x, f, flag__)
integer *n;
doublereal *x, *f;
integer *flag__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;

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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
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
    i__1 = probdata_1.npun;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = x[(i__ << 1) - 1] - probdata_1.put[i__ - 1];
/* Computing 2nd power */
	d__2 = x[i__ * 2] - probdata_1.put[i__ + 131];
	*f = *f + d__1 * d__1 + d__2 * d__2;
    }
    *f *= .5;
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALF. */
/*     ****************************************************************** */
} /* america_evalf__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int america_evalg__(n, x, g, flag__)
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
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
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
    i__1 = probdata_1.npun;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[(i__ << 1) - 1] = x[(i__ << 1) - 1] - probdata_1.put[i__ - 1];
	g[i__ * 2] = x[i__ * 2] - probdata_1.put[i__ + 131];
    }
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALG. */
/*     ****************************************************************** */
} /* america_evalg__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int america_evalh__(n, x, hlin, hcol, hval, nnzh, flag__)
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
} /* america_evalh__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int america_evalc__(n, x, ind, c__, flag__)
integer *n;
doublereal *x;
integer *ind;
doublereal *c__;
integer *flag__;
{
    extern doublereal america_comparea__();

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
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     COMMON BLOCKS */
/*     EXTERNAL FUNCTIONS */
/*     Constraints */
/*     ****************************************************************** */
/*     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR */
/*     CONSTRAINTS: */
/*     ****************************************************************** */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *flag__ = 0;
    *c__ = america_comparea__(n, &x[1], ind) - probdata_1.area[*ind - 1];
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALC. */
/*     ****************************************************************** */
} /* america_evalc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int america_evaljac__(n, x, ind, indjac, valjac, nnzjac, 
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
    static integer i__, i1, i2, ip;

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
/*     PARAMETERS */
/*     COMMON SCALARS */
/*     COMMON ARRAYS */
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
/*     LOCAL SCALARS */
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
    *nnzjac = probdata_1.nbord[*ind - 1] << 1;
    i__1 = *nnzjac;
    for (i__ = 1; i__ <= i__1; ++i__) {
	valjac[i__] = 0.;
    }
    i__1 = probdata_1.nbord[*ind - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ < probdata_1.nbord[*ind - 1]) {
	    ip = i__ + 1;
	} else {
	    ip = 1;
	}
	i1 = probdata_1.nplin[*ind + i__ * 17 - 18];
	i2 = probdata_1.nplin[*ind + ip * 17 - 18];
	indjac[(i__ << 1) - 1] = (i1 << 1) - 1;
	valjac[(i__ << 1) - 1] += x[i2 * 2];
	indjac[i__ * 2] = i1 << 1;
	valjac[i__ * 2] -= x[(i2 << 1) - 1];
	indjac[(ip << 1) - 1] = (i2 << 1) - 1;
	valjac[(ip << 1) - 1] -= x[i1 * 2];
	indjac[ip * 2] = i2 << 1;
	valjac[ip * 2] += x[(i1 << 1) - 1];
    }
    i__1 = *nnzjac;
    for (i__ = 1; i__ <= i__1; ++i__) {
	valjac[i__] *= .5;
    }
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALJAC. */
/*     ****************************************************************** */
} /* america_evaljac__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
doublereal america_comparea__(n, x, ind)
integer *n;
doublereal *x;
integer *ind;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, i1, i2, ip;

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
    i__1 = probdata_1.nbord[*ind - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ < probdata_1.nbord[*ind - 1]) {
	    ip = i__ + 1;
	} else {
	    ip = 1;
	}
	i1 = probdata_1.nplin[*ind + i__ * 17 - 18];
	i2 = probdata_1.nplin[*ind + ip * 17 - 18];
	ret_val = ret_val + x[i2 * 2] * x[(i1 << 1) - 1] - x[i1 * 2] * x[(i2 
		<< 1) - 1];
    }
    ret_val *= .5;
    return ret_val;
} /* america_comparea__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int america_evalhc__(n, x, ind, hclin, hccol, hcval, nnzhc, 
	flag__)
integer *n;
doublereal *x;
integer *ind, *hclin, *hccol;
doublereal *hcval;
integer *nnzhc, *flag__;
{
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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
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
} /* america_evalhc__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int america_evalhlp__(n, x, m, lambda, p, hp, goth, flag__)
integer *n;
doublereal *x;
integer *m;
doublereal *lambda, *p, *hp;
logical *goth;
integer *flag__;
{
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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
    /* Parameter adjustments */
    --hp;
    --p;
    --x;
    --lambda;

    /* Function Body */
    *flag__ = -1;
} /* america_evalhlp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int america_endp__(n, x, l, u, m, lambda, rho, equatn, 
	linear)
integer *n;
doublereal *x, *l, *u;
integer *m;
doublereal *lambda, *rho;
logical *equatn, *linear;
{
    extern /* Subroutine */ int america_drawsol__();

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
/*     SCALAR ARGUMENTS */
/*     ARRAY ARGUMENTS */
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
/*    america_drawsol__(n, &x[1], &c_false, "america-endsol.tex", (ftnlen)18);*/
/*     ****************************************************************** */
/*     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE ENDP. */
/*     ****************************************************************** */
} /* america_endp__ */

/*     ****************************************************************** */
/*     ****************************************************************** */
/* Subroutine */ int america_drawsol__(n, x, printnum, filename, filename_len)
integer *n;
doublereal *x;
logical *printnum;
char *filename;
ftnlen filename_len;
{
    /* Format strings */
    static char fmt_100[] = "(\002 \\put(\002,f10.3,\002,\002,f10.3,\002){\\\
circle*{\002,f9.3,\002}}\002)";
    static char fmt_110[] = "(\002 \\put(\002,f10.3,\002,\002,f10.3,\002)\
{\002,i3,\002}\002)";
    static char fmt_120[] = "(\002 \\put(\002,f10.3,\002,\002,f10.3,\002)\
{\002,\002$\\cdot$\002,\002}\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_wsle(), do_lio(), e_wsle(), s_wsfe(), do_fio(), 
	    e_wsfe();
    double sqrt();
    integer f_clos();

    /* Local variables */
    static doublereal diam, xpic, ypic, zpic;
    static integer i__, j;
    static doublereal ancho, p1, p2, q1, q2, v1, v2;
    static integer jp;
    static doublereal xx, yy, altura, amb;
    static integer ind1, ind2;

    /* Fortran I/O blocks */
    static cilist io___16 = { 0, 10, 0, 0, 0 };
    static cilist io___17 = { 0, 10, 0, 0, 0 };
    static cilist io___18 = { 0, 10, 0, 0, 0 };
    static cilist io___19 = { 0, 10, 0, 0, 0 };
    static cilist io___20 = { 0, 10, 0, 0, 0 };
    static cilist io___21 = { 0, 10, 0, 0, 0 };
    static cilist io___34 = { 0, 10, 0, fmt_100, 0 };
    static cilist io___35 = { 0, 10, 0, fmt_110, 0 };
    static cilist io___42 = { 0, 10, 0, fmt_120, 0 };
    static cilist io___43 = { 0, 10, 0, 0, 0 };
    static cilist io___44 = { 0, 10, 0, 0, 0 };
    static cilist io___45 = { 0, 10, 0, 0, 0 };


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
    ancho = 10.;
    altura = 19.;
    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = 18;
    o__1.ofnm = filename;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_wsle(&io___16);
    do_lio(&c__9, &c__1, "\\documentstyle[11pt]{article}", (ftnlen)29);
    e_wsle();
    s_wsle(&io___17);
    do_lio(&c__9, &c__1, "\\setlength{\\unitlength}{0.9cm}", (ftnlen)30);
    e_wsle();
    s_wsle(&io___18);
    do_lio(&c__9, &c__1, "\\pagestyle{empty}", (ftnlen)17);
    e_wsle();
    s_wsle(&io___19);
    do_lio(&c__9, &c__1, "\\begin{document}", (ftnlen)16);
    e_wsle();
    s_wsle(&io___20);
    do_lio(&c__9, &c__1, "\\begin{center}", (ftnlen)14);
    e_wsle();
    s_wsle(&io___21);
    do_lio(&c__9, &c__1, "\\begin{picture}(", (ftnlen)16);
    do_lio(&c__5, &c__1, (char *)&ancho, (ftnlen)sizeof(doublereal));
    do_lio(&c__9, &c__1, ",", (ftnlen)1);
    do_lio(&c__5, &c__1, (char *)&altura, (ftnlen)sizeof(doublereal));
    do_lio(&c__9, &c__1, ")", (ftnlen)1);
    e_wsle();
    i__1 = probdata_1.nstates;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = probdata_1.nbord[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    if (j < probdata_1.nbord[i__ - 1]) {
		jp = j + 1;
	    } else {
		jp = 1;
	    }
	    ind1 = probdata_1.nplin[i__ + j * 17 - 18];
	    ind2 = probdata_1.nplin[i__ + jp * 17 - 18];
	    if (! (ind1 == 2 && ind2 == 73) && ! (ind1 == 75 && ind2 == 74)) {
		p1 = x[(ind1 << 1) - 1] * .75;
		p2 = x[ind1 * 2] * .75;
		q1 = x[(ind2 << 1) - 1] * .75;
		q2 = x[ind2 * 2] * .75;
		if (*printnum) {
		    xpic = p1;
		    ypic = p2;
		    diam = .05;
		    s_wsfe(&io___34);
		    do_fio(&c__1, (char *)&xpic, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&ypic, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&diam, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		    s_wsfe(&io___35);
		    do_fio(&c__1, (char *)&xpic, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&ypic, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&ind1, (ftnlen)sizeof(integer));
		    e_wsfe();
		}
/* Computing 2nd power */
		d__1 = p1 - q1;
/* Computing 2nd power */
		d__2 = p2 - q2;
		zpic = sqrt(d__1 * d__1 + d__2 * d__2);
		if (zpic != 0.) {
		    v1 = (q1 - p1) / zpic;
		    v2 = (q2 - p2) / zpic;
		    amb = 0.;
L10:
		    if (amb <= zpic) {
			xx = p1 + amb * v1;
			yy = p2 + amb * v2;
			amb += .08;
			s_wsfe(&io___42);
			d__1 = xx - .05;
			do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				doublereal));
			d__2 = yy - .1;
			do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(
				doublereal));
			e_wsfe();
			goto L10;
		    }
		}
	    }
	}
    }
    s_wsle(&io___43);
    do_lio(&c__9, &c__1, "\\end{picture}", (ftnlen)13);
    e_wsle();
    s_wsle(&io___44);
    do_lio(&c__9, &c__1, "\\end{center}", (ftnlen)12);
    e_wsle();
    s_wsle(&io___45);
    do_lio(&c__9, &c__1, "\\end{document}", (ftnlen)14);
    e_wsle();
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*     NON-EXECUTABLE STATEMENTS */
} /* america_drawsol__ */

