C     =================================================================
C     File: america.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     May 20th, 2005.

C     Users are encouraged to download periodically updated versions of 
C     this code at the COLLECTION home page:
C
C     www.ime.usp.br/~egbirgin/collection/
C
C     and periodically updated versions of the TANGO Project solvers at
C     the TANGO home page:
C 
C     www.ime.usp.br/~egbirgin/tango/ 

C     =================================================================

C     America problem
C     ---------------

C     The objective is to draw a map of America in which the countries 
C     appear with areas that are proportional to their real values. The 
C     unknowns are $132$ points in $\R^2$, which define the boundaries
C     of $17$ countries. Each point is assigned to the boundary of one 
C     or more countries. The computed area of each country is calculated
C     as a function of its boundary points using Green's formula. The
C     constraints of the problem are:
C
C                   Computed area = True area 
C
C     for each country. As initial approximation we took the coordinates
C     of the~$132$ points in the New York Times map of America which, of
C     course, do not fit with true areas. The objective function is 
C     $\frac{1}{2}\sum_{j=1}^{132} \|P_j - Q_j\|_2^2$, where 
C     $P_1,\ldots,P_{132}$ are the unknowns and $Q_1,\ldots,Q_{132}$ are 
C     the initial approximation. 
C
C     Considered "countries":
C
C      1 Cuba
C      2 Canada
C      3 USA
C      4 Mexico
C      5 America Central
C      6 Colombia
C      7 Venezuela
C      8 Guianas
C      9 Brasil
C     10 Ecuador
C     11 Peru
C     12 Bolivia
C     13 Paraguay
C     14 Chile
C     15 Argentina
C     16 Alaska
C     17 Uruguay

C     ******************************************************************
C     ******************************************************************

      subroutine america_inip(n,x,l,u,m,lambda,rho,equatn,linear)

      implicit none

C     This subroutine must set some problem data. For achieving this 
C     objective YOU MUST MODIFY it according to your problem. See below 
C     where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     This subroutine has no input parameters.
C
C     On Return
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              initial point,
C
C     l        double precision l(n),
C              lower bounds on x,
C
C     u        double precision u(n),
C              upper bounds on x,
C
C     m        integer,
C              number of constraints (excluding the bounds),
C
C     lambda   double precision lambda(m),
C              initial estimation of the Lagrange multipliers,
C
C     rho      double precision rho(m),
C              initial penalty parameters.
C
C     equatn   logical equatn(m)
C              for each constraint j, set equatn(j) = .true. if it is an 
C              equality constraint of the form c_j(x) = 0, and set 
C              equatn(j) = .false. if it is an inequality constraint of 
C              the form c_j(x) <= 0,
C
C     linear   logical linear(m)
C              for each constraint j, set linear(j) = .true. if it is a 
C              linear constraint, and set linear(j) = .false. if it is a
C              nonlinear constraint.

C     PARAMETERS
      integer nnstates,nnpun
      parameter ( nnstates =  17 )
      parameter ( nnpun    = 132 )

C     COMMON SCALARS
      integer nstates,npun

C     COMMON ARRAYS
      integer nbord(nnstates),nplin(nnstates,nnpun)
      double precision area(nnstates),put(nnpun,2)

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(*),linear(*)
      double precision l(*),lambda(*),rho(*),u(*),x(*)

C     LOCAL SCALARS
      integer i
      double precision totarea, falsearea

C     COMMON BLOCKS
      common /probdata/ put,area,nbord,nplin,nstates,npun

C     EXTERNAL FUNCTIONS
      double precision america_comparea

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR PROBLEM
C     DATA:
C     ******************************************************************

C     Set problem data

C     Number of states and number of points

      nstates =  17
      npun    = 132

C     Number of boundary points for each state and indices of the 
C     boundary points

      nbord(1) = 4

      nplin( 1, 1) =  49
      nplin( 1, 2) =  50
      nplin( 1, 3) =  51
      nplin( 1, 4) =  52


      nbord(2) = 21 

      nplin( 2, 1) =  83 
      nplin( 2, 2) =  35 
      nplin( 2, 3) =  64 
      nplin( 2, 4) =  63 
      nplin( 2, 5) =  84 
      nplin( 2, 6) =  36 
      nplin( 2, 7) =  99 
      nplin( 2, 8) = 100 
      nplin( 2, 9) =  37 
      nplin( 2,10) =  38 
      nplin( 2,11) =  39 
      nplin( 2,12) = 101 
      nplin( 2,13) = 102 
      nplin( 2,14) =  40 
      nplin( 2,15) =  41 
      nplin( 2,16) =  42 
      nplin( 2,17) = 127 
      nplin( 2,18) =  81 
      nplin( 2,19) =  82 
      nplin( 2,20) =  43 
      nplin( 2,21) =  46 

      nbord(3) = 22 

      nplin( 3, 1) =  33 
      nplin( 3, 2) = 108 
      nplin( 3, 3) = 109 
      nplin( 3, 4) = 110 
      nplin( 3, 5) = 111 
      nplin( 3, 6) = 112 
      nplin( 3, 7) =  34 
      nplin( 3, 8) =  94 
      nplin( 3, 9) =  62 
      nplin( 3,10) =  98 
      nplin( 3,11) =  61 
      nplin( 3,12) =  60 
      nplin( 3,13) =  96 
      nplin( 3,14) =  97 
      nplin( 3,15) =  36 
      nplin( 3,16) =  84 
      nplin( 3,17) =  63 
      nplin( 3,18) =  64 
      nplin( 3,19) =  35 
      nplin( 3,20) = 120 
      nplin( 3,21) =  56 
      nplin( 3,22) = 121 

      nbord(4) = 26 

      nplin( 4, 1) =  86 
      nplin( 4, 2) = 113 
      nplin( 4, 3) = 114 
      nplin( 4, 4) =  31 
      nplin( 4, 5) =  32 
      nplin( 4, 6) =  87 
      nplin( 4, 7) = 115 
      nplin( 4, 8) =  88 
      nplin( 4, 9) = 103 
      nplin( 4,10) =  89 
      nplin( 4,11) =  85 
      nplin( 4,12) =  34 
      nplin( 4,13) = 112 
      nplin( 4,14) = 111 
      nplin( 4,15) = 110 
      nplin( 4,16) = 109 
      nplin( 4,17) = 108 
      nplin( 4,18) =  33 
      nplin( 4,19) = 118 
      nplin( 4,20) = 104 
      nplin( 4,21) = 116 
      nplin( 4,22) = 119 
      nplin( 4,23) = 105 
      nplin( 4,24) = 106 
      nplin( 4,25) = 117  
      nplin( 4,26) = 107 

      nbord(5) = 8 

      nplin( 5, 1) =  30 
      nplin( 5, 2) =  29 
      nplin( 5, 3) =  91 
      nplin( 5, 4) =  90 
      nplin( 5, 5) =  32 
      nplin( 5, 6) =  31 
      nplin( 5, 7) =  92 
      nplin( 5, 8) =  93  

      nbord(6) = 9 

      nplin( 6, 1) =  28 
      nplin( 6, 2) =  25 
      nplin( 6, 3) =  23 
      nplin( 6, 4) =  22 
      nplin( 6, 5) =  55 
      nplin( 6, 6) =  95 
      nplin( 6, 7) =  24 
      nplin( 6, 8) =  29 
      nplin( 6, 9) =  30 

      nbord(7) = 9 

      nplin( 7, 1) =  22 
      nplin( 7, 2) =  21 
      nplin( 7, 3) =  20 
      nplin( 7, 4) =  19 
      nplin( 7, 5) = 122 
      nplin( 7, 6) = 123 
      nplin( 7, 7) =  24 
      nplin( 7, 8) =  95 
      nplin( 7, 9) =  55 

      nbord(8) = 5 

      nplin( 8, 1) =  20 
      nplin( 8, 2) = 126 
      nplin( 8, 3) =  18 
      nplin( 8, 4) =  17 
      nplin( 8, 5) =  19 

      nbord(9) = 22 

      nplin( 9, 1) =  53 
      nplin( 9, 2) =  16 
      nplin( 9, 3) =  15 
      nplin( 9, 4) =  10 
      nplin( 9, 5) =   8 
      nplin( 9, 6) =   7 
      nplin( 9, 7) =   6 
      nplin( 9, 8) =  65 
      nplin( 9, 9) =   5 
      nplin( 9,10) =  68 
      nplin( 9,11) =  76 
      nplin( 9,12) =  48 
      nplin( 9,13) = 129 
      nplin( 9,14) =  47 
      nplin( 9,15) = 125 
      nplin( 9,16) =  17 
      nplin( 9,17) =  18 
      nplin( 9,18) = 126 
      nplin( 9,19) =  20 
      nplin( 9,20) =  21 
      nplin( 9,21) =  22 
      nplin( 9,22) =  23 

      nbord(10) = 4 

      nplin(10, 1) =  27 
      nplin(10, 2) =  26 
      nplin(10, 3) =  25 
      nplin(10, 4) =  28 

      nbord(11) = 10 

      nplin(11, 1) = 128 
      nplin(11, 2) =  59 
      nplin(11, 3) =  12  
      nplin(11, 4) =  13 
      nplin(11, 5) =  16 
      nplin(11, 6) =  53 
      nplin(11, 7) =  23 
      nplin(11, 8) =  25 
      nplin(11, 9) =  26 
      nplin(11,10) =  27 

      nbord(12) = 7 

      nplin(12, 1) =  14 
      nplin(12, 2) =   9 
      nplin(12, 3) =  10 
      nplin(12, 4) =  15 
      nplin(12, 5) =  16 
      nplin(12, 6) =  13 
      nplin(12, 7) =  12 

      nbord(13) = 6 

      nplin(13, 1) =  54 
      nplin(13, 2) =   7 
      nplin(13, 3) =   8 
      nplin(13, 4) =  10 
      nplin(13, 5) =   9 
      nplin(13, 6) =  70 

      nbord(14) = 16 

      nplin(14, 1) =   1 
      nplin(14, 2) =  11  
      nplin(14, 3) = 132 
      nplin(14, 4) = 131 
      nplin(14, 5) =  14 
      nplin(14, 6) =  12 
      nplin(14, 7) =  59 
      nplin(14, 8) = 130 
      nplin(14, 9) =  57 
      nplin(14,10) =  58 
      nplin(14,11) =   2 
      nplin(14,12) =  73 
      nplin(14,13) =  77 
      nplin(14,14) =  71 
      nplin(14,15) =  75 
      nplin(14,16) =  74 

      nbord(15) = 20 

      nplin(15, 1) =   1 
      nplin(15, 2) =   2 
      nplin(15, 3) =  73 
      nplin(15, 4) =  71 
      nplin(15, 5) =  72 
      nplin(15, 6) =  75 
      nplin(15, 7) =  74 
      nplin(15, 8) =  69  
      nplin(15, 9) =  67 
      nplin(15,10) =  66 
      nplin(15,11) =   4 
      nplin(15,12) =   6 
      nplin(15,13) =   7 
      nplin(15,14) =  54 
      nplin(15,15) =  70 
      nplin(15,16) =   9 
      nplin(15,17) =  14 
      nplin(15,18) = 131 
      nplin(15,19) = 132 
      nplin(15,20) =  11 

      nbord(16) = 8 

      nplin(16, 1) =  43 
      nplin(16, 2) =  80 
      nplin(16, 3) =  44 
      nplin(16, 4) =  78 
      nplin(16, 5) = 124 
      nplin(16, 6) =  45 
      nplin(16, 7) =  79 
      nplin(16, 8) =  46 

      nbord(17) = 5 

      nplin(17, 1) =   4 
      nplin(17, 2) =   3 
      nplin(17, 3) =   5 
      nplin(17, 4) =  65 
      nplin(17, 5) =   6 

C     Real area of each state

      area( 1) =  115.0d0
      area( 2) = 9922.0d0
      area( 3) = 7885.0d0
      area( 4) = 1973.0d0
      area( 5) =  543.0d0
      area( 6) = 1139.0d0
      area( 7) =  912.0d0
      area( 8) =  470.0d0
      area( 9) = 8512.0d0
      area(10) =  461.0d0
      area(11) = 1285.0d0
      area(12) = 1099.0d0
      area(13) =  407.0d0
      area(14) =  752.0d0
      area(15) = 2767.0d0
      area(16) = 1478.0d0
      area(17) =  187.0d0

C     Putative coordinates of each boundary point (initial approximation)

      put(  1,1) = 11.70d0
      put(  1,2) =  5.70d0
      put(  2,1) = 12.10d0
      put(  2,2) =  5.60d0
      put(  3,1) = 13.50d0
      put(  3,2) =  8.30d0
      put(  4,1) = 13.10d0
      put(  4,2) =  8.40d0
      put(  5,1) = 13.70d0
      put(  5,2) =  8.60d0
      put(  6,1) = 13.20d0
      put(  6,2) =  9.00d0
      put(  7,1) = 13.60d0
      put(  7,2) =  9.40d0
      put(  8,1) = 13.40d0
      put(  8,2) = 10.00d0
      put(  9,1) = 12.50d0
      put(  9,2) = 10.00d0
      put( 10,1) = 13.10d0
      put( 10,2) = 10.40d0
      put( 11,1) = 11.70d0
      put( 11,2) =  6.40d0
      put( 12,1) = 11.60d0
      put( 12,2) = 10.60d0
      put( 13,1) = 11.60d0
      put( 13,2) = 11.20d0
      put( 14,1) = 12.00d0
      put( 14,2) =  9.90d0
      put( 15,1) = 12.70d0
      put( 15,2) = 11.20d0
      put( 16,1) = 11.80d0
      put( 16,2) = 11.60d0
      put( 17,1) = 13.80d0
      put( 17,2) = 13.50d0
      put( 18,1) = 13.50d0
      put( 18,2) = 13.20d0
      put( 19,1) = 12.90d0
      put( 19,2) = 14.10d0
      put( 20,1) = 12.80d0
      put( 20,2) = 13.70d0
      put( 21,1) = 12.40d0
      put( 21,2) = 13.30d0
      put( 22,1) = 11.90d0
      put( 22,2) = 13.30d0
      put( 23,1) = 11.70d0
      put( 23,2) = 12.70d0
      put( 24,1) = 11.60d0
      put( 24,2) = 14.60d0
      put( 25,1) = 11.10d0
      put( 25,2) = 12.80d0
      put( 26,1) = 10.80d0
      put( 26,2) = 12.50d0
      put( 27,1) = 10.50d0
      put( 27,2) = 12.50d0
      put( 28,1) = 10.70d0
      put( 28,2) = 13.20d0
      put( 29,1) = 10.90d0
      put( 29,2) = 14.20d0
      put( 30,1) = 10.80d0
      put( 30,2) = 14.10d0
      put( 31,1) =  9.20d0
      put( 31,2) = 15.00d0
      put( 32,1) =  9.40d0
      put( 32,2) = 15.30d0
      put( 33,1) =  6.60d0
      put( 33,2) = 17.60d0
      put( 34,1) =  8.40d0
      put( 34,2) = 16.60d0
      put( 35,1) =  6.40d0
      put( 35,2) = 20.00d0
      put( 36,1) = 12.50d0
      put( 36,2) = 19.50d0
      put( 37,1) = 13.80d0
      put( 37,2) = 20.60d0
      put( 38,1) = 11.80d0
      put( 38,2) = 22.40d0
      put( 39,1) = 11.20d0
      put( 39,2) = 20.50d0
      put( 40,1) = 10.00d0
      put( 40,2) = 22.20d0
      put( 41,1) = 11.50d0
      put( 41,2) = 23.00d0
      put( 42,1) = 11.00d0
      put( 42,2) = 24.60d0
      put( 43,1) =  5.80d0
      put( 43,2) = 23.80d0
      put( 44,1) =  3.50d0
      put( 44,2) = 23.80d0
      put( 45,1) =  3.00d0
      put( 45,2) = 21.50d0
      put( 46,1) =  5.50d0
      put( 46,2) = 22.00d0
      put( 47,1) = 15.70d0
      put( 47,2) = 12.20d0
      put( 48,1) = 15.30d0
      put( 48,2) = 10.80d0
      put( 49,1) = 10.40d0
      put( 49,2) = 16.20d0
      put( 50,1) = 11.00d0
      put( 50,2) = 15.80d0
      put( 51,1) = 11.30d0
      put( 51,2) = 16.00d0
      put( 52,1) = 10.40d0
      put( 52,2) = 16.30d0
      put( 53,1) = 11.60d0
      put( 53,2) = 12.00d0
      put( 54,1) = 13.00d0
      put( 54,2) =  9.30d0
      put( 55,1) = 12.00d0
      put( 55,2) = 13.80d0
      put( 56,1) =  6.00d0
      put( 56,2) = 18.70d0
      put( 57,1) = 11.50d0
      put( 57,2) =  5.70d0
      put( 58,1) = 11.80d0
      put( 58,2) =  5.50d0
      put( 59,1) = 11.40d0
      put( 59,2) = 10.60d0
      put( 60,1) = 10.70d0
      put( 60,2) = 17.40d0
      put( 61,1) = 10.60d0
      put( 61,2) = 16.60d0
      put( 62,1) = 10.40d0
      put( 62,2) = 17.20d0
      put( 63,1) = 11.00d0
      put( 63,2) = 19.10d0
      put( 64,1) = 10.00d0
      put( 64,2) = 20.00d0
      put( 65,1) = 13.50d0
      put( 65,2) =  8.80d0
      put( 66,1) = 13.20d0
      put( 66,2) =  7.70d0
      put( 67,1) = 12.60d0
      put( 67,2) =  7.70d0
      put( 68,1) = 14.20d0
      put( 68,2) =  9.70d0
      put( 69,1) = 12.30d0
      put( 69,2) =  6.70d0
      put( 70,1) = 13.10d0
      put( 70,2) =  9.60d0
      put( 71,1) = 12.20d0
      put( 71,2) =  5.20d0
      put( 72,1) = 12.50d0
      put( 72,2) =  5.30d0
      put( 73,1) = 12.20d0
      put( 73,2) =  5.40d0
      put( 74,1) = 12.10d0
      put( 74,2) =  5.60d0
      put( 75,1) = 12.20d0
      put( 75,2) =  5.40d0
      put( 76,1) = 15.00d0
      put( 76,2) = 10.00d0
      put( 77,1) = 11.90d0
      put( 77,2) =  5.30d0
      put( 78,1) =  2.70d0
      put( 78,2) = 23.00d0
      put( 79,1) =  4.00d0
      put( 79,2) = 22.00d0
      put( 80,1) =  4.50d0
      put( 80,2) = 24.00d0
      put( 81,1) =  9.00d0
      put( 81,2) = 23.50d0
      put( 82,1) =  7.30d0
      put( 82,2) = 23.70d0
      put( 83,1) =  6.00d0
      put( 83,2) = 21.00d0
      put( 84,1) = 12.00d0
      put( 84,2) = 19.50d0
      put( 85,1) =  8.50d0
      put( 85,2) = 16.00d0
      put( 86,1) =  7.70d0
      put( 86,2) = 15.80d0
      put( 87,1) =  9.30d0
      put( 87,2) = 15.50d0
      put( 88,1) =  9.90d0
      put( 88,2) = 16.00d0
      put( 89,1) =  9.00d0
      put( 89,2) = 15.60d0
      put( 90,1) = 10.20d0
      put( 90,2) = 15.30d0
      put( 91,1) = 10.20d0
      put( 91,2) = 14.50d0
      put( 92,1) =  9.70d0
      put( 92,2) = 14.70d0
      put( 93,1) = 10.00d0
      put( 93,2) = 14.40d0
      put( 94,1) =  9.30d0
      put( 94,2) = 17.20d0
      put( 95,1) = 11.50d0
      put( 95,2) = 14.00d0
      put( 96,1) = 11.30d0
      put( 96,2) = 17.80d0
      put( 97,1) = 11.40d0
      put( 97,2) = 18.40d0
      put( 98,1) = 10.50d0
      put( 98,2) = 16.60d0
      put( 99,1) = 13.10d0
      put( 99,2) = 19.50d0
      put(100,1) = 12.40d0
      put(100,2) = 19.90d0
      put(101,1) = 11.20d0
      put(101,2) = 20.20d0
      put(102,1) = 10.00d0
      put(102,2) = 21.50d0
      put(103,1) =  9.50d0
      put(103,2) = 16.00d0
      put(104,1) =  7.10d0
      put(104,2) = 16.30d0
      put(105,1) =  6.80d0
      put(105,2) = 17.50d0
      put(106,1) =  6.90d0
      put(106,2) = 17.50d0
      put(107,1) =  7.80d0
      put(107,2) = 15.90d0
      put(108,1) =  6.80d0
      put(108,2) = 17.60d0
      put(109,1) =  7.40d0
      put(109,2) = 17.40d0
      put(110,1) =  7.70d0
      put(110,2) = 17.50d0
      put(111,1) =  8.10d0
      put(111,2) = 17.10d0
      put(112,1) =  8.30d0
      put(112,2) = 17.30d0
      put(113,1) =  8.50d0
      put(113,2) = 15.30d0
      put(114,1) =  8.90d0
      put(114,2) = 15.30d0
      put(115,1) =  9.70d0
      put(115,2) = 15.60d0
      put(116,1) =  7.20d0
      put(116,2) = 16.30d0
      put(117,1) =  7.40d0
      put(117,2) = 16.60d0
      put(118,1) =  6.80d0
      put(118,2) = 16.80d0
      put(119,1) =  7.00d0
      put(119,2) = 16.80d0
      put(120,1) =  6.10d0
      put(120,2) = 19.50d0
      put(121,1) =  6.20d0
      put(121,2) = 18.20d0
      put(122,1) = 12.50d0
      put(122,2) = 14.50d0
      put(123,1) = 12.10d0
      put(123,2) = 14.50d0
      put(124,1) =  2.75d0
      put(124,2) = 22.15d0
      put(125,1) = 14.60d0
      put(125,2) = 12.80d0
      put(126,1) = 13.00d0
      put(126,2) = 13.20d0
      put(127,1) = 10.20d0
      put(127,2) = 23.30d0
      put(128,1) = 10.70d0
      put(128,2) = 11.25d0
      put(129,1) = 15.30d0
      put(129,2) = 11.30d0
      put(130,1) = 11.40d0
      put(130,2) =  7.50d0
      put(131,1) = 11.65d0
      put(131,2) =  9.00d0
      put(132,1) = 11.60d0
      put(132,2) =  7.25d0

C     Number of variables

      n = 2 * npun

C     Initial point
      do i = 1,npun
          x(2*i-1) = put(i,1)
          x(2*i  ) = put(i,2)
      end do
 
      call america_drawsol(n,x,.false.,'america-inisol.tex')

C     Scale areas

      totarea = 0.0d0
      do i = 1,nstates
          totarea = totarea + area(i)
      end do

      falsearea = 0.0d0
      do i = 1,nstates
         falsearea = falsearea + america_comparea(n,x,i)
      end do

      do i = 1,nstates
          area(i) = area(i) * falsearea / totarea
      end do

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 1.0d+03
          u(i) =   1.0d+03
      end do

C     Number of constraints (equalities plus inequalities)

      m = nstates

C     Lagrange multipliers approximation. Most users prefer to use the 
C     null initial Lagrange multipliers estimates. However, if the 
C     problem that you are solving is "slightly different" from a 
C     previously solved problem of which you know the correct Lagrange 
C     multipliers, we encourage you to set these multipliers as initial 
C     estimates. Of course, in this case you are also encouraged to use 
C     the solution of the previous problem as initial estimate of the 
C     solution. Similarly, most users prefer to use rho = 10 as initial 
C     penalty parameters. But in the case mentioned above (good 
C     estimates of solution and Lagrange multipliers) larger values of 
C     the penalty parameters (say, rho = 1000) may be more useful. More 
C     warm-start procedures are being elaborated.

      do i = 1,m
          lambda(i) =  0.0d0
      end do

C     Initial penalty parameters

      do i = 1,m
          rho(i) = 10.0d0
      end do

C     For each constraint i, set equatn(i) = .true. if it is an equality
C     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if 
C     it is an inequality constraint of the form c_i(x) <= 0.

      do i = 1,m
          equatn(i) = .true.
      end do

C     For each constraint i, set linear(i) = .true. if it is a linear
C     constraint, otherwise set linear(i) = .false.

      do i = 1,m
          linear(i) = .false.
      end do

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE INIP.
C     ******************************************************************

      end

C     ******************************************************************
C     ******************************************************************

      subroutine america_evalf(n,x,f,flag)

      implicit none

C     This subroutine must compute the objective function. For achieving 
C     this objective YOU MUST MODIFY it according to your problem. See 
C     below where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     On Return
C
C     f        double precision,
C              objective function value at x,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the objective 
C              function. (For example, trying to compute the square root 
C              of a negative number, dividing by zero or a very small 
C              number, etc.) If everything was o.k. you must set it 
C              equal to zero.

C     PARAMETERS
      integer nnstates,nnpun
      parameter ( nnstates =  17 )
      parameter ( nnpun    = 132 )

C     COMMON SCALARS
      integer nstates,npun

C     COMMON ARRAYS
      integer nbord(nnstates),nplin(nnstates,nnpun)
      double precision area(nnstates),put(nnpun,2)

C     SCALAR ARGUMENTS
      integer flag,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /probdata/ put,area,nbord,nplin,nstates,npun

C     Objective function

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR OBJECTIVE
C     FUNCTION:
C     ******************************************************************

      flag = 0

      f = 0.0d0
      do i = 1,npun
          f = f + ( x(2*i-1) - put(i,1) ) ** 2 + 
     +            ( x(2*i  ) - put(i,2) ) ** 2
      end do

      f = 0.5d0 * f

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALF.
C     ******************************************************************

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine america_evalg(n,x,g,flag)

      implicit none

C     This subroutine must compute the gradient vector of the objective 
C     function. For achieving these objective YOU MUST MODIFY it in the 
C     way specified below. However, if you decide to use numerical 
C     derivatives (we dont encourage this option at all!) you dont need
C     to modify evalg.
C
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     On Return
C
C     g        double precision g(n),
C              gradient vector of the objective function evaluated at x,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of any component 
C              of the gradient vector. (For example, trying to compute 
C              the square root of a negative number, dividing by zero or 
C              a very small number, etc.) If everything was o.k. you 
C              must set it equal to zero.

C     PARAMETERS
      integer nnstates,nnpun
      parameter ( nnstates =  17 )
      parameter ( nnpun    = 132 )

C     COMMON SCALARS
      integer nstates,npun

C     COMMON ARRAYS
      integer nbord(nnstates),nplin(nnstates,nnpun)
      double precision area(nnstates),put(nnpun,2)

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /probdata/ put,area,nbord,nplin,nstates,npun

C     Gradient vector

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENT
C     VECTOR OF YOUR OBJECTIVE FUNCTION: 
C     ******************************************************************

      flag = 0

      do i = 1,npun
          g(2*i-1) = x(2*i-1) - put(i,1)
          g(2*i)   = x(2*i)   - put(i,2)
      end do

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALG. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine america_evalh(n,x,hlin,hcol,hval,nnzh,flag)

      implicit none

C     This subroutine might compute the Hessian matrix of the objective 
C     function. For achieving this objective YOU MAY MODIFY it according 
C     to your problem. To modify this subroutine IS NOT MANDATORY. See 
C     below where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     On Return
C
C     nnzh     integer,
C              number of perhaps-non-null elements of the computed 
C              Hessian,
C
C     hlin     integer hlin(nnzh),
C              see below,
C
C     hcol     integer hcol(nnzh),
C              see below,
C
C     hval     double precision hval(nnzh),
C              the non-null value of the (hlin(k),hcol(k)) position 
C              of the Hessian matrix of the objective function must 
C              be saved at hval(k). Just the lower triangular part of
C              Hessian matrix must be computed,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the Hessian
C              matrix of the objective funtion. (For example, trying 
C              to compute the square root of a negative number, 
C              dividing by zero or a very small number, etc.) If 
C              everything was o.k. you must set it equal to zero.

C     SCALAR ARGUMENTS
      integer flag,n,nnzh

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     ******************************************************************
C     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE 
C     HESSIAN MATRIX OF YOUR OBJECTIVE FUNCTION: 
C     ******************************************************************

      flag = - 1

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALH. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine america_evalc(n,x,ind,c,flag)

      implicit none

C     This subroutine must compute the ind-th constraint of your problem. 
C     For achieving this objective YOU MUST MOFIFY it according to your 
C     problem. See below the places where your modifications must be 
C     inserted.
C
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     ind      integer,
C              index of the constraint to be computed,
C
C     On Return
C
C     c        double precision,
C              ind-th constraint evaluated at x,
C
C     flag     integer
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the 
C              constraint. (For example, trying to compute the square 
C              root of a negative number, dividing by zero or a very 
C              small number, etc.) If everything was o.k. you must set 
C              it equal to zero.
 
C     PARAMETERS
      integer nnstates,nnpun
      parameter ( nnstates =  17 )
      parameter ( nnpun    = 132 )

C     COMMON SCALARS
      integer nstates,npun

C     COMMON ARRAYS
      integer nbord(nnstates),nplin(nnstates,nnpun)
      double precision area(nnstates),put(nnpun,2)

C     SCALAR ARGUMENTS
      integer ind,flag,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

C     COMMON BLOCKS
      common /probdata/ put,area,nbord,nplin,nstates,npun

C     EXTERNAL FUNCTIONS
      double precision america_comparea

C     Constraints

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR
C     CONSTRAINTS: 
C     ******************************************************************

      flag = 0

      c = america_comparea(n,x,ind) - area(ind)

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALC. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine america_evaljac(n,x,ind,indjac,valjac,nnzjac,flag)

      implicit none

C     This subroutine must compute the gradient of the constraint i. For 
C     achieving these objective YOU MUST MODIFY it in the way specified 
C     below.
C
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     ind      integer,
C              index of the constraint whose gradient will be computed,
C
C     On Return
C
C     nnzjac   integer,
C              number of perhaps-non-null elements of the computed 
C              gradient,
C
C     indjac   integer indjac(nnzjac),
C              see below,
C
C     valjac   double precision valjac(nnzjac),
C              the non-null value of the partial derivative of the ind-th
C              constraint with respect to the indjac(k)-th variable must 
C              be saved at valjac(k).
C
C     flag     integer
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the 
C              constraint. (For example, trying to compute the square 
C              root of a negative number, dividing by zero or a very 
C              small number, etc.) If everything was o.k. you must set 
C              it equal to zero.

C     PARAMETERS
      integer nnstates,nnpun
      parameter ( nnstates =  17 )
      parameter ( nnpun    = 132 )

C     COMMON SCALARS
      integer nstates,npun

C     COMMON ARRAYS
      integer nbord(nnstates),nplin(nnstates,nnpun)
      double precision area(nnstates),put(nnpun,2)

C     SCALAR ARGUMENTS
      integer flag,ind,n,nnzjac

C     ARRAY ARGUMENTS
      integer indjac(n)
      double precision x(n),valjac(n)

C     LOCAL SCALARS
      integer i,i1,i2,ip

C     COMMON BLOCKS
      common /probdata/ put,area,nbord,nplin,nstates,npun

C     Sparse gradient vector of the ind-th constraint

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENTS  
C     OF YOUR CONSTRAINTS: 
C     ******************************************************************

      flag = 0

      nnzjac = 2 * nbord(ind)

      do i = 1,nnzjac
          valjac(i) = 0.0d0
      end do

      do i = 1,nbord(ind)

          if ( i .lt. nbord(ind) ) then
              ip = i + 1
          else
              ip = 1
          end if

          i1 = nplin(ind,i)
          i2 = nplin(ind,ip)

          indjac(2*i-1)  = 2 * i1 - 1
          valjac(2*i-1)  = valjac(2*i-1)  + x(2*i2)

          indjac(2*i)    = 2 * i1
          valjac(2*i)    = valjac(2*i)    - x(2*i2-1)

          indjac(2*ip-1) = 2 * i2 - 1
          valjac(2*ip-1) = valjac(2*ip-1) - x(2*i1)

          indjac(2*ip)   = 2 * i2
          valjac(2*ip)   = valjac(2*ip)   + x(2*i1-1)

      end do

      do i = 1,nnzjac
          valjac(i) = 0.5d0 * valjac(i)
      end do      

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALJAC. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************
 
      double precision function america_comparea(n,x,ind)

      implicit none

C     PARAMETERS
      integer nnstates,nnpun
      parameter ( nnstates =  17 )
      parameter ( nnpun    = 132 )

C     COMMON SCALARS
      integer nstates,npun

C     COMMON ARRAYS
      integer nbord(nnstates),nplin(nnstates,nnpun)
      double precision area(nnstates),put(nnpun,2)

C     SCALAR ARGUMENTS
      integer ind,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer i,i1,i2,ip

C     COMMON BLOCKS
      common /probdata/ put,area,nbord,nplin,nstates,npun

      america_comparea = 0.0d0

      do i = 1,nbord(ind)

          if ( i .lt. nbord(ind) ) then
              ip = i + 1
          else
              ip = 1
          end if

          i1 = nplin(ind,i)
          i2 = nplin(ind,ip)

          america_comparea = america_comparea + 
     +               x(2*i2) * x(2*i1-1) - x(2*i1) * x(2*i2-1)

      end do

      america_comparea = 0.5d0 * america_comparea

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine america_evalhc(n,x,ind,hclin,hccol,hcval,nnzhc,flag)

      implicit none

C     This subroutine might compute the Hessian matrix of the ind-th
C     constraint. For achieving this objective YOU MAY MODIFY it 
C     according to your problem. To modify this subroutine IS NOT 
C     MANDATORY. See below where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     ind      integer,
C              index of the constraint whose Hessian will be computed,
C
C     On Return
C
C     nnzhc    integer,
C              number of perhaps-non-null elements of the computed 
C              Hessian,
C
C     hclin    integer hclin(nnzhc),
C              see below,
C
C     hccol    integer hccol(nnzhc),
C              see below,
C
C     hcval    double precision hcval(nnzhc),
C              the non-null value of the (hclin(k),hccol(k)) position 
C              of the Hessian matrix of the ind-th constraint must 
C              be saved at hcval(k). Just the lower triangular part of
C              Hessian matrix must be computed,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the Hessian
C              matrix of the ind-th constraint. (For example, trying 
C              to compute the square root of a negative number, 
C              dividing by zero or a very small number, etc.) If 
C              everything was o.k. you must set it equal to zero.

C     SCALAR ARGUMENTS
      integer flag,ind,n,nnzhc

C     ARRAY ARGUMENTS
      integer hccol(*),hclin(*)
      double precision hcval(*),x(n)

C     ******************************************************************
C     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE 
C     HESSIANS OF YOUR CONSTRAINTS: 
C     ******************************************************************

      flag = - 1

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALHC. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************

      subroutine america_evalhlp(n,x,m,lambda,p,hp,goth,flag)

      implicit none

C     This subroutine might compute the product of the Hessian of the
C     Lagrangian times vector p (just the Hessian of the objective 
C     function in the unconstrained or bound-constrained case). 
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     m        integer,
C              number of constraints,
C
C     lambda   double precision lambda(m),
C              vector of Lagrange multipliers,
C
C     p        double precision p(n),
C              vector of the matrix-vector product,
C
C     goth     logical,
C              can be used to indicate if the Hessian matrices were
C              computed at the current point. It is set to .false.
C              by the optimization method every time the current
C              point is modified. Sugestion: if its value is .false. 
C              then compute the Hessians, save them in a common 
C              structure and set goth to .true.. Otherwise, just use 
C              the Hessians saved in the common block structure,
C
C     On Return
C
C     hp       double precision hp(n),
C              Hessian-vector product,
C
C     goth     logical,
C              see above,
C              
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the 
C              Hessian-vector product. (For example, trying to compute 
C              the square root of a negative number, dividing by zero 
C              or a very small number, etc.) If everything was o.k. you 
C              must set it equal to zero.

C     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine america_endp(n,x,l,u,m,lambda,rho,equatn,linear)

      implicit none

C     This subroutine can be used to do some extra job after the solver
C     has found the solution,like some extra statistics, or to save the
C     solution in some special format or to draw some graphical
C     representation of the solution. If the information given by the
C     solver is enough for you then leave the body of this subroutine
C     empty.
C     
C     Parameters of the subroutine:
C
C     The paraemters of this subroutine are the same parameters of
C     subroutine inip. But in this subroutine there are not output
C     parameter. All the parameters are input parameters.

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),rho(m),u(n),x(n)

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO COMPLEMENT THE
C     INFORMATION RELATED TO THE SOLUTION GIVEN BY THE SOLVER
C     ******************************************************************

      call america_drawsol(n,x,.false.,'america-endsol.tex')

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE ENDP.
C     ******************************************************************

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine america_drawsol(n,x,printnum,filename)

      implicit none

C     PARAMETERS
      integer nnstates,nnpun
      parameter ( nnstates =  17 )
      parameter ( nnpun    = 132 )

C     COMMON SCALARS
      integer nstates,npun

C     COMMON ARRAYS
      integer nbord(nnstates),nplin(nnstates,nnpun)
      double precision area(nnstates),put(nnpun,2)

C     SCALAR ARGUMENTS
      character * 18 filename
      logical printnum
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer i,ind1,ind2,j,jp
      double precision amb,ancho,altura,diam,p1,p2,q1,q2,v1,v2,xpic,xx,
     +        ypic,yy,zpic

C     COMMON BLOCKS
      common /probdata/ put,area,nbord,nplin,nstates,npun

      ancho  = 10.0d0
      altura = 19.0d0

      open(unit=10,file=filename)

      write(10,*) '\\documentstyle[11pt]{article}'
      write(10,*) '\\setlength{\\unitlength}{0.9cm}'
      write(10,*) '\\pagestyle{empty}'
      write(10,*) '\\begin{document}'
      write(10,*) '\\begin{center}'
      write(10,*) '\\begin{picture}(',ancho,',', altura,')'

      do i = 1,nstates

          do j = 1,nbord(i)

              if ( j .lt. nbord(i) ) then
                  jp = j + 1
              else
                  jp = 1
              end if

              ind1 = nplin(i,j)
              ind2 = nplin(i,jp)

              if ( .not. ( ind1 .eq.  2 .and. ind2 .eq. 73 ) .and.
     +             .not. ( ind1 .eq. 75 .and. ind2 .eq. 74 ) ) then

                  p1 = x(2*ind1-1) * 0.75d0
                  p2 = x(2*ind1  ) * 0.75d0
                  q1 = x(2*ind2-1) * 0.75d0
                  q2 = x(2*ind2  ) * 0.75d0

                  if ( printnum ) then
                      xpic = p1
                      ypic = p2
                      diam = 0.05d0
                      write(10,100) xpic,ypic,diam
                      write(10,110) xpic,ypic,ind1
                  end if

                  zpic = sqrt ( (p1 - q1) ** 2 + (p2 - q2) ** 2 )

                  if ( zpic .ne. 0.0d0 ) then

                      v1 = ( q1 - p1 ) / zpic
                      v2 = ( q2 - p2 ) / zpic

                      amb = 0.0d0

 10                   continue

                      if ( amb .le. zpic ) then

                          xx = p1 + amb * v1
                          yy = p2 + amb * v2

                          amb = amb + 0.08d0

                          write(10,120) xx - 0.05d0,yy - 0.1d0
	
                          go to 10

                      end if

                  end if

              end if

          end do

      end do

      write(10,*) '\\end{picture}'
      write(10,*) '\\end{center}'
      write(10,*) '\\end{document}'

      close(10)

C     NON-EXECUTABLE STATEMENTS

 100  format(' \\put(',f10.3,',',f10.3,'){\\circle*{',f9.3,'}}')
 110  format(' \\put(',f10.3,',',f10.3,'){',i3,'}')
 120  format(' \\put(',f10.3,',',f10.3,'){','$\\cdot$','}')

      end
