C     =================================================================
C     File: genpack-csq-mina.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     July 3rd, 2006.

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

C     PACKING CIRCLES INTO A SQUARE
C     -----------------------------

C     Given a fixed number of fixed-size identical circles (called
C     items), the objective is to find the smallest square (called 
C     object) that fits the  items without overlapping. The model 
C     follows:
C
C     Mininimize L
C     subject to ci in [r,L-r] x [r,L-r], for all i
C                sum_{i < j} max{ 0, (2r)^2 - d(ci,cj)^2 }^2 = 0
C
C     (The variables are the centers ci in R^2 of the circular items 
C     and the side L of the squared object that fits the items.)
C
C     [1] E. G. Birgin and F. N. C. Sobral, "A tool based on 
C     nonlinear programming models for packing 2D and 3D circular 
C     items", Computers and Operations Research, 2006 
C     (DOI: 10.1016/j.cor.2006.11.002).
C
C     [2] J. M. Martínez and D. P. Ronconi, "Optimizing the 
C     Packing of Cylinders into a Rectangular Container: A Nonlinear 
C     Approach", European Journal of Operational Research 160, pp. 
C     19-33, 2005.

C     ******************************************************************
C     ******************************************************************

      subroutine genpack_csq_mina_inip(n,x,l,u,m,lambda,rho,equatn,
     +                                 linear,iterad_in,nite_in,seed_in)

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
      integer nitemax,nregmax
      double precision pi
      parameter ( nitemax =        10000 )
      parameter ( nregmax =        10000 )
      parameter ( pi      =  3.1415927d0 )

C     COMMON SCALARS
      integer nbr,ndim,nite,nreg,nite_in
      double precision iterad,iterad_in,seed_in

C     COMMON ARRAYS
      integer start(0:nregmax+1,0:nregmax+1),next(nitemax),br(2,nitemax)

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(*),linear(*)
      double precision l(*),lambda(*),rho(*),u(*),x(*)

C     LOCAL SCALARS
      integer i,ind,j
      double precision objsizeub,seed

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad
      common /partdata/ start,next,br,nbr,nreg

C     FUNCTIONS
      integer genpack_csq_mina_ceil
      double precision drand

      ndim = 2

      write(*,FMT=7000)

C     Read problem data

C     Enter the circular items radius (iterad)
C     and the number of items (nite)
      iterad = iterad_in
      nite   = nite_in

C     Seed for the initial point random generation
      seed = seed_in

C     Initialize regions

      objsizeub = sqrt( 2.0d0 * pi * iterad ** 2 * nite )

      nreg = genpack_csq_mina_ceil( objsizeub / iterad )
      nreg = min( nreg, nregmax )

      do i = 0,nreg + 1
          do j = 0,nreg + 1
              start(i,j) = 0
          end do 
      end do

C     Number of variables

      n = ndim * nite + 1

C     Initial point

      x(n) = sqrt( iterad ** 2 * pi * nite * ( 1.0d0 + drand(seed) ) )

      i = 1
      do i = 1,nite
          do j = 1,ndim
              ind = (i - 1) * ndim + j
              x(ind) = iterad + drand(seed) * ( x(n) - 2.0d0 * iterad )
          end do
      end do

C     Lower and upper bounds
      do i = 1,n
          l(i) = iterad
          u(i) = 1.0d+20
      end do

      l(n) = 2.0d0 * iterad

C     Number of constraints (equalities plus inequalities)

      m = ndim * nite + 1

C     Lagrange multipliers approximation 

      do i = 1,m
          lambda(i) = 0.0d0
      end do

C     Initial penalty parameters

      do i = 1,m
          rho(i) = 10.0d0
      end do

C     For each constraint i, set equatn(i) = .true. if it is an equality
C     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if 
C     it is an inequality constraint of the form c_i(x) <= 0.

      do i = 1,m - 1
          equatn(i) = .false.
      end do
      equatn(m) = .true.

C     For each constraint i, set linear(i) = .true. if it is a linear
C     constraint, otherwise set linear(i) = .false.

      do i = 1,m - 1
          linear(i) = .true.
      end do
      linear(m) = .false.

      return

 7000 format(/,' Welcome to GENPACK !!!',/
     +       /,' This program uses ALGENCAN and the strategy',
     +         ' described in:',/
     +       /,' E. G. Birgin and F. N. C. Sobral, A tool based on',
     +         ' nonlinear programming',/,' models for packing 2D and',
     +         ' 3D circular items, Computers and Operations',
     +       /,' Research, 2006 (DOI: 10.1016/j.cor.2006.11.002).',/,
     +       /,' (See also: E. G. Birgin, J. M. Martínez and',
     +         ' D. P. Ronconi, Optimizing the',/,' Packing of',
     +         ' Cylinders into a Circular Container: A Nonlinear',
     +         ' Approach,',/,' European Journal of Operational',
     +         ' Research 160, pp. 19-33, 2005.)',/,
     +       /,' to solve the problem of packing a fixed number of',
     +         ' (fixed-size) identical',
     +       /,' circles into a squared region, minimizing the size',
     +         ' of the region.',/)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine genpack_csq_mina_evalf(n,x,f,flag)

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

C     SCALAR ARGUMENTS
      integer flag,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

      flag = 0

      f = x(n)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_evalg(n,x,g,flag)

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

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     LOCAL SCALARS
      integer i

      flag = 0

      do i = 1,n - 1
          g(i) = 0.0d0
      end do

      g(n) = 1.0d0

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_evalh(n,x,hlin,hcol,hval,nnzh,flag)

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

      flag = -1

      nnzh = 0

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_evalc(n,x,ind,c,flag)

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

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad

C     SCALAR ARGUMENTS
      integer ind,flag,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer m

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad

C     EXTERNAL FUNCTIONS
      double precision genpack_csq_mina_belong,genpack_csq_mina_overlap

      flag = 0
      m = ndim * nite + 1

      if ( 1 .le. ind .and. ind .le. m - 1 ) then
          c = genpack_csq_mina_belong(n,x,ind)
      else if ( ind .eq. m ) then
          c = genpack_csq_mina_overlap(n,x)
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_evaljac(n,x,ind,indjac,valjac,nnzjac,
     +                                    flag)

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

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad

C     SCALAR ARGUMENTS
      integer flag,ind,n,nnzjac

C     ARRAY ARGUMENTS
      integer indjac(n)
      double precision x(n),valjac(n)

C     LOCAL SCALARS
      integer i,m

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad

      flag = 0
      m = ndim * nite + 1

      if ( 1 .le. ind .and. ind .le. m - 1 ) then
          call genpack_csq_mina_gbelong(n,x,ind,indjac,valjac,nnzjac)
      else if ( ind .eq. m ) then
          call genpack_csq_mina_goverlap(n,x,valjac)
          nnzjac = n - 1
          do i = 1,nnzjac
              indjac(i) = i
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_evalhc(n,x,ind,hclin,hccol,hcval,
     +                                   nnzhc,flag)

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

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad

C     SCALAR ARGUMENTS
      integer flag,ind,n,nnzhc

C     ARRAY ARGUMENTS
      integer hccol(*),hclin(*)
      double precision hcval(*),x(n)

C     LOCAL SCALARS
      integer m

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad

      flag = 0
      m = ndim * nite + 1

      if ( 1 .le. ind .and. ind .le. m - 1 ) then
          call genpack_csq_mina_hbelong(n,x,ind,hclin,hccol,hcval,nnzhc)
      else if ( ind .eq. m ) then
          call genpack_csq_mina_hoverlap(n,x,hclin,hccol,hcval,nnzhc)
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine genpack_csq_mina_evalhlp(n,x,m,lambda,p,hp,goth,flag)

C     This subroutine might compute the product of the Hessian of the
C     Lagrangian times vector p. For achieving this objective YOU MAY 
C     MODIFY it according to your problem. To modify this subroutine IS 
C     NOT MANDATORY. See below where your modifications must be 
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

C     ******************************************************************
C     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO COMPUTE 
C     THE PRODUCT OF THE HESSIAN OF THE LAGRANGIAN TIMES A VECTOR 
C     ******************************************************************

      flag = -1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine genpack_csq_mina_endp(n,x,l,u,m,lambda,rho,equatn,
     +                                 linear)

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

C     COMMON SCALARS
      integer ndim,nite
      double precision iterad

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),rho(m),u(n),x(n)

C     LOCAL SCALARS
      double precision vover,valoc

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad

C     FUNCTIONS
      double precision genpack_csq_mina_maxover,genpack_csq_mina_maxaloc

      vover = genpack_csq_mina_maxover(n,x)
      valoc = genpack_csq_mina_maxaloc(n,x)

      write(*,FMT=0001)x(n),vover,valoc

      call genpack_csq_mina_drawsol(nite,iterad,n,x,'solution.mp')
      call genpack_csq_mina_savesol(nite,iterad,ndim,n,x,vover,valoc,
     +                              'solution.dat')

 0001 format(/,' Side                          = ',F12.8,/,
     +         ' Maximum overlapping violation = ',3X,1P,D9.1,/,
     +         ' Maximum allocation violation  = ',3X,1P,D9.1)

      end

C     ******************************************************************
C     ******************************************************************

      double precision function genpack_csq_mina_maxover(n,x)

C     PARAMETERS
      integer nitemax,nregmax
      parameter ( nitemax =        10000 )
      parameter ( nregmax =        10000 )

C     COMMON SCALARS
      integer nbr,nreg

C     COMMON ARRAYS
      integer start(0:nregmax+1,0:nregmax+1),next(nitemax),br(2,nitemax)

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer i,j,k,l
      double precision genpack_csq_mina_mnrover

C     COMMON BLOCKS
      common /partdata/ start,next,br,nbr,nreg

      call genpack_csq_mina_classify(n,x)
 
      genpack_csq_mina_maxover = 0.0d0
      do k = 1,nbr
          i = br(1,k)
          j = br(2,k)
          l = start(i,j)
 10       if ( l .ne. 0 ) then
              genpack_csq_mina_maxover = max( genpack_csq_mina_maxover,
     +                 genpack_csq_mina_mnrover(n,x,l,next(l)) )
              genpack_csq_mina_maxover = max( genpack_csq_mina_maxover,
     +                 genpack_csq_mina_mnrover(n,x,l,start(i-1,j+1)) )
              genpack_csq_mina_maxover = max( genpack_csq_mina_maxover,
     +                 genpack_csq_mina_mnrover(n,x,l,start(i,  j+1)) )
              genpack_csq_mina_maxover = max( genpack_csq_mina_maxover,
     +                 genpack_csq_mina_mnrover(n,x,l,start(i+1,j+1)) )
              genpack_csq_mina_maxover = max( genpack_csq_mina_maxover,
     +                 genpack_csq_mina_mnrover(n,x,l,start(i+1,j))   )
              l = next(l)
              go to 10
          end if
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      double precision function genpack_csq_mina_mnrover(n,x,i,lstart)

C     PARAMETERS
      integer nitemax,nregmax
      parameter ( nitemax =        10000 )
      parameter ( nregmax =        10000 )

C     COMMON SCALARS
      integer nbr,ndim,nite,nreg
      double precision iterad

C     COMMON ARRAYS
      integer start(0:nregmax+1,0:nregmax+1),next(nitemax),br(2,nitemax)

C     SCALAR ARGUMENTS
      integer i,lstart,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer ind1,ind2,j,k
      double precision dist

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad
      common /partdata/ start,next,br,nbr,nreg

      genpack_csq_mina_mnrover = 0.0d0

      j = lstart

 10   if ( j .ne. 0 ) then

          dist = 0.0d0
          do k = 1,ndim
              ind1 = ( i - 1 ) * ndim + k
              ind2 = ( j - 1 ) * ndim + k
              dist = dist + ( x(ind1) - x(ind2) ) ** 2
          end do
          dist = sqrt( dist )

          genpack_csq_mina_mnrover = max( genpack_csq_mina_mnrover,
     +                                    2.0d0 * iterad - dist )

          j = next(j)
          go to 10

      end if

      end

C     ******************************************************************
C     ******************************************************************

      double precision function genpack_csq_mina_maxaloc(n,x)

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer i,ind,k

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad

      genpack_csq_mina_maxaloc = 0.0d0

      do i = 1,nite

          do k = 1,ndim
              ind = ( i - 1 ) * ndim + k
              genpack_csq_mina_maxaloc = max( genpack_csq_mina_maxaloc,
     +                                        iterad     -    x(ind) )
              genpack_csq_mina_maxaloc = max( genpack_csq_mina_maxaloc,
     +                                        x(ind) - x(n) + iterad )
          end do

      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine genpack_csq_mina_drawsol(nite,iterad,n,x,solfile)

C     This subroutine generate a metapost file with the
C     graphical representation of the problem.

C     SCALAR ARGUMENTS
      integer n,nite
      double precision iterad

C     ARRAY ARGUMENTS
      double precision x(n)
      character * 11 solfile

C     LOCAL SCALARS
      integer i
      double precision scale

      open(unit=10,file=solfile)

C     SCALING
      scale = min( 0.5d0, 15.0d0 / x(n) )
      write(10,10) scale

C     CIRCULAR OBJECT
      write(10,40)x(n),x(n),x(n),x(n)

C     CIRCULAR ITEMS
      do i = 1,nite
          write(10,20) 2.0d0*iterad,2.0d0*iterad,x(2*i-1),x(2*i)
      end do

      write(10,30)

      close(10)

C     NON-EXECUTABLE STATEMENTS

  10  format('beginfig(1);'/,'u = ',f20.10,' cm;') 
  20  format('draw fullcircle',
     +       /,'     xscaled  ',f20.10,'u',
     +       /,'     yscaled  ',f20.10,'u',
     +       /,'     shifted (',f20.10,'u,',f20.10,'u);')
  30  format('endfig;',/,'end;') 
  40  format('draw (0,0)--(0,',F20.10,'u)--(',F20.10,'u,',F20.10,'u)',
     +       '--(',F20.10,'u,0)--cycle;',/)
      end

C     ******************************************************************
C     ******************************************************************

      subroutine genpack_csq_mina_savesol(nite,iterad,ndim,n,x,vover,
     +                                    valoc,solfile)

C     SCALAR ARGUMENTS
      integer n,ndim,nite
      double precision iterad,vover,valoc

C     ARRAY ARGUMENTS
      double precision x(n)
      character * 12 solfile

C     LOCAL SCALARS
      integer i,j

      open(unit=10,file=solfile)

C     PROBLEM DATA
      write(10,10) x(n),iterad,nite,vover,valoc

C     ITEMS LOCATION
      do i = 1,nite
          write(10,20) i,(x((i-1)*ndim+j),j=1,ndim)
      end do

      close(10)

C     NON-EXECUTABLE STATEMENTS

  10  format(/,'Size of the squared object    = ',F16.10,
     +       /,'Circular items radius         = ',F16.10,
     +       /,'Number of packed items        = ',8X,I8,
     +       /,'Maximum overlapping violation = ',7X,1P,D9.1,
     +       /,'Maximum allocation violation  = ',7X,1P,D9.1,/
     +       /,'Items location:',/) 
  20  format(I8,3(1X,1P,D16.8))

      end

C     ******************************************************************
C     ******************************************************************

      double precision function genpack_csq_mina_overlap(n,x)

C     PARAMETERS
      integer nitemax,nregmax
      parameter ( nitemax =        10000 )
      parameter ( nregmax =        10000 )

C     COMMON SCALARS
      integer nbr,nreg

C     COMMON ARRAYS
      integer start(0:nregmax+1,0:nregmax+1),next(nitemax),br(2,nitemax)

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer i,j,k,l
      double precision nrdist

C     COMMON BLOCKS
      common /partdata/ start,next,br,nbr,nreg

C     CLASSIFY ITEMS BY REGIONS
      call genpack_csq_mina_classify(n,x)
 
C     COMPUTE SPARSE OVERLAPPING
      genpack_csq_mina_overlap = 0.0d0
      do k = 1,nbr
          i = br(1,k)
          j = br(2,k)
          l = start(i,j)
 10       if ( l .ne. 0 ) then
              genpack_csq_mina_overlap = 
     +        genpack_csq_mina_overlap + nrdist(n,x,l,next(l))        + 
     +                                   nrdist(n,x,l,start(i-1,j+1)) + 
     +                                   nrdist(n,x,l,start(i,  j+1)) + 
     +                                   nrdist(n,x,l,start(i+1,j+1)) + 
     +                                   nrdist(n,x,l,start(i+1,j))
              l = next(l)
              go to 10
          end if
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine genpack_csq_mina_goverlap(n,x,g)

C     PARAMETERS
      integer nitemax,nregmax
      parameter ( nitemax =        10000 )
      parameter ( nregmax =        10000 )

C     COMMON SCALARS
      integer nbr,nreg

C     COMMON ARRAYS
      integer start(0:nregmax+1,0:nregmax+1),next(nitemax),br(2,nitemax)

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     LOCAL SCALARS
      integer i,j,k,l

C     COMMON BLOCKS
      common /partdata/ start,next,br,nbr,nreg

C     CLASSIFY ITEMS BY REGIONS
      call genpack_csq_mina_classify(n,x)
 
C     COMPUTE SPARSE OVERLAPPING

      do i = 1,n
          g(i) = 0.0d0
      end do

      do k = 1,nbr
          i = br(1,k)
          j = br(2,k)
          l = start(i,j)
 10       if ( l .ne. 0 ) then
              call genpack_csq_mina_gnrdist(n,x,g,l,next(l)) 
              call genpack_csq_mina_gnrdist(n,x,g,l,start(i-1,j+1))
              call genpack_csq_mina_gnrdist(n,x,g,l,start(i,  j+1))  
              call genpack_csq_mina_gnrdist(n,x,g,l,start(i+1,j+1))  
              call genpack_csq_mina_gnrdist(n,x,g,l,start(i+1,j  ))
              l = next(l)
              go to 10
          end if
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_hoverlap(n,x,hlin,hcol,hval,nnzh)

C     PARAMETERS
      integer nitemax,nregmax,ndimmax
      parameter ( nitemax =        10000 )
      parameter ( nregmax =        10000 )
      parameter ( ndimmax =            2 )

C     COMMON SCALARS
      integer nbr,nreg,nite,ndim
      double precision iterad

C     COMMON ARRAYS
      integer start(0:nregmax+1,0:nregmax+1),next(nitemax),br(2,nitemax)
      double precision diagb(ndimmax,ndimmax,nitemax)

C     SCALAR ARGUMENTS
      integer n,nnzh

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     LOCAL SCALARS
      integer i,j,k,l

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad
      common /partdata/ start,next,br,nbr,nreg
      common /diblocks/ diagb

C     CLASSIFY ITEMS BY REGIONS
      call genpack_csq_mina_classify(n,x)

C     INITALIZE DIAGONAL BLOCKS
      do k = 1,nite
          do j = 1,ndim
              do i = j,ndim
                  diagb(i,j,k) = 0.0d0
              end do
          end do
      end do

C     COMPUTE SPARSE OVERLAPPING SECOND DERIVATIVES

      nnzh = 0

      do k = 1,nbr
        i = br(1,k)
        j = br(2,k)
        l = start(i,j)
 10     if ( l .ne. 0 ) then
         call genpack_csq_mina_hnrdist(n,x,hlin,hcol,hval,nnzh,
     +                                 l,next(l)) 
         call genpack_csq_mina_hnrdist(n,x,hlin,hcol,hval,nnzh,
     +                                 l,start(i-1,j+1))
         call genpack_csq_mina_hnrdist(n,x,hlin,hcol,hval,nnzh,
     +                                 l,start(i,  j+1))  
         call genpack_csq_mina_hnrdist(n,x,hlin,hcol,hval,nnzh,
     +                                 l,start(i+1,j+1))  
         call genpack_csq_mina_hnrdist(n,x,hlin,hcol,hval,nnzh,
     +                                 l,start(i+1,j  ))
         l = next(l)
         go to 10
        end if
      end do

C     ADD DIAGONAL BLOCKS TO THE SPARSE STRUCTURE
      do k = 1,nite
          do j = 1,ndim
              do i = j,ndim
                  if ( diagb(i,j,k) .ne. 0.0d0 ) then
                      nnzh = nnzh + 1
                      hlin(nnzh) = ( k - 1 ) * ndim + i
                      hcol(nnzh) = ( k - 1 ) * ndim + j
                      hval(nnzh) = diagb(i,j,k)
                  end if
              end do
          end do
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      double precision function genpack_csq_mina_nrdist(n,x,i,lstart)

C     PARAMETERS
      integer nitemax,nregmax
      parameter ( nitemax =        10000 )
      parameter ( nregmax =        10000 )

C     COMMON SCALARS
      integer nbr,ndim,nite,nreg
      double precision iterad

C     COMMON ARRAYS
      integer start(0:nregmax+1,0:nregmax+1),next(nitemax),br(2,nitemax)

C     SCALAR ARGUMENTS
      integer i,lstart,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer ind1,ind2,j,k
      double precision dist

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad
      common /partdata/ start,next,br,nbr,nreg

      genpack_csq_mina_nrdist = 0.0d0

      j = lstart

 10   if ( j .ne. 0 ) then

          dist = 0.0d0
          do k = 1,ndim
              ind1 = ( i - 1 ) * ndim + k
              ind2 = ( j - 1 ) * ndim + k
              dist = dist + ( x(ind1) - x(ind2) ) ** 2
          end do

          dist = max( 0.0d0, ( 2.0d0 * iterad ) ** 2 - dist )
          genpack_csq_mina_nrdist = genpack_csq_mina_nrdist + dist ** 2

          j = next(j)
          go to 10

      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_gnrdist(n,x,g,i,lstart)

C     PARAMETERS
      integer nitemax,nregmax
      parameter ( nitemax =        10000 )
      parameter ( nregmax =        10000 )

C     COMMON SCALARS
      integer nbr,ndim,nite,nreg
      double precision iterad

C     COMMON ARRAYS
      integer start(0:nregmax+1,0:nregmax+1),next(nitemax),br(2,nitemax)

C     SCALAR ARGUMENTS
      integer i,lstart,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     LOCAL SCALARS
      integer ind1,ind2,j,k
      double precision dist,tmp

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad
      common /partdata/ start,next,br,nbr,nreg

      j = lstart

 10   if ( j .ne. 0 ) then

          dist = 0.0d0
          do k = 1,ndim
              ind1 = ( i - 1 ) * ndim + k
              ind2 = ( j - 1 ) * ndim + k
              dist = dist + ( x(ind1) - x(ind2) ) ** 2
          end do

          dist = max( 0.0d0, ( 2.0d0 * iterad ) ** 2 - dist )

          if ( dist .ne. 0.d0 ) then
              do k = 1,ndim
                  ind1 = ( i - 1 ) * ndim + k
                  ind2 = ( j - 1 ) * ndim + k
                  tmp = 4.0d0 * dist * ( x(ind1) - x(ind2) )
                  g(ind1) = g(ind1) - tmp 
                  g(ind2) = g(ind2) + tmp
              end do
          end if

          j = next(j)
          go to 10

      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_hnrdist(n,x,hlin,hcol,hval,nnzh,
     +                                    i,lstart)

C     PARAMETERS
      integer nitemax,nregmax,ndimmax
      parameter ( nitemax =        10000 )
      parameter ( nregmax =        10000 )
      parameter ( ndimmax =            2 )

C     COMMON SCALARS
      integer nbr,ndim,nite,nreg
      double precision iterad

C     COMMON ARRAYS
      integer start(0:nregmax+1,0:nregmax+1),next(nitemax),br(2,nitemax)
      double precision diagb(ndimmax,ndimmax,nitemax)

C     SCALAR ARGUMENTS
      integer i,lstart,n,nnzh

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     LOCAL SCALARS
      integer col,ind1,ind2,ind3,ind4,j,k,l,lin
      double precision dist,tmp

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad
      common /partdata/ start,next,br,nbr,nreg
      common /diblocks/ diagb

      j = lstart

 10   if ( j .ne. 0 ) then

          if ( j .gt. i ) then
              col = i
              lin = j
          else
              col = j
              lin = i
          end if

          dist = 0.0d0
          do k = 1,ndim
              ind1 = ( col - 1 ) * ndim + k
              ind2 = ( lin - 1 ) * ndim + k
              dist = dist + ( x(ind1) - x(ind2) ) ** 2
          end do

          if ( dist .le. ( 2.0d0 * iterad ) ** 2 ) then
              do k = 1,ndim
                  ind1 = ( col - 1 ) * ndim + k
                  ind2 = ( lin - 1 ) * ndim + k
                  tmp = 8.0d0 * ( x(ind1) - x(ind2) ) ** 2
     +                - 4.0d0 * ( ( 2.0d0 * iterad ) ** 2 - dist )
                  if ( tmp .ne. 0.0d0 ) then
C                     H(ind1,ind1) = H(ind1,ind1) + tmp
                      diagb(k,k,col) = diagb(k,k,col) + tmp
C                     H(ind2,ind2) = H(ind2,ind2) + tmp
                      diagb(k,k,lin) = diagb(k,k,lin) + tmp
C                     H(ind2,ind1) = H(ind2,ind1) - tmp
                      nnzh = nnzh + 1
                      hlin(nnzh) = ind2
                      hcol(nnzh) = ind1
                      hval(nnzh) = - tmp
                  end if
                  do l = 1,k - 1
                      ind3 = ( col - 1 ) * ndim + l
                      ind4 = ( lin - 1 ) * ndim + l
                      tmp = 8.0d0 * ( x(ind3) - x(ind4) )
     +                            * ( x(ind1) - x(ind2) )
                      if ( tmp .ne. 0.0d0 ) then
C                         H(ind1,ind3) = H(ind1,ind3) + tmp
                          diagb(k,l,col) = diagb(k,l,col) + tmp
C                         H(ind2,ind4) = H(ind2,ind4) + tmp
                          diagb(k,l,lin) = diagb(k,l,lin) + tmp
C                         H(ind2,ind3) = H(ind2,ind3) - tmp
                          nnzh = nnzh + 1
                          hlin(nnzh) = ind2
                          hcol(nnzh) = ind3
                          hval(nnzh) = - tmp
                      end if
                  end do
                  do l = k + 1,ndim
                      ind3 = ( col - 1 ) * ndim + l
                      ind4 = ( lin - 1 ) * ndim + l
                      tmp = 8.0d0 * ( x(ind3) - x(ind4) )
     +                            * ( x(ind1) - x(ind2) )
                      if ( tmp .ne. 0.0d0 ) then
C                         H(ind2,ind3) = H(ind2,ind3) - tmp
                          nnzh = nnzh + 1
                          hlin(nnzh) = ind2
                          hcol(nnzh) = ind3
                          hval(nnzh) = - tmp
                      end if
                  end do
              end do
          end if

          j = next(j)
          go to 10

      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_classify(n,x)

C     PARAMETERS
      integer nitemax,nregmax
      parameter ( nitemax =        10000 )
      parameter ( nregmax =        10000 )

C     COMMON SCALARS
      integer nbr,ndim,nite,nreg
      double precision iterad

C     COMMON ARRAYS
      integer start(0:nregmax+1,0:nregmax+1),next(nitemax),br(2,nitemax)

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer i,j,k

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad
      common /partdata/ start,next,br,nbr,nreg

C     CLEAN-UP THE START STRUCTURE

      do k = 1,nbr
          i = br(1,k)
          j = br(2,k)
          start(i,j) = 0
      end do

C     FILL-IN THE START STRUCTURE AGAIN

      nbr = 0
      do k = 1,nite
          call genpack_csq_mina_region(nreg,x((k-1)*ndim+1),
     +                                 x((k-1)*ndim+2),i,j)
          if ( start(i,j) .eq. 0 ) then
              nbr = nbr + 1
              br(1,nbr) = i
              br(2,nbr) = j
          end if
          next(k) = start(i,j)
          start(i,j) = k
      end do
 
      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_region(nreg,x,y,i,j)

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad

C     SCALAR ARGUMENTS
      integer i,j,nreg
      double precision x,y

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad

      i = 1 + int( x / ( 2.0d0 * iterad ) )
      i = min( max( 1, i ), nreg )

      j = 1 + int( y / ( 2.0d0 * iterad ) )
      j = min( max( 1, j ), nreg )

      end

C     ******************************************************************
C     ******************************************************************
 
      integer function genpack_csq_mina_ceil(x)

C     SCALAR ARGUMENTS
      double precision x

      if ( x .eq. int( x ) ) then
          genpack_csq_mina_ceil = int( x )
      else
          genpack_csq_mina_ceil = int( x + 1 )
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      double precision function genpack_csq_mina_belong(n,x,ind)

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad

C     SCALAR ARGUMENTS
      integer ind,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad

      genpack_csq_mina_belong = x(ind) + iterad - x(n)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_gbelong(n,x,ind,gind,gval,nnzg)

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad

C     SCALAR ARGUMENTS
      integer ind,n,nnzg

C     ARRAY ARGUMENTS
      integer gind(n)
      double precision x(n),gval(n)

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad

      nnzg = 2

      gind(1)    =     ind
      gval(1)    =   1.0d0

      gind(nnzg) =       n
      gval(nnzg) = - 1.0d0

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine genpack_csq_mina_hbelong(n,x,ind,hlin,hcol,hval,nnzh)

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad

C     SCALAR ARGUMENTS
      integer ind,n,nnzh

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad

      nnzh = 0

      end
