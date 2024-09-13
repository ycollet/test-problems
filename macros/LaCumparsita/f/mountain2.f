C     =================================================================
C     File: mountain2.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     April 10, 2006.

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

C     Mountain Pass problem
C     ---------------------

C     Given a surface S(x,y) and initial and final points pi, pf in R^2 
C     the problem consists on finding a path pi, p1, p2, ..., pN, pf 
C     from pi to pf such that max_{1 <= k <= N} S(pk) is as small as 
C     possible. Moreover, the distance between consecutive points in the 
C     path must be less than or equal to a prescribed tolerance.
C
C     The nonlinear programming formulation of the problem follows:
C
C     min z
C
C     subject to
C
C     d(pi     ,p1)^2 <= dmax^2
C     d(p_{k-1},pk)^2 <= dmax^2, k = 2, ..., N
C     d(pN     ,pf)^2 <= dmax^2
C
C     S(pk) <= z, k = 1, 2, ..., N
C
C     where dmax is the prescribed maximum distance and d(.,.) is the
C     Euclidian distance.
C
C     The problem has n = 2 N + 1 variables and m = 2 N + 1 inequality 
C     constraints, where N is the number of intermediate points in the
C     path.
C
C     In the implementation, pi = (-10,-10), pf = (10,10). N is a user
C     defined parameter. dmax is also defined by the user and must
C     satisfy dmax >= d(pi,pf) / (N + 1). The initial values for p1, p2,
C     ..., pN correspond to a perturbation of the equally spaced points 
C     in the segment [pi,pf]. The intial value of z is z = max { S(pi), 
C     S(p1), S(p2), ..., S(pN), S(pf) }. The considered surfaces are
C
C     Mountain 1: S(x,y) = sin( x * y ) + sin( x + y )
C
C     Mountain 2: S(x,y) = sin(x) + cos(y)

C     References:
C
C     [1] J. Horák, Constrained mountain pass algorithm for the 
C     numerical solution of semilinear elliptic problems, Numerische 
C     Mathematik 98, pp. 251-276, 2004. 
C
C     [2] J. J. Moré and T. S. Munson, Computing mountain passes and 
C     transition states, Mathematical Programming 100, pp. 151-182, 2004.

C     ******************************************************************
C     ******************************************************************

      subroutine mountain2_inip(n,x,l,u,m,lambda,rho,equatn,linear,
     &                          np_in,dmax2_in,seed_in)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n,np_in

C     ARRAY ARGUMENTS
      logical equatn(*),linear(*)
      double precision l(*),lambda(*),rho(*),u(*),x(*)
      double precision dmax2_in,seed_in

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
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     LOCAL SCALARS
      integer b,i,j
      double precision dist,drand,mountain2_fmoun,seed,tmp

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

C     Set problem data

      pi(1) = - 10.0d0
      pi(2) = - 10.0d0

      pf(1) =   10.0d0
      pf(2) =   10.0d0

      dist = sqrt( ( pi(1) - pf(1) ) ** 2 + ( pi(2) - pf(2) ) ** 2 )

C     The distance between the extremes of the path is 'dist'

C     Number of points in the path excluding extremes (maximum 100,000)
      np = np_in

      dist = dist / ( np + 1 ) 

C     The minimum distance between consecutive points in the path is 'dist'

      dmax2 = dmax2_in

      if ( dmax2 .lt. 0.0d0 ) then
          dmax2 = 2.0d0 * dist
      end if

      dmax2 = dmax2 ** 2

C     Seed for the initial point random generation
      seed = seed_in

C     Number of variables

      n = 2 * np + 1

C     Initial point
      x(n) = max( mountain2_fmoun(pi(1),pi(2)),
     +            mountain2_fmoun(pf(1),pf(2)) )
      do i = 1,np
          b = 2 * ( i - 1 )
          do j = 1,2
              tmp    = pi(j) + i * ( pf(j) - pi(j) ) / ( np + 1 )
              x(b+j) = tmp * ( 1.0d0 + 0.1d0 * drand(seed) )  
          end do
          x(n) = max( x(n), mountain2_fmoun(x(b+1),x(b+2)) )
      end do

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = np + 1 + np

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
          equatn(i) = .false.
      end do

C     For each constraint i, set linear(i) = .true. if it is a linear
C     constraint, otherwise set linear(i) = .false.

      do i = 1,m
          linear(i) = .false.
      end do

 100  format(/,1X,'Enter the maximum allowed distance between ',
     +         'consecutive points in the path',/,1X,'(It must be ',
     +         'greater than the minimum distance. Enter a negative ',
     +         'number to',/,1X,'use twice the minimum distance): ')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine mountain2_evalf(n,x,f,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

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
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

      flag = 0

      f = x(n)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine mountain2_evalg(n,x,g,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

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
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

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
 
      subroutine mountain2_evalh(n,x,hlin,hcol,hval,nnzh,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n,nnzh

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

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

      flag = -1

      nnzh = 0

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine mountain2_evalc(n,x,ind,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,flag,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subroutine must compute the ind-th constraint of your 
C     problem. For achieving this objective YOU MUST MOFIFY it 
C     according to your problem. See below the places where your 
C     modifications must be inserted.
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
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     LOCAL SCALARS
      integer b,b1,b2
      double precision mountain2_fmoun

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

      flag = 0

      if ( ind .eq. 1 ) then

          b1 = 2 * ( ind - 1 )
          c  = (x(b1+1) - pi(1))**2   + (x(b1+2) - pi(2))**2   - dmax2

      else if ( ind .le. np ) then

          b1 = 2 * ( ind - 1 )
          b2 = 2 * ( ind - 2 )
          c  = (x(b1+1) - x(b2+1))**2 + (x(b1+2) - x(b2+2))**2 - dmax2

      else if ( ind .eq. np + 1 ) then

          b2 =  2 * ( ind - 2 )
          c  = (pf(1) - x(b2+1))**2   + (pf(2) - x(b2+2))**2   - dmax2

      else if ( ind .le. np + 1 + np ) then

          b = 2 * ( ind - np - 2 )
          c = mountain2_fmoun(x(b+1),x(b+2)) - x(n)

      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine mountain2_evaljac(n,x,ind,indjac,valjac,nnzjac,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,n,nnzjac

C     ARRAY ARGUMENTS
      integer indjac(n)
      double precision x(n),valjac(n)

C     This subroutine must compute the gradient of the constraint ind. 
C     For achieving these objective YOU MUST MODIFY it in the way 
C     specified below.
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
C              the non-null value of the partial derivative of the 
C              ind-th constraint with respect to the indjac(k)-th 
C              variable must be saved at valjac(k).
C
C     flag     integer
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the 
C              constraint. (For example, trying to compute the square 
C              root of a negative number, dividing by zero or a very 
C              small number, etc.) If everything was o.k. you must set 
C              it equal to zero.

C     PARAMETERS
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     LOCAL SCALARS
      integer b,b1,b2,j

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

      flag = 0

      if ( ind .eq. 1 ) then

          nnzjac = 2

          b1 = 2 * ( ind - 1 )

          do j = 1,2
              indjac(j) = b1 + j
              valjac(j) = 2.0d0 * ( x(b1 + j) - pi(j) )
          end do

      else if ( ind .le. np ) then

          nnzjac = 4

          b1 = 2 * ( ind - 1 )
          b2 = 2 * ( ind - 2 )

          do j = 1,2
              indjac(j  ) = b1 + j
              valjac(j  ) =   2.0d0 * ( x(b1 + j) - x(b2 + j) )
              indjac(2+j) = b2 + j
              valjac(2+j) = - 2.0d0 * ( x(b1 + j) - x(b2 + j) )
          end do

      else if ( ind .eq. np + 1 ) then

          nnzjac = 2

          b2 = 2 * ( ind - 2 )

          do j = 1,2
              indjac(j) = b2 + j
              valjac(j) = - 2.0d0 * ( pf(j) - x(b2 + j) )
          end do

      else if ( ind .le. np + 1 + np ) then

          nnzjac = 3

          b = 2 * ( ind - np - 2 )

          call mountain2_gmoun(x(b+1),x(b+2),valjac(1),valjac(2))

          indjac(1) = b + 1
          indjac(2) = b + 2

          indjac(3) = n
          valjac(3) = - 1.0d0

      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine mountain2_evalhc(n,x,ind,hclin,hccol,hcval,nnzhc,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,n,nnzhc

C     ARRAY ARGUMENTS
      integer hccol(*),hclin(*)
      double precision hcval(*),x(n)

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

C     PARAMETERS
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     LOCAL SCALARS
      integer b,b1,b2,j

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

      flag = 0

      if ( ind .eq. 1 ) then

          nnzhc = 2

          b1 = 2 * ( ind - 1 )

          do j = 1,2
              hclin(j) = b1 + j
              hccol(j) = b1 + j
              hcval(j) = 2.0d0
          end do

      else if ( ind .le. np ) then

          nnzhc = 6

          b1 = 2 * ( ind - 1 )
          b2 = 2 * ( ind - 2 )

          do j = 1,2
              hclin(j)   = b1 + j
              hccol(j)   = b1 + j
              hcval(j)   = 2.0d0

              hclin(2+j) = b2 + j
              hccol(2+j) = b2 + j
              hcval(2+j) = 2.0d0

              hclin(4+j) = b1 + j
              hccol(4+j) = b2 + j
              hcval(4+j) = - 2.0d0
          end do

      else if ( ind .eq. np + 1 ) then

          nnzhc = 2

          b2 = 2 * ( ind - 2 )

          do j = 1,2
              hclin(j) = b2 + j
              hccol(j) = b2 + j
              hcval(j) = 2.0d0
          end do

      else if ( ind .le. np + 1 + np ) then

          nnzhc = 3

          b = 2 * ( ind - np - 2 )

          call mountain2_hmoun(x(b+1),x(b+2),hcval(1),hcval(2),hcval(3))

          hclin(1) = b + 1
          hccol(1) = b + 1

          hclin(2) = b + 2
          hccol(2) = b + 2

          hclin(3) = b + 2
          hccol(3) = b + 1

      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine mountain2_evalhlp(n,x,m,lambda,p,hp,goth,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),x(n)

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

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine mountain2_endp(n,x,l,u,m,lambda,rho,equatn,linear)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),rho(m),u(n),x(n)

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

C     PARAMETERS
      integer npmax
      parameter ( npmax = 100000 )

C     COMMON SCALARS
      integer np
      double precision dmax2

C     COMON ARRAYS
      double precision pi(npmax),pf(npmax)

C     LOCAL SCALARS
      integer b,i,imax
      double precision f,fmax,mountain2_fmoun

C     LOCAL ARRAYS
      double precision g(2),h(2,2)

C     COMMON BLOCKS
      common /probdata/ pi,pf,dmax2,np

      write(*,*)

      fmax = - 1.0d+99
      do i = 1,np
          b = 2 * ( i - 1 )
          f = mountain2_fmoun(x(b+1),x(b+2))
          if ( f .gt. fmax ) then
              fmax = f
              imax = i
          end if
      end do

      b = 2 * ( imax - 1 )

      write(*,*) 'Maximizer on optimal path: ',x(b+1),x(b+2)
      write(*,*) 'Maximum: ',fmax

      call mountain2_gmoun(x(b+1),x(b+2),g(1),g(2))
      write(*,*) 'Gradient of the energy: ',g(1),g(2)

      call mountain2_hmoun(x(b+1),x(b+2),h(1,1),h(2,2),h(1,2))
      write(*,*) 'Hessian of the energy: '

      write(*,*) 'Determinant: ',h(1,1) * h(2,2) - h(1,2) ** 2

      if ( imax .ne. 1 ) then
          b = 2 * ( imax - 2 )
          write(*,*) 'Previous point in the path: ',x(b+1),x(b+2)
          write(*,*) 'Energy at previous point in the path: ',
     +                mountain2_fmoun(x(b+1),x(b+2))
      end if

      if ( imax .ne. np ) then
          b = 2 * imax
          write(*,*) 'Following point in the path: ',x(b+1),x(b+2)
          write(*,*) 'Energy at following point in the path: ',
     +                mountain2_fmoun(x(b+1),x(b+2))
      end if

      end

C     ******************************************************************
C     ******************************************************************

      double precision function mountain2_fmoun(x,y)

      implicit none

C     SCALAR ARGUMENTS
      double precision x,y

      mountain2_fmoun = sin( x ) + cos( y )

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine mountain2_gmoun(x,y,dfdx,dfdy)

      implicit none

C     SCALAR ARGUMENTS
      double precision x,y,dfdx,dfdy

      dfdx =   cos( x )
      dfdy = - sin( y )

      end

C     ******************************************************************
C     ******************************************************************

      subroutine mountain2_hmoun(x,y,dfdxx,dfdyy,dfdxy)

      implicit none

C     SCALAR ARGUMENTS
      double precision x,y,dfdxx,dfdyy,dfdxy

      dfdxx = - sin( x )
      dfdyy = - cos( y )
      dfdxy = 0.0d0

      end

