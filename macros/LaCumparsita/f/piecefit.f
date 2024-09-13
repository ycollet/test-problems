C     =================================================================
C     File: piecefit.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     July 10th, 2006.

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

C     Piecefit problem
C     ----------------

C     We defined usol, a piecewise linear function defined on [0,4] 
C     with kinks in 1 and 3 and discotinuity in 2. We discretized [0,4] 
C     with n points, including extremes. We defined uobs, a 10% 
C     perturbation of usol. The objective function is the sum of the p 
C     smaller squares of the discretized second derivatives of x. p, 
C     the OVO parameter, is given by the user. The restriction is that 
C     the root mean squared deviation of the solution found with 
C     respect to uobs is smaller than 0.2. In other words, using LOVO, 
C     this code aims to find a piecewise harmonic one-dimensional 
C     function. 

C     ******************************************************************
C     ******************************************************************

      subroutine piecefit_inip(n,x,l,u,m,lambda,rho,equatn,linear,
     &                         n_in,p_in)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n,n_in,p_in

C     ARRAY ARGUMENTS
      logical equatn(*),linear(*)
      double precision l(*),lambda(*),rho(*),u(*),x(*)

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

C     COMMON SCALARS
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i
      double precision drand,seed,t

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
C     Number of variables (points in [0,4])
      n = n_in

C     p (integer close to n-2), perhaps n-3
      p = p_in

C     Generation of uobs

      seed = 17172937.0d0

      h = 4.0d0 / dfloat(n - 1)

      do i = 1,n
          t = dfloat(i-1) * h

          if ( t .le. 1.0d0 ) then
              uobs(i) = 1.0d0 + t

          else if ( t .le. 2.0d0 ) then
              uobs(i) = 1.0d0 - ( t - 2.0d0 )

          else if ( t .le. 3.0d0 ) then
              uobs(i) = 2.0d0

          else 
              uobs(i) = 2.0d0 - ( t - 3.0d0 )
          end if

          usol(i) = uobs(i)
          uobs(i) = uobs(i) + 
     +              uobs(i) * 0.1d0 * ( 2.0d0 * drand(seed) - 1.0d0 )
      end do

C     write(*,*) 'uobs: ',(uobs(i),i=1,n)

C     Initial point

      do i = 1,n
          x(i) = uobs(i)
      end do
      
C     Lower and upper bounds

      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = 1

C     Lagrange multipliers approximation.

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

      end

C     ******************************************************************
C     ******************************************************************

      subroutine piecefit_evalf(n,x,f,flag)

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

C     COMMON SCALARS
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i,imax,j
      double precision zmax

C     LOCAL ARRAYS
      double precision resid(2000)

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
      flag = 0

      f = 0.0d0

      do i = 1,n - 2
          resid(i) = ( x(i+2) - 2.0d0 * x(i+1) + x(i) ) / h ** 2
      end do

      do j = 1,p
          zmax = 1.0d+99
          do i = 1,n - 2
              if ( abs( resid(i) ) .lt. zmax ) then
                  zmax = abs( resid(i) )
                  imax = i
              end if
          end do

          f = f + resid(imax) ** 2
          resid(imax) = 1.0d+99
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine piecefit_evalg(n,x,g,flag)

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

C     COMMON SCALARS
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i,imax,j
      double precision zmax

C     LOCAL ARRAYS
      double precision resid(2000)

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
      flag = 0

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,n - 2
          resid(i) = ( x(i+2) - 2.0d0 * x(i+1) + x(i) ) / h ** 2
      end do

      do j = 1,p
          zmax = 1.0d+99
          do i = 1,n - 2
              if ( abs( resid(i) ) .lt. zmax ) then
                  zmax = abs( resid(i) )
                  imax = i
              end if
          end do

          g(imax+2) = g(imax+2) + 2.0d0 * resid(imax) / h ** 2
          g(imax+1) = g(imax+1) - 4.0d0 * resid(imax) / h ** 2
          g(imax)   = g(imax)   + 2.0d0 * resid(imax) / h ** 2   

          resid(imax) = 1.0d+99
      end do
 
      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine piecefit_evalh(n,x,hlin,hcol,hval,nnzh,flag)

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

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine piecefit_evalc(n,x,ind,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,flag,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subroutine computes the ind-th constraint.
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
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i
      double precision z

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
      flag = 0

      z = 0.0d0
      do i = 1,n
          z = z + ( x(i) - uobs(i) ) ** 2
      end do

      z = z / dfloat(n)

      c = z - 0.01d0 

      end
 
C     ******************************************************************
C     ******************************************************************
 
      subroutine piecefit_evaljac(n,x,ind,indjac,valjac,nnzjac,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,n,nnzjac

C     ARRAY ARGUMENTS
      integer indjac(n)
      double precision x(n),valjac(n)

C     This subroutine computes the gradient of the ind-th constraint. 
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

C     COMMON SCALARS
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
      flag = 0

      nnzjac = n

      do i = 1,n
          indjac(i) = i
          valjac(i) = 2.0d0 * ( x(i) - uobs(i) ) / dfloat(n)
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine piecefit_evalhc(n,x,ind,hclin,hccol,hcval,nnzhc,flag)

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

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine piecefit_evalhlp(n,x,m,lambda,p,hp,goth,flag)

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

      subroutine piecefit_endp(n,x,l,u,m,lambda,rho,equatn,linear)

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

C     COMMON SCALARS
      integer p
      double precision h

C     COMMON ARRAYS
      double precision uobs(2000),usol(2000)

C     LOCAL SCALARS
      integer i
      double precision z,zmax

C     COMMON BLOCKS
      common /ovo/ uobs,usol,h,p
      
      zmax = 0.0d0

      write(*,*) 'i, t, sol-found, true-sol, second-derivative'

      do i = 1,n
          if ( i .eq. 1 .or. i .eq. n ) then
              z = 0.0d0
          else
              z = ( x(i+1) - 2.0d0 * x(i) + x(i-1) ) / h ** 2
          end if

          write(*,*) i,(i-1)*h,x(i),usol(i),z

          zmax = max( zmax, z ) 
      end do

      write(*, *)' Maximum second derivative: ',zmax

      end
