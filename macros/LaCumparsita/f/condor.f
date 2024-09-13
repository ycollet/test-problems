C     =================================================================
C     File: condor.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     July 7th, 2006.

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

C     Condor problem
C     --------------

C     We generate (using subroutine evalu) a function called usol as a 
C     Fourier polynomial with up to q = 20 terms. A grid in [0, 2 pi] 
C     with up to 1000 points is established. We generate function fi(t) 
C     in such a way that usol is a solution of the discretization of 
C     the differential equation u''(t) + [u'(t)]^2 + u(t) = fi(t). We
C     select up to 1000 observation points equally spaced in [-2, 18]
C     and observe the solution usol at these observation points, 
C     obtaining the observed solution uobs. Finally, uint, the 
C     interpolation of uobs in the grid points is computed. The 
C     optimization problem is to obtain the solution of the 
C     discretized differential equation which is closest, in the 
C     least-squares sense, to uint. Instances with q = 20, 500 
C     observation points and 1000 grid points are challenging.

C     ******************************************************************
C     ******************************************************************

      subroutine condor_inip(n,x,l,u,m,lambda,rho,equatn,linear,n_in,
     &                       p_in,q_in)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n,n_in,p_in,q_in

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
      double precision h

C     COMMON ARRAYS
      double precision uint(1000),fon(1000),usol(1000)

C     LOCAL SCALARS
      integer i,j,p,q
      double precision a0,z,pen,aux,u1,u2

C     LOCAL ARRAYS
      double precision a(20),b(20),t(1000),uobs(1000)

C     COMMON BLOCKS
      common /objective/ uint
      common /constraints/ fon,h
      common /solution/ usol
      
C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR PROBLEM
C     DATA:
C     ******************************************************************

C     Number of points in the grid (nmax=1000)
      n = n_in

C     Number of observed point (pmax<=nmax)
      p = p_in

C     Number of terms of the Fourier polynomial (qmax=20)
      q = q_in

C     Fourier polynomial coefficients

      do i = 1,q
          a(i) = cos( dfloat(i) )      / dfloat(i)
          b(i) = sin( dfloat(i) ) ** 2 / dfloat(i)
      end do

      a0 = 2.0d0
      z  = 0.0d0

      do i = 1,p
          t(i) = -2.0d0 + ( 20.0d0 / dfloat(p) ) * i
      end do

      do i = 1,p
          call condor_evalu(a0,a,b,q,t(i),uobs(i),u1,u2,z)
      end do

C     write(*,*) 't = ',(t(i),i=1,p)
C     write(*,*) 'uobs = ',(uobs(i),i=1,p)

C     Independent term of the differential equation (fon)

      h = atan(1.0d0) * 8.0d0 / dfloat(n - 1)

      do i = 1,n - 2
          z = i * h
          call condor_evalu(a0,a,b,q,z,aux,u1,u2,fon(i))
      end do

C     Analytic solution at each grid point (usol)

      do i = 1,n
          z = ( i - 1 ) * h
          call condor_evalu(a0,a,b,q,z,usol(i),u1,u2,aux)
      end do
 
C     Interpolation

      do i = 1,n
          z = dfloat(i - 1) * h

          if( z .le. t(2) ) then
              pen = ( uobs(2) - uobs(1) ) / ( t(2) - t(1) )
              uint(i) = uobs(1) + pen * ( z - t(1) )

          else if( z .ge. t(p-1) )then
              pen = ( uobs(p) - uobs(p-1) ) / ( t(p) - t(p-1) )
              uint(i) = uobs(p-1) + pen * ( z - t(p-1) )

          else
              do j = 2,p - 2
                  if ( z .ge. t(j) .and. z .le. t(j+1) ) then
                      pen = ( uobs(j+1) - uobs(j) ) / ( t(j+1) - t(j) )
                      uint(i) = uobs(j) + pen * ( z - t(j) )
                  end if
              end do
          end if
      end do      

C     write(*,*) 'Interpolated points: '
C     do i = 1,n
C         write(*,*) i,(i-1)*h,uint(i)
C     end do
      
C     Initial point

      do i = 1,n
          x(i) = uint(i)
      end do
      
C     Lower and upper bounds

      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = n - 2

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

      subroutine condor_evalu(a0,a,b,m,x,u,u1,u2,fi)

      implicit none

C     SCALAR ARGUMENTS
      integer m
      double precision a0,fi,u,u1,u2,x

C     ARRAY ARGUMENTS
      double precision a(20),b(20)

C     LOCAL SCALARS
      integer k

      u  = a0
      u1 = 0.0d0
      u2 = 0.0d0

      do k = 1,m
          u  = u  + a(k) * cos(k*x)         + b(k) * sin(k*x)
          u1 = u1 - a(k) * k * sin(k*x)     + b(k) * k * cos(k*x)
          u2 = u2 - a(k) * k * k * cos(k*x) - b(k) * k * k * sin(k*x)
      end do

      fi = u2 + u1 ** 2 + u

      end
      
C     ******************************************************************
C     ******************************************************************

      subroutine condor_evalf(n,x,f,flag)

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

C     COMMON ARRAYS
      double precision uint(1000)

C     LOCAL SCALAR
      integer i

C     COMMON BLOCKS
      common /objective/ uint

C     Objective function

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR OBJECTIVE
C     FUNCTION:
C     ******************************************************************

      flag = 0

      f = 0.0d0
      do i = 1,n
          f = f + ( x(i) - uint(i) ) ** 2
      end do

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALF.
C     ******************************************************************

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine condor_evalg(n,x,g,flag)

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

C     COMMON ARRAYS
      double precision uint(1000)

C     LOCAL SCALAR
      integer i

C     COMMON BLOCKS
      common /objective/ uint

C     Gradient vector

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENT
C     VECTOR OF YOUR OBJECTIVE FUNCTION: 
C     ******************************************************************

      flag = 0

      do i = 1,n 
          g(i) = 2.0d0 * ( x(i) - uint(i) ) 
      end do

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALG. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine condor_evalh(n,x,hlin,hcol,hval,nnzh,flag)

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

      subroutine condor_evalc(n,x,ind,c,flag)

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
 
C     COMMON SCALARS
      double precision h

C     COMMON ARRAYS
      double precision fon(1000)

C     COMMON BLOCKS
      common /constraints/ fon,h
      
C     Constraints

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR
C     CONSTRAINTS: 
C     ******************************************************************

      flag = 0

      c = ( x(ind+2) - 2.0d0 * x(ind+1) + x(ind) ) / h ** 2 + 
     +    ( ( x(ind+2) - x(ind) ) / ( 2.0d0 * h ) ) ** 2 + 
     +    x(ind+1) - fon(ind)

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALC. 
C     ******************************************************************
 
      end
 
C     ******************************************************************
C     ******************************************************************
 
      subroutine condor_evaljac(n,x,ind,indjac,valjac,nnzjac,flag)

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

C     COMMON SCALARS
      double precision h

C     COMMON ARRAYS
      double precision fon(1000)

C     COMMON BLOCKS
      common /constraints/ fon,h
      
C     Sparse gradient vector of the ind-th constraint

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENTS  
C     OF YOUR CONSTRAINTS: 
C     ******************************************************************

      flag = 0

      nnzjac = 3

      indjac(1) = ind
      valjac(1) = 1.0d0 / h ** 2 - 
     +            2.0d0 * ( x(ind+2) - x(ind) ) / ( 2.0d0 * h ) ** 2 

      indjac(2) = ind + 1
      valjac(2) = - 2.0d0 / h ** 2 + 1.0d0

      indjac(3) = ind + 2

      valjac(3) = 1.0d0 / h ** 2 + 
     +            2.0d0 * ( x(ind+2) - x(ind) ) / ( 2.0d0 * h ) ** 2 

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALJAC. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine condor_evalhc(n,x,ind,hclin,hccol,hcval,nnzhc,flag)

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

      subroutine condor_evalhlp(n,x,m,lambda,p,hp,goth,flag)

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

C     ******************************************************************
C     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO CUMPUTE 
C     THE PRODUCT OF THE HESSIANS OF THE LAGRANGIAN TIME VECTOR P: 
C     ******************************************************************

      flag = - 1

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALHC. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************

      subroutine condor_endp(n,x,l,u,m,lambda,rho,equatn,linear)

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

C     COMMON ARRAYS
      double precision usol(1000)

C     LOCAL SCALARS
      integer i
      double precision zmax

C     COMMON BLOCKS
      common /solution/ usol
      
C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO COMPLEMENT THE
C     INFORMATION RELATED TO THE SOLUTION GIVEN BY THE SOLVER
C     ******************************************************************

      zmax = 0.0d0
      do i = 1,n
          zmax = max( zmax, abs( x(i) - usol(i) ) )
      end do

      write(*,*)
      write(*,*) 'Comparison between analytic and computed solution: '

      write(*,*) 'Analytic Computed Difference'
      do i = 1,n
          write(*,*) usol(i),x(i),usol(i) - x(i)
      end do

      write(*, *) 'Maximum error: ',zmax        

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE ENDP
C     ******************************************************************

      end
