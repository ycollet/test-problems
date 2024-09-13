C     =================================================================
C     File: pedazos4.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     July 20th, 2006.

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

C     Pedazos4 problem: Fit a continuous piecewise polynomial
C     -------------------------------------------------------

C     We generate a table (t_i, y_i) with 30 equally spaced points 
C     between 1 and 30. The first 10 points are on a line, the second 
C     10 on a parabola and the third 10 on a cubic. The overall 
C     function is continuous at 10 and 20. We fit a line to the first 
C     10, a quadratic to the second 10 and a cubic to the third ten, 
C     with the conditions that we preserve continuity at 10 and 20. So,
C     we have 9 unknowns and 2 equality constraints. The objective 
C     function is quadratic and the constraints are linear. 

C     ******************************************************************
C     ******************************************************************

      subroutine pedazos4_inip(n,x,l,u,m,lambda,rho,equatn,linear)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(*),linear(*)
      double precision l(*),lambda(*),rho(*),u(*),x(*)

C     This subroutine sets some problem data. 
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

C     COMMON ARRAYS
      double precision t(30),y(30)

C     LOCAL SCALARS
      integer i
      double precision drand,seed

C     COMMON BLOCKS
      common /probdata/ t,y
      
C     Number of variables

      n = 9

C     Adjust a line, a parabola and a cubic

      do i = 1,10
          t(i) = dfloat(i)
          y(i) = t(i)
      end do

      do i = 11,20
          t(i) = dfloat(i)
          y(i) = 0.4d0 * ( t(i) - 15.0d0 ) ** 2
      end do

      do i = 21,30
          t(i) = dfloat(i)
          y(i) = - 10.0d0 * ( t(i) - 25.0d0 ) ** 3 / 125.0d0
      end do

C     Perturbation

      seed = 17172937.0d0

      do i = 1,30
          y(i) = y(i) +
     +           y(i) * 0.20d0 * ( drand(seed) * 2.0d0 - 1.0d0 )
      end do


C     Initial point

      do i = 1,n
          x(i) = 0.0d0
      end do
      
C     Lower and upper bounds

      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = 2

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
          linear(i) = .true.
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine pedazos4_evalf(n,x,f,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subroutine computes the objective function.
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
      double precision t(30),y(30)

C     LOCAL SCALARS
      integer i
      double precision r

C     COMMON BLOCKS
      common /probdata/ t,y
      
      flag = 0

      f = 0.0d0

      do i = 1,30

          if ( t(i) .le. 10.0d0 ) then
              r = ( x(1) + x(2) * t(i) - y(i) ) ** 2

          else if ( t(i) .le. 20.0d0 ) then
              r = ( x(3) + x(4) * t(i) + x(5) * t(i) ** 2 - y(i) ) ** 2

          else
              r = ( x(6) + x(7) * t(i) + x(8) * t(i) ** 2 + 
     +              x(9) * t(i) ** 3 - y(i) ) ** 2
          end if

          f = f + r

      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine pedazos4_evalg(n,x,g,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     This subroutine computes the gradient vector of the objective 
C     function. 
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
      double precision t(30),y(30)

C     LOCAL SCALARS
      integer i
      double precision r

C     COMMON BLOCKS
      common /probdata/ t,y
      
      flag = 0

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,30

          if ( t(i) .le. 10.0d0 ) then
              r    = x(1) + x(2) * t(i) - y(i)
              g(1) = g(1) + 2.0d0 * r
              g(2) = g(2) + 2.0d0 * r * t(i)

          else if ( t(i) .le. 20.0d0 ) then
              r    = x(3) + x(4) * t(i) + x(5) * t(i) ** 2 - y(i)
              g(3) = g(3) + 2.0d0 * r
              g(4) = g(4) + 2.0d0 * r * t(i)
              g(5) = g(5) + 2.0d0 * r * t(i) ** 2

          else
              r    = x(6) + x(7) * t(i) + x(8) * t(i) ** 2 + 
     +               x(9) * t(i) ** 3 - y(i)
              g(6) = g(6) + 2.0d0 * r
              g(7) = g(7) + 2.0d0 * r * t(i)
              g(8) = g(8) + 2.0d0 * r * t(i) ** 2
              g(9) = g(9) + 2.0d0 * r * t(i) ** 3
          end if

      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine pedazos4_evalh(n,x,hlin,hcol,hval,nnzh,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n,nnzh

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     This subroutine computes the Hessian matrix of the objective 
C     function. 
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

      flag = 0

      nnzh = 19 

      hlin(1)  = 1
      hcol(1)  = 1
      hval(1)  = 20.0d0  

      hlin(2)  = 2
      hcol(2)  = 1
      hval(2)  = 110.0d0 

      hlin(3)  = 2
      hcol(3)  = 2
      hval(3)  = 770.0d0  

      hlin(4)  = 3
      hcol(4)  = 3
      hval(4)  = 20.0d0  

      hlin(5)  = 4
      hcol(5)  = 3
      hval(5)  = 310.0d0

      hlin(6)  = 4    
      hcol(6)  = 4
      hval(6)  = 4970.0d0  

      hlin(7)  = 5
      hcol(7)  = 3
      hval(7)  = 4970.0d0

      hlin(8)  = 5     
      hcol(8)  = 4
      hval(8)  = 82150.0d0  

      hlin(9)  = 5
      hcol(9)  = 5
      hval(9)  = 1394666.0d0 

      hlin(10) = 6
      hcol(10) = 6 
      hval(10) = 20.0d0  

      hlin(11) = 7    
      hcol(11) = 6
      hval(11) = 510.0d0  

      hlin(12) = 7 
      hcol(12) = 7 
      hval(12) = 13170.0d0

      hlin(13) = 8
      hcol(13) = 6
      hval(13) = 13170.0d0  

      hlin(14) = 8
      hcol(14) = 7
      hval(14) = 344250.0d0  

      hlin(15) = 8
      hcol(15) = 8
      hval(15) = 9102666.0d0 

      hlin(16) = 9
      hcol(16) = 6
      hval(16) = 344250.0d0  

      hlin(17) = 9
      hcol(17) = 7
      hval(17) = 9102666.0d0  

      hlin(18) = 9
      hcol(18) = 8
      hval(18) = 243308250.0d0  

      hlin(19) = 9
      hcol(19) = 9
      hval(19) = 6.56895081d+09

      end

C     ******************************************************************
C     ******************************************************************

      subroutine pedazos4_evalc(n,x,ind,c,flag)

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
 
      flag = 0

      if ( ind .eq. 1 ) then
          c = ( x(1) + x(2) * 10.d0 ) - 
     +        ( x(3) + x(4) * 10.d0 + x(5) * 100.d0 )

      else if( ind .eq. 2 ) then
          c = ( x(3) + x(4) * 20.d0 + x(5) * 400.d0 ) - 
     +        ( x(6) + x(7) * 20.d0 + x(8) * 400.d0 + x(9) * 8000.d0 )
      end if 

      end
 
C     ******************************************************************
C     ******************************************************************
 
      subroutine pedazos4_evaljac(n,x,ind,indjac,valjac,nnzjac,flag)

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

      flag = 0

      if ( ind .eq. 1 ) then

          nnzjac = 5

          indjac(1) =        1
          valjac(1) =     1.d0

          indjac(2) =        2
          valjac(2) =    10.d0

          indjac(3) =        3
          valjac(3) = -   1.d0

          indjac(4) =        4
          valjac(4) = -  10.d0

          indjac(5) =        5
          valjac(5) = - 100.d0

      else if( ind .eq. 2 ) then

          nnzjac = 7

          indjac(1) =         3
          valjac(1) =      1.d0

          indjac(2) =         4
          valjac(2) =     20.d0

          indjac(3) =         5
          valjac(3) =    400.d0

          indjac(4) =         6
          valjac(4) = -    1.d0

          indjac(5) =         7
          valjac(5) = -   20.d0

          indjac(6) =         8
          valjac(6) = -  400.d0

          indjac(7) =         9
          valjac(7) = - 8000.d0

      end if 

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine pedazos4_evalhc(n,x,ind,hclin,hccol,hcval,nnzhc,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,n,nnzhc

C     ARRAY ARGUMENTS
      integer hccol(*),hclin(*)
      double precision hcval(*),x(n)

C     This subroutine computes the Hessian matrix of the ind-th
C     constraint. 
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

      flag = -1

      nnzhc = 0

      end

C     ******************************************************************
C     ******************************************************************

      subroutine pedazos4_evalhlp(n,x,m,lambda,p,hp,goth,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),x(n)

C     This subroutine computes the product of the Hessian of the
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

      subroutine pedazos4_endp(n,x,l,u,m,lambda,rho,equatn,linear)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),rho(m),u(n),x(n)

C     This subroutine can be used to do some extra job after the solver
C     has found the solution, like some extra statistics, or to save the
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
      double precision t(30),y(30)

C     LOCAL SCALARS
      integer i
      double precision z,tt

C     COMMON BLOCKS
      common /probdata/ t,y
      
      write(*,*) 't, observed y, computed y'

      do i = 1,30

          if ( t(i) .le. 10.0d0 ) then
              z = x(1) + x(2) * t(i)

          else if ( t(i) .le. 20.0d0 ) then
              z = x(3) + x(4) * t(i) + x(5) * t(i) ** 2

          else
              z = x(6) + x(7) * t(i) + x(8) * t(i) ** 2 + 
     +            x(9) * t(i) ** 3
          end if

          write(*,*) t(i),y(i),z

      end do

      tt = 1.0d0
 10   if ( tt .le. 30.01d0 ) then
          if ( tt .le. 10.0d0 ) then
              z = x(1) + x(2) * tt

          else if ( tt .le. 20.0d0 ) then
              z = x(3) + x(4) * tt + x(5) * tt ** 2

          else
              z = x(6) + x(7) * tt + x(8) * tt ** 2 + 
     +            x(9) * tt ** 3
          end if

          write(*,*) tt,z

          tt = tt + 0.01d0
          go to 10

      end if

      end
