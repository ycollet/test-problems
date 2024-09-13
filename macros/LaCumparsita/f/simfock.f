C     =================================================================
C     File: simfock.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     May 12, 2005.

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

C     Simulation of Hartree-Fock problem
C     ----------------------------------

C     Given kdim, ndim in 1, 2, 3, ..., ndim <= kdim, n = kdim * kdim, 
C     x in R^n, we define mat(x) in R^{kdim x kdim} as the matrix whose 
C     columnwise representation is x. The problem is:
C     
C                  Minimize 0.5 x^T A x + b^T x
C
C                  subject to
C
C                          mat(x) mat(x) = mat(x) 
C
C                          trace[mat(x)] = ndim 
C
C                          mat(x) = mat(x)^T
C
C     Therefore, the problem has kdim^2 variables, kdim^2 nonlinear 
C     equality constraints and kdim * ( kdim - 1 ) / 2 + 1 linear 
C     equality constraints.

C     ******************************************************************
C     ******************************************************************

      subroutine simfock_inip(n,x,l,u,m,lambda,rho,equatn,linear,
     &                        n_in,k_in)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

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
      integer kdim,ndim,k_in,n_in

C     LOCAL SCALARS
      integer i,j

C     COMMON BLOCKS
      common /probdata/ kdim,ndim

C     Set problem data

C     Dimension of the space (kdim) maximum 100
      kdim = k_in

C     Rank (ndim) of the desired idempotent solution: (ndim must be less than or equal to kdim)
      ndim = n_in

C     Number of variables

      n = kdim ** 2

C     Initial point
      do i = 1,kdim
          do j = 1,kdim
              x((j-1)*kdim+i) = 1.0d0
          end do
      end do

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 1.0d+05
          u(i) =   1.0d+05
      end do

C     Number of constraints (equalities plus inequalities)

      m = 2 * kdim ** 2 + 1

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

      end

C     ******************************************************************
C     ******************************************************************

      subroutine simfock_evalf(n,x,f,flag)

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
      integer kdim,ndim

C     LOCAL SCALARS
      integer i,j

C     LOCAL ARRAYS
      double precision g(10000000)

C     COMMON BLOCKS
      common /probdata/ kdim,ndim

      flag = 0

      f = 0.0d0
 
      do i = 2,kdim - 1
          f = f - x((i-2)*kdim+i) + 2.0d0 * x((i-1)*kdim+i) -
     +        x(i*kdim+i)
      end do

      f = f + 2.0d0 * x(1) - x(kdim+1) - x((kdim-2)*kdim+kdim) +
     +    2.0d0 * x((kdim-1)*kdim+kdim)

      f = 2.0d0 * f 
 
      call simfock_gedede(n,x,g) 

      do i = 1,kdim
          do j = 1,kdim
              f = f + g((i-1)*kdim+j) * x((i-1)*kdim+j)
          end do
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine simfock_evalg(n,x,g,flag)

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
      integer kdim,ndim

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /probdata/ kdim,ndim

      flag = 0

      call simfock_gedede (n,x,g) 

      do i = 2,kdim - 1
          g((i-2)*kdim+i) = g((i-2)*kdim+i) - 1.0d0
          g((i-1)*kdim+i) = g((i-1)*kdim+i) + 2.0d0
          g(i*kdim+i)     = g(i*kdim+i)     - 1.0d0
      end do

      g(1) = g(1) + 2.0d0

      g((kdim-2)*kdim+kdim) = g((kdim-2)*kdim+kdim) - 1.0d0
      g((kdim-1)*kdim+kdim) = g((kdim-1)*kdim+kdim) + 2.0d0
      g(kdim+1)             = g(kdim+1)             - 1.0d0
 
      do i = 1,n
          g(i) = 2.0d0 * g(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine simfock_evalh(n,x,hlin,hcol,hval,nnzh,flag)

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
 
      subroutine simfock_evalc(n,x,ind,c,flag)

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
      integer kdim,ndim

C     LOCAL SCALARS
      integer i,j,k,iminus

C     COMMON BLOCKS
      common /probdata/ kdim,ndim

      flag = 0

      if ( ind .le. kdim ** 2 ) then

          j = mod(ind,kdim)
          if ( j .eq. 0 ) then
              j = kdim
          end if

          i = ( ind - j ) / kdim + 1

          c = - x((j-1)*kdim+i)
          do k = 1,kdim
              c = c + x((k-1)*kdim+i) * x((j-1)*kdim+k)
          end do
     
          return

      end if

      if ( ind .le. 2 * kdim ** 2 ) then

          iminus = ind - kdim ** 2
 
          j = mod(iminus,kdim)
          if( j .eq. 0 ) then
              j = kdim
          end if

          i = ( iminus - j ) / kdim + 1

          c =  x((j-1)*kdim+i) - x((i-1)*kdim+j)
          
          return
      
      end if

      c = - dfloat( ndim )

      do i = 1,kdim
          c = c + x((i-1)*kdim+i)
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine simfock_evaljac(n,x,ind,indjac,valjac,nnzjac,flag)

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
      integer kdim,ndim

C     LOCAL SCALARS
      integer i,j,k,iminus

C     COMMON BLOCKS
      common /probdata/ kdim,ndim

      flag = 0

      if ( ind .le. kdim ** 2 ) then

          nnzjac = 2 * kdim + 1 

          j = mod(ind,kdim)
          if ( j .eq. 0 ) then
              j = kdim
          end if

          i = ( ind - j ) / kdim + 1

          do k = 1,kdim
              indjac(k) = ( j - 1 ) * kdim + k
              valjac(k) = x((k-1)*kdim+i)
              indjac(kdim+k) = ( k - 1 ) * kdim + i
              valjac(kdim+k) = x((j-1)*kdim+k)
          end do

          indjac(2*kdim+1) = ( j - 1 ) * kdim + i
          valjac(2*kdim+1) = - 1.0d0
 
          return

      end if

      if ( ind .le. 2 * kdim ** 2 ) then

          nnzjac = 2

          iminus = ind - kdim ** 2
 
          j = mod(iminus,kdim)
          if ( j .eq. 0 ) then
              j = kdim
          end if

          i = ( iminus - j ) / kdim + 1

          indjac(1) = ( j - 1 ) * kdim + i
          valjac(1) =   1.0d0
          indjac(2) = ( i - 1 ) * kdim + j
          valjac(2) = - 1.0d0     
    
          return

      end  if

      nnzjac = kdim

      do i = 1,kdim
          indjac(i) = ( i - 1 ) * kdim + i
          valjac(i) = 1.0d0
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine simfock_gedede(n,x,g)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n),g(n)

C     LOCAL SCALARS
      integer i,j,jj

      do i = 1,n

          g(i) = x(i)

          do jj = 1,10

              j = i + jj

              if ( j .le. n ) then
                  g(i) = g(i) + 1.0d0 / dfloat(i+j) * x(j)
              end if

              j = i - jj

              if ( j .ge. 1 ) then
                  g(i) = g(i) + 1.0d0 / dfloat(i+j) * x(j)
              end if

          end do

      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine simfock_evalhc(n,x,ind,hclin,hccol,hcval,nnzhc,flag)

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

      subroutine simfock_evalhlp(n,x,m,lambda,p,hp,goth,flag)

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

      subroutine simfock_endp(n,x,l,u,m,lambda,rho,equatn,linear)

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

      end
