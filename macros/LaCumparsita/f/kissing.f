C     =================================================================
C     File: kissing.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     October 30th, 2006.

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

C     Kissing Number problem
C     ----------------------

C     Given npun and ndim, try to find npun points p_1, p_2, ..., 
C     p_npun in the unitary sphere in R^ndim such that 
C
C          dist(pi,pj) >= 1, for all i < j.

C     ******************************************************************
C     ******************************************************************

      subroutine kissing_inip(n,x,l,u,m,lambda,rho,equatn,linear,
     &                        nd_in,np_in,seed_in)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n,nd_in,np_in

C     ARRAY ARGUMENTS
      logical equatn(*),linear(*)
      double precision l(*),lambda(*),rho(*),u(*),x(*)
      double precision seed_in

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
      integer npumax
      parameter ( npumax = 1000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer pair(2,npumax * ( npumax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,j,k
      double precision seed

C     COMMON BLOCKS
      common /probdata/ pair,ndim,npun

C     EXTERNAL FUNCTIONS
      double precision drand

C     Set problem data

C     Dimension of the space (ndim) and number of points (npun):
      ndim = nd_in
      npun = np_in

C     Seed for the initial point random generation
      seed = seed_in

C     Set pairs
      k = 0
      do i = 1,npun
          do j = i + 1,npun
              k = k + 1
              pair(1,k) = i
              pair(2,k) = j
          end do
      end do

C     Number of variables
      n = ndim * npun

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 2.0d0
          u(i) =   2.0d0
      end do

C     Initial point
      do i = 1,n
          x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
      end do

C     Number of constraints (equalities plus inequalities)

      m = npun + npun * ( npun - 1 ) / 2

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
          rho(i) = 1.0d0
      end do

C     For each constraint i, set equatn(i) = .true. if it is an equality
C     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if 
C     it is an inequality constraint of the form c_i(x) <= 0.

      do i = 1,npun
          equatn(i) = .true.
      end do

      do i = 1,npun + 1, npun + npun * ( npun - 1 ) / 2
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

      subroutine kissing_evalf(n,x,f,flag)

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

      flag = 0

      f = 0.0d0

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine kissing_evalg(n,x,g,flag)

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

C     LOCAL SCALARS
      integer i

      flag = 0

      do i = 1,n
          g(i) = 0.0d0
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine kissing_evalh(n,x,hlin,hcol,hval,hnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n,hnnz

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
C     hnnz     integer,
C              number of perhaps-non-null elements of the computed 
C              Hessian,
C
C     hlin     integer hlin(hnnz),
C              see below,
C
C     hcol     integer hcol(hnnz),
C              see below,
C
C     hval     double precision hval(hnnz),
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

      hnnz = 0

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine kissing_evalc(n,x,ind,c,flag)

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
      integer npumax
      parameter ( npumax = 1000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer pair(2,npumax * ( npumax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,i1,i2,ind1,ind2

C     COMMON BLOCKS
      common /probdata/ pair,ndim,npun

      flag = 0

      if ( 1 .le. ind .and. ind .le. npun ) then

          c = - 1.0d0

          do i = 1,ndim
              c = c + x(ndim * ( ind - 1 ) + i) ** 2
          end do

      else if ( npun + 1 .le. ind .and. 
     +          ind .le. npun + npun * (npun - 1) / 2 ) then

          i1 = pair(1,ind - npun)
          i2 = pair(2,ind - npun)

          c = 1.0d0
          do i = 1,ndim
              ind1 = ndim * (i1 - 1) + i
              ind2 = ndim * (i2 - 1) + i

              c = c - ( x(ind1) - x(ind2) ) ** 2
          end do

      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine kissing_evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

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
C     jcnnz   integer,
C              number of perhaps-non-null elements of the computed 
C              gradient,
C
C     jcvar   integer jcvar(jcnnz),
C              see below,
C
C     jcval   double precision jcval(jcnnz),
C              the non-null value of the partial derivative of the 
C              ind-th constraint with respect to the jcvar(k)-th 
C              variable must be saved at jcval(k).
C
C     flag     integer
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the 
C              constraint. (For example, trying to compute the square 
C              root of a negative number, dividing by zero or a very 
C              small number, etc.) If everything was o.k. you must set 
C              it equal to zero.

C     PARAMETERS
      integer npumax
      parameter ( npumax = 1000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer pair(2,npumax * ( npumax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,i1,i2,ind1,ind2

C     COMMON BLOCKS
      common /probdata/ pair,ndim,npun

      flag = 0

      if ( 1 .le. ind .and. ind .le. npun ) then

          jcnnz = ndim

          do i = 1,ndim
              jcvar(i) = ndim * ( ind - 1 ) + i
              jcval(i) = 2.0d0 * x(ndim * ( ind - 1 ) + i)
          end do

      else if ( npun + 1 .le. ind .and. 
     +          ind .le. npun + npun * (npun - 1) / 2 ) then

          i1 = pair(1,ind - npun)
          i2 = pair(2,ind - npun)

          jcnnz = 2 * ndim
            
          do i = 1,ndim
              ind1 = ndim * (i1 - 1) + i
              ind2 = ndim * (i2 - 1) + i

              jcvar(i) = ind1
              jcval(i) = - 2.0d0 * ( x(ind1) - x(ind2) )

              jcvar(ndim+i) = ind2
              jcval(ndim+i) = 2.0d0 * ( x(ind1) - x(ind2) )
          end do

      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine kissing_evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,n,hcnnz

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
C     hcnnz    integer,
C              number of perhaps-non-null elements of the computed 
C              Hessian,
C
C     hclin    integer hclin(hcnnz),
C              see below,
C
C     hccol    integer hccol(hcnnz),
C              see below,
C
C     hcval    double precision hcval(hcnnz),
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
      integer npumax
      parameter ( npumax = 1000 )

C     COMMON SCALARS
      integer ndim,npun

C     COMMON ARRAYS
      integer pair(2,npumax * ( npumax - 1 ) / 2)

C     LOCAL SCALARS
      integer i,i1,i2,ind1,ind2

C     COMMON BLOCKS
      common /probdata/ pair,ndim,npun

      flag = 0

      if ( 1 .le. ind .and. ind .le. npun ) then

          hcnnz = ndim

          do i = 1,ndim
              hclin(i) = ndim * ( ind - 1 ) + i
              hccol(i) = ndim * ( ind - 1 ) + i
              hcval(i) = 2.0d0
          end do

      else if ( npun + 1 .le. ind .and. 
     +          ind .le. npun + npun * (npun - 1) / 2 ) then

          i1 = pair(1,ind - npun)
          i2 = pair(2,ind - npun)

          hcnnz = 3 * ndim
            
          do i = 1,ndim
              ind1 = ndim * ( i1 - 1 ) + i
              ind2 = ndim * ( i2 - 1 ) + i

              hclin(i) = ind1
              hccol(i) = ind1 
              hcval(i) = - 2.0d0

              hclin(ndim+i) = ind2
              hccol(ndim+i) = ind2 
              hcval(ndim+i) = - 2.0d0

              hclin(2*ndim+i) = ind2
              hccol(2*ndim+i) = ind1 
              hcval(2*ndim+i) = 2.0d0
          end do

      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine kissing_evalhlp(n,x,m,lambda,p,hp,goth,flag)

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

      subroutine kissing_endp(n,x,l,u,m,lambda,rho,equatn,linear)

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
