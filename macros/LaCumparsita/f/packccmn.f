C     =================================================================
C     File: packccmn.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     April 7, 2006.

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

C     Packing fixed-dimension circular items within a fixed-dimension
C     --------------------------------------------------------------- 
C     circular object maximizing the number of items
C     ----------------------------------------------
C
C     We wish to place k circles of radii ri (i=1, ..., k) into a 
C     circular object with radius R in such a way that the circles are 
C     not overlapped. Therefore, given k, the radii ri (i=1, ..., k) 
C     and R, the goal is to solve the problem:
C
C     Minimize \sum_{i not equal j} max(0, (ri+rj)^2 - d(pi,pj)^2 )^2
C
C     subject to d(0,pi)^2 <= (R - ri)^2, i = 1, ..., k,
C
C     where d(.,.) is the Euclidian distance. If the objective function 
C     value at the global minimizer of this problem is zero then the 
C     answer of the decision problem: "Given k circular items of radii 
C     ri (i=1, ..., k) and a circular object of radius R, whether is it 
C     possible to locate all the items within the object or not." is YES, 
C     otherwise, the answer is NO.

C     ******************************************************************
C     ******************************************************************

      subroutine packccmn_inip(n,x,l,u,m,lambda,rho,equatn,linear)

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
      integer nite,ndim
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i
      double precision drand,seed,r,t

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

C     SET PROBLEM DATA

C     Dimension of the space
      ndim   =      2

C     Number of identical circular items to be packed
      nite   =     80

C     Radius of the circular items to be packed
      iterad =  1.0d0

C     Radius of the circular object within which the items will 
C     be packed
      objrad = 10.0d0

C     Number of variables

      n = 2 * nite

C     Initial point

      seed = 12337.0d0

      do i = 1,nite
          r = ( objrad - iterad ) * drand(seed)
          t = 2.0d0 * 3.14159d0 * drand(seed)
          x(2*i-1) = r * cos( t )
          x(2*i  ) = r * sin( t )
      end do

C     Lower and upper bounds
      do i = 1,n
          l(i) = - 1.0d+20
          u(i) =   1.0d+20
      end do

C     Number of constraints (equalities plus inequalities)

      m = nite

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

      subroutine packccmn_evalf(n,x,f,flag)

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
      integer ndim,nite
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i,j,k,ind1,ind2
      double precision fparc,dist

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0
 
C     COMPUTE DENSE OVERLAPPING

      f = 0.0d0
      do i = 1,nite
          do j = i + 1,nite
              dist = 0.0d0
              do k = 1,ndim
                  ind1 = ( i - 1 ) * ndim + k
                  ind2 = ( j - 1 ) * ndim + k
                  dist = dist + ( x(ind1) - x(ind2) ) ** 2
              end do
              fparc = max( 0.0d0, ( 2.0d0 * iterad ) ** 2 - dist )
              f = f + fparc ** 2
          end do
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine packccmn_evalg(n,x,g,flag)

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
      integer ndim,nite
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i,ind1,ind2,j,k
      double precision dist,fparc,tmp

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0

C     COMPUTE DENSE OVERLAPPING

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,nite
          do j = i + 1,nite

              dist = 0.0d0
              do k = 1,ndim
                  ind1 = ( i - 1 ) * ndim + k
                  ind2 = ( j - 1 ) * ndim + k
                  dist = dist + ( x(ind1) - x(ind2) ) ** 2
              end do
              fparc = max( 0.0d0, ( 2.0d0 * iterad ) ** 2 - dist )

              if ( fparc .ne. 0.d0 ) then
                  do k = 1,ndim
                      ind1 = ( i - 1 ) * ndim + k
                      ind2 = ( j - 1 ) * ndim + k
                      tmp = 4.0d0 * fparc * ( x(ind1) - x(ind2) )
                      g(ind1) = g(ind1) - tmp 
                      g(ind2) = g(ind2) + tmp
                  end do
              end if

          end do
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine packccmn_evalh(n,x,hlin,hcol,hval,nnzh,flag)

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

C     PARAMETERS
      integer nitemax,ndimmax
      parameter ( nitemax =        10000 )
      parameter ( ndimmax =            2 )

C     COMMON SCALARS
      integer ndim,nite
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i,j,k,l,ind1,ind2,ind3,ind4
      double precision dist,tmp

C     LOCAL ARRAYS
      double precision diagb(ndimmax,ndimmax,nitemax)

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0

C     INITALIZE DIAGONAL BLOCKS

      do k = 1,nite
          do j = 1,ndim
              do i = j,ndim
                  diagb(i,j,k) = 0.0d0
              end do
          end do
      end do

C     COMPUTE DENSE OVERLAPPING SECOND DERIVATIVES

      nnzh = 0

      do i = 1,nite
          do j = i + 1,nite

              dist = 0.0d0
              do k = 1,ndim
                  ind1 = ( i - 1 ) * ndim + k
                  ind2 = ( j - 1 ) * ndim + k
                  dist = dist + ( x(ind1) - x(ind2) ) ** 2
              end do

              if ( dist .le. ( 2.0d0 * iterad ) ** 2 ) then
                  do k = 1,ndim
                      ind1 = ( i - 1 ) * ndim + k
                      ind2 = ( j - 1 ) * ndim + k
                      tmp = 8.0d0 * ( x(ind1) - x(ind2) ) ** 2
     +                    - 4.0d0 * ( ( 2.0d0 * iterad ) ** 2 - dist )
                      if ( tmp .ne. 0.0d0 ) then
C                         H(ind1,ind1) = H(ind1,ind1) + tmp
                          diagb(k,k,i) = diagb(k,k,i) + tmp
C                         H(ind2,ind2) = H(ind2,ind2) + tmp
                          diagb(k,k,j) = diagb(k,k,j) + tmp
C                         H(ind2,ind1) = H(ind2,ind1) - tmp
                          nnzh = nnzh + 1
                          hlin(nnzh) = ind2
                          hcol(nnzh) = ind1
                          hval(nnzh) = - tmp
                      end if
                      do l = 1,k - 1
                          ind3 = ( i - 1 ) * ndim + l
                          ind4 = ( j - 1 ) * ndim + l
                          tmp = 8.0d0 * ( x(ind3) - x(ind4) )
     +                                * ( x(ind1) - x(ind2) )
                          if ( tmp .ne. 0.0d0 ) then
C                             H(ind1,ind3) = H(ind1,ind3) + tmp
                              diagb(k,l,i) = diagb(k,l,i) + tmp
C                             H(ind2,ind4) = H(ind2,ind4) + tmp
                              diagb(k,l,j) = diagb(k,l,j) + tmp
C                             H(ind2,ind3) = H(ind2,ind3) - tmp
                              nnzh = nnzh + 1
                              hlin(nnzh) = ind2
                              hcol(nnzh) = ind3
                              hval(nnzh) = - tmp
                          end if
                      end do
                      do l = k + 1,ndim
                          ind3 = ( i - 1 ) * ndim + l
                          ind4 = ( j - 1 ) * ndim + l
                          tmp = 8.0d0 * ( x(ind3) - x(ind4) )
     +                                * ( x(ind1) - x(ind2) )
                          if ( tmp .ne. 0.0d0 ) then
C                             H(ind2,ind3) = H(ind2,ind3) - tmp
                              nnzh = nnzh + 1
                              hlin(nnzh) = ind2
                              hcol(nnzh) = ind3
                              hval(nnzh) = - tmp
                          end if
                      end do
                  end do
              end if

          end do
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
 
      subroutine packccmn_evalc(n,x,ind,c,flag)

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
      integer nite,ndim
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0

      if ( 0 .lt. ind .and. ind .le. nite ) then

          c = - ( objrad -iterad ) ** 2

          do i = 1,ndim
              c = c + x(ndim * ( ind - 1 ) + i) ** 2
          end do

          return

      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine packccmn_evaljac(n,x,ind,indjac,valjac,nnzjac,flag)

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
      integer nite,ndim
      double precision iterad,objrad

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0

      if ( 0 .lt. ind .and. ind .le. nite ) then

          nnzjac = ndim

          do i = 1,ndim
              indjac(i) = ndim * ( ind - 1 ) + i
              valjac(i) = 2.0d0 * x(ndim * ( ind - 1 ) + i)
          end do

          return

      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine packccmn_evalhc(n,x,ind,hclin,hccol,hcval,nnzhc,flag)

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

C     COMMON SCALARS
      integer nite,ndim
      double precision iterad,objrad

C     LOCAL SCALARS
      integer k

C     COMMON BLOCKS
      common /packdata/ nite,ndim,iterad,objrad

      flag = 0

      if ( 0 .lt. ind .and. ind .le. nite ) then

          nnzhc = ndim

          do k = 1,ndim
              hclin(k) = ndim * ( ind - 1 ) + k
              hccol(k) = ndim * ( ind - 1 ) + k
              hcval(k) = 2.0d0
          end do

          return

      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine packccmn_evalhlp(n,x,m,lambda,p,hp,goth,flag)

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

      subroutine packccmn_endp(n,x,l,u,m,lambda,rho,equatn,linear)

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

