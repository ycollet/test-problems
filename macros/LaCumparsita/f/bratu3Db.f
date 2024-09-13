C     =================================================================
C     File: bratu3Db.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     November 30, 2006.

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

C     Discretized 3D Bratu-based problem
C     ----------------------------------

C     Description: 

C     The discretized three-dimensional Bratu-based optimization 
C     problem is:
C
C     Minimize \sum_{(i,j,k) \in S} ( u(i,j,k) - u_*(i,j,k) )^2
C
C     subject to 
C
C     \phi_\theta(u,i,j,k) = \phi_\theta(u_*,i,j,k), i,j,k=2,...,np-1,
C
C     where u_* was choosed as
C
C     u_*(i,j,k) = 10 q(i) q(j) q(k) 
C                     (1-q(i)) (1-q(j)) (1-q(k)) e^{q(k)^{4.5}}
C
C     with 
C
C     q(\ell) = \frac{ np - \ell }{ np - 1 } for i,j,k=1,...,np 
C
C     and 
C
C     \phi_\theta(v,i,j,k) = - \Delta v(i,j,k) + \theta e^{v(i,j,k)}, 
C
C     \Delta v(i,j,k) =  \frac{v(i \pm 1,j,k) + 
C                              v(i,j \pm 1,k) + 
C                              v(i,j,k \pm 1) - 6 v(i,j,k)} {h^2},
C
C     for i,j,k=2,...,np-1. 
C
C     The number of variables is n=np^3 and the number of (equality) 
C     constraints is m=(np-2)^3. 
C
C     We set \theta=-100, h=1/(np-1), |S|=7 and the 3-uples of 
C     indices in S are randomly selected in [1,np]^3.
C
C     The initial point is randomly generated in [0,1]^n.
C
C     Several problems can be considered varying the value of np, for
C     example, np = 5, 6, ..., 20.

C     A much more friendly description of the problem can be found in:
C
C     R. Andreani, E.G. Birgin, J.M. Martinez and M.L. Schuverdt, On 
C     Augmented Lagrangian Methods with general lower-level constraints, 
C     Technical Report MCDO-051015 (see www.ime.usp.br/~egbirgin/), 
C     Department of Applied Mathematics, UNICAMP, Brazil, 2005.

C     ******************************************************************
C     ******************************************************************

      subroutine bratu3Db_inip(n,x,l,u,m,lambda,rho,equatn,linear,
     +                         nppar,seed)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n,nppar

C     ARRAY ARGUMENTS
      logical equatn(*),linear(*)
      double precision l(*),lambda(*),rho(*),u(*),x(*),seed

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
      integer npmax,tmax
      parameter ( npmax =  100 )
      parameter ( tmax  = 1000 )

C     COMMON SCALARS
      integer np,t
      double precision h,theta

C     COMON ARRAYS
      integer vind(tmax,3),cind(3,(npmax-2)**3)
      double precision phi(2:npmax-1,2:npmax-1,2:npmax-1),vval(tmax)

C     LOCAL SCALARS
      integer i,ind,j,k
      double precision a,b,c

C     LOCAL ARRAYS
      double precision uopt(npmax,npmax,npmax)

C     COMMON BLOCKS
      common /probdata/ phi,vind,vval,cind,np,t,h,theta

C     EXTERNAL FUNCTIONS
      integer bratu3Db_xind
      double precision drand

C     READ PROBLEM VARIABLE DATA

      theta = - 100.0d0
      t     =   7

C     write(*,*) 'Number of points (np): '
C     read(*,*) np

C     write(*,*) 'Seed for the initial point random generation: '
C     read(*,*) seed

C     read(*,*) np,theta,t,seed
C     seed = 123456.0d0 * seed

C     DEFINE SOME PROBLEM PARAMETERS

C     Discretization step
      np = nppar
      h = 1.0d0 / ( np - 1 )

C     "Known solution" used to generate the problem data

      do k = 1,np
          do j = 1,np
              do i = 1,np
                  a = float( np - k ) / ( np - 1 )
                  b = float( np - j ) / ( np - 1 )
                  c = float( np - i ) / ( np - 1 )
                  uopt(i,j,k) = 10.0d0 * a * b * c
     +                        * (1.0d0 - a) * (1.0d0 - b) * (1.0d0 - c)
     +                        * exp( a ** 4.5d0 )
              end do
          end do
      end do

C     Right hand side of the system

      do k = 2,np - 1
          do j = 2,np - 1
              do i = 2,np - 1
                  phi(i,j,k) = ( uopt(i,j,k+1) + uopt(i,j,k-1) +
     +                           uopt(i,j+1,k) + uopt(i,j-1,k) +
     +                           uopt(i+1,j,k) + uopt(i-1,j,k) -
     +                           6.0d0 * uopt(i,j,k) ) / h ** 2
     +                         - theta * exp( uopt(i,j,k) )
              end do
          end do
      end do

C     Points of the grid at which the solution is known

      do i = 1,t
          do j = 1,3
              vind(i,j) = 1 + int( drand(seed) * np )
          end do
          vval(i) = uopt(vind(i,1),vind(i,2),vind(i,3))
      end do

C     Indices i, j and k related to ind-th constraint

      ind = 0
      do i = 2,np - 1
          do j = 2,np - 1
              do k = 2,np - 1
                  ind = ind + 1
                  cind(1,ind) = i
                  cind(2,ind) = j
                  cind(3,ind) = k
              end do
          end do
      end do

C     Number of variables (elements of a lower-triangular ndim x ndim matrix)

      n = np ** 3

C     Lower and upper bounds

      do k = 1,np
          do j = 1,np
              do i = 1,np
                  ind = bratu3Db_xind(i,j,k)
                  l(ind) = - 1.0d+20
                  u(ind) =   1.0d+20
              end do
          end do
      end do

C     Initial guess

      do k = 1,np
          do j = 1,np
              do i = 1,np
                  ind = bratu3Db_xind(i,j,k)
                  x(ind) = drand(seed)
              end do
          end do
      end do

C     Number of constraints (equalities plus inequalities)

      m = ( np - 2 ) ** 3

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

      end

C     ******************************************************************
C     ******************************************************************

      subroutine bratu3Db_evalf(n,x,f,flag)

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
      integer npmax,tmax
      parameter ( npmax =  100 )
      parameter ( tmax  = 1000 )

C     COMMON SCALARS
      integer np,t
      double precision h,theta

C     COMON ARRAYS
      integer vind(tmax,3),cind(3,(npmax-2)**3)
      double precision phi(2:npmax-1,2:npmax-1,2:npmax-1),vval(tmax)

C     LOCAL SCALARS
      integer i,ind

C     COMMON BLOCKS
      common /probdata/ phi,vind,vval,cind,np,t,h,theta

C     EXTERNAL FUNCTIONS
      integer bratu3Db_xind

      flag = 0

      f = 0.0d0
      do i = 1,t
          ind = bratu3Db_xind(vind(i,1),vind(i,2),vind(i,3))
          f = f + ( x(ind) - vval(i) ) ** 2
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine bratu3Db_evalg(n,x,g,flag)

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
      integer npmax,tmax
      parameter ( npmax =  100 )
      parameter ( tmax  = 1000 )

C     COMMON SCALARS
      integer np,t
      double precision h,theta

C     COMON ARRAYS
      integer vind(tmax,3),cind(3,(npmax-2)**3)
      double precision phi(2:npmax-1,2:npmax-1,2:npmax-1),vval(tmax)

C     LOCAL SCALARS
      integer i,ind

C     COMMON BLOCKS
      common /probdata/ phi,vind,vval,cind,np,t,h,theta

C     EXTERNAL FUNCTIONS
      integer bratu3Db_xind

      flag = 0

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,t
          ind = bratu3Db_xind(vind(i,1),vind(i,2),vind(i,3))
          g(ind) = g(ind) + 2.0d0 * ( x(ind) - vval(i) )
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine bratu3Db_evalh(n,x,hlin,hcol,hval,hnnz,flag)

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

C     PARAMETERS
      integer npmax,tmax
      parameter ( npmax =  100 )
      parameter ( tmax  = 1000 )

C     COMMON SCALARS
      integer np,t
      double precision h,theta

C     COMON ARRAYS
      integer vind(tmax,3),cind(3,(npmax-2)**3)
      double precision phi(2:npmax-1,2:npmax-1,2:npmax-1),vval(tmax)

C     LOCAL SCALARS
      integer i,ind

C     COMMON BLOCKS
      common /probdata/ phi,vind,vval,cind,np,t,h,theta

C     EXTERNAL FUNCTIONS
      integer bratu3Db_xind

      flag = 0

      hnnz = t

      do i = 1,t
          ind = bratu3Db_xind(vind(i,1),vind(i,2),vind(i,3))
          hlin(i) = ind
          hcol(i) = ind
          hval(i) = 2.0d0
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine bratu3Db_evalc(n,x,ind,c,flag)

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

C     PARAMETERS
      integer npmax,tmax
      parameter ( npmax =  100 )
      parameter ( tmax  = 1000 )

C     COMMON SCALARS
      integer np,t
      double precision h,theta

C     COMON ARRAYS
      integer vind(tmax,3),cind(3,(npmax-2)**3)
      double precision phi(2:npmax-1,2:npmax-1,2:npmax-1),vval(tmax)

C     LOCAL SCALARS
      integer i,i0,i1,i2,i3,i4,i5,i6,j,k

C     COMMON BLOCKS
      common /probdata/ phi,vind,vval,cind,np,t,h,theta

C     EXTERNAL FUNCTIONS
      integer bratu3Db_xind

      flag = 0

      i = cind(1,ind)
      j = cind(2,ind)
      k = cind(3,ind)

      i0 = bratu3Db_xind(i,j,k)
      i1 = bratu3Db_xind(i,j,k+1)
      i2 = bratu3Db_xind(i,j,k-1)
      i3 = bratu3Db_xind(i,j+1,k)
      i4 = bratu3Db_xind(i,j-1,k)
      i5 = bratu3Db_xind(i+1,j,k)
      i6 = bratu3Db_xind(i-1,j,k)

      c = ( x(i1) + x(i2) + x(i3) + x(i4) + x(i5) + x(i6) - 
     +    6.0d0 * x(i0) ) / h ** 2 - theta * exp( x(i0) ) - phi(i,j,k)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine bratu3Db_evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

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
C     jcnnz    integer,
C              number of perhaps-non-null elements of the computed 
C              gradient,
C
C     jcvar    integer jcvar(jcnnz),
C              see below,
C
C     jcval    double precision jcval(jcnnz),
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
      integer npmax,tmax
      parameter ( npmax =  100 )
      parameter ( tmax  = 1000 )

C     COMMON SCALARS
      integer np,t
      double precision h,theta

C     COMON ARRAYS
      integer vind(tmax,3),cind(3,(npmax-2)**3)
      double precision phi(2:npmax-1,2:npmax-1,2:npmax-1),vval(tmax)

C     LOCAL SCALARS
      integer i,i0,i1,i2,i3,i4,i5,i6,j,k

C     COMMON BLOCKS
      common /probdata/ phi,vind,vval,cind,np,t,h,theta

C     EXTERNAL FUNCTIONS
      integer bratu3Db_xind

      flag = 0

      i = cind(1,ind)
      j = cind(2,ind)
      k = cind(3,ind)

      i0 = bratu3Db_xind(i,j,k)
      i1 = bratu3Db_xind(i,j,k+1)
      i2 = bratu3Db_xind(i,j,k-1)
      i3 = bratu3Db_xind(i,j+1,k)
      i4 = bratu3Db_xind(i,j-1,k)
      i5 = bratu3Db_xind(i+1,j,k)
      i6 = bratu3Db_xind(i-1,j,k)

      jcnnz = 7

      jcvar(1) = i1
      jcval(1) = 1.0d0 / h ** 2

      jcvar(2) = i2
      jcval(2) = 1.0d0 / h ** 2

      jcvar(3) = i3
      jcval(3) = 1.0d0 / h ** 2

      jcvar(4) = i4
      jcval(4) = 1.0d0 / h ** 2

      jcvar(5) = i5
      jcval(5) = 1.0d0 / h ** 2

      jcvar(6) = i6
      jcval(6) = 1.0d0 / h ** 2

      jcvar(7) = i0
      jcval(7) = - 6.0d0 / h ** 2 - theta * exp( x(i0) )

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine bratu3Db_evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag)

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
      integer npmax,tmax
      parameter ( npmax =  100 )
      parameter ( tmax  = 1000 )

C     COMMON SCALARS
      integer np,t
      double precision h,theta

C     COMON ARRAYS
      integer vind(tmax,3),cind(3,(npmax-2)**3)
      double precision phi(2:npmax-1,2:npmax-1,2:npmax-1),vval(tmax)

C     LOCAL SCALARS
      integer i,i0,j,k

C     COMMON BLOCKS
      common /probdata/ phi,vind,vval,cind,np,t,h,theta

C     EXTERNAL FUNCTIONS
      integer bratu3Db_xind

      flag = 0

      i = cind(1,ind)
      j = cind(2,ind)
      k = cind(3,ind)

      i0 = bratu3Db_xind(i,j,k)

      hcnnz = 1

      hclin(1) = i0
      hccol(1) = i0
      hcval(1) = - theta * exp( x(i0) )

      end

C     ******************************************************************
C     ******************************************************************

      subroutine bratu3Db_evalhlp(n,x,m,lambda,p,hp,goth,flag)

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

      subroutine bratu3Db_endp(n,x,l,u,m,lambda,rho,equatn,linear)

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

C     ******************************************************************
C     ******************************************************************

      integer function bratu3Db_xind(i,j,k)

C     SCALAR ARGUMENTS
      integer i,j,k

C     PARAMETERS
      integer npmax,tmax
      parameter ( npmax =  100 )
      parameter ( tmax  = 1000 )

C     COMMON SCALARS
      integer np,t
      double precision h,theta

C     COMON ARRAYS
      integer vind(tmax,3),cind(3,(npmax-2)**3)
      double precision phi(2:npmax-1,2:npmax-1,2:npmax-1),vval(tmax)

C     COMMON BLOCKS
      common /probdata/ phi,vind,vval,cind,np,t,h,theta

      bratu3Db_xind = np ** 2 * ( i - 1 ) + np * ( j - 1 ) + k

      return

      end
