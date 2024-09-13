      subroutine ssfifj(nx,ny,x,fvec,fjac,ldfjac,task,lambda)
      character*(*) task
      integer nx,ny,ldfjac
      real lambda
      real x(nx*ny),fvec(nx*ny),fjac(ldfjac,nx*ny)
c     ********
c
c     Subroutine ssfifj
c
c     This subroutine computes the function and Jacobian matrix of
c     the Solid Fuel Ignition problem. This problem arises in the
c     analysis of steady state solid fuel ignition.
c
c     The solid fuel ignition problem is modeled by the partial
c     differential equation
c
c          -Laplacian u - lambda*exp(u) = 0,  0 < x,y < 1 ,
c
c     with boundary conditions
c
c           u = 0  for  x = 0, x = 1, y = 0, y = 1.
c
c     A finite difference approximation is used to discretize the
c     boundary value problem to obtain a nonlinear system of equations.
c
c     The subroutine statement is:
c
c      subroutine ssfifj(nx,ny,x,fvec,fjac,ldfjac,task,lambda)
c
c     where
c
c       nx is an integer variable.
c         On entry nx is the number of grid points in the first
c            coordinate direction.
c         On exit nx is unchanged.
c
c       ny is an integer variable.
c         On entry ny is the number of grid points in the second
c            coordinate direction.
c         On exit ny is unchanged.
c
c       x is a real array of dimension n.
c         On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
c            x is set according to task.
c
c       fvec is a real array of dimension n.
c         On entry fvec need not be specified.
c         On exit fvec contains the function evaluated at x if
c            task = 'F' or 'FJ'.
c
c       fjac is a real array of dimension (ldfjac,n).
c         On entry fjac need not be specified.
c         On exit fjac contains the Jacobian matrix evaluated at x if
c            task = 'J' or 'FJ'.
c
c       ldfjac is an integer variable.
c          On entry ldfjac is the leading dimension of fjac.
c          On exit ldfjac is unchanged.
c
c       task is a character variable.
c         On entry task specifies the action of the subroutine:
c
c            task               action
c            ----               ------
c             'F'     Evaluate the function at x.
c             'J'     Evaluate the Jacobian matrix at x.
c             'FJ'    Evaluate the function and the Jacobian at x.
c             'XS'    Set x to the standard starting point xs.
c
c         On exit task is unchanged.
c
c       lambda is a real variable.
c         On entry lambda is the Frank-Kamenetski parameter.
c         On exit lambda is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     R.S. Maier and R.G. Carter
c
c     **********
      real zero,one,two,four
      parameter (zero=0.0,one=1.0,two=2.0,four=4.0)

      integer i,j,k,n
      real hx,hy,hxhy,hxdhy,hydhx,ub,ut,ul,ur,u,uxx,uyy,temp,temp1

      n = nx*ny
      hx = one/real(nx+1)
      hy = one/real(ny+1)
      hydhx = hy/hx
      hxdhy = hx/hy
      hxhy = hx*hy

c     Compute the standard starting point if task = 'XS'

      if (task.eq.'XS') then
         temp1 = lambda/ (lambda+one)
         do 20 j = 1,ny
            temp = real(min(j,ny-j+1))*hy
            do 10 i = 1,nx
               k = nx* (j-1) + i
               x(k) = temp1*sqrt(min(real(min(i,nx-i+1))*hx,temp))
   10       continue
   20    continue

         return

      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task.eq.'F' .or. task.eq.'FJ') then
         do 40 j = 1,ny
            do 30 i = 1,nx
               k = (j-1)*nx + i
               ut = zero
               ub = zero
               ul = zero
               ur = zero
               u = x(k)
               if (i.ne.1) ul = x(k-1)
               if (i.ne.nx) ur = x(k+1)
               if (j.ne.1) ub = x(k-nx)
               if (j.ne.ny) ut = x(k+nx)
               uxx = (-ur+two*u-ul)*hydhx
               uyy = (-ut+two*u-ub)*hxdhy
               fvec(k) = uxx + uyy - hxhy*lambda*exp(x(k))
   30       continue
   40    continue

         if (task.eq.'F') return

      end if

      if (task.eq.'J' .or. task.eq.'FJ') then

         do 60 j = 1,n
            do 50 i = 1,n
               fjac(i,j) = zero
   50       continue
   60    continue

c        Evaluate the Jacobian at x.

         do 80 j = 1,ny
            do 70 i = 1,nx
               k = (j-1)*nx + i
               if (i.ne.1) fjac(k,k-1) = -hydhx
               if (i.ne.nx) fjac(k,k+1) = -hydhx
               if (j.ne.1) fjac(k,k-nx) = -hxdhy
               if (j.ne.ny) fjac(k,k+nx) = -hxdhy
               fjac(k,k) = two* (hydhx+hxdhy) - hxhy*lambda*exp(x(k))
   70       continue
   80    continue

         return

      end if

      end
