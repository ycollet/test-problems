      subroutine dsysjs(n,nx,ny,x,s,y,prob,par)
      character*6 prob
      integer n, nx, ny
      double precision par
      double precision x(n), s(n), y(n)
c     *********
c
c     Subroutine dsysjs
c
c     This subroutine computes the product f'(x)*s = y, where
c     f'(x) is the Jacobian matrix for the nonlinear equations 
c     problem from the MINPACK-2 test problem collection specified
c     by the character variable prob.
c
c     The subroutine statement is
c
c       subroutine dsysjs(n,nx,ny,x,y,prob,par)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       nx is an integer variable.
c         On entry nx is the number of grid points in the first
c            coordinate direction.
c         On exit nx is unchanged.
c
c       ny is an integer variable.
c         On entry ny is the number of grid points in the second
c            coordinate direction. If the problem is formulated in
c            one spatial dimension, ny = 1.
c         On exit ny is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the vector x.
c         On exit x is unchanged.
c
c       s is a double precision array of dimension n.
c         On entry s contains the vector s.
c         On exit s is unchanged.
c
c       y is a double precision array of dimension n.
c         On entry out need not be specified.
c         On exit y contains the product f'(x)*s.
c
c       prob is a character*6 variable.
c         On entry prob specifies the problem.
c         On exit prob is set to 'ERROR' if prob is not an
c            acceptable problem name. Otherwise prob is unchanged.
c
c       par is a double precision variable.
c         On entry par specifies a probem-dependent parameter.
c         On exit par is unchanged.
c
c     Subprograms called
c
c       MINPACK-2  ...  dfdcjs, dficjs, dierjs, dsfdjs, dsfijs
c
c     MINPACK-2 Project. April 2000.
c     Argonne National Laboratory.
c     Jorge J. More'.
c
c     **********
      integer nint
      double precision a, b
      parameter (a=0.2d0,b=0.0d0)

c     Select a problem.

      if (prob .eq. 'DFDC') then
         call dfdcjs(nx,ny,x,s,y,par)
      else if (prob .eq. 'DFIC') then
         nint = n/8
         call dficjs(n,x,s,y,par,nint)
       else if (prob .eq. 'DIER') then
         nint = n/15
         call dierjs(n,x,s,y,nint)
      else if (prob .eq. 'DSFD') then
         nint = n/14
         call dsfdjs(n,x,s,y,par,nint)
      else if (prob .eq. 'DSFI') then
         call dsfijs(nx,ny,x,s,y,par)
      else
         prob = 'ERROR'
      endif 

      end

