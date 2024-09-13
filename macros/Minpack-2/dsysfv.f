      subroutine dsysfv(n,nx,ny,x,fvec,task,prob,par)
      character*60 task
      character*6 prob
      integer n, nx, ny
      double precision par
      double precision x(n), fvec(n)
c     *********
c
c     Subroutine dsysfv
c
c     This subroutine selects a nonlinear equations problem from
c     the MINPACK-2 test problem collection specified by the
c     character*6 variable prob.  Appropriate values for the problem
c     parameters can be found in the corresponding subroutines.
c
c     The subroutine statement is:
c
c       subroutine dsysfv(n,x,fvec,task,prob,par)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the vector x if task = 'F'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F'.
c            Otherwise x is set according to task.
c
c       fvec is a double precision array of dimension n.
c         On entry fvec need not be specified.
c         On exit fvec contains the function evaluated
c            at x if task = 'F'.
c
c       task is a character variable.
c         On entry task specifies the action of the subroutine:
c
c            task               action
c            ----               ------
c             'F'     Evaluate the function at x.
c             'XS'    Set x to the standard starting point xs.
c
c         On exit task is unchanged.
c
c       prob is a character*6 variable.
c         On entry prob specifies which problem is being solved.
c         On exit prob is set to 'ERROR' if prob is not an
c            acceptable problem name. Otherwise prob is unchanged.
c
c       par is a double precision variable.
c         On entry par sepcifies a problem-dependent parameter.
c         On exit x is unchanged.
c
c     Subprograms called
c
c       MINPACK-2 ...  dficfj, dsfdfj, dierfj, dsfifj,
c                      dfdcfj, dhhdfj, dspffj, dcprfj
c
c     MINPACK-2 Project. March 2000
c     Argonne National Laboratory.
c     Brett M. Averick and Jorge J. More'.
c
c     **********
      integer nint
      double precision a, b
      parameter (a=0.2d0,b=0.0d0)

      integer ldfjac
      parameter (ldfjac=1)
      double precision fjac(1,1)

c     Select a problem.

      if (prob .eq. 'DFIC') then
         nint = n/8
         call dficfj(n,x,fvec,fjac,ldfjac,task,par,nint)
      else if (prob .eq. 'DSFD') then
         nint = n/14
         call dsfdfj(n,x,fvec,fjac,ldfjac,task,par,nint)
      else if (prob .eq. 'DIER') then
         nint = n/15
         call dierfj(n,x,fvec,fjac,ldfjac,task,a,b,par,nint)
      else if (prob .eq. 'DSFI') then
         call dsfifj(nx,ny,x,fvec,fjac,ldfjac,task,par)
      else if (prob .eq. 'DFDC') then
         call dfdcfj(nx,ny,x,fvec,fjac,ldfjac,task,par)
      else if (prob(1:4) .eq. 'DHHD') then
         call dhhdfj(n,x,fvec,fjac,ldfjac,task,prob)
      else if (prob .eq. 'DCPF') then
         call dcpffj(n,x,fvec,fjac,ldfjac,task)
      else if (prob .eq. 'DCPR') then
         call dcprfj(n,x,fvec,fjac,ldfjac,task)
      else
         prob = 'ERROR'
      end if

      end
