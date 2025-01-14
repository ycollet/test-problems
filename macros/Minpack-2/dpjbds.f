      subroutine dpjbds(nx,ny,w,ecc,b)
      integer nx, ny
      double precision ecc, b
      double precision w(nx*ny)
c     **********
c
c     Subroutine dpjbds
c
c     This subroutine computes a scaling for the pressure distribution
c     in a journal bearing problem as the square roots of the diagonal
c     elements of the Hessian matrix.
c
c     The subroutine statement is
c
c       subroutine dpjbds(nx,ny,w,ecc,b)
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
c       w is a double precision array of dimension nx*ny.
c         On entry w need not be specified.
c         On exit w contains the square roots of the diagonal elements
c            of the Hessian matrix.
c
c       ecc is a double precision variable.
c         On entry ecc is the eccentricity in (0,1).
c         On exit ecc is unchanged.
c
c       b is a double precision variable.
c         On entry b defines the domain as D = (0,2*pi) X (0,2*b).
c         On exit b is unchanged.
c
c     MINPACK-2 Project. March 1999.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick and Jorge J. More'.
c
c     **********
      double precision four, one, six, two
      parameter (one=1.0d0,two=2.0d0,four=4.0d0,six=6.0d0)

      integer i, j, k
      double precision hx, hxhx, hxhy, hy, hyhy, pc, pi, pl, pr, ruled,
     +                 rulel, ruler, ruleu, xi

      double precision p

      p(xi) = (1+ecc*cos(xi))**3

      pi = four*atan(one)
      hx = two*pi/dble(nx+1)
      hy = two*b/dble(ny+1)
      hxhx = hx*hx
      hyhy = hy*hy
      hxhy = hx*hy

c     Computation of the scaling.

      do 20 i = 1, nx
         xi = dble(i)*hx
         do 10 j = 1, ny
            k = nx*(j-1) + i
            pl = p(xi-hx)
            pc = p(xi)
            pr = p(xi+hx)
            ruler = hxhy*(pc+pr)/two
            rulel = hxhy*(pc+pl)/two
            ruleu = hxhy*(four*pc+pr+pl)/six
            ruled = ruleu
            w(k) = sqrt(ruler/hxhx+ruleu/hyhy+rulel/hxhx+ruled/hyhy)
   10    continue
   20 continue

      end
