      subroutine septfg(nx,ny,x,f,fgrad,task,c)
      character*(*) task
      integer nx,ny
      real f,c
      real x(nx*ny),fgrad(nx*ny)
c     **********
c
c     Subroutine septfg
c
c     This subroutine computes the function and gradient of the
c     Elastic-Plastic Torsion problem as formulated by, eg., Glowinski.
c     This problem arises in the determination of the stress field on an
c     infinitely long cylindrical bar subject to a torsional couple
c     characterized by an angle of twist per unit length c. A finite
c     element approximation to the torsion problem is obtained by using
c     piecewise linear basis functions over a triangulation of the unit
c     square.
c
c     The subroutine statement is:
c
c       subroutine septfg(nx,ny,x,f,fgrad,task,c)
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
c       x is a real array of dimension nx*ny.
c         On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
c            x is set according to task.
c
c       f is a real variable.
c         On entry f need not be specified.
c         On exit f is set to the function evaluated at x if task = 'F'
c            or 'FG'.
c
c       fgrad is a real array of dimension nx*ny.
c         On entry fgrad need not be specified.
c         On exit fgrad contains the gradient evaluated at x if
c            task = 'G' or 'FG'.
c
c       task is a character variable.
c         On entry task specifies the action of the subroutine:
c
c            task               action
c            ----               ------
c             'F'     Evaluate the function at x.
c             'G'     Evaluate the gradient at x.
c             'FG'    Evaluate the function and the gradient at x.
c             'XS'    Set x to the standard starting point xs.
c             'XL'    Set x to the lower bound xl.
c             'XU'    Set x to the upper bound xu.
c
c         On exit task is unchanged.
c
c       c is a real variable.
c         On entry c is the angle of twist per unit length.
c         On exit c is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick and Jorge J. More'.
c
c     **********
      real zero,p5,one,three
      parameter (zero=0.0,p5=0.5,one=1.0,three=3.0)

      integer i,j,k
      real flin,fquad,dvdx,dvdy,hx,hy,v,vb,vl,vr,vt,temp,temp1,area,
     +     cdiv3

      hx = one/real(nx+1)
      hy = one/real(ny+1)
      area = p5*hx*hy
      cdiv3 = c/three

c     Compute a lower bound for x if task = 'XL' or an upper bound if
c     task = 'XU'.

      if (task.eq.'XL' .or. task.eq.'XU') then
         if (task.eq.'XL') temp1 = -one
         if (task.eq.'XU') temp1 = one
         do 20 j = 1,ny
            temp = real(min(j,ny-j+1))*hy
            do 10 i = 1,nx
               k = nx* (j-1) + i
               x(k) = sign(min(real(min(i,nx-i+1))*hx,temp),temp1)
   10       continue
   20    continue

         return

      end if

c     Compute the standard starting point if task = 'XS'.

      if (task.eq.'XS') then
         do 40 j = 1,ny
            temp = real(min(j,ny-j+1))*hy
            do 30 i = 1,nx
               k = nx* (j-1) + i
               x(k) = min(real(min(i,nx-i+1))*hx,temp)
   30       continue
   40    continue

         return

      end if

c     Compute the function if task = 'F', the gradient if task = 'G', or
c     both if task = 'FG'.

      if (task.eq.'F' .or. task.eq.'FG') then
         fquad = zero
         flin = zero
      end if

      if (task.eq.'G' .or. task.eq.'FG') then
         do 50 k = 1,nx*ny
            fgrad(k) = zero
   50    continue
      end if

c     Computation of the function and the gradient over the lower
c     triangular elements.

      do 70 j = 0,ny
         do 60 i = 0,nx
            k = nx* (j-1) + i
            v = zero
            vr = zero
            vt = zero
            if (i.ge.1 .and. j.ge.1) v = x(k)
            if (i.lt.nx .and. j.gt.0) vr = x(k+1)
            if (i.gt.0 .and. j.lt.ny) vt = x(k+nx)
            dvdx = (vr-v)/hx
            dvdy = (vt-v)/hy
            if (task.eq.'F' .or. task.eq.'FG') then
               fquad = fquad + dvdx**2 + dvdy**2
               flin = flin - cdiv3* (v+vr+vt)
            end if

            if (task.eq.'G' .or. task.eq.'FG') then
               if (i.ne.0 .and. j.ne.0) fgrad(k) = fgrad(k) - dvdx/hx -
     +         dvdy/hy - cdiv3
               if (i.ne.nx .and. j.ne.0) fgrad(k+1) = fgrad(k+1) +
     +         dvdx/hx - cdiv3
               if (i.ne.0 .and. j.ne.ny) fgrad(k+nx) = fgrad(k+nx) +
     +         dvdy/hy - cdiv3
            end if

   60    continue
   70 continue

c     Computation of the function and the gradient over the upper
c     triangular elements.

      do 90 j = 1,ny + 1
         do 80 i = 1,nx + 1
            k = nx* (j-1) + i
            vb = zero
            vl = zero
            v = zero
            if (i.le.nx .and. j.gt.1) vb = x(k-nx)
            if (i.gt.1 .and. j.le.ny) vl = x(k-1)
            if (i.le.nx .and. j.le.ny) v = x(k)
            dvdx = (v-vl)/hx
            dvdy = (v-vb)/hy
            if (task.eq.'F' .or. task.eq.'FG') then
               fquad = fquad + dvdx**2 + dvdy**2
               flin = flin - cdiv3* (vb+vl+v)
            end if

            if (task.eq.'G' .or. task.eq.'FG') then
               if (i.ne.nx+1 .and. j.ne.1) fgrad(k-nx) = fgrad(k-nx) -
     +         dvdy/hy - cdiv3
               if (i.ne.1 .and. j.ne.ny+1) fgrad(k-1) = fgrad(k-1) -
     +         dvdx/hx - cdiv3
               if (i.ne.nx+1 .and. j.ne.ny+1) fgrad(k) = fgrad(k) +
     +         dvdx/hx + dvdy/hy - cdiv3
            end if

   80    continue
   90 continue

c     Scale the result.

      if (task.eq.'F' .or. task.eq.'FG') then
         f = area* (p5*fquad+flin)
      end if

      if (task.eq.'G' .or. task.eq.'FG') then
         do 100 k = 1,nx*ny
            fgrad(k) = area*fgrad(k)
  100    continue
      end if

      return

      end
