      subroutine dfdcfj(nx,ny,x,fvec,fjac,ldfjac,task,r)
      character*60 task
      integer nx, ny, ldfjac
      double precision r
      double precision x(nx*ny), fvec(nx*ny), fjac(ldfjac,nx*ny)
c     **********
c
c     Subroutine dfdcfj
c
c     This subroutine computes the function and Jacobian matrix of the
c     flow in a driven cavity problem.
c
c     The subroutine statement is
c
c       subroutine dfdcfj(nx,ny,x,fvec,fjac,ldfjac,task,r)
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
c       x is a double precision array of dimension nx*ny.
c         On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
c            x is set according to task.
c
c       fvec is a double precision array of dimension nx*ny.
c         On entry fvec need not be specified.
c         On exit fvec contains the function evaluated at x if
c            task = 'F' or 'FJ'.
c
c       fjac is a double precision array of dimension (ldfjac,nx*ny).
c         On entry fjac need not be specified.
c         On exit fjac contains the Jacobian matrix evaluated at x if
c            task = 'J' or 'FJ'.
c
c       ldfjac is an integer variable.
c          On entry ldfjac is the leading dimension of fjac.
c          On exit ldfjac is unchanged.
c
c       task is a character*60 variable.
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
c       r is a double precision variable.
c         On entry r is the Reynolds number.
c         On exit r is unchanged.
c
c     MINPACK-2 Project. March 1999.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      double precision four, one, three, two, zero
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,four=4.0d0)

      logical feval, jeval
      integer i, j, k, n
      double precision dpdx, dpdy, hx, hx2, hx4, hx2hy2, hy, hy2, hy4,
     +                 p, pb, pbb, pbl, pblap, pbr, pl, plap, pll,
     +                 pllap, pr, prlap, prr, pt, ptl, ptlap, ptr, ptt,
     +                 xx, yy

      n = nx*ny
      hx = one/dble(nx+1)
      hy = one/dble(ny+1)
      hx2 = hx*hx
      hy2 = hy*hy
      hx4 = hx2*hx2
      hy4 = hy2*hy2
      hx2hy2 = hx2*hy2

      if (task(1:2) .eq. 'XS') then
         yy = hy
         do 20 i = 1, ny
            xx = hx
            do 10 j = 1, nx
               k = (i-1)*nx + j
               x(k) = -xx*(one-xx)*yy*(one-yy)
               xx = xx + hx
   10       continue
            xx = hx
            yy = yy + hy
   20    continue

         return

      end if

      if (task(1:1) .eq. 'F' .or. task(1:2) .eq. 'FJ') then
         feval = .true.
      else
         feval = .false.
      end if
      if (task(1:1) .eq. 'J' .or. task(1:2) .eq. 'FJ') then
         jeval = .true.
      else
         jeval = .false.
      end if

      if (jeval) then
         do 40 j = 1, n
            do 30 i = 1, n
               fjac(i,j) = zero
   30       continue
   40    continue
      end if

      do 60 i = 1, ny
         do 50 j = 1, nx
            k = (i-1)*nx + j
            if (i .eq. 1 .or. j .eq. 1) then
               pbl = zero
            else
               pbl = x(k-nx-1)
            end if
            if (i .eq. 1) then
               pb = zero
               pbb = x(k)
            else if (i .eq. 2) then
               pb = x(k-nx)
               pbb = zero
            else
               pb = x(k-nx)
               pbb = x(k-2*nx)
            end if
            if (i .eq. 1 .or. j .eq. nx) then
               pbr = zero
            else
               pbr = x(k-nx+1)
            end if
            if (j .eq. 1) then
               pl = zero
               pll = x(k)
            else if (j .eq. 2) then
               pl = x(k-1)
               pll = zero
            else
               pl = x(k-1)
               pll = x(k-2)
            end if
            p = x(k)
            if (j .eq. nx-1) then
               pr = x(k+1)
               prr = zero
            else if (j .eq. nx) then
               pr = zero
               prr = x(k)
            else
               pr = x(k+1)
               prr = x(k+2)
            end if
            if (i .eq. ny .or. j .eq. 1) then
               ptl = zero
            else
               ptl = x(k+nx-1)
            end if
            if (i .eq. ny-1) then
               pt = x(k+nx)
               ptt = zero
            else if (i .eq. ny) then
               pt = zero
               ptt = x(k) + two*hy
            else
               pt = x(k+nx)
               ptt = x(k+2*nx)
            end if
            if (i .eq. ny .or. j .eq. nx) then
               ptr = zero
            else
               ptr = x(k+nx+1)
            end if

            dpdy = (pt-pb)/(two*hy)
            dpdx = (pr-pl)/(two*hx)

c           Laplacians at each point in the 5 point stencil.

            pblap = (pbr-two*pb+pbl)/hx2 + (p-two*pb+pbb)/hy2
            pllap = (p-two*pl+pll)/hx2 + (ptl-two*pl+pbl)/hy2
            plap = (pr-two*p+pl)/hx2 + (pt-two*p+pb)/hy2
            prlap = (prr-two*pr+p)/hx2 + (ptr-two*pr+pbr)/hy2
            ptlap = (ptr-two*pt+ptl)/hx2 + (ptt-two*pt+p)/hy2

            if (feval) then
               fvec(k) = (prlap-two*plap+pllap)/hx2 +
     +                   (ptlap-two*plap+pblap)/hy2 -
     +                   r*(dpdy*(prlap-pllap)/(two*hx)-
     +                   dpdx*(ptlap-pblap)/(two*hy))
            end if

            if (jeval) then
               if (i .gt. 2) fjac(k,k-2*nx) = one/hy4 -
     +             r*dpdx/hy2/(two*hy)
               if (j .gt. 2) fjac(k,k-2) = one/hx4 + r*dpdy/hx2/(two*hx)
               if (j .lt. nx-1) fjac(k,k+2) = one/hx4 -
     +             r*dpdy/hx2/(two*hx)
               if (i .lt. ny-1) fjac(k,k+2*nx) = one/hy4 +
     +             r*dpdx/hy2/(two*hy)
               if (i .ne. 1 .and. j .ne. 1) fjac(k,k-nx-1) = two/
     +             hx2hy2 + r*(dpdy/hy2/(two*hx)-dpdx/hx2/(two*hy))
               if (i .ne. 1 .and. j .ne. nx) fjac(k,k-nx+1) = two/
     +             hx2hy2 - r*(dpdy/hy2/(two*hx)+dpdx/hx2/(two*hy))
               if (i .ne. ny .and. j .ne. 1) fjac(k,k+nx-1) = two/
     +             hx2hy2 + r*(dpdy/hy2/(two*hx)+dpdx/hx2/(two*hy))
               if (i .ne. ny .and. j .ne. nx) fjac(k,k+nx+1) = two/
     +             hx2hy2 - r*(dpdy/hy2/(two*hx)-dpdx/hx2/(two*hy))
               if (i .ne. 1) fjac(k,k-nx) = -four*(one/hx2hy2+one/hy4) +
     +             r*((prlap-pllap)/(two*hx)/(two*hy)+
     +             dpdx*(one/hx2+one/hy2)/hy)
               if (j .ne. 1) fjac(k,k-1) = -four*(one/hx2hy2+one/hx4) -
     +                                     r*((ptlap-pblap)/(two*hx)/
     +                                     (two*hy)+dpdy*
     +                                     (one/hx2+one/hy2)/hx)
               if (j .ne. nx) fjac(k,k+1) = -four*(one/hx2hy2+one/hx4) +
     +             r*((ptlap-pblap)/(two*hx)/(two*hy)+
     +             dpdy*(one/hx2+one/hy2)/hx)
               if (i .ne. ny) fjac(k,k+nx) = -four*
     +             (one/hx2hy2+one/hy4) - r*((prlap-pllap)/(two*hx)/
     +             (two*hy)+dpdx*(one/hx2+one/hy2)/hy)
               fjac(k,k) = two*(four/hx2hy2+three/hx4+three/hy4)
               if (i .eq. 1) fjac(k,k) = fjac(k,k) + one/hy4 -
     +                                   r*(dpdx/hy2/(two*hy))
               if (j .eq. 1) fjac(k,k) = fjac(k,k) + one/hx4 +
     +                                   r*(dpdy/hx2/(two*hx))
               if (j .eq. nx) fjac(k,k) = fjac(k,k) + one/hx4 -
     +                                    r*(dpdy/hx2/(two*hx))
               if (i .eq. ny) fjac(k,k) = fjac(k,k) + one/hy4 +
     +                                    r*(dpdx/hy2/(two*hy))
            end if
   50    continue
   60 continue

c     Scale the Result.  This is not desired if preconditioning.

      if (task(1:1) .eq. 'F' .or. task(1:2) .eq. 'FJ') then
         do 70 k = 1, n
            fvec(k) = hx2hy2*fvec(k)
   70    continue
      end if

      if (task(1:1) .eq. 'J' .or. task(1:2) .eq. 'FJ') then
         do 90 j = 1, n
            do 80 i = 1, n
               fjac(i,j) = hx2hy2*fjac(i,j)
   80       continue
   90    continue
      end if

      end
