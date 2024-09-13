      subroutine dchqfj(m,n,x,fvec,fjac,ldfjac,task)
      character*(*) task
      integer m,n,ldfjac
      double precision x(n),fvec(m),fjac(ldfjac,n)
c     **********
c
c     Subroutine dchqfj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the Chebyshev Quadrature problem. This problem arises in the 
c     determination of the nodes of a quadrature formula with equal 
c     weights.
c
c     The subroutine statement is:
c
c       subroutine dchqfj(m,n,x,fvec,fjac,ldfjac,task)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions.
c            For the Chebyshev quadrature problem, m >= n.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables.
c            For the Chebyshev quadrature problem, n <= m.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
c            x is set according to task.
c
c       fvec is a double precision array of dimension m.
c         On entry fvec need not be specified.
c         On exit fvec contains the function evaluated at x if
c            task = 'F' or 'FJ'.
c
c       fjac is a double precision array of dimension (ldfjac,n).
c         On entry fjac need not be specified.
c         On exit fjac contains the Jacobian matrix evaluated at x if
c            task = 'J' or 'FJ'.
c
c       ldfjac is an integer variable.
c          On entry is specifies the leading dimension of fjac.
c          On exit ldfjac is unchanged.
c
c            task               action
c            ----               ------
c             'F'     Evaluate the function at x.
c             'J'     Evaluate the Jacobian matrix at x.
c             'FJ'    Evaluate the function and the Jacobian at x.
c             'XS'    Set x to the standard starting point xs.
c             'XL'    Set x to the lower bound xl.
c             'XU'    Set x to the upper bound xu.
c
c         On exit task is unchanged.
c
c     MINPACK-2 Project. October 1992.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      double precision zero,one,two,four
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0)

      integer i,iev,j
      double precision dx,temp,temp1,temp2,temp3,temp4,ti

c     Check input arguments for errors.

      if (m .lt. n) then
         task = 'ERROR: M MUST BE .GE. N'
         return
      endif

c     Compute a lower bound for x if task = 'XL'.

      if (task .eq. 'XL') then
         do 10 i = 1, n
            x(i) = zero
   10    continue

         return

      endif

c     Compute an upper bound for x if task = 'XU'.

      if (task .eq. 'XU') then
         do 20 i = 1, n
            x(i) = one
   20    continue

         return

      endif

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         dx = one/dble(n + 1)
         do 30 j = 1, n
            x(j) = dble(j)*dx
   30    continue

         return

      endif

c     Evaluate the function if task = 'F', the Jacobian matrix if 
c     task = 'J', or both if task = 'FJ'.

      dx = one/dble(n)

      if (task .eq. 'F' .or. task .eq. 'FJ') then
         do 40 i = 1, m
            fvec(i) = zero
   40    continue
         do 60 j = 1, n
            temp1 = one
            temp2 = two*x(j) - one
            temp = two*temp2
            do 50 i = 1, m
               fvec(i) = fvec(i) + temp2
               ti = temp*temp2 - temp1
               temp1 = temp2
               temp2 = ti
   50       continue
   60    continue
         iev = -1
         do 70 i = 1, m
            fvec(i) = dx*fvec(i)
            if (iev .gt. 0) fvec(i) = fvec(i) + one/dble(i**2-1)
            iev = -iev
   70    continue

         if (task .eq. 'F') return

      endif

      if (task .eq. 'J' .or. task .eq. 'FJ') then
         do 90 j = 1, n
            temp1 = one
            temp2 = two*x(j) - one
            temp = two*temp2
            temp3 = zero
            temp4 = two
            do 80 i = 1, m
               fjac(i,j) = dx*temp4
               ti = four*temp2 + temp*temp4 - temp3
               temp3 = temp4
               temp4 = ti
               ti = temp*temp2 - temp1
               temp1 = temp2
               temp2 = ti
   80       continue
   90    continue
 
         return

      endif

      end
