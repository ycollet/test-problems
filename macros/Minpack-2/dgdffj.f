      subroutine dgdffj(m,n,x,fvec,fjac,ldfjac,task)
      character*60 task
      integer m, n, ldfjac
      double precision x(n), fvec(m), fjac(ldfjac,n)
c     **********
c
c     Subroutine dgdffj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the Gaussian data fitting problem.
c
c     The subroutine statement is
c
c       subroutine dgdffj(m,n,x,fvec,fjac,ldfjac,task)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions.
c            For the exponential data fitting II problem m = 65.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables.
c            For the exponential data fitting II problem n = 11.
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
c             'XL'    Set x to the lower bound xl.
c             'XU'    Set x to the upper bound xu.
c
c         On exit task is unchanged.
c
c     MINPACK-2 Project. March 1999.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      double precision ten, two, zero
      parameter (zero=0.0d0,two=2.0d0,ten=1.0d1)

      integer i
      double precision temp, temp1, temp2, temp3, temp4
      double precision y(65)

      data y/1.366d0, 1.191d0, 1.112d0, 1.013d0, 9.91d-1, 8.85d-1,
     +     8.31d-1, 8.47d-1, 7.86d-1, 7.25d-1, 7.46d-1, 6.79d-1,
     +     6.08d-1, 6.55d-1, 6.16d-1, 6.06d-1, 6.02d-1, 6.26d-1,
     +     6.51d-1, 7.24d-1, 6.49d-1, 6.49d-1, 6.94d-1, 6.44d-1,
     +     6.24d-1, 6.61d-1, 6.12d-1, 5.58d-1, 5.33d-1, 4.95d-1, 5.0d-1,
     +     4.23d-1, 3.95d-1, 3.75d-1, 3.72d-1, 3.91d-1, 3.96d-1,
     +     4.05d-1, 4.28d-1, 4.29d-1, 5.23d-1, 5.62d-1, 6.07d-1,
     +     6.53d-1, 6.72d-1, 7.08d-1, 6.33d-1, 6.68d-1, 6.45d-1,
     +     6.32d-1, 5.91d-1, 5.59d-1, 5.97d-1, 6.25d-1, 7.39d-1, 7.1d-1,
     +     7.29d-1, 7.2d-1, 6.36d-1, 5.81d-1, 4.28d-1, 2.92d-1, 1.62d-1,
     +     9.8d-2, 5.4d-2/

c     Check input arguments for errors.

      if (m .ne. 65 .or. n .ne. 11) then
         task = 'ERROR: M .NE. 65 OR N .NE. 11 IN DGDFFJ'

         return

      end if

c     Compute a lower bound for x if task = 'XL'.

      if (task(1:2) .eq. 'XL') then
         do 10 i = 1, n
            x(i) = zero
   10    continue

         return

      end if

c     Compute an upper bound for x if task = 'XU'.

      if (task(1:2) .eq. 'XU') then
         do 20 i = 1, n
            x(i) = ten
   20    continue

         return

      end if

c     Compute the standard starting point if task = 'XS'.

      if (task(1:2) .eq. 'XS') then
         x(1) = 1.3d0
         x(2) = 6.5d-1
         x(3) = 6.5d-1
         x(4) = 7.0d-1
         x(5) = 6.0d-1
         x(6) = 3.0d0
         x(7) = 5.0d0
         x(8) = 7.0d0
         x(9) = 2.0d0
         x(10) = 4.5d0
         x(11) = 5.5d0

         return

      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      do 30 i = 1, m
         temp = dble(i-1)/ten
         temp1 = exp(-x(5)*temp)
         temp2 = exp(-x(6)*(temp-x(9))**2)
         temp3 = exp(-x(7)*(temp-x(10))**2)
         temp4 = exp(-x(8)*(temp-x(11))**2)
         if (task(1:1) .eq. 'F' .or. task(1:2) .eq. 'FJ') fvec(i) = 
     +       y(i) - (x(1)*temp1+x(2)*temp2+x(3)*temp3+x(4)*temp4)
         if (task(1:1) .eq. 'J' .or. task(1:2) .eq. 'FJ') then
            fjac(i,1) = -temp1
            fjac(i,2) = -temp2
            fjac(i,3) = -temp3
            fjac(i,4) = -temp4
            fjac(i,5) = temp*x(1)*temp1
            fjac(i,6) = x(2)*(temp-x(9))**2*temp2
            fjac(i,7) = x(3)*(temp-x(10))**2*temp3
            fjac(i,8) = x(4)*(temp-x(11))**2*temp4
            fjac(i,9) = -two*x(2)*x(6)*(temp-x(9))*temp2
            fjac(i,10) = -two*x(3)*x(7)*(temp-x(10))*temp3
            fjac(i,11) = -two*x(4)*x(8)*(temp-x(11))*temp4
         end if
   30 continue

      end
