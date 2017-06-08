C     ***QUARTIC************************************************25.03.98
C     Solution of a quartic equation
C     ref.: J. E. Hacke, Amer. Math. Monthly, Vol. 48, 327-328, (1941)
C     NO WARRANTY, ALWAYS TEST THIS SUBROUTINE AFTER DOWNLOADING
C     ******************************************************************
C     dd(0:4)     (i)  vector containing the polynomial coefficients
C     sol(1:4)    (o)  results, real part
C     soli(1:4)   (o)  results, imaginary part
C     Nsol        (o)  number of real solutions 
C     ==================================================================
      subroutine quartic(dd,sol,soli,Nsol)
      implicit double precision (a-h,o-z)
      dimension dd(0:4),sol(4),soli(4)
      dimension AA(0:3),z(3)
C
      Nsol = 0
      a = dd(4)
      b = dd(3)
      c = dd(2)
      d = dd(1)
      e = dd(0)
C
      if (dd(4).eq.0.d+0) then
	write(6,*)'ERROR: NOT A QUARTIC EQUATION'
	return
      endif
C
      p = (-3.d+0*b**2 + 8.d+0*a*c)/(8.d+0*a**2)
      q = (b**3 - 4.d+0*a*b*c + 8.d+0*d*a**2)/(8.d+0*a**3)
      r = (-3.d+0*b**4 + 16.d+0*a*b**2*c - 64.d+0*a**2*b*d + 
     &      256.d+0*a**3*e)/(256.d+0*a**4)
C
C     solve cubic resolvent
      AA(3) =  8.d+0
      AA(2) = -4.d+0*p 
      AA(1) = -8.d+0*r
      AA(0) =  4.d+0*p*r - q**2
      call cubic(AA,z,ncube)
C      
      zsol = -1.d+99
      do 5 i=1,ncube
 5      zsol = max(zsol,z(i))
      z(1) = zsol
      xK2 = 2.d+0 * z(1) - p
      xK  = sqrt(xK2)
C-----------------------------------------------
      if (xK.eq.0.d+0) then
        xL2 = z(1)**2 - r
	if (xL2.lt.0.d+0) then
	  write(6,*)'Sorry, no solution'
	  return
        endif
	xL  = sqrt(xL2)
      else
        xL = q/(2.d+0 * xK)
      endif
C-----------------------------------------------
      sqp = xK2 - 4.d+0*(z(1) + xL)
      sqm = xK2 - 4.d+0*(z(1) - xL)
C
      do 10 i=1,4
 10     soli(i) = 0.d+0
      if       (sqp.ge.0.d+0 .and. sqm.ge.0.d+0) then
	sol(1) = 0.5d+0*( xK + sqrt(sqp))
	sol(2) = 0.5d+0*( xK - sqrt(sqp))
	sol(3) = 0.5d+0*(-xK + sqrt(sqm))
	sol(4) = 0.5d+0*(-xK - sqrt(sqm))
	Nsol = 4
      else if  (sqp.ge.0.d+0 .and. sqm.lt.0.d+0) then
	sol(1) =  0.5d+0*(xK + sqrt(sqp))
	sol(2) =  0.5d+0*(xK - sqrt(sqp))
	sol(3) = -0.5d+0*xK 
	sol(4) = -0.5d+0*xK 
	soli(3) =  sqrt(-0.25d+0 * sqm)
	soli(4) = -sqrt(-0.25d+0 * sqm)
	Nsol = 2
      else if  (sqp.lt.0.d+0 .and. sqm.ge.0.d+0) then
	sol(1) = 0.5d+0*(-xK + sqrt(sqm))
	sol(2) = 0.5d+0*(-xK - sqrt(sqm))
	sol(3) =  0.5d+0*xK 
	sol(4) =  0.5d+0*xK 
	soli(3) =  sqrt(-0.25d+0 * sqp)
	soli(4) = -sqrt(-0.25d+0 * sqp)
	Nsol = 2
      else if  (sqp.lt.0.d+0 .and. sqm.lt.0.d+0) then
	sol(1) = -0.5d+0*xK 
	sol(2) = -0.5d+0*xK 
	soli(1) =  sqrt(-0.25d+0 * sqm)
	soli(2) = -sqrt(-0.25d+0 * sqm)
	sol(3) =  0.5d+0*xK 
	sol(4) =  0.5d+0*xK 
	soli(3) =  sqrt(-0.25d+0 * sqp)
	soli(4) = -sqrt(-0.25d+0 * sqp)
	Nsol = 0
      endif
      do 20 i=1,4
 20     sol(i) = sol(i) - b/(4.d+0*a)
C
      return
      END


