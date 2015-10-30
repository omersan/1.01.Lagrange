!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Lagrange Polinomial Interpolation
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Aug. 12, 2015
!-----------------------------------------------------------------------------!

program lagrange_test
implicit none
integer::j,nd,np,io
real*8 ::x,f,fexact,xp,fp,pi
real*8,allocatable::xd(:),fd(:)


!---------------------------------------------!
!Given data set on equally spaced grid points
!---------------------------------------------!
nd = 10
allocate(xd(0:nd),fd(0:nd))
do j=0,nd
xd(j) =-1.0d0 + dfloat(j)*2.0d0/(dfloat(nd))
fd(j) = 1.0d0/(1.0d0+25.0d0*xd(j)**2)
end do

!Clustered towards to ends
!Chebyshev-Lobatto points (cosine distribution)
io=0 !set io=1 for clustered points
if (io.eq.1) then
	pi = 4.0d0*datan(1.0d0)
	do j=0,nd
	xd(j) =-dcos((j)*pi/dfloat(nd))
	fd(j) = 1.0d0/(1.0d0+25.0d0*xd(j)**2)
	end do
end if


! Writing data to text file
open(11, file="dataset.plt")
write(11,*)'variables ="x","f"'
do j=0,nd
write(11,*) xd(j),fd(j)
end do
close(11)


! Writing exact solution using np points
np = 2000 !use 2000 points to plot curve
open(12, file="exact_curve.plt")
write(12,*)'variables ="x","f"'
	do j=0,np
		xp =-1.0d0 + dfloat(j)*2.0d0/(dfloat(np))
		fp = 1.0d0/(1.0d0+25.0d0*xp**2)
		write(12,*) xp,fp
	end do
close(12)

!compare solutions at a desired point
x = 0.7d0
fexact = 1.0d0/(1.0d0+25.0d0*x**2)
call lagrange(nd,xd,fd,x,f)
write(*,100)"numerical:",f
write(*,100)"exact    :",fexact

100 format (A20,F12.4)

! Writing numerical solution using np points 
open(13, file="numerical_curve.plt")
write(13,*)'variables ="x","f"'
	do j=0,np
		xp =-1.0d0 + dfloat(j)*2.0d0/(dfloat(np))
		call lagrange(nd,xd,fd,xp,f)
		write(13,*) xp,f
	end do
close(13)


       
end

!-----------------------------------------------------------------------------!
!Construct Lagrange interpolation from a given data set fd(xd) 
!fd(i): given data values 
!xd(i): given points
!     where i=0,1,2,...nd 
!     nd: order of Lagrange polynomial (equal to the number of data point - 1)
!
!
!f: interpolated value
!x: interpolation location 
!-----------------------------------------------------------------------------!
subroutine lagrange(nd,xd,fd,x,f)
implicit none
integer::nd,i,j
real*8 ::xd(0:nd),fd(0:nd)
real*8 ::x,f,b

!compute Lagrange basis functions b (multiplication operation) 
!compute Lagrange interpolated value f at point x (summation operation)
f = 0.0d0
do j=0,nd
b = 1.0d0
	do i = 0, nd
		if (i.ne.j) then	
			b = b*(x - xd(i))/(xd(j)-xd(i))		
		end if
	end do
f = f + fd(j)*b
end do

end
