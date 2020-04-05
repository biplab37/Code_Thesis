!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                                   !!!
!!!   A program which solves the intego-differential coupled equations given in the   !!!
!!!   paper by Carsten Bauer, Andreas Rückriegel, Anand Sharma, and Peter Kopietz     !!!
!!!   DOI: 10.1103/PhysRevB.92.121409                                                 !!!
!!!                                                                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program bauer

	real :: phi,momentum,cutoff,dcutoff
	integer,parameter  :: n=501,m=297,num=2000
	real :: pi,cphi1,cphi2,integral_vel,integral_eps,eps(n,m),vel(n,m)
	integer :: index1,index2

	pi = 3.14159
	dcutoff = 1.0/n

   ! Initialise velocity and dielectric function to unity at all momenta at the highest 
   ! cutoff
	do j=1,m
		eps(n,j) = 1.0
		vel(n,j) = 1.0
	end do

	do i=1,n-1

		cutoff = 1.0*(n-i)/n ! Start from the highest cutoff 1 and gradually decrease it

		do j=1,m

			momentum = 1.0*j/m

			integral_vel = 0.0
			integral_eps = 0.0

			do k=1,num
				! phi integration

				phi = k*pi/num

				! Calculating cos(phi) to reduce number of function calls.
				! For the dielctric function integral limits are from 0 to pi/2 hence
				! phi/2 is used.

				cphi1 = COS(phi)
				cphi2 = COS(phi/2.)

				! calculation of the velocities and dielectric functions to be used to
				! calculate velocity and epsilon integrand at the next step.

				! epsilon to be calculated at momentum |k-q|. if it exceeds the highest
				! cutoff then 1 is used.

				index1 = NINT(SQRT(cutoff**2 + momentum**2 - 2*momentum*cutoff*cphi1)*m)+1
				if (index1<m) then
					epsilon = eps(n-i+1,index1)
				else
					epsilon = 1.
				end if

				! velocity at momentum = cutoff
				vel_1 = vel(n-i+1,NINT(cutoff*m)+1)

				! other velocity to be calcualted at mometa = cutoff + momentum*cos(phi).
				! If it exceeds the highest cutoff 1, which is the v_F, is used.
				index2 = NINT(m*(cutoff+momentum*cphi2))+1

				if (index2<m) then
					vel_2 = vel(n-i+1,index2)
				else
					vel_2 = 1.
				end if

				! Theta function implementation
				if (cphi2>1 - 2*cutoff/momentum) then
					integral_eps = integral_eps + dielectric(phi/2.,momentum,cutoff,vel_1,vel_2)
				end if

				integral_vel = integral_vel + velocity(cphi1,momentum,cutoff,epsilon)

			end do
			! update for the next cutoff
			eps(n-i,j) = eps(n-i+1,j) + dcutoff*integral_eps*pi/(num*2.)
			vel(n-i,j) = vel(n-i+1,j) + dcutoff*integral_vel*pi/num
		end do
	end do

	open(2,file="vel.dat")

	do i=1,n
		do j=1,m
			write(2,*) vel(i,j)
		end do
		write(2,*)
	end do
	close(2)

	STOP

	contains

	function velocity(cphi,momentum,cutoff,epsilon) result(value)
	! Returns the integrand for the velocity renormalisation
		real :: momentum,cutoff,value,cphi,pi,epsilon

		pi = 3.14159

		value = 2.2*cphi*cutoff/(2*pi*epsilon*momentum*SQRT(cutoff**2 - 2*cphi*momentum*cutoff + momentum**2))

	end function velocity

	function dielectric(phi,momentum,cutoff,vel_1,vel_2) result(value)
	! Returns the integrand for dielectric function renormalisation

		real :: phi,momentum,cutoff,value,cphi,sphi,pi,vel_1,vel_2
		real,parameter :: highcutoff = .2

		cphi = COS(phi)
		sphi = SIN(phi)
		pi = 3.14159

		value = 4.4*momentum*sphi**2/(pi*(cutoff*vel_1+(cutoff+momentum*cphi)*vel_2)* &
			SQRT((2*cutoff + momentum*cphi)**2 - momentum**2))

	end function dielectric

end program bauer
