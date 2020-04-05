!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                                   !!!
!!!   A program which solves the intego-differential coupled equations given in the   !!!
!!!   paper by Carsten Bauer, Andreas Rückriegel, Anand Sharma, and Peter Kopietz     !!!
!!!   DOI: 10.1103/PhysRevB.92.121409. This program in addition also takes into       !!!
!!!   account the effect of frequency. Velocity is calculated using earlier methods   !!!
!!!   and used for the calculation of dielectric function here.                       !!!
!!!                                                                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program frequency_dependence

	real :: phi,momentum,cutoff,dcutoff
	integer,parameter :: n=501,m=297,num=2000
	real :: pi,cphi,integral_vel,integral_eps,eps(n,m),vel(n,m)
	real :: epsilon_1,epsilon_2,vel_1,vel_2,omega
	integer :: index1,index2

	pi = 3.14159
	dcutoff = 1.0/n
	dphi = pi/(2.0*num)
	omega = 1.0

	! Reading data of the velocity from file 
	open(1,file='vel.dat',status='old')

	do i=1,n
		do j=1,m
			read(1,*) vel(i,j)
		end do
	end do

	close(1)

   ! Initialise dielectric function to unity at all momenta at the highest cutoff
	do j=1,m
		eps(n,j) = 1.0
	end do

	do i=1,n-1

		cutoff = 1.0*(n-i)/n ! Start from the highest cutoff 1 and gradually decrease it

		do j=1,m

			momentum = 1.0*j/m

			integral_eps = 0.0

			do k=1,num

				phi = k*dphi

				cphi = COS(phi)

				! calculating velocities to be used as the inputs of dielectric function.
				index1 = NINT(cutoff*m)+1

				vel_1 = vel(n-i+1,index1)

				index2 = NINT(m*(cutoff+momentum*cphi))+1

				if (index2<m) then
					vel_2 = vel(n-i+1,index2)
				else
					vel_2 = 1.
				end if

				! Theta function implementation
					if (cphi>1 - 2*cutoff/momentum) then
						integral_eps = integral_eps + dielectric(phi,momentum,cutoff,vel_1,vel_2,omega)
					end if

			end do
			eps(n-i,j) = eps(n-i+1,j) + dcutoff*dphi*integral_eps
		end do
	end do

	open(10,file="epson0_1.0.dat")

	! do i=1,n
		do j=1,m
			write(10,*) 1.0*j/m, eps(2,j)
		end do
		! write(10,*)
	! end do
	close(10)

	STOP

	contains

	function dielectric(phi,momentum,cutoff,vel_1,vel_2,omega) result(value)
	! Returns the integrand for dielectric function renormalisation

		real :: phi,momentum,cutoff,value,cphi,sphi,pi,vel_1,vel_2,omega

		cphi = COS(phi)
		sphi = SIN(phi)
		pi = 3.14159

		value = 4.4*momentum*(sphi**2)*(cutoff*vel_1+(cutoff+momentum*cphi)*vel_2) &
			/(pi*((cutoff*vel_1+(cutoff+momentum*cphi)*vel_2)**2 + omega**2)* &
			SQRT((2*cutoff + momentum*cphi)**2 - momentum**2))

	end function dielectric

end program frequency_dependence
