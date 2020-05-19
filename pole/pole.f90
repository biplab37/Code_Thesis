!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                                   !!!
!!!   A program which solves the intrego-differential coupled equations given in the   !!!
!!!   paper by Carsten Bauer, Andreas Rückriegel, Anand Sharma, and Peter Kopietz     !!!
!!!   DOI: 10.1103/PhysRevB.92.121409. Here we have kept the bosonic mometa inside    !!!
!!!   the mometum shell to be integrated out.                                         !!!
!!!                                                                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program pole

	real :: phi,momentum,cutoff,dcutoff,dphi
	integer,parameter :: n=501,m=297,num=2000
	real :: pi,cphi,integral_vel,integral_eps,eps(n,m),vel(n,m),eps2(n,m),wavefunc(n,m)
	real :: epsilon_1,epsilon_2,vel_1,vel_2,integral_eps2,Y_1,Y_2,integral_wavefunc
	integer :: index1,index2

	pi = 3.14159

	dcutoff = 1.0/n
	dphi = pi/(2.*num)

   ! Initialise velocity and dielectric function to unity at all momenta at the highest 
   ! cutoff
	do j=1,m
		eps(n,j) = 1.0
		vel(n,j) = 1.0
		eps2(n,j) = 0.0
		wavefunc(n,j) = 0.0
	end do

	do i=1,n-1

		cutoff = 1.0*(n-i)/n ! Start from the highest cutoff 1 and gradually decrease it

		do j=1,m

			momentum = 1.0*j/m

			integral_vel = 0.0
			integral_eps = 0.0
			integral_eps2 = 0.0
			integral_wavefunc = 0.0

			do k=1,num

				phi = k*pi/(2.0*num)

				! Calculating cos(phi) to reduce number of function calls.
				! For the dielctric function integral limits are from 0 to pi/2 hence
				! phi/2 is used.

				cphi = COS(phi)

				! calculating velocities and dielectric functions to be used as the inputs
				! of dielectric function and velocity integrand respectively. Velocity and
				! dielectric function lying outside the highest cutoff is taken to be unity.
				index1 = NINT(cutoff*m)+1

				vel_1 = vel(n-i+1,index1)
				epsilon_1 = eps(n-i+1,index1)
				Y_1 = eps2(n-i+1,index1)

				index2 = NINT(m*(cutoff+momentum*cphi))+1

				if (index2<m) then
					vel_2 = vel(n-i+1,index2)
					epsilon_2 = eps(n-i+1,index2)
					Y_2 = eps2(n-i+1,index2)
				else
					vel_2 = 1.
					epsilon_2 = 1.
					Y_2 = 0.
				end if

				! Theta function implementation
					if (cphi>1.0 - 2.0*cutoff/momentum) then
						integral_eps = integral_eps + dielectric(phi,momentum,cutoff,vel_1,vel_2)
						integral_vel = integral_vel + velocity(cphi,momentum,cutoff,epsilon_1,epsilon_2, &
						vel_1,vel_2,0.0,0.0)
						integral_eps2 = integral_eps2 + dielectric2(phi,momentum,cutoff,vel_1,vel_2)
						integral_wavefunc = integral_wavefunc + wavefunc_renorm(cphi,momentum,cutoff, &
							epsilon_1,epsilon_2,vel_1,vel_2,Y_1,Y_2)
					end if

			end do
			vel(n-i,j) = vel(n-i+1,j) + dcutoff*dphi*integral_vel
			eps(n-i,j) = eps(n-i+1,j) + dcutoff*dphi*integral_eps
			eps2(n-i,j) = eps2(n-i+1,j) + dcutoff*dphi*integral_eps2
			wavefunc(n-i,j) = wavefunc(n-i+1,j) + dcutoff*dphi*integral_wavefunc
		end do
	end do

	open(1,file='epscon0.dat')
	open(2,file="velcon0.dat")
	open(3,file="epsc2on0.dat")
	open(4,file="wavefunc.dat")

	! do i=1,n
		do j=1,m
			write(1,*) 1.0*j/m, eps(2,j)
			write(2,*) 1.0*j/m, vel(2,j)
			write(3,*) 1.0*j/m, eps2(2,j)
			write(4,*) 1.0*j/m, wavefunc(2,j)
		end do
		! write(1,*)
		! write(2,*)
	! end do
	close(1)
	close(2)
	close(3)
	close(4)

	STOP

	contains

	function velocity(cphi,momentum,cutoff,epsilon_1,epsilon_2,vel_1,vel_2,Y_1,Y_2) result(value)
	! Returns the integrand for the velocity renormalisation
		real :: momentum,cutoff,value,cphi,pi,epsilon_1,epsilon_2
		real :: vel_1,vel_2,Y_1,Y_2,k,m1,m2

		pi = 3.14159

		m1 = cutoff
		m2 = cutoff + momentum*cphi
		k = momentum
		
		value = 2.2*((k**2 - m1**2 + m2**2)/(epsilon_1 + m2*vel_2*SQRT(epsilon_1*Y_1)) - &
		(-k**2-m1**2 + m2**2)/(epsilon_2 + m1*vel_1*SQRT(epsilon_2*Y_2))) &
		/(2.*pi*k**2*SQRT((m1+m2)**2 - k**2))

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

	function dielectric2(phi,momentum,cutoff,vel_1,vel_2) result(value)

		real :: phi,momentum,cutoff,value,cphi,sphi,pi,vel_1,vel_2
		real,parameter :: highcutoff = .2

		cphi = COS(phi)
		sphi = SIN(phi)
		pi = 3.14159

		value = 4.4*momentum*sphi**2/(pi*(cutoff*vel_1+(cutoff+momentum*cphi)*vel_2)**3* &
			SQRT((2*cutoff + momentum*cphi)**2 - momentum**2))		
	end function dielectric2

	function wavefunc_renorm(cphi,momentum,cutoff,epsilon_1,epsilon_2,vel_1,vel_2,Y_1,Y_2) result(value)
		
		real :: momentum,cutoff,value,cphi,pi,epsilon_1,epsilon_2
		real :: vel_1,vel_2,Y_1,Y_2,k,m1,m2,sqrtep1,sqrtep2,sqrty1,sqrty2

		pi = 3.14159

		m1 = cutoff
		m2 = cutoff + momentum*cphi
		k = momentum
		sqrtep1 = SQRT(epsilon_1)
		sqrtep2 = SQRT(epsilon_2)
		sqrty1 = SQRT(Y_1)
		sqrty2 = SQRT(Y_2)

		value = 2.2*((m2*sqrty1)/(sqrtep1*(sqrtep1 + m2*vel_2*sqrty1)**2) + &
			(m1*sqrty2)/(sqrtep2*(sqrtep2 + m1*vel_1*sqrty2)**2))/(pi*SQRT((m1+m2)**2 - k**2))

	end function wavefunc_renorm

end program pole