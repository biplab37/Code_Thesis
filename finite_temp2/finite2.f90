!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                                   !!!
!!!   A program which solves the intego-differential coupled equations given in the   !!!
!!!   paper by Carsten Bauer, Andreas Rückriegel, Anand Sharma, and Peter Kopietz     !!!
!!!   DOI: 10.1103/PhysRevB.92.121409. Here we have added temperature dependence.     !!!
!!!                                                                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program finite2 

	real :: phi,momentum,cutoff,dcutoff,dphi
	integer,parameter :: n=501,m=297,num=2000
	real :: pi,cphi,integral_vel,integral_eps,eps(n,m),vel(n,m)
	real :: epsilon_1,epsilon_2,vel_1,vel_2
	integer :: index1,index2

	pi = 3.14159

	dcutoff = 1.0/n
	dphi = pi/(2.*num)
	temp = 0.10

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

				phi = k*pi/(2.0*(num))

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

				index2 = NINT(m*(cutoff+momentum*cphi))+1

				if (index2<m) then
					vel_2 = vel(n-i+1,index2)
					epsilon_2 = eps(n-i+1,index2)
				else
					vel_2 = 1.
					epsilon_2 = 1.
				end if
				! write(*,*) epsilon_1,epsilon_2
				! Theta function implementation
				! if (cphi<(1 -cutoff)/momentum) then
					if (cphi>1 - 2*cutoff/momentum) then
						integral_eps = integral_eps + dielectric_finite(phi,momentum,cutoff,vel_1,vel_2,temp)
						integral_vel = integral_vel + velocity_finite(cphi,momentum,cutoff,epsilon_1,epsilon_2,vel_1,vel_2,temp)
						if (integral_eps==NaN) then
							write(*,*) momentum,cutoff
						end if
					end if
				! end if

			end do
			vel(n-i,j) = vel(n-i+1,j) + dcutoff*dphi*integral_vel
			eps(n-i,j) = eps(n-i+1,j) + dcutoff*dphi*integral_eps
		end do
	end do

	! open(1,file='eps.dat')
	! open(2,file="vel.dat")

	! do i=1,n
	! 	do j=1,m
	! 		write(1,*) 1.0*i/n,1.0*j/m, eps(i,j)
	! 		write(2,*) 1.0*i/n,1.0*j/m, vel(i,j)
	! 	end do
	! 	write(1,*)
	! 	write(2,*)
	! end do
	! close(1)
	! close(2)

	open(3,file="epson0_0.10.dat")
	open(4,file="velon0_0.10.dat")

	write(3,*) 'This file contain the data of renormalised dielectric constant at finite temperature'
	write(4,*) 'This file contain the data of renormalised velocity at finite temperature'

	write(3,*) 'Momentum      ', 'Dielectric'
	write(4,*) 'Momentum      ', 'Velocity'

	do j=1,m
		write(3,*) 1.0*j/m, eps(2,j)
		write(4,*) 1.0*j/m, vel(2,j)
	end do

	close(3)
	close(4)

	STOP

	contains

	function velocity_finite(cphi,momentum,cutoff,epsilon_1,epsilon_2,vel_1,vel_2,temp) result(value)
	! Returns the integrand for the velocity_finite renormalisation
		real :: momentum,cutoff,value,cphi,pi,epsilon_1,epsilon_2,val1
		real :: vel_1,vel_2,temp

		pi = 3.14159

		val1 = 2*cutoff + momentum*cphi
		value = 2.2*((val1*cphi + k)*(TANH(vel_1*momentum/(2.*temp)))/(k*epsilon_1) & 
		- (val1*cphi - k)*(TANH(vel_2*(cutoff+momentum*cphi)/(2.*temp)))/(k*epsilon_2)) &
		/(2.*pi*SQRT(val1**2 - momentum**2))
		! write(*,*) value

	end function velocity_finite


	function velocity(cphi,momentum,cutoff,epsilon_1,epsilon_2,vel_1,vel_2,temp) result(value)
	! Returns the integrand for the velocity renormalisation
		real :: momentum,cutoff,value,cphi,pi,epsilon_1,epsilon_2,val1
		real :: vel_1,vel_2,temp

		pi = 3.14159

		val1 = 2*cutoff + momentum*cphi
		value = 2.2*((val1*cphi + k)/(k*epsilon_1) - (val1*cphi - k)/(k*epsilon_2))/(2.*pi*SQRT(val1**2 - momentum**2))
		! write(*,*) value

	end function velocity


	function dielectric_finite(phi,momentum,cutoff,vel_1,vel_2,temp) result(value)
	! Returns the integrand for dielectric function renormalisation

		real :: phi,momentum,cutoff,value,cphi,sphi,pi,vel_1,vel_2,temp
		real :: momentum1,momentum2,momsum,exp1,exp2,expA,expB,expC

		cphi = COS(phi)
		! sphi = SIN(phi)
		pi = 3.14159

		momentum1 = cutoff
		momentum2 = cutoff + momentum*cphi
		tanh1 = TANH((vel_1*momentum1)/(2.0*temp))
		tanh2 = TANH((vel_2*momentum2)/(2.0*temp))
		denomin1 = SQRT((momentum1+momentum2)**2-momentum**2)
		! expA = vel_1*(momentum1**2 + momentum2**2 - momentum**2) + 2*vel_2*momentum2**2
		! expB = 2*vel_1*momentum1**2 + vel_2*(momentum1**2 - momentum**2 + momentum2**2)
		! expC = (vel_2**2)*(momentum2**2) - (vel_1**2)*(momentum1**2)
		if (abs((momentum1*vel_1 - momentum2*vel_2))<0.0001) then
			value = 4.4*(-denomin1+ 2*momentum1*momentum2*(tanh1 + tanh2)/denomin1)/ &
					(pi*momentum*(momentum1*vel_1 + momentum2*vel_2))
		else	
		value = 4.4*(denomin1*(momentum2*vel_2*tanh1 - momentum1*vel_1*tanh2)/(momentum1*vel_1 - momentum2*vel_2) +&
				2*momentum1*momentum2*(tanh1 + tanh2)/denomin1)/(pi*momentum*(momentum1*vel_1 + momentum2*vel_2))
		end if

	end function dielectric_finite

	function dielectric(phi,momentum,cutoff,vel_1,vel_2,temp) result(value)
	! Returns the integrand for dielectric function renormalisation

		real :: phi,momentum,cutoff,value,cphi,sphi,pi,vel_1,vel_2
		real,parameter :: highcutoff = .2

		cphi = COS(phi)
		sphi = SIN(phi)
		pi = 3.14159

		value = 4.4*momentum*sphi**2/(pi*(cutoff*vel_1+(cutoff+momentum*cphi)*vel_2)* &
			SQRT((2*cutoff + momentum*cphi)**2 - momentum**2))

	end function dielectric

end program finite2 