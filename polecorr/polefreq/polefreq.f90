!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                                   !!!
!!!   A program which solves the intrego-differential coupled equations given in the  !!!
!!!   paper by Carsten Bauer, Andreas Rückriegel, Anand Sharma, and Peter Kopietz     !!!
!!!   DOI: 10.1103/PhysRevB.92.121409. Here we have kept the bosonic momenta inside   !!!
!!!   the mometum shell to be integrated out.                                         !!!
!!!                                                                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program pole

	real :: phi,momentum,cutoff,dcutoff,dphi
	integer,parameter :: n=101,m=97,num=500,fr=20
	real :: pi,cphi,integral_vel,integral_eps,eps(m),vel(m),eps2(m),wavefunc(m)
	real :: selfreal(fr,m),selfimag(fr,m),integral_real,integral_imag,integral_wavefunc
	real :: epsilon_1,epsilon_2,vel_1,vel_2,integral_eps2,Y_1,Y_2,omega
	integer :: index1,index2

	pi = 3.14159

	dcutoff = 1.0/n
	dphi = pi/(2.*num)

	open(2,file="velcon0w.dat")
	! open(3,file="epsc2on0.dat")
	! open(4,file="wavefunc.dat")
	open(5,file='self_realw.dat')
	open(6,file='self_imagw.dat')
   ! Initialise velocity and dielectric function to unity at all momenta at the highest 
   ! cutoff
do l=1,fr+1
	omega = 1.0*(l-1)/fr
	do j=1,m
		eps(j) = 1.0
		vel(j) = 1.0
		eps2(j) = 0.0
		wavefunc(j) = 0.0
		selfreal(l,j) = 1.0
		selfimag(l,j) = 0.0
	end do

	do i=1,n-1

		cutoff = 1.0*(n-i)/n ! Start from the highest cutoff 1 and gradually decrease it

		do j=1,m

			momentum = 1.0*j/m

			integral_vel = 0.0
			integral_eps = 0.0
			integral_eps2 = 0.0
			integral_wavefunc = 0.0
			integral_real = 0.0
			integral_imag = 0.0

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

				vel_1 = vel(index1)
				epsilon_1 = eps(index1)
				Y_1 = eps2(index1)

				index2 = NINT(m*(cutoff+momentum*cphi))+1

				if (index2<m) then
					vel_2 = vel(index2)
					epsilon_2 = eps(index2)
					Y_2 = eps2(index2)
				else
					vel_2 = 1.
					epsilon_2 = 1.
					Y_2 = 0.
				end if

				! Theta function implementation
					if (cphi>1.0 - 2.0*cutoff/momentum) then
						integral_eps = integral_eps + dielectric(phi,momentum,cutoff,vel_1,vel_2)
						integral_vel = integral_vel + velocity(cphi,momentum,cutoff,epsilon_1,epsilon_2, &
						vel_1,vel_2,Y_1,Y_2)
						integral_eps2 = integral_eps2 + dielectric2(phi,momentum,cutoff,vel_1,vel_2)
						integral_wavefunc = integral_wavefunc + wavefunc_renorm(cphi,momentum,cutoff, &
							epsilon_1,epsilon_2,vel_1,vel_2,Y_1,Y_2)
						integral_real = integral_real + selfenergy_real(cphi,momentum,cutoff,epsilon_1, &
							epsilon_2,vel_1,vel_2,Y_1,Y_2,omega)
						integral_imag = integral_imag + selfenergy_imag(cphi,momentum,cutoff,epsilon_1, &
							epsilon_2,vel_1,vel_2,Y_1,Y_2,omega)
					end if

			end do
			vel(j) = vel(j) + dcutoff*dphi*integral_vel/(1.0 - wavefunc(j))
			eps(j) = eps(j) + dcutoff*dphi*integral_eps
			eps2(j) = eps2(j) + dcutoff*dphi*integral_eps2
			wavefunc(j) = wavefunc(j) + dcutoff*dphi*integral_wavefunc
			selfreal(l,j) = selfreal(l,j) + dcutoff*dphi*integral_real
			selfimag(l,j) = selfimag(l,j) + dcutoff*dphi*integral_imag
		end do
	end do
			write(2,*) 1.0*(l-1)/fr, vel(10)
			! write(3,*) 1.0*j/m, eps2(j)
			! write(4,*) 1.0*j/m, wavefunc(j)
			write(5,*) 1.0*(l-1)/fr, selfreal(l,10)
			write(6,*) 1.0*(l-1)/fr, selfimag(l,10)
end do

	! open(1,file='epscon0.dat')

	! do i=1,n
	! do l=1,fr+1
		! do j=1,m
			! write(1,*) 1.0*j/m, eps(j)
		! end do
		! write(1,*)
		! write(2,*)
		! write(5,*) 
		! write(6,*)
	! end do
	! end do
	! close(1)
	close(2)
	! close(3)
	! close(4)
	close(5)
	close(6)

	STOP

	contains

	function velocity(cphi,momentum,cutoff,epsilon_1,epsilon_2,vel_1,vel_2,Y_1,Y_2) result(value)
	! Returns the integrand for the velocity renormalisation
		real :: momentum,cutoff,value,cphi,pi,epsilon_1,epsilon_2
		real :: vel_1,vel_2,Y_1,Y_2,k,m1,m2 ,first,second

		pi = 3.14159

		m1 = cutoff
		m2 = cutoff + momentum*cphi
		k = momentum
		if ( Y_1.eq.0.0 .and. epsilon_1.eq.1.0 ) then
			first = (k**2 - m1**2 + m2**2)
		else
			first = (k**2 - m1**2 + m2**2)*(1 + (epsilon_1 - 1.0)*SQRT(epsilon_1 - 1.0)/&
				(epsilon_1*SQRT(epsilon_1 - 1.0) + m2*vel_2*SQRT(epsilon_1*Y_1)))
		end if

		if ( Y_2.eq.0.0 .and. epsilon_2.eq.1.0 ) then
			second = (k**2 - m2**2 + m1**2)
		else
			second = (k**2 - m2**2 + m1**2)*(1 - (epsilon_2 - 1.0)*SQRT(epsilon_2 - 1.0)/&
				(epsilon_2*SQRT(epsilon_2 - 1.0) + m1*vel_1*SQRT(epsilon_2*Y_2)))
		end if

		value = 2.2*(first + second)/(2.0*pi*k**2*SQRT((m1+m2)**2 - k**2))

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
		real :: vel_1,vel_2,Y_1,Y_2,k,m1,m2,first,second

		pi = 3.14159

		m1 = cutoff
		m2 = cutoff + momentum*cphi
		k = momentum

		if ( Y_1.eq.0.0 .and. epsilon_1.eq.1.0 ) then
			first = 0.0
		else
			first = -(m2*(epsilon_1 - 1.0)*(-(epsilon_1 - 1.0)*SQRT(Y_1) + m2*vel_2*Y_1*&
				SQRT((epsilon_1 - 1.0)/epsilon_1)))/(SQRT(epsilon_1*(epsilon_1 - 1.0)) + &
				m2*vel_2*SQRT(Y_1))
		end if

		if ( Y_2.eq.0.0 .and. epsilon_2.eq.1.0 ) then
			second = 0.0
		else
			second = -(m1*(epsilon_2 - 1.0)*(-(epsilon_2 - 1.0)*SQRT(Y_2) + m1*vel_1*Y_2*&
				SQRT((epsilon_2 - 1.0)/epsilon_2)))/(SQRT(epsilon_2*(epsilon_2 - 1.0)) &
				+ m1*vel_1*SQRT(Y_2))
		end if

		value = 2.2*(first + second)/(pi*k**2*SQRT((m1+m2)**2 - k**2))

	end function wavefunc_renorm

	function selfenergy_real(cphi,momentum,cutoff,epsilon_1,epsilon_2,vel_1,vel_2,Y_1,Y_2,omega) result(value)
		
		real :: momentum,cutoff,value,cphi,pi,epsilon_1,epsilon_2,omega
		real :: vel_1,vel_2,Y_1,Y_2,k,m1,m2

		pi = 3.14159

		m1 = cutoff
		m2 = cutoff + momentum*cphi
		k = momentum

		if ( Y_1.eq.0.0 .and. epsilon_1.eq.1.0 ) then
			first = 0.0
		else
			first = (k**2-m1**2+m2**2)*((epsilon_1 - 1.0)**2 + m2*vel_2*SQRT((epsilon_1-1.0)*Y_1/epsilon_1)) &
				/(epsilon_1*(epsilon_1-1.0) + Y_1*(omega**2 + m2**2*vel_2**2) + 2.0*m2*vel_2*SQRT(Y_1*epsilon_1 &
				*(epsilon_1 - 1.0)))
		end if

		if ( Y_2.eq.0.0 .and. epsilon_2.eq.1.0 ) then
			second = 0.0
		else
			second = (k**2-m2**2+m1**2)*((epsilon_2 - 1.0)**2 + m1*vel_1*SQRT((epsilon_2-1.0)*Y_2/epsilon_2)) &
				/(epsilon_2*(epsilon_2-1.0) + Y_2*(omega**2 + m1**2*vel_1**2) + 2.0*m1*vel_1*SQRT(Y_2*epsilon_2 &
				*(epsilon_2 - 1.0)))
		end if

		value = 2.2*(first + second)/(2.0*pi*k**2*SQRT((m1+m2)**2 - k**2))

	end function selfenergy_real

	function selfenergy_imag(cphi,momentum,cutoff,epsilon_1,epsilon_2,vel_1,vel_2,Y_1,Y_2,omega) result(value)
		
		real :: momentum,cutoff,value,cphi,pi,epsilon_1,epsilon_2,omega
		real :: vel_1,vel_2,Y_1,Y_2,k,m1,m2

		pi = 3.14159

		m1 = cutoff
		m2 = cutoff + momentum*cphi
		k = momentum

		if ( Y_1.eq.0.0 .and. epsilon_1.eq.1.0 ) then
			first = 0.0
		else
			first = (m2*(epsilon_1 - 1.0)*SQRT((epsilon_1 - 1.0)*Y_1/epsilon_1))/(epsilon_1*(epsilon_1-1.0) &
				+ Y_1*(omega**2 + m2**2*vel_2**2) + 2.0*m2*vel_2*SQRT(Y_1*epsilon_1*(epsilon_1 - 1.0)))
		end if

		if ( Y_2.eq.0.0 .and. epsilon_2.eq.1.0 ) then
			second = 0.0
		else
			second = (m1*(epsilon_2 - 1.0)*SQRT((epsilon_2 - 1.0)*Y_2/epsilon_2))/(epsilon_2*(epsilon_2-1.0) + &
				Y_2*(omega**2 + m1**2*vel_1**2) + 2.0*m1*vel_1*SQRT(Y_2*epsilon_2*(epsilon_2 - 1.0)))
		end if

		value = 2.2*omega*(first + second)/(pi*SQRT((m1+m2)**2 - k**2))

	end function selfenergy_imag

end program pole