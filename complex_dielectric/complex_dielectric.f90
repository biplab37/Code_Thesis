!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                                   !!!
!!!   A program which solves the intrego-differential coupled equations given in the  !!!
!!!   paper by Carsten Bauer, Andreas Rückriegel, Anand Sharma, and Peter Kopietz     !!!
!!!   DOI: 10.1103/PhysRevB.92.121409. Here we have kept the bosonic momenta inside   !!!
!!!   the mometum shell to be integrated out.                                         !!!
!!!                                                                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program pole

   real :: phi, momentum, cutoff, dcutoff, dphi
   integer, parameter :: n = 101, m = 197, num = 800, fr = 50
   real :: pi, cphi, integral_vel, integral_eps, eps(m), vel(m)
   real :: vel_1, vel_2, int_eps_im, omega, eps_imag(fr, m)
   integer :: index1, index2

   pi = 3.14159

   dcutoff = 1.0/n
   dphi = pi/(2.*num)

   ! omega = 1.0

   ! open(2,file="velcon0w.dat")
   ! open(3,file="epsc2on0.dat")
   ! open(4,file="wavefunc.dat")
   ! open(5,file='self_realw.dat')
   ! open(6,file='self_imagw.dat')
   ! Initialise velocity and dielectric function to unity at all momenta at the highest
   ! cutoff
   do l = 1, fr + 1
      omega = 1.0*(l - 1)/fr
      do j = 1, m
         eps(j) = 1.0
         vel(j) = 1.0
         eps_imag(l, j) = 0.0
         ! eps2(j) = 0.0
         ! wavefunc(j) = 0.0
         ! selfreal(l,j) = 1.0
         ! selfimag(l,j) = 0.0
      end do

      do i = 1, n - 1

         cutoff = 1.0*(n - i)/n ! Start from the highest cutoff 1 and gradually decrease it

         do j = 1, m

            momentum = 1.0*j/m

            integral_vel = 0.0
            integral_eps = 0.0
            int_eps_im = 0.0
            ! integral_eps2 = 0.0
            ! integral_wavefunc = 0.0
            ! integral_real = 0.0
            ! integral_imag = 0.0

            do k = 1, num

               phi = k*pi/(2.0*num)

               ! Calculating cos(phi) to reduce number of function calls.
               ! For the dielctric function integral limits are from 0 to pi/2 hence
               ! phi/2 is used.

               cphi = COS(phi)

               ! calculating velocities and dielectric functions to be used as the inputs
               ! of dielectric function and velocity integrand respectively. Velocity and
               ! dielectric function lying outside the highest cutoff is taken to be unity.
               index1 = NINT(cutoff*m) + 1

               vel_1 = vel(index1)
               epsilon_1 = eps(index1)
               ! Y_1 = eps2(index1)

               index2 = NINT(m*(cutoff + momentum*cphi)) + 1

               if (index2 < m) then
                  vel_2 = vel(index2)
                  epsilon_2 = eps(index2)
                  ! Y_2 = eps2(index2)
               else
                  vel_2 = 1.
                  epsilon_2 = 1.
                  ! Y_2 = 0.
               end if

               ! Theta function implementation
               if (cphi > 1.0 - 2.0*cutoff/momentum) then
                  integral_eps = integral_eps + dielectric(phi, momentum, cutoff, vel_1, vel_2)
                  integral_vel = integral_vel + velocity(cphi, momentum, cutoff, epsilon_1, epsilon_2)
                  int_eps_im = int_eps_im + dielectric_imag(phi, momentum, cutoff, vel_1, vel_2, omega)
               end if

            end do
            vel(j) = vel(j) + dcutoff*dphi*integral_vel
            eps(j) = eps(j) + dcutoff*dphi*integral_eps
            eps_imag(l, j) = eps_imag(l, j) + dcutoff*dphi*int_eps_im
            ! eps2(j) = eps2(j) + dcutoff*dphi*integral_eps2
            ! wavefunc(j) = wavefunc(j) + dcutoff*dphi*integral_wavefunc
            ! selfreal(l,j) = selfreal(l,j) + dcutoff*dphi*integral_real
            ! selfimag(l,j) = selfimag(l,j) + dcutoff*dphi*integral_imag
         end do
      end do
      ! write(2,*) 1.0*(l-1)/fr, vel(10)
      ! write(3,*) 1.0*j/m, eps2(j)
      ! write(4,*) 1.0*j/m, wavefunc(j)
      ! write(5,*) 1.0*(l-1)/fr, selfreal(l,10)*10/m
      ! write(6,*) 1.0*(l-1)/fr, selfimag(l,10)*(l-1)/fr
! end do

      ! open(1,file='epscon0.dat')
      open (1, file='vel.dat')
      open (2, file='eps.dat')
      open (3, file='eps_imag.dat')

      ! do i=1,n
      ! do l=1,fr+1
      do j = 1, m
         ! write(1,*) 1.0*j/m, vel(j)
         ! write(2,*) 1.0*j/m, eps(j)
         write (3, *) 1.0*(l - 1)/fr, 1.0*j/m, eps_imag(l, j)

      end do
      ! write(1,*)
      ! write(2,*)
      write (3, *)
      ! write(5,*)
      ! write(6,*)
   end do
   ! end do
   close (1)
   close (2)
   close (3)
   ! close(4)
   ! close(5)
   ! close(6)

   STOP

contains

   function velocity(cphi, momentum, cutoff, epsilon_1, epsilon_2) result(value)
      ! Returns the integrand for the velocity renormalisation
      real :: momentum, cutoff, value, cphi, pi, epsilon_1, epsilon_2
      real :: k, m1, m2

      pi = 3.14159

      m1 = cutoff
      m2 = cutoff + momentum*cphi
      k = momentum

      value = 2.2*((k**2 - m1**2 + m2**2)/(epsilon_1) - (-k**2 - m1**2 + m2**2) &
                   /(epsilon_2))/(2.*pi*k**2*SQRT((m1 + m2)**2 - k**2))

   end function velocity

   function dielectric(phi, momentum, cutoff, vel_1, vel_2) result(value)
      ! Returns the integrand for dielectric function renormalisation

      real :: phi, momentum, cutoff, value, cphi, sphi, pi, vel_1, vel_2
      real, parameter :: highcutoff = .2

      cphi = COS(phi)
      sphi = SIN(phi)
      pi = 3.14159

      value = 4.4*momentum*sphi**2/(pi*(cutoff*vel_1 + (cutoff + momentum*cphi)*vel_2)* &
                                    SQRT((2*cutoff + momentum*cphi)**2 - momentum**2))

   end function dielectric

   function dielectric_imag(phi, momentum, cutoff, vel_1, vel_2, omega) result(value)

      real :: phi, momentum, cutoff, value, cphi, sphi, pi, vel_1, vel_2
      real :: m1, m2, k, omega
      real, parameter :: del = 0.001

      sphi = SIN(phi)
      pi = 3.14159
      m1 = cutoff
      m2 = cutoff + momentum*COS(phi)
      k = momentum

      value = 2.2*k*del**2*sphi**2/(pi*(del**2 + (omega - m1*vel_1 - m2*vel_2)**2)*m1*m2 &
                                    *SQRT((m1 + m2)**2 - k**2))
   end function dielectric_imag

   ! function wavefunc_renorm(cphi,momentum,cutoff,epsilon_1,epsilon_2,vel_1,vel_2,Y_1,Y_2) result(value)

   !         real :: momentum,cutoff,value,cphi,pi,epsilon_1,epsilon_2
   !         real :: vel_1,vel_2,Y_1,Y_2,k,m1,m2,sqrtep1,sqrtep2,sqrty1,sqrty2

   !         pi = 3.14159

   !         m1 = cutoff
   !         m2 = cutoff + momentum*cphi
   !         k = momentum
   !         sqrtep1 = SQRT(epsilon_1)
   !         sqrtep2 = SQRT(epsilon_2)
   !         sqrty1 = SQRT(Y_1)
   !         sqrty2 = SQRT(Y_2)

   !         value = 2.2*((m2*sqrty1)/(sqrtep1*(sqrtep1 + m2*vel_2*sqrty1)**2) + &
   !                 (m1*sqrty2)/(sqrtep2*(sqrtep2 + m1*vel_1*sqrty2)**2))/(pi*SQRT((m1+m2)**2 - k**2))

   ! end function wavefunc_renorm

   ! function selfenergy_real(cphi,momentum,cutoff,epsilon_1,epsilon_2,vel_1,vel_2,Y_1,Y_2,omega) result(value)

   !         real :: momentum,cutoff,value,cphi,pi,epsilon_1,epsilon_2,omega
   !         real :: vel_1,vel_2,Y_1,Y_2,k,m1,m2,sqrtep1,sqrtep2,sqrty1,sqrty2

   !         pi = 3.14159

   !         m1 = cutoff
   !         m2 = cutoff + momentum*cphi
   !         k = momentum
   !         sqrtep1 = SQRT(epsilon_1)
   !         sqrtep2 = SQRT(epsilon_2)
   !         sqrty1 = SQRT(Y_1)
   !         sqrty2 = SQRT(Y_2)

   !         value = 2.2*(((sqrtep1 + m2*vel_2*sqrty1)*(k**2 - m1**2 + m2**2))/(sqrtep1*(epsilon_1 + &
   !                 2*m2*vel_2*sqrtep1*sqrty1 + Y_1*(omega**2 + m2**2*vel_2**2))) + ((sqrtep2 + m1*vel_1*sqrty2)* &
   !                 (k**2 + m1**2 - m2**2))/(sqrtep2*(epsilon_2 + 2*m1*vel_1*sqrtep2*sqrty2 + Y_2*(omega**2 + &
   !                 m1**2*vel_1**2))))/(2.0*pi*k**2*SQRT((m1+m2)**2-k**2))

   ! end function selfenergy_real

   ! function selfenergy_imag(cphi,momentum,cutoff,epsilon_1,epsilon_2,vel_1,vel_2,Y_1,Y_2,omega) result(value)

   !         real :: momentum,cutoff,value,cphi,pi,epsilon_1,epsilon_2,omega
   !         real :: vel_1,vel_2,Y_1,Y_2,k,m1,m2,sqrtep1,sqrtep2,sqrty1,sqrty2

   !         pi = 3.14159

   !         m1 = cutoff
   !         m2 = cutoff + momentum*cphi
   !         k = momentum
   !         sqrtep1 = SQRT(epsilon_1)
   !         sqrtep2 = SQRT(epsilon_2)
   !         sqrty1 = SQRT(Y_1)
   !         sqrty2 = SQRT(Y_2)

   !         value = 2.2*((m2*sqrty1)/(sqrtep1*(epsilon_1 + 2*m2*vel_2*sqrtep1*sqrty1 + Y_1*(omega**2 + &
   !                 m2**2*vel_2**2))) + (m1*sqrty2)/(sqrtep2*(epsilon_2 + 2*m1*vel_1*sqrtep2*sqrty2 + &
   !                 Y_2*(omega**2 + m1**2*vel_1**2))))/(pi*SQRT((m1 + m2)**2 - k**2))

   ! end function selfenergy_imag

end program pole
