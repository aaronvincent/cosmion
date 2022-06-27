! This is where the action lives
! Units are cgs, except mass
! Follow the recipe: https://arxiv.org/abs/2111.06895


subroutine spawn(x,v)
  use star
  implicit none
  double precision :: x(3),v(3),T
  logical accept
  double precision :: r, N,pdf,a,b,ctheta,phi
  double precision random_normal

  ! print*,"Spawining a walker"

  !Basic rejection method to get r
  accept = .false.
  do while (accept .eqv. .false.)
  call random_number(a)
  a = a*Rsun
  pdf = a**2*exp(-(a/rchi)**2)
  call random_number(b)
  b= b*rchi**2*exp(-1.d0)
  ! print*," Rsun ", rsun, " rchi ", rchi,  " a ", a, ' b ', b
  if (b .lt. pdf) then
    accept = .true.
  end if
  end do
  r = a
  ! get theta, recycle a and b
  call random_number(a)
  ctheta = 2.*a-1.
  a = acos(ctheta)
  call random_number(b)
  phi = 2.*pi*b
  x(1) = r*sin(a)*sin(phi)
  x(2) = r*sin(a)*cos(phi)
  x(3) = r*ctheta
  T = temperature(r)

  ! print*,"Spawining velocities"
  v(1) = Sqrt(kB*T/mdm)*random_normal()
  v(2) = Sqrt(kB*T/mdm)*random_normal()
  v(3) = Sqrt(kB*T/mdm)*random_normal()

  ! Throws the particle out of the star in one step, for testing
  !x = (/37d9,37d9,37d9/)
  !v = (/0d7,54d7,32d7/)

! print*,'A random number ', random_normal()

end subroutine spawn

subroutine propagate(xin,vin,xout,vout,tout)
use star
use init_conds
implicit none
interface
  subroutine step(t,y,yp)
    !should be rk kind but whatever
    double precision, intent(in) :: t,y
    double precision, intent(out) :: yp
  end subroutine
  subroutine pets(y,t,yp)
    !should be rk kind but whatever
    double precision, intent(in) :: t,y
    double precision, intent(out) :: yp
  end subroutine
  subroutine pets_sph(t,y,yp)     !integrand for numerical potential
    double precision, intent(in) :: t,y(3)
    double precision, intent(out) :: yp(3)
  end subroutine
  subroutine pets_to_surf(t,y,yp) !integrand for propagation to surface without collision
    double precision, intent(in) :: t,y
    double precision, intent(out) :: yp
  end subroutine
  subroutine pets_to_surf2d(t,y,yp) !integrand for propagation to surface without collision
    double precision, intent(in) :: t,y(3)
    double precision, intent(out) :: yp(3)
  end subroutine
  ! subroutine pets_to_turning(t,y,yp) !integrand for propagation to surface without collision
  !   double precision, intent(in) :: t,y(3)
  !   double precision, intent(out) :: yp(3)
  ! end subroutine
  function turnaroundEnergy(rin) result(eout)
  ! used to find turnaround for inward-bound orbits
    implicit none
    double precision, intent(in) :: rin
    double precision eout
   end function

  function turnaroundEnergyPrime(rin) result(eout)
    ! used to find turnaround for inward-bound orbits
    implicit none
    double precision, intent(in) :: rin
    double precision eout
  end function
end interface

double precision, intent(in) :: xin(3),vin(3)
double precision, intent(out) :: xout(3),vout(3)
double precision :: a, tau,r,vx
double precision :: T,n,mfp,tvec(3000),yarrout(3000,3),dt


double precision :: phase_r,amplitude_r,cosine
! double precision :: phase_i(3), amplitude_i(3) !initial conditions now in module
!for the rkf45
integer :: neqn, flag,counter,fcounter
double precision :: y, yp, time, tout, relerr, abserr
double precision :: yarr(3), yparr(3) !these are used in the numerical potential case
double precision ::  vr, vtheta,taustart,Rbar
double precision :: ellvec(3) !angular momentum over m ( = r x v) and its magnitude
integer :: intcounter,iters

time = 0.d0
tout = 0.d0
relerr = 1.d-5
abserr = 1.d-10
flag = 1
counter = 0
!used for tracking position in integrator
tvec(:) = tvec(:)*0.d0
yarrout(:,:)= yarrout(:,:)*0.d0

!optical depth that we travel before collision
!1) draw an optical depth
call random_number(a)
tau = -log(1.d0-a)
! tau = 5.
! print*,'WARNING Optical depth is hardcoded in; Walk.f90 L95'

! 2) Integrate the time it takes to reach that optical depth


r  = sqrt(sum(xin**2))
vx = sqrt(sum(vin**2))

!
! ! I'm gonna do what you might call a pro-gamer move
! !this integrates t until we reach optical depth tau
!
! ! print*,'time ', tout
! end do
! print*," tau = ", tau, " y = ", y, ' tol ', abs(y-tau)/tau, ' tout ', time
time = 0.d0
flag = 1
y = 0.d0

if (anPot) then
  !first get the initial conditions for this step
  phase_i = atan(-vin,OmegaSHO*xin)
  amplitude_i = xin/cos(phase_i)

  call rkf45 (pets, time, yp, y,tau, relerr, abserr, flag )
  tout = time

  do while (flag .eq. 4)
    intcounter = intcounter + 1
  call rkf45 (pets, time, yp, y,tau, relerr, abserr, flag )
  tout = time
  if (intcounter .eq. 1000) then
    print*,"You might be stuck in an infinite loop?"
  end if
  ! print*,'time = ', tout, 'tau after integration: ', taustart, 'flag: ', flag
  end do

  ! print*, "took ", counter, ' tries ', ' guess ', tau*mfp/4., 'actual ', tout
  ! print*,"using pro move: "
  ! print*," tau = ", tau, " y = ", y, ' tol ', abs(y-tau)/tau, ' tout ', t
  ! print*," Did the loop, tout = ", tout


  xout(1) =  amplitude_i(1)*cos(OmegaSHO*tout+phase_i(1))
  xout(2) =  amplitude_i(2)*cos(OmegaSHO*tout+phase_i(2))
  xout(3) =  amplitude_i(3)*cos(OmegaSHO*tout+phase_i(3))

  ! vout(1) = -amplitude_i(1)*OmegaSHO*sin(OmegaSHO*tout+phase_i(1))
  ! vout(2) = -amplitude_i(2)*OmegaSHO*sin(OmegaSHO*tout+phase_i(2))
  ! vout(3) = -amplitude_i(3)*OmegaSHO*sin(OmegaSHO*tout+phase_i(3))
  ! print*,'xout ', xout


  r = sqrt(sum(xout**2))

  ! This checks if the particle left the star during the integration.
  if (outside_flag .ne. 0 .or. r>Rsun) then

    if (sqrt(sum(vout**2)) .ge. vescape(r)) then
      outside_flag =2
      print*, "particle has evaporated"
      return
    end if

    ! We now have to determine the path from the particle's initial position to the solar radius.

    ! The Keplerian orbit places the particle at the boundary of the star.
    ! So, if the particle is already outside the star, this moves it inside the boudary.
    ! This should only take one step, which is a negligible time that is later ignored.
    tout = 0.d0
    dt = 2.d-4 ! Determines precision
    xout = xin
    r = sqrt(sum(xout**2))
    do while (r/Rsun>=0.99999d0)
      tout = tout+dt
      xout(1) =  amplitude_i(1)*cos(OmegaSHO*tout+phase_i(1))
      xout(2) =  amplitude_i(2)*cos(OmegaSHO*tout+phase_i(2))
      xout(3) =  amplitude_i(3)*cos(OmegaSHO*tout+phase_i(3))
      r = sqrt(sum(xout**2))
    end do
    vout(1) = -amplitude_i(1)*OmegaSHO*sin(OmegaSHO*tout+phase_i(1))
    vout(2) = -amplitude_i(2)*OmegaSHO*sin(OmegaSHO*tout+phase_i(2))
    vout(3) = -amplitude_i(3)*OmegaSHO*sin(OmegaSHO*tout+phase_i(3))
    if (tout/dt > 1.) then
      print*,"Multiple steps for particle to re-enter star:",tout/dt
    end if
    phase_i = atan(-vout,OmegaSHO*xout)
    amplitude_i = xout/cos(phase_i)



    ! This determines how long it takes the particle to reach the boundary of the star.
    phase_r = atan(-(sum(xout*vout))/(omegaSHO*(sum(xout**2)-(sum(amplitude_i**2)/2))))
    amplitude_r = (sum(xout**2)-sum(amplitude_i**2)/2)/cos(phase_r)
    cosine = acos(((Rsun*0.99999005d0)**2-sum(amplitude_i**2)/2)/amplitude_r)
    if (phase_r>0 .and. cosine<abs(phase_r)) then
      cosine = 2*pi - cosine
    else if (phase_r<0 .and. cosine<abs(phase_r)) then
      cosine = -cosine
    end if
    tout = (cosine-phase_r)/(2.*OmegaSHO)
  end if
  xout(1) =  amplitude_i(1)*cos(OmegaSHO*tout+phase_i(1))
  xout(2) =  amplitude_i(2)*cos(OmegaSHO*tout+phase_i(2))
  xout(3) =  amplitude_i(3)*cos(OmegaSHO*tout+phase_i(3))

  vout(1) = -amplitude_i(1)*OmegaSHO*sin(OmegaSHO*tout+phase_i(1))
  vout(2) = -amplitude_i(2)*OmegaSHO*sin(OmegaSHO*tout+phase_i(2))
  vout(3) = -amplitude_i(3)*OmegaSHO*sin(OmegaSHO*tout+phase_i(3))

else ! not anPot: numerically integrate potential
  ! Here we are integrating the EoM for R; thanks to spherical
  !symmetry and conservation of angular momentum, the other coordinates do not matter.
  ! So we'll track velocities before and after propagation,
  !but we won't actually move the angular position. In this way we only need to solve 3 equations of motion:
  !t(tau), r(tau) and rdot(tau). Optical depth tau is still our dependent variable
  !tangential velocity recovered from angular momentum.


!angular momentum/m = R x V
call cross(xin,vin,ellvec)

  ! ellvec(3) = xin(1)*vin(2)-xin(2)*vin(1)
  ! ellvec(2) = xin(3)*vin(1)-xin(1)*vin(3)
  ! ellvec(1) = xin(2)*vin(3)-xin(3)*vin(2)



  taustart = 0.d0
!magnitude of angular momentum.
!This is stored in initial conditions module since it is required in the eom
  ell = sqrt(sum(ellvec**2))

  ! print*, "L = ", ellvec, " ell = ", ell
  yarr(1) = 0.d0 !initial time
  yarr(2) = r

  r  = sqrt(sum(xin**2))
  vx = sqrt(sum(vin**2))
  vr = sum(vin*xin)/r
  yarr(3) = vr !velocity in R direction


  !energy over m is conserved
  !This is stored in initial conditions module since it can be required in the eom
   eoverm = .5*vx**2  + potential(r)

intcounter = 0
if (fullHistory) then
  print*,"fullHistory Doesn't work yet!!"
  stop
!   ! call rkf45full (pets_sph,3, yarr, yparr, taustart,tau, relerr, abserr, flag )
!   call rkf45fullhistory (pets_sph,3, yarr, yparr, taustart,tau, relerr, abserr, flag,tvec,yarrout )
! !arguments from the original function ( f, neqn, y, yp, t, tout, relerr, abserr, flag )
! ! print*,"called"
! tout = yarr(1) !time
! print*,"tvec: ", tvec(1:10)
! print*,"yarr: ", yarrout(1:5,:)
!
! !write all computed trajectory positions. first column is tau, then time, r, vr
! ! open(99,file='rep_pos.dat',access='append')
! ! fcounter = 0
! ! do while (tvec(fcounter+1) .ne. 0.d0)
! !   fcounter = fcounter+1
! ! write(99,*) tvec(fcounter), yarrout(fcounter,1),yarrout(fcounter,2),yarrout(fcounter,3)
! ! end do
! ! close(99)
! ! tvec = tvec*0.d0
!
! ! print*,'time = ', tout, 'tau after integration: ', taustart, 'flag: ', flag
! !If the integrator reaches a max number of steps it aborts with flag 4
! do while (flag .eq. 4)
!   intcounter = intcounter + 1
! ! call rkf45full (pets_sph,3, yarr, yparr, taustart,tau, relerr, abserr, flag )
! call rkf45fullhistory (pets_sph,3, yarr, yparr, taustart,tau, relerr, abserr, flag,tvec,yarrout )
! tout = yarr(1)
!
! !!!!write trajectory positions again
! open(99,file='rep_pos.dat',access='append')
! fcounter = 0
! do while (tvec(fcounter+1) .ne. 0.d0)
!   fcounter = fcounter+1
! write(99,*) tvec(fcounter), yarrout(fcounter,:)
! end do
! close(99)
! tvec = tvec*0.d0
! if (intcounter .eq. 1000) then
!   print*,"You might be stuck in an infinite loop?"
! end if
! ! print*,'time = ', tout, 'tau after integration: ', taustart, 'flag: ', flag
! end do

else !not full history

!!! Main propagation
  ! call pets_sph(taustart,yarr,yparr)
  ! print*,"flag before rkf", flag
  call rkf45full (pets_sph,3, yarr, yparr, taustart,tau, relerr, abserr, flag )
  tout = yarr(1)
  do while (flag .eq. 4)
    intcounter = intcounter + 1
  ! call rkf45full (pets_sph,3, yarr, yparr, taustart,tau, relerr, abserr, flag )
  call rkf45full (pets_sph,3, yarr, yparr, taustart,tau, relerr, abserr, flag )
  tout = yarr(1)
  if (intcounter .eq. 1000) then
    print*,"You might be stuck in an infinite loop?"
    print*, "at time ", tout, 'tau', taustart,'going to ',tau,'yarr', yarr
    ! print*, yarr
  end if
  ! print*,'time = ', tout, 'tau after integration: ', taustart, 'flag: ', flag
  end do
!!! Main propagation done

  !if at any point the integrator realized it had left the star, we ditch any
  !work it did and figure out the keplerian bit
  if (outside_flag .ne. 0 .or. r>Rsun) then
    print*,"outside"
  if (vx .ge. vescape(r)) then
    outside_flag = 2
    print*, "Escaped"
    return
  end if !escape
  ! print*,'eoverm ', eoverm, 'r = ', r, 'ell = ', ell, 'v = ', vx, 'vesc ', vescape(r), 'potential ', potential(r)
  time = 0.d0
  !integrate to just under Rsun, otherwise we get problems
  ! call pets_to_surf(r,0.d0,yp)
  !1d integrator *should* work, but it goes berzerk. I don't get it
  ! call rkf45 (pets_to_surf,time, yp, r,Rsun*(1.-1.d-10), relerr, abserr, flag )
  vr = sum(vin*xin)/r

  ! call pets_to_surf2d(r,yarr,yparr)
 flag = 1 !reset integrator flag
  if (vr .lt. 0.d0) then
    !we need to integrate until vr changes sign, then integrate up to r
    !vr = 0 when we hit the angular momentum barrier
    !solve for rdot = 0 using energy conservation
    yarr(1) = time
    yarr(2) = vr
    yarr(3) = 0.d0 !along for the ride, does nothing
    call newton(turnaroundEnergy, turnaroundEnergyPrime, r/10., Rbar, iters, .false.)
    if (rbar .gt. Rsun) then
      print*, "Major problem, particle turnaround point is outside the star"
      stop
    end if
    ! print*,'Found reversal rbar = ',Rbar
    ! integrate vr from vr to zero
    call rkf45full (pets_to_surf2d,3,yarr, yparr, r,Rbar, relerr, abserr, flag )
    time = yarr(1)
    vr = yarr(2)
    ! print*,"time to rbar = ", time, " reached with velocity ",vr

    ! vr = 1.d-10! start it off positive
    !reverse course!
    ! yarr(2) = abs(yarr(2))
    flag = 1
    ! stop
  end if

  yarr(1) = time
  yarr(2) = dabs(vr)
  yarr(3) = 0.d0 !along for the ride, does nothing
  call pets_to_surf2d(Rbar,yarr,yparr)
  ! print*,"intergrating from ", Rbar, " to ", Rsun, "yarr ", yarr
  !place it just
  call rkf45full (pets_to_surf2d,3,yarr, yparr, Rbar,Rsun*(1.-1.d-4), relerr, abserr, flag )
  xout(1) = Rsun*(1.-1.d-4)
  xout(2) = 0.d0
  xout(3) = 0.d0
  ! print*,'result of rkf time,', yarr(1), 'r ', r, 'vesc at rsun ', vescape(Rsun),'flag = ',flag
  ! vout(1) = sqrt(2.*eoverm - ell**2/Rsun**2 - 2.*potential(Rsun)) !get radial velocity from position and conservation of energy
  time = yarr(1)
  vout(1) = yarr(2)
  ! print*, "eoverm", eoverm, "phi", potential(Rsun),"ell ", ell
  ! print*, "ell, ", ell, ", r ", xout(1), "vout: ", ell/xout(1)
  vout(2) = ell/dble(xout(1)) ! stick tangential velocity in the y direction
  vout(3) = 0.d0

  tout = time
  ! print*,'vr is now ', dble(vout(1)), 'Vtheta is ', vout(2)

else !outside flag
  xout(1) = yarr(2)
  xout(2) = 0.d0
  xout(3) = 0.d0

  vout(1) = yarr(3) !sqrt(2.*eoverm - ell**2/xout(1)**2 - 2.*potential(r)) !get radial velocity from position and conservation of energy
  ! print*, "ell, ", ell, ", r ", xout(1), "vout: ", ell/xout(1)
  vout(2) = ell/xout(1) ! stick tangential velocity in the y direction
  vout(3) = 0.d0
end if !outside flag != 0



end if !fullhistory
!because we're not tracking angles here, we will just reset the particle position
!at x = r, other coordinates zero.
!We assign radial and tangential velocity in the xy plane, with vz = 0


end if !numerical solution


end subroutine propagate



subroutine collide(x,vin,vout)
  !turns old velocity into new velocity
  use star
  implicit none
  integer niso
  double precision, intent(in) :: x(3),vin(3)
  double precision, intent(out) :: vout(3)
  double precision :: v(3),vnuc(3),unuc,s(3),T,r,vcm,a,b
  double precision :: ctheta, phi,random_normal !outgoing angles in COM frame
!) select a species to collide with. Hardcoded for now


  niso = 1
  !a little different from Hannah's method: we draw 3 nuclear velocities from a local MB distribution
  v = vin
  r = sqrt(sum(x**2))
  T = temperature(r)

  vnuc(1) = Sqrt(kB*T/(AtomicNumber(niso)*mnucg))*random_normal()
  vnuc(2) = Sqrt(kB*T/(AtomicNumber(niso)*mnucg))*random_normal()
  vnuc(3) = Sqrt(kB*T/(AtomicNumber(niso)*mnucg))*random_normal()

! 2) Boost to CM frame
s = (mdm*v + AtomicNumber(niso)*mnucg*vnuc)/(mdm+AtomicNumber(niso)*mnucg)
vcm = sqrt(sum((v-s)**2)) !velocity in CM frame

!Scattering is isotropic, so the new angle does not depend on the old one
call random_number(a)
ctheta = 2.*a-1.
a = acos(ctheta)
call random_number(b)
phi = 2.*pi*b
vout(1) = vcm*sin(a)*sin(phi)
vout(2) = vcm*sin(a)*cos(phi)
vout(3) = vcm*ctheta
!boost back to star frame

vout = vout + s



end subroutine collide

! this reprocesses the run to provide equal-time positions and velocities. Required to
!adequately represent the phase space distribution
!it doesn't work atm
! subroutine reprocess(Nsteps,times,ratio,infile,outfile)
!   use star
!   implicit none
!   interface
!     subroutine step_sph(t,y,yp)     !integrand for numerical potential
!       double precision, intent(in) :: t,y(2)
!       double precision, intent(out) :: yp(2)
!     end subroutine
!   end interface
!
!   integer, intent(in):: Nsteps, ratio
!   double precision, intent(in) :: times(Nsteps)
!   double precision, parameter :: relerr = 1.d-5, abserr = 1.d-8
!   character*100 :: infile,outfile
!   ! double precision :: x1(Nsteps),x2(Nsteps),x3(Nsteps),r(Nsteps)
!   double precision :: x(3), v(3), vout(3), pot, r, vx
!   double precision :: dt, time, tcoll, yarr(2),yparr(2),tvar
!   double precision :: ellvec(3), ell,vtheta,vr,tend,toffset, ttot
!   integer :: i, intcounter, flag
!
!   ! compute time steps
!   ttot = sum(times)
!   dt = ttot/dble(nsteps*ratio)
!   toffset = 0.d0
!   flag = 1
!
!   open(45,file = infile)
!   open(46,file = outfile)
!   do i = 1,Nsteps
!     print*,'reprocessing, collision number ', i
!      !ignoring potential, vout should do nothing
!     read(45,*)  x(1),x(2),x(3), v(1),v(2),v(3), vout(1),vout(2),vout(3),time
!
!     ellvec(3) = x(1)*v(2)-x(2)*v(1)
!     ellvec(2) = x(3)*v(1)-x(1)*v(3)
!     ellvec(1) = x(2)*v(3)-x(3)*v(2)
!     ell = sqrt(sum(ellvec**2))
!     r  = sqrt(sum(x**2))
!     vx = sqrt(sum(v**2))
!     vr = sum(v*x)/r
!
!
!     tvar = toffset
!     print*,'toffset = ', toffset
!     ! tend = 0;
!
!     !step through trajectory between collisions, integrating each time
!     do while (tvar+dt .lt. times(i+1))
!       !offset time to account for the difference between the end of the
!       !last step and subsequent collision
!       !computed before step is updated
!       toffset = times(i+1) - tvar
!       print*," time: ", tvar, "next time: ", times(i+1), "dt: ", dt
!       yarr(1) = r
!       yarr(2) = vr!velocity in R direction
!       tend = tend+dt !so many time variables
!       print*, "tvar ", tvar, "tend ", tend
!       call rkf45full_n2 (step_sph,2, yarr, yparr, tvar,tend, relerr, abserr, flag )
!       !catch incomplete integrations
!         do while (flag .eq. 4)
!           intcounter = intcounter + 1
!         call rkf45full_n2 (step_sph,2, yarr, yparr, tvar,tend, relerr, abserr, flag )
!
!         if (intcounter .eq. 1000) then
!           print*,"You might be stuck in an infinite loop?"
!         end if
!         end do
!         print*, "tvar ", tvar, "tend ", tend, "flag ", flag
!         tvar = tend
!         r = yarr(1)
!         vr = yarr(2)
!         vtheta = ell/r
!         !you've integrated one time step, write to file and move on to the next One
!         write(46,*) r, vr, vtheta, tvar
!
!     end do
!
!   end do
!
!   close(45)
!   close(46)
!
! end subroutine reprocess


  !this goes into the RK solver
  !y in this subroutine is the function we're integrating (i.e. the optical depth, tau)
subroutine step(t,y,yprime)
  use init_conds
  use star
  implicit none
  double precision, intent(in) :: t,y
  double precision, intent(out) ::  yprime
  double precision :: ri(3), vi(3)
  double precision :: x(3),vx(3) !the integrator does not need to know these
  x = amplitude_i*cos(OmegaSHO*t+phase_i)
  vx = -amplitude_i*OmegaSHO*sin(OmegaSHO*t+phase_i)
  ! print*,'calling step, x = ', x
!this needs to be a loop if you have multiple species
  !y is not used
  call omega(x,vx,yprime)

end subroutine step

!inverse of steps(). dt/dtau = yprimeinv = 1/w
subroutine pets(y,t,yprimeinv)
  use init_conds
  use star
  implicit none
  double precision, intent(in) :: t,y
  double precision, intent(out) ::  yprimeinv
  double precision :: ri(3), vi(3),yprime
  double precision :: x(3),vx(3) !the integrator does not need to know these
  x = amplitude_i*cos(OmegaSHO*t+phase_i)
  vx = -amplitude_i*OmegaSHO*sin(OmegaSHO*t+phase_i)


  ! open(92,file = "intvals_SHO.dat",access='append')
  ! write(92,*) y, t, x, vx
  ! close(92)
  ! print*,'calling step, x = ', x
!this needs to be a loop if you have multiple species
  !y is not used
  call omega(x,vx,yprime)
  yprimeinv = 1.d0/yprime

end subroutine pets

!inverse of steps(). dt/dtau = yprimeinv = 1/w
!pets, but integrating the r-dependent potential
!arranged the arguments so that they make sense
! y(1) is time, y(2) is r, y(3) is vr
subroutine pets_sph(tau,y,yprime)
  use init_conds
  use star
  implicit none
  double precision, intent(in) :: tau,y(3)
  double precision, intent(out) ::  yprime(3)
  double precision :: omega_i, time, r, vr,vdot
  ! double precision ::

  time = y(1)
  ! Zero crossing
  r = (y(2))
  vr = y(3)
  if (r .lt. 0.d0) then
    r = -r
    vr = -vr
  end if


  !!commented out attempt to integrate single equation only
  ! things become problematic when the particle turns around.
  !vdot is used to get the sign of v_r

  ! vdot = gofr(r) + ell**2/r**3
  ! print*, "g", gofr(r)
  !this is a fudge. Check later that energy is conserved.
!   if (2.*eoverm - ell**2/r**2 - 2*potential(r) .lt. 0.d0) then
!     vr = 0.d0
!   else
!   vr = sqrt(2.*eoverm - ell**2/r**2 - 2*potential(r))
! end if

  ! vr = -1.*vdot/abs(vdot)*vr


! open(92,file = "intvals.dat",access='append')
! write(92,*)tau, time, r, potential(r), vr, eoverm, ell,.5*(vr**2+ell**2/r**2)
! close(92)
! print*,'calling step, r = ', r, "vr = ", vr, "potential = ", potential(r), "eoverm = ", eoverm,"ell = ", ell
!this needs to be a loop if you have multiple species
  !y is not used
  ! print*,'calling Omega with ell = ', ell
  call omegaofR(r,vr,ell,omega_i)

  ! print*, "Called"
  yprime(1) = 1.d0/omega_i != dt/dtau

!the yprime(1) here is to go from d/dt to d/dtau using the chain rule
  yprime(2) = yprime(1)*vr !eom for r

  yprime(3) = yprime(1)*(gofr(r) + ell**2/r**3) !EOM for rdot

! print*,'Arrays assigned: tau = ', tau, ' y = ', y, 'yprime = ', yprime
  ! print*, "g of r ", gofr(r)

end subroutine pets_sph

!used to find time to surface when the particle leaves the star
  subroutine pets_to_surf(t,y,yprime)
    use init_conds
    use star
    implicit none
    double precision, intent(in) :: t,y
    double precision, intent(out) ::  yprime
    double precision :: time, r
    ! double precision ::

    time = y
    ! Zero crossing
    r = t
    yprime = 1./sqrt(abs(2.*eoverm-ell**2/r**2-2.*potential(r)))

    ! print*,"E is ", .5/yprime**2 + ell**2/r**2/2.d0+potential(r),' and should be ', eoverm
    ! open(92,file = "inching.dat",access='append')
    ! write(92,*) r/rsun, aux, potential(r),1./yprime,ell, eoverm,time
    ! close(92)
    ! print*, "r/rsun", aux2, 'aux ', aux
  end subroutine pets_to_surf

!the above subroutine should work, but the integrator goes berzerk when it takes backwards steps
  subroutine pets_to_surf2d(r,y,yprime)
    use init_conds
    use star
    implicit none
    double precision, intent(in) :: r,y(3)
    double precision, intent(out) ::  yprime(3)
    ! double precision ::  aux, aux2
    ! double precision ::
    !only 2 elements of the y vector are used. 3rd is along for the ride
    !y(1) is time, y(2) is v
    yprime(1) = 1./y(2) !dt/dr
    yprime(2) = 1./y(2)*(gofr(r) + ell**2/r**3) ! dv/dr = dt/dr dv/dt
    yprime(3) = 0.d0
    ! print*, 'r = ', r, 'time = ', y(1),'vr = ', y(2),'gofr,', gofr(r)
    ! ! print*,"E is ", .5/yprime**2 + ell**2/r**2/2.d0+potential(r),' and should be ', eoverm
    ! open(92,file = "inching.dat",access='append')
    ! write(92,*) r, potential(r),y(1),y(2),ell, eoverm
    ! close(92)
    ! print*, "r/rsun", aux2, 'aux ', aux
  end subroutine pets_to_surf2d

!This subroutine takes the radial velocity v as the independent variable, and integrates r and t
!to find the turning point. This is to track particles with v_tot > vesc, but are on an inward trajectory
! !this is stupid and doesn't work (v not monotonic) stage for deletion
!   subroutine pets_to_turning(v,y,yprime)
!     use init_conds
!     use star
!     implicit none
!     double precision, intent(in) :: v,y(3)
!     double precision, intent(out) ::  yprime(3)
!     double precision ::  aux, aux2, r
!     ! double precision ::
!     !only 2 elements of the y vector are used. 3rd is along for the ride
!     ! time = y(1)
!     r = y(2)
!     yprime(1) = 1.d0/(gofr(r) + ell**2/r**3) !dt/dv
!     yprime(2) = v/(gofr(r) + ell**2/r**3) ! dr/dv = dr/dt dt/dv. Minus sign necessary here
!     yprime(3) = 0.d0
!     print*, 'r = ', r, 'vr = ', v,'a = ', (gofr(r) + ell**2/r**3), "time", y(1)
!     ! print*,"E is ", .5/yprime**2 + ell**2/r**2/2.d0+potential(r),' and should be ', eoverm
!     open(92,file = "inching.dat",access='append')
!     write(92,*) r, potential(r),y(1),v,ell, eoverm,gofr(r)
!     close(92)
!     ! print*, "r/rsun", aux2, 'aux ', aux
!   end subroutine pets_to_turning

!This integrand is used for reprocessing. It integrates
!orbits in time (not optical depth)
subroutine step_sph(time,y,yprime)
  use init_conds
  use star
  implicit none
  double precision, intent(in) :: time, y(2)
  double precision, intent(out) ::  yprime(2)
  double precision :: omega_i, r, vr,vdot

  ! time
  ! Zero crossing
  r = y(1)
  vr = y(2)
  if (r .lt. 0.d0) then
    r = -r
    vr = -vr
  end if

  yprime(1) = vr != dr/dt
  yprime(2) = (gofr(r) + ell**2/r**3) !dv/dt

end subroutine step_sph



subroutine omega(xin,vin,omega_out) !,omegaprime)
  !compute omega and its derivative given the position and velocity of a particle
  !only scattering with a single species (hydrogen)
  use star
  implicit none
  double precision, intent(in) :: xin(3),vin(3)
  double precision :: vT,r,v2,y,omega_out,omegaprime,yprime,accel(3),wprefactor
  ! The following conditional checks if the particle is inside the star.
  ! If it left the star, it raises a flag to be detected after the integration is complete.
  if (sqrt(sum(xin**2)) .ge. Rsun) then
    !omega_out = 0.0
    !omega_out = omega_out/omega_out
    outside_flag = 1
    return
  end if
  r = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
  vT = sqrt(2.*kB*temperature(r)/mdm)
  v2 = vin(1)**2 + vin(2)**2 + vin(3)**2
  y = sqrt(v2/mu/vT**2)

! print*, "Omega: R ", r, " vT ", vT, " v ", sqrt(v2), " y ", y

  !niso = 1 = hydrogen hardcoded
  wprefactor = 2.*sigSD*ndensity(r,1)*vT*sqrt(mu)
  omega_out = wprefactor*((y+.5/y)*erf(y)+1./sqrt(pi)*exp(-y**2))
  !the next bit is me not understanding wtf is going on. I made derivatives yay. Ignore
!   accel = -OmegaSHO**2*xin
!
!   yprime = 2.d0*sum(accel*vin)/mdm !y^2'
!   yprime = .5/y
!   omegaprime = dndr(r,1)/ndensity(r,1)*omega_out + yprime*wprefactor*(erf(y)*(1.-1./y**2) + exp(-y**2)/sqrt(pi)/y)
! print*,'y ', y, 'yprime ', yprime
end subroutine omega


subroutine omegaofR(rin,vr,Lin,omega_out) !,omegaprime)
  !compute omega and its derivative given the position and velocity of a particle
  !different signature from omega used in the analytic potential version
  !since we are alread in spherical coordinates
  !input: rin (radial coord), vr (velocity in r direction), Lin (angular momentum divided by m)
  !only scattering with a single species (hydrogen)
  use star
  implicit none
  double precision, intent(in) :: rin,Lin,vr
  double precision :: vT,r,v2,y,omega_out,omegaprime,yprime,accel(3),wprefactor
  r = rin
  if (r .ge. Rsun) then
    outside_flag = 1
    return
  end if
  vT = sqrt(2.*kB*temperature(r)/mdm)
  v2 = vr**2 + (Lin/r)**2 !v_theta = L/r
  y = sqrt(v2/mu/vT**2)

! print*, "Omega: R ", r, " vT ", vT, " v ", sqrt(v2), " y ", y

  !niso = 1 = hydrogen hardcoded
  wprefactor = 2.*sigSD*ndensity(r,1)*vT*sqrt(mu)
  omega_out = wprefactor*((y+.5/y)*erf(y)+1./sqrt(pi)*exp(-y**2))

end subroutine omegaofR



!probably could be somewhere else

subroutine cross(x,y,z)
  implicit none
  double precision, intent(in) :: x(3),y(3)
  double precision, intent(out) :: z(3)
  z(1) = x(2)*y(3)-x(3)*y(2)
  z(2) = -x(1)*y(3)+x(3)*y(1)
  z(3) = x(1)*y(2)-x(2)*y(1)
end subroutine

subroutine keplerian(xin,vin,xout,vout,tout)
  use star
  implicit none
  ! interface
  !   subroutine cross(x,y,z)
  !     double precision, intent(in) :: x(3),y(3)
  !     double precision, intent(out) :: z(3)
  !   end subroutine
  !   subroutine keplerian_rad(xin,vin,xout,vout,tout)
  !     double precision, intent(in) :: xin(3),vin(3)
  !     double precision, intent(out) :: xout(3),vout(3),tout
  !   end subroutine
  ! end interface
  double precision, intent(in) :: xin(3),vin(3)
  double precision, intent(out) :: xout(3),vout(3),tout
  double precision :: r,vr,vtot,vesc
  double precision :: h(3),vh(3),smaj,e,theta,norm(3),plvec(3),ang
  double precision :: area,areatot,period,c
  ! These are for Hannah's method
  !double precision :: e1(3),e2(3),bigtheta,anomaly,timediff,vt,theta0,newh,newe,M

  ! Here we must include code for the Keplerian orbit of the particle.
  ! Using the initial position and velocity, we find the re-entry parameters.
  r = sqrt(sum(xin**2))
  vr = sum(xin*vin) / r
  vtot = sqrt(sum(vin**2))
  !check now done before this sub is called
  ! vesc = 2.*Rsun*sqrt(2.*pi*GN*rhoSHO/3.)
  ! Check if the particle exceeds the escape velocity.
  ! if (vtot >= vesc) then
  !   print*,"The particle has escaped!"
  !   ! We'll have to stop and respawn the particle in this case.
  !   outside_flag = 2
  if (vr/vtot > 1.d0-1.d-10) then
    ! Use the Keplerian function for a radial emission.
    call keplerian_rad(xin,vin,xout,vout,tout)
  else
    !print*,"It's coming back"
    ! Do Keplerian stuff
    smaj = 4.*pi*GN*Rsun**3*rhoSHO*r / (8.*pi*GN*Rsun**3*rhoSHO-3.*vtot**2*r)
    call cross(xin,vin,h)
    e = sqrt(1.-sum(h**2)/(smaj*4.*pi*GN*Rsun**3*rhoSHO/3.))
    theta = acos((smaj-e**2*smaj-r)/(e*r))
    call cross(xin,vin,norm)
    norm = norm/sqrt(sum(norm**2))
    call cross(norm,xin,plvec)
    xout = cos(2.*pi-theta*2.)*xin+sin(2.*pi-theta*2.)*plvec
    ang = acos(vr/vtot)
    call cross(norm,xout,plvec)
    if (cos(ang)>0) then
      vout = cos(pi-ang)*xout+sin(pi-ang)*plvec
    else
      vout = cos(ang)*xout+sin(ang)*plvec
    end if
    vout = (vout/r)*vtot
    period = sqrt(3.*pi*smaj**3/(GN*Rsun**3*rhoSHO))
    areatot = pi*smaj**2*sqrt(1-e**2)
    c = 2.*atanh((e-1.)*tan(theta/2.)/sqrt(cmplx(e**2-1.)))/sqrt(cmplx(e**2-1.))
    area = smaj**2*(e**2-1.)*(e*sin(theta)/(e*cos(theta)+1.)-c)
    tout = (1.-area/areatot)*period
    !print*,"a=",smaj
    !print*,"e=",e
    !print*,"theta=",theta
    !print*,"norm=",norm
    !print*,"xin=",xin
    !print*,"vin=",vin
    !print*,"xout=",xout
    !print*,"vout=",vout
    !print*,"area",area
    !print*,"areatot",areatot
    !print*,"tout=",tout
    !print*,"tfrac=",tout/period
    !print*,""

    ! The below is Hannah's method of finding the position and time of re-entry.
    ! It was used to verify that the above method is correct.
    !M = 4./3. * pi * Rsun**3 * rhoSHO
    !e1 = xin
    !e1 = e1 / sqrt(sum(e1**2))
    !e2 = vin - (sum(vin*e1))*e1
    !e2 = e2 / sqrt(sum(e2**2))
    !vr = sum(vin*e1)
    !vt = sum(vin*e2)
    !newh = Rsun*vt
    !theta0 = atan(vr/(GN*M/newh-vt))
    !if (theta0>0) then
    !  theta0 = theta0+pi
    !else if (theta0<0) then
    !  theta0 = theta0+2.*pi
    !end if
    !newe = -(newh*vr)/(GN*M*sin(theta0))
    !bigtheta = acos((newh**2-Rsun*GN*M)/(Rsun*GN*M*newe))+theta0
    !anomaly = 2.*atan(sqrt((1.-e)/(1.+e))*tan((bigtheta-theta0)/2.))
    !do while (abs(anomaly-(bigtheta-theta0))>pi .and. anomaly>(bigtheta-theta0))
    !  anomaly = anomaly - 2.*pi
    !end do
    !do while (abs(anomaly-(bigtheta-theta0))>pi .and. anomaly<(bigtheta-theta0))
    !  anomaly = anomaly + 2.*pi
    !end do
    !timediff = smaj*(smaj*sqrt(1.-e**2))*(anomaly-e*sin(anomaly))/newh
    !print*,"original outside time:",tout
    !print*,"new time difference:",timediff
    !print*,"original period:",period
    !print*,""
  end if
end subroutine

subroutine keplerian_rad(xin,vin,xout,vout,tout)
  use star
  implicit none
  ! interface
  !   subroutine cross(x,y,z)
  !     double precision, intent(in) :: x(3),y(3)
  !     double precision, intent(out) :: z(3)
  !   end subroutine
  ! end interface
  double precision, intent(in) :: xin(3),vin(3)
  double precision, intent(out) :: xout(3),vout(3),tout
  double precision :: r,vtot,vesc,vr(3)
  double precision :: h(3),vh(3),smaj,e,theta,norm(3),plvec(3),ang
  double precision :: area,areatot,period,c
  ! Here we must include code for the Keplerian orbit of the particle.
  ! Using the initial position and velocity, we find the re-entry parameters.
  ! In this version, the particle re-enters at the same position.
  r = sqrt(sum(xin**2))
  vtot = sqrt(sum(vin**2))
  ! this is now done in propagate, before walk is called.
  ! vesc = 2.*Rsun*sqrt(2.*pi*GN*rhoSHO/3.)
  ! Check if the particle exceeds the escape velocity
  ! if (vtot >= vesc) then
  !   print*,"The particle has escaped!"
  !   ! We'll have to stop and respawn the particle in this case.
  !   outside_flag = 2
  ! else
    !print*,"It's coming back"
    ! Do Keplerian stuff
    xout = xin
    vr = sum(xin*vin)/r**2 * xin
    vout = vin-2.*vr
    call cross(xin,vin,h)
    smaj = 4.*pi*GN*Rsun**3*rhoSHO*r / (8.*pi*GN*Rsun**3*rhoSHO-3.*vtot**2*r)
    e = sqrt(1.-sum(h**2)/(smaj*4.*pi*GN*Rsun**3*rhoSHO/3.))
    theta = acos((smaj-e**2*smaj-r)/(e*r))
    period = sqrt(3.*pi*smaj**3/(GN*Rsun**3*rhoSHO))
    areatot = pi*smaj**2*sqrt(1.-e**2)
    c = 2.*atanh((e-1.)*tan(theta/2.)/sqrt(cmplx(e**2-1.)))/sqrt(cmplx(e**2-1.))
    area = smaj**2*(e**2-1.)*(e*sin(theta)/(e*cos(theta)+1.)-c)
    tout = (1.-area/areatot)*period
  ! end if
end subroutine

function turnaroundEnergy(rin) result(eout)

  ! used to find turnaround for inward-bound orbits
  use star
  use init_conds
  implicit none
  double precision, intent(in) :: rin
  double precision eout,r
  r = dabs(rin) !avoid funny business
  eout = 2.d0*eoverm - ell**2/r**2 - 2.d0*potential(r)
end function

function turnaroundEnergyPrime(rin) result(eout)

  ! used to find turnaround for inward-bound orbits
  use star
  use init_conds
  implicit none
  double precision, intent(in) :: rin
  double precision eout
  double precision r
  r = dabs(rin)
  eout = 2.d0*ell**2/r**3 +2.d0*gofr(r)
end function
! function spawnCDF(r,r0)
!   double precision :: r,r0, CDF,N
!   N = sqrt(pi)*r0**3/4.
!   CDF =
!   return CDF
! end function spawnCDF
