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
  end interface

double precision, intent(in) :: xin(3),vin(3)
double precision, intent(out) :: xout(3),vout(3)
double precision :: a, tau,r,vx
double precision :: T,n,mfp
! double precision :: phase_i(3), amplitude_i(3) !initial conditions now in module
!for the rkf45
integer :: neqn, flag,counter
double precision :: y, yp, time, tout, relerr, abserr
double precision :: yarr(3), yparr(3) !these are used in the numerical potential case
double precision ::  vr, vtheta,taustart
double precision :: ellvec(3) !angular momentum over m ( = r x v) and its magnitude
integer :: intcounter

time = 0.d0
tout = 0.d0
relerr = 1.d-4
abserr = 1.d-7
flag = 1
counter = 0

!optical depth that we travel before collision
!1) draw an optical depth
call random_number(a)
tau = -log(1.d0-a)
! tau = 5.
! print*,'WARNING Optical depth is hardcoded in; Walk.f90 L95'

! 2) Integrate the time it takes to reach that optical depth


!this is where we call RKF45, passin the step sub

!Replaced the external block with an interface
! Instead of integrating until we hit the right optical depth, we'll guess then adjust
!a little janky
! 1) Guess a time to get to tau. Hardcoded to H only

r  = sqrt(sum(xin**2))
vx = sqrt(sum(vin**2))

! n = ndensity(R,1)
!initial guess for how long we should go
! mfp = 1./(n*sigSD*vx) !a mean free path in seconds
! print*,"mean free path in cm ", mfp*vx, "and in s ", mfp
! tout = tau*mfp
! tout = tout/4. !this speeds things up
!tout is no longer used
!old iterative method:
! do while ((abs(y-tau)/tau .gt. 1.d-6 ) .and. (counter .lt. 10000))
! counter = counter + 1
! time = 0.d0
! flag = 1
! y = 0.d0 !initial column density
! call rkf45 (step, y, yp, time, tout, relerr, abserr, flag )
! ! print*, "flag is", flag
!  ! Newton's method to find the stopping point-- this does very well for constant density/small cross section,
!  ! not well at all if n(r) varies a lot over the trajectory
! ! print*," tau = ", tau, " y = ", y, ' tol ', abs(y-tau)/tau, ' tout ', tout
! tout = tout - (y-tau)/yp
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

  vout(1) = -amplitude_i(1)*OmegaSHO*sin(OmegaSHO*tout+phase_i(1))
  vout(2) = -amplitude_i(2)*OmegaSHO*sin(OmegaSHO*tout+phase_i(2))
  vout(3) = -amplitude_i(3)*OmegaSHO*sin(OmegaSHO*tout+phase_i(3))
  ! print*,'xout ', xout

else !numerically integrate potential
  ! Here we start to cheat : we are integrating the EoM for R; thanks to spherical
  !symmetry and conservation of angular momentum, the other coordinates do not matter.
  ! So we'll track velocities before and after propagation,
  !but we won't actually move the angular position. In this way we only need to solve 3 equations of motion:
  !t(tau), r(tau) and rdot(tau). Optical depth tau is still our dependent variable
  !tangential velocity recovered from angular momentum.


!angular momentum/m = R x V
  ellvec(3) = xin(1)*vin(2)-xin(2)*vin(1)
  ellvec(2) = xin(3)*vin(1)-xin(1)*vin(3)
  ellvec(1) = xin(2)*vin(3)-xin(3)*vin(2)

  taustart = 0.d0
!magnitude of angular momentum.
!This is stored in initial conditions module since it is required in the eom
  ell = sqrt(sum(ellvec**2))

  ! print*, "L = ", ellvec, " ell = ", ell
  yarr(1) = 0.d0 !initial time
  yarr(2) = r

  r  = sqrt(sum(xin**2))
  vx = sqrt(sum(vin**2))
  yarr(3) = sum(vin*xin)/r !velocity in R direction

  !energy over m is conserved
  !This is stored in initial conditions module since it is required in the eom
  ! eoverm = .5*vx**2 + .5*ell**2/r**2 + potential(r)
  eoverm = .5*vx**2  + potential(r) !can be used to track error
! print*, "E before propagation: ", eoverm
  ! print*,"callking rkf, eoverm = ", eoverm

intcounter = 0
  call rkf45full (pets_sph,3, yarr, yparr, taustart,tau, relerr, abserr, flag )
!arguments from the original function ( f, neqn, y, yp, t, tout, relerr, abserr, flag )
! print*,"called"
tout = yarr(1) !time

! print*,'time = ', tout, 'tau after integration: ', taustart, 'flag: ', flag
!If the integrator reaches a max number of steps it aborts with flag 4
do while (flag .eq. 4)
  intcounter = intcounter + 1
call rkf45full (pets_sph,3, yarr, yparr, taustart,tau, relerr, abserr, flag )
tout = yarr(1)
if (intcounter .eq. 1000) then
  print*,"You might be stuck in an infinite loop?"
end if
! print*,'time = ', tout, 'tau after integration: ', taustart, 'flag: ', flag
end do
!because we're not tracking angles here, we will just reset the particle position
!at x = r, other coordinates zero.
!We assign radial and tangential velocity in the xy plane, with vz = 0
xout(1) = yarr(2)
xout(2) = 0.d0
xout(3) = 0.d0





! !we have chosen r to be parallel to x
vout(1) = yarr(3) !sqrt(2.*eoverm - ell**2/xout(1)**2 - 2.*potential(r)) !get radial velocity from position and conservation of energy
! print*, "ell, ", ell, ", r ", xout(1), "vout: ", ell/xout(1)
vout(2) = ell/xout(1) ! stick tangential velocity in the y direction
vout(3) = 0.d0
! print*, "walked"
! eoverm = .5*(sum(vout**2)) + potential(sqrt(sum(xout**2)))
! print*, "E after propagation: ", eoverm

end if

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

!this reprocesses
! subroutine reprocess(Nsteps,ratio,infile)
!   use star
!   implicit none
!   integer, intent(in):: Nsteps, ratio
!   character*100 :: infile
!   double precision :: x1(Nsteps),x2(Nsteps),x3(Nsteps),r(Nsteps)
!
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
  call omega(x,vx,phase_i,amplitude_i,yprime)

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
  call omega(x,vx,phase_i,amplitude_i,yprime)
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


open(92,file = "intvals.dat",access='append')
write(92,*)tau, time, r, potential(r), vr, eoverm, ell,.5*(vr**2+ell**2/r**2)
close(92)
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



subroutine omega(xin,vin,phase_i,amplitude_i,omega_out) !,omegaprime)
  !compute omega and its derivative given the position and velocity of a particle
  !only scattering with a single species (hydrogen)
  use star
  implicit none
  double precision, intent(in) :: xin(3),vin(3),phase_i(3),amplitude_i(3)
  double precision :: vT,r,v2,y,omega_out,omegaprime,yprime,accel(3),wprefactor
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
  vT = sqrt(2.*kB*temperature(r)/mdm)
  v2 = vr**2 + (Lin/r)**2 !v_theta = L/r
  y = sqrt(v2/mu/vT**2)

! print*, "Omega: R ", r, " vT ", vT, " v ", sqrt(v2), " y ", y

  !niso = 1 = hydrogen hardcoded
  wprefactor = 2.*sigSD*ndensity(r,1)*vT*sqrt(mu)
  omega_out = wprefactor*((y+.5/y)*erf(y)+1./sqrt(pi)*exp(-y**2))

end subroutine omegaofR



!probably could be somewhere else
subroutine isinside(x,isinside_flag)
  use star
  implicit none
  double precision, intent(in) :: x(3)
  logical :: isinside_flag
  ! print*, "R/Rsun ", sqrt(sum(x**2))/Rsun
  if (sqrt(sum(x**2)) .ge. rsun) then
    isinside_flag = .false.
    ! print*,"flagged"
  end if
end subroutine

! function spawnCDF(r,r0)
!   double precision :: r,r0, CDF,N
!   N = sqrt(pi)*r0**3/4.
!   CDF =
!   return CDF
! end function spawnCDF
