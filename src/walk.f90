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
end interface

double precision, intent(in) :: xin(3),vin(3)
double precision, intent(out) :: xout(3),vout(3)
double precision :: a,tau,r,vx
double precision :: T,n,mfp,dt
double precision :: phase_r,amplitude_r,cosine
! double precision :: phase_i(3), amplitude_i(3) !initial conditions now in module
!for the rkf45
integer :: neqn, flag,counter
double precision :: y, yp, time, tout, relerr, abserr
integer :: intcounter

time = 0.d0
tout = 0.d0
relerr = 1.d-10
abserr = 1.d-10
flag = 1
counter = 0

!optical depth that we travel before collision
!1) draw an optical depth
call random_number(a)
tau = -log(1.d0-a)

! 2) Integrate the time it takes to reach that optical depth
!first get the initial conditions for this step
phase_i = atan(-vin,OmegaSHO*xin)
amplitude_i = xin/cos(phase_i)

!this is where we call RKF45, passin the step sub

!Replaced the external block with an interface
! Instead of integrating until we hit the right optical depth, we'll guess then adjust
!a little janky
! 1) Guess a time to get to tau. Hardcoded to H only

r  = sqrt(sum(xin**2))
vx = sqrt(sum(vin**2))

n = ndensity(R,1)
!initial guess for how long we should go
mfp = 1./(n*sigSD*vx) !a mean free path in seconds
! print*,"mean free path in cm ", mfp*vx, "and in s ", mfp
tout = tau*mfp
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
call rkf45 (pets, time, yp, y,tau, relerr, abserr, flag )
tout = time

do while (flag .eq. 4)
  print*,"Too many steps for integration."
  intcounter = intcounter + 1
  call rkf45 (pets, time, yp, y,tau, relerr, abserr, flag )
  tout = time
  if (intcounter .eq. 1000) then
    print*,"You might be stuck in an infinite loop?"
  end if
  ! print*,'time = ', tout, 'tau after integration: ', taustart, 'flag: ', flag
end do

xout(1) =  amplitude_i(1)*cos(OmegaSHO*tout+phase_i(1))
xout(2) =  amplitude_i(2)*cos(OmegaSHO*tout+phase_i(2))
xout(3) =  amplitude_i(3)*cos(OmegaSHO*tout+phase_i(3))
r = sqrt(sum(xout**2))
! This checks if the particle left the star during the integration.
if (outside_flag==1 .or. r>Rsun) then
  ! print*,"Elvis has left the building"
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

  ! Loop and keep checking if it's at the boundary yet.
  ! This worked, but we found a more efficient analytic solution.
  !vr = sum(xin*vin) / r
  !if (vr > 0.) then
  !  going_away = .true.
  !else
  !  going_away = .false.
  !end if
  !! This moves the particle to the edge of the star.
  !dt = 10.d0
  !do while (r<Rsun*0.99999)
  !  tout = tout+dt
  !  xout(1) =  amplitude_i(1)*cos(OmegaSHO*tout+phase_i(1))
  !  xout(2) =  amplitude_i(2)*cos(OmegaSHO*tout+phase_i(2))
  !  xout(3) =  amplitude_i(3)*cos(OmegaSHO*tout+phase_i(3))
  !  vout(1) = -amplitude_i(1)*OmegaSHO*sin(OmegaSHO*tout+phase_i(1))
  !  vout(2) = -amplitude_i(2)*OmegaSHO*sin(OmegaSHO*tout+phase_i(2))
  !  vout(3) = -amplitude_i(3)*OmegaSHO*sin(OmegaSHO*tout+phase_i(3))
  !  r = sqrt(sum(xout**2))
  !  vr = sum(xout*vout) / r
  !  if (r>Rsun*0.9999901 .or. vr<0 .and. going_away) then
  !  !if (r>Rsun*0.9999901) then
  !    tout = tout-dt
  !    xout(1) =  amplitude_i(1)*cos(OmegaSHO*tout+phase_i(1))
  !    xout(2) =  amplitude_i(2)*cos(OmegaSHO*tout+phase_i(2))
  !    xout(3) =  amplitude_i(3)*cos(OmegaSHO*tout+phase_i(3))
  !    r = sqrt(sum(xout**2))
  !    dt = dt/10.
  !  end if
  !end do

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

end subroutine propagate



subroutine collide(x,vin,vout)
  !turns old velocity into new velocity
  use star
  implicit none
  interface
    subroutine omega(xin,vin,omega_out,niso)
      double precision, intent(in) :: xin(3),vin(3)
      double precision, intent(out) :: omega_out
      integer, optional :: niso
    end subroutine
  end interface
  integer niso
  double precision, intent(in) :: x(3),vin(3)
  double precision, intent(out) :: vout(3)
  double precision :: v(3),vnuc(3),unuc,s(3),T,r,vcm,a,b
  double precision :: ctheta, phi,random_normal !outgoing angles in COM frame
  double precision :: tot,omegas(29),ratio(29),cumsum(29)
  integer :: i
!) select a species to collide with
  niso = 1
  if (.not. spinDep) then
    ! Randomly select what species to collide with based on their interaction rates.
    do i=1,29
      call omega(x,vin,omegas(i),i)
    end do
    tot = sum(omegas)
    ratio = omegas / tot
    cumsum = [(sum(ratio(1:i)), i=1, size(ratio))]
  
    call random_number(a)
    a = a * cumsum(29)
    do while (a>cumsum(niso))
      niso = niso + 1
    end do
    !print*,"element:",niso
    !print*,"probability:",ratio(niso)
    !print*,"radius:",r
  end if
  
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
species = niso


end subroutine collide

  !this goes into the RK solver
  !y in this subroutine is the function we're integrating (i.e. the optical depth, tau)
subroutine step(t,y,yprime)
  use init_conds
  use star
  implicit none
  interface
    subroutine omega(xin,vin,omega_out,niso)
      double precision, intent(in) :: xin(3),vin(3)
      double precision, intent(out) :: omega_out
      integer, optional :: niso
    end subroutine
  end interface
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
  interface
    subroutine omega(xin,vin,omega_out,niso)
      double precision, intent(in) :: xin(3),vin(3)
      double precision, intent(out) :: omega_out
      integer, optional :: niso
    end subroutine
  end interface
  double precision, intent(in) :: t,y
  double precision, intent(out) ::  yprimeinv
  double precision :: ri(3), vi(3),yprime
  double precision :: x(3),vx(3) !the integrator does not need to know these
  x = amplitude_i*cos(OmegaSHO*t+phase_i)
  vx = -amplitude_i*OmegaSHO*sin(OmegaSHO*t+phase_i)
  ! print*,'calling step, x = ', x
!this needs to be a loop if you have multiple species
  !y is not used
  call omega(x,vx,yprime)
  yprimeinv = 1.d0/yprime

end subroutine pets


subroutine omega(xin,vin,omega_out,niso) !,omegaprime)
  !compute omega and its derivative given the position and velocity of a particle
  !only scattering with a single species (hydrogen)
  use star
  implicit none
  double precision, intent(in) :: xin(3),vin(3)
  double precision, intent(out) :: omega_out
  double precision :: vT,r,v2,y,omegaprime,yprime,accel(3),wprefactor,sigma
  integer, optional :: niso
  integer :: i
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
  if (present(niso)) then
    y = sqrt(v2*AtomicNumber(niso)/mu/vT**2)
    sigma = sigSD * AtomicNumber(niso)**4 * ((mdm+mnucg)/(mdm+AtomicNumber(niso)*mnucg))**2
    wprefactor = 2.*sigma*ndensity(r,niso)*vT*sqrt(mu/AtomicNumber(niso))
    omega_out = wprefactor*((y+.5/y)*erf(y)+1./sqrt(pi)*exp(-y**2))
  else
    if (spinDep) then
      y = sqrt(v2/mu/vT**2)

      ! print*, "Omega: R ", r, " vT ", vT, " v ", sqrt(v2), " y ", y

      !niso = 1 = hydrogen hardcoded
      wprefactor = 2.*sigSD*ndensity(r,1)*vT*sqrt(mu)
      omega_out = wprefactor*((y+.5/y)*erf(y)+1./sqrt(pi)*exp(-y**2))
      !the next bit is me not understanding wtf is going on. I made derivatives yay. Ignore
      ! accel = -OmegaSHO**2*xin
      !
      ! yprime = 2.d0*sum(accel*vin)/mdm !y^2'
      ! yprime = .5/y
      ! omegaprime = dndr(r,1)/ndensity(r,1)*omega_out + yprime*wprefactor*(erf(y)*(1.-1./y**2) + exp(-y**2)/sqrt(pi)/y)
      ! print*,'y ', y, 'yprime ', yprime
    else
      omega_out = 0.d0
      do i=1,29
        y = sqrt(v2*AtomicNumber(i)/mu/vT**2)
        sigma = sigSD * AtomicNumber(i)**4 * ((mdm+mnucg)/(mdm+AtomicNumber(i)*mnucg))**2
        wprefactor = 2.*sigma*ndensity(r,i)*vT*sqrt(mu/AtomicNumber(i))
        omega_out = omega_out + wprefactor*((y+.5/y)*erf(y)+1./sqrt(pi)*exp(-y**2))
      end do
    end if
  end if
end subroutine omega

subroutine cross(x,y,z)
  ! Compute the 3-component vector product of x cross y.
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
  interface
    subroutine cross(x,y,z)
      double precision, intent(in) :: x(3),y(3)
      double precision, intent(out) :: z(3)
    end subroutine
    subroutine keplerian_rad(xin,vin,xout,vout,tout)
      double precision, intent(in) :: xin(3),vin(3)
      double precision, intent(out) :: xout(3),vout(3),tout
    end subroutine
  end interface
  double precision, intent(in) :: xin(3),vin(3)
  double precision, intent(out) :: xout(3),vout(3),tout
  double precision :: r,vr,vtot,vesc,Mstar
  double precision :: h(3),smaj,e,theta,norm(3),plvec(3),ang
  double precision :: area,areatot,period,c
  ! These are for Hannah's method
  !double precision :: e1(3),e2(3),bigtheta,anomaly,timediff,vt,theta0,newh,newe,M

  ! Here we must include code for the Keplerian orbit of the particle.
  ! Using the initial position and velocity, we find the re-entry parameters.
  Mstar = 4./3.*pi*Rsun**3*rhoSHO
  r = sqrt(sum(xin**2))
  vr = sum(xin*vin) / r
  vtot = sqrt(sum(vin**2))
  vesc = sqrt(2.*GN*Mstar/Rsun)
  ! Check if the particle exceeds the escape velocity.
  if (vtot >= vesc) then
    print*,"The particle has escaped!"
    ! We'll have to stop and respawn the particle in this case.
    outside_flag = 2
  else if (vr/vtot > 1.d0-1.d-10) then
    ! Use the Keplerian function for a radial emission.
    call keplerian_rad(xin,vin,xout,vout,tout)
  else
    !print*,"It's coming back"
    ! Do Keplerian stuff
    smaj = GN*Mstar*r / (2.*GN*Mstar-vtot**2*r)
    call cross(xin,vin,h)
    e = sqrt(1.-sum(h**2)/(smaj*GN*Mstar))
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
    period = sqrt(4.*pi**2*smaj**3/(GN*Mstar))
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
  interface
    subroutine cross(x,y,z)
      double precision, intent(in) :: x(3),y(3)
      double precision, intent(out) :: z(3)
    end subroutine
  end interface
  double precision, intent(in) :: xin(3),vin(3)
  double precision, intent(out) :: xout(3),vout(3),tout
  double precision :: r,vtot,vesc,vr(3),Mstar
  double precision :: h(3),smaj,e,theta
  double precision :: area,areatot,period,c
  ! Here we must include code for the Keplerian orbit of the particle.
  ! Using the initial position and velocity, we find the re-entry parameters.
  ! In this version, the particle re-enters at the same position.
  Mstar = 4./3.*pi*Rsun**3*rhoSHO
  r = sqrt(sum(xin**2))
  vtot = sqrt(sum(vin**2))
  vesc = sqrt(2.*GN*Mstar/Rsun)
  ! Check if the particle exceeds the escape velocity
  if (vtot >= vesc) then
    print*,"The particle has escaped!"
    ! We'll have to stop and respawn the particle in this case.
    outside_flag = 2
  else
    !print*,"It's coming back"
    ! Do Keplerian stuff
    xout = xin
    vr = sum(xin*vin)/r**2 * xin
    vout = vin-2.*vr
    smaj = GN*Mstar*r / (2.*GN*Mstar-vtot**2*r)
    if (sqrt(sum(vr**2))/vtot > 1.d0-1.d-10) then
      c = 2.*smaj
      tout = 2.*sqrt(c**3/(2.*GN*Mstar))*(sqrt(r/c*(1.-r/c))+acos(sqrt(r/c)))
    else
      call cross(xin,vin,h)
      e = sqrt(1.-sum(h**2)/(smaj*GN*Mstar))
      theta = acos((smaj-e**2*smaj-r)/(e*r))
      period = sqrt(4.*pi**2*smaj**3/(GN*Mstar))
      areatot = pi*smaj**2*sqrt(1.-e**2)
      c = 2.*atanh((e-1.)*tan(theta/2.)/sqrt(cmplx(e**2-1.)))/sqrt(cmplx(e**2-1.))
      area = smaj**2*(e**2-1.)*(e*sin(theta)/(e*cos(theta)+1.)-c)
      tout = (1.-area/areatot)*period
    end if
  end if
end subroutine

! function spawnCDF(r,r0)
!   double precision :: r,r0, CDF,N
!   N = sqrt(pi)*r0**3/4.
!   CDF =
!   return CDF
! end function spawnCDF
