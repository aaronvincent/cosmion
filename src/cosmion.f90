program cosmion
use star
implicit none

double precision :: xi(3),vi(3),x(3),v(3),vout(3),r,time !the DM coordinates
double precision, parameter :: GeV = 1.78266d-24
logical isinside_flag

character*100 :: outfile
! logical antemp, andens, anpot
! double precision mdm
integer Nsteps, i

!set parameters


call random_seed
!set up the star


!masses in g
mdm = 2.*GeV
sigsd = 1.d-37 !cm^2
anTemp = .false.
anDens = .false.
anPot = .false.
isinside_flag = .true.

Nsteps =1e5

outfile = 'positions.dat'

call init_star(anTemp,anDens,anPot,mdm,sigSD)

open(94,file = outfile)


! call spawn(x,y,z,vx,vy,vz)
call spawn(xi,vi)
print*,xi, vi

! big loop
call timestamp
do i = 1,Nsteps
call propagate(xi,vi,x,v,time)

!this check doesn't actually work, since the RKF solver will keep trying to integrate past rsun
!it dies without closing the file, I think we lose everything
call isinside(x,isinside_flag)
if (isinside_flag .eqv. .false.) then
print*, "Elvis has left the building"
return
end if
call collide(x,v,vout)
! print*, "vin ", v, "vout ", vout
write(94,*) x(1),x(2),x(3), v(1),v(2),v(3), vout(1),vout(2),vout(3),time
xi = x
vi = vout





! print*,x,v
end do
close(94)
print*,"Simulation complete"
call timestamp
! call propagate()
! call collide()
!
! end do



end program
