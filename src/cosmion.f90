program cosmion
use star
implicit none

double precision :: xi(3),vi(3),x(3),v(3),vout(3),r,time !the DM coordinates
double precision, allocatable :: times(:) !makes reprocessing a bit easier to keep this in memory
double precision, parameter :: GeV = 1.78266d-24
logical isinside_flag

character*100 :: outfile, reprofile
! logical antemp, andens, anpot
! double precision mdm
integer Nsteps, i,ratio



outfile = 'positions.dat'
reprofile = 'rep_pos.dat' !only used if fullHistory = true

!set parameters


call random_seed
!set up the star


!masses in g
mdm = 10.d0*GeV
sigsd = 1.d-36 !cm^2
Nsteps =1d6


!anXXX flags: if false, interpolate from a stellar model. If true, use analytic
!functions defined in star.f90
!For testing only.
anTemp = .false.
anDens = .false.
!treat potential as SHO and use analytic trajectory expressions for x, v
anPot = .false.
!turn this on if you want the full trajectory history, not just timestamps at every collision
!Note this will take a ludicrous amount of HD space.
fullHistory = .false.

!SHO_debug overrides the tabulated phi(r) to provide SHO potential
!for testing of phi(r) and comparison with anPot. Don't set to true for realistic sims
!this flag does nothing if anPot is already true.
SHO_debug = .true.

if (anPot .or. SHO_debug) then
  print*, "Watch out, you are using a SHO potential"
end if


isinside_flag = .true.



allocate(times(Nsteps))


call init_star(anTemp,anDens,anPot,mdm,sigSD)

! print*,"potential at 0.5 ", potential(0.5d0)

!wipe the trajectory history
open(99,file=reprofile,status='replace')
close(99)

open(94,file = outfile)



call spawn(xi,vi)
! spawining at a specific place, for testing
! xi = (/2705710525.4906921,       -3873534938.3634562 ,      -2681433813.0402393 /)
! vi = (/-11372871.430080282,        73591.840957018765,       -16765518.336228890     /)
print*,xi, vi

! big loop
call timestamp
do i = 1,Nsteps
call propagate(xi,vi,x,v,time)
! print*,"time: ", time, 'r: ', sqrt(sum(x**2))
!this check doesn't actually work, since the RKF solver will keep trying to integrate past rsun
!it dies without closing the file, I think we lose everything
call isinside(x,isinside_flag)
if (isinside_flag .eqv. .false.) then
print*, "Elvis has left the building"
return
end if
call collide(x,v,vout)
! print*, "vin ", v, "vout ", vout
r = sqrt(sum(x**2))
write(94,*) x(1),x(2),x(3), v(1),v(2),v(3), vout(1),vout(2),vout(3),time, potential(r)
xi = x
vi = vout
times(i) = time

! print*,x,v
end do
close(94)
print*,"Simulation complete"
call timestamp
! call propagate()
! call collide()
!
! end do
! Now we reprocess our file
! ratio = 10
! call reprocess(Nsteps,times,ratio,outfile,reprofile)


end program
