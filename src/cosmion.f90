program cosmion
use star
use init_conds ! used to check energy conservation. Can remove if not needed
implicit none

double precision :: xi(3),vi(3),x(3),v(3),vout(3),r,time,start,finish !the DM coordinates
double precision, allocatable :: times(:) !makes reprocessing a bit easier to keep this in memory
double precision, parameter :: GeV = 1.78266d-24

character*100 :: outfile, reprofile
! logical antemp, andens, anpot
! double precision mdm
integer Nsteps, i,ratio
debug_flag = .false. !used for debugging (duh)



outfile = 'positions.dat'
reprofile = 'rep_pos.dat' !only used if fullHistory = true

!set parameters


call random_seed
!set up the star


!masses in g
mdm = 5.d0*GeV
sigsd = 1.d-37 !cm^2
Nsteps =1e5

!FOR A REALISTIC STAR SET EVERY FLAG TO FALSE
!anXXX flags: if false, interpolate from a stellar model. If true, use analytic
!functions defined in star.f90
!For testing only.
anTemp = .false.
anDens = .false.
!treat potential as SHO and use analytic trajectory expressions for x, v
anPot = .false.
!Spin-dependent? (true = only hydrogen)
spinDep = .true.


!SHO_debug overrides the tabulated phi(r) to provide SHO potential
!for testing of phi(r) and comparison with anPot. Don't set to true for realistic sims
!this flag does nothing if anPot is already true.
SHO_debug = .false.
!turn this on if you want the full trajectory history, not just timestamps at every collision
!Note this will take a ludicrous amount of HD space. It also doesn't work right now
fullHistory = .false.

if (anPot .or. SHO_debug) then
  print*, "Watch out, you are using a SHO potential"
end if

if ((.not. anPot) .and. (.not. spinDep)) then
  print*, "Spin-independent collisions do not yet work for non-analytic potential"
  stop
end if


allocate(times(Nsteps))


call init_star(anTemp,anDens,anPot,mdm,sigSD)

! open(94,file = "potential.dat")
! do i = 1,100
!   r = Rsun*dble(i-1)/100.
!   write(94,*) r, potential(r)
! end do
! close(94)



!wipe the trajectory history
open(99,file=reprofile,status='replace')
close(99)

open(94,file = outfile)



call spawn(xi,vi)
! spawining at a specific place, for testing
! xi = (/2705710525.4906921,       -3873534938.3634562 ,      -2681433813.0402393 /)
! vi = (/-11372871.430080282,        73591.840957018765,       -16765518.336228890     /)
! xi = (/5.4906921d10,       0.12d0 ,      -2.0402393d0 /)
! vi = (/500.d5,        100.d5,       0.d0    /)
print*,xi, vi
vout = vi
time = 0.d0
species = 0
write(94,*) xi(1),xi(2),xi(3), vi(1),vi(2),vi(3), vout(1),vout(2),vout(3), time, eoverm,outside_flag,species

! big loop
call timestamp
call cpu_time(start)
do i = 1,Nsteps
    call propagate(xi,vi,x,v,time)

    if (outside_flag == 0) then
        call collide(x,v,vout)
        xi = x
        vi = vout
    else if (outside_flag == 1) then
        vout = v
        xi = x
        vi = v
        outside_flag = 0
        species = 0
        write(94,*) x(1),x(2),x(3), v(1),v(2),v(3), vout(1),vout(2),vout(3), time, eoverm,outside_flag,species
        ! print*,"r before keplerian ",sqrt(sum(x**2)), "v ", v
        call keplerian(xi,vi,x,v,time)
        ! print*,"r after keplerian ",sqrt(sum(x**2)), "v ", v
        !call keplerian_rad(xi,vi,x,v,time)
        outside_flag = 1
        vout = v
        xi = x
        vi = v


    else if (outside_flag == 2) then
        call spawn(x,v)
        vout = v
        time = 0.d0
        species = 0
        xi = x
        vi = v
    end if
    write(94,*) x(1),x(2),x(3), v(1),v(2),v(3), vout(1),vout(2),vout(3), time, eoverm, outside_flag,species
    outside_flag = 0

end do
close(94)
print*,"Simulation complete"
call timestamp
call cpu_time(finish)
Print*, Nsteps, "Collisions perfomed in: ",finish-start, " seconds (", (finish-start)/dble(Nsteps), " sec/coll.)"

! call propagate()
! call collide()
!
! end do
! Now we reprocess our file
! ratio = 10
! call reprocess(Nsteps,times,ratio,outfile,reprofile)


end program
