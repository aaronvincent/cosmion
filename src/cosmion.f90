!Cosmion takes 4 arguments and produces a formatted output text file


program cosmion
use star
use init_conds ! used to check energy conservation. Can remove if not needed
implicit none
interface
  subroutine omega(xin,vin,omega_out,niso)
    double precision, intent(in) :: xin(3),vin(3)
    double precision, intent(out) :: omega_out
    integer, optional :: niso
  end subroutine
end interface
character*100 :: massin, sigmain, Nstepsin, FileNameIn
double precision :: xi(3),vi(3),x(3),v(3),vout(3),xsamp(3),vsamp(3)
double precision :: r,time,start,finish,weight !the DM coordinates
double precision :: species_precision
double precision, parameter :: GeV = 1.78266d-24

character*100 :: outfile, reprofile
! logical antemp, andens, anpot
! double precision mdm
integer Nsteps, i,ratio
debug_flag = .false. !used for debugging (duh)

IF(COMMAND_ARGUMENT_COUNT().NE.4)THEN
  WRITE(*,*)'ERROR, FOUR COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
  STOP
ENDIF


!set parameters


call random_seed

CALL GET_COMMAND_ARGUMENT(1,massin)   !first, read in the two values
CALL GET_COMMAND_ARGUMENT(2,sigmain)
CALL GET_COMMAND_ARGUMENT(3,Nstepsin)   !first, read in the two values
CALL GET_COMMAND_ARGUMENT(4,FileNameIn)
print*, "m = ", trim(massin), " GeV"
! massin = trim(massin)


read(massin,*)mdm
read(sigmain,* )sigsd
read(Nstepsin,*)Nsteps
! read(FileNameIn,*)outfile
outfile = FileNameIn

print*, "initializing with "
print*, "m = ", mdm, " GeV"
print*, "sigma = ", sigsd, " cm^2"
print*, Nsteps, " Collisions"
print*, "Printing to ", outfile

mdm = mdm*GeV



!masses in g
! mdm = 5.d0*GeV
! sigsd = 1.d-37 !cm^2
! Nsteps =5e5
species_precision = 1.d-2 ! A precision of 0 uses all species for collisions.


!set up the star
!FOR A REALISTIC STAR SET EVERY FLAG TO FALSE
!anXXX flags: if false, interpolate from a stellar model. If true, use analytic
!functions defined in star.f90
!For testing only.
anTemp = .false.
anDens = .false.
!treat potential as SHO and use analytic trajectory expressions for x, v
!This is an implementation of the original Banks et al Simulations
anPot = .false.
!Spin-dependent? (true = only hydrogen)
!if false, species_precision above will determine how many isotopes to use
spinDep = .false.


!SHO_debug overrides the tabulated phi(r) to provide SHO potential
!for testing of phi(r) and comparison with anPot. Don't set to true for realistic sims
!this flag does nothing if anPot is already true.
SHO_debug = .false.


if (anPot .or. SHO_debug) then
  print*, "Watch out, you are using a SHO potential"
end if



! Set the elements that the particle will collide with,
! based on collision probabilities above the specified precision.
if (.not. spinDep) then
  call select_species(mdm/GeV,species_precision)
  print '("Elements:"(29I4.2))', elements
end if

!initialize the star
call init_star(anTemp,anDens,anPot,mdm,sigSD)

!wipe the trajectory history
! open(99,file=reprofile,status='replace')
! close(99)

open(94,file = outfile)



call spawn(xi,vi)
! spawining at a specific place, for testing
! xi = (/4576709851.6707411d0,       -0.d0 ,      -0.d0 /)
! vi = (/-6068728.6153145507d0,        96408852.454135373d0,       0.d0     /)
! xi = (/50.d9,       0.12d0 ,      -2.0402393d0 /)
! vi = (/-60.d5,        000.d5,       0.d0    /)

print*,xi/1d5, "km"
print*, vi/1d5, "km/s"
vout = vi
time = 0.d0
species = 0
weight = 1.d10 !this ensures the initial point is not really counted
write(94,*) xi(1),xi(2),xi(3), vi(1),vi(2),vi(3), vout(1),vout(2),vout(3), time, eoverm,outside_flag,weight,species

! big loop
call timestamp
call cpu_time(start)
do i = 1,Nsteps
    call propagate(xi,vi,x,v,time)

    if (outside_flag == 0) then
        call collide(x,v,vout)
        xi = x
        vi = vout
        !get the weight of this position
        call omega(xi,vi,weight) !note this doesn't work for SI!!!
    else if (outside_flag == 1) then

        outside_flag = 3 !this indicates that the weights need to be time/time_total. You need to include this weighting in your analysis script since it can't be done on the fly
        species = 0
!travel to the surface, making friends along the way
        call propagate_to_surface(xi,vi,x,v,time,xsamp,vsamp)
! We've found that we are leaving the star. Current position is now at the surface. We need to record the time that took though
! it would also help to take a sample from that path
! call omega(xi,vi,weight)
        vout = v
!this counts as in the star, but we'll write some random position sampled from this last trajectory
        weight = 1. ! this is wrong
        write(94,*) xsamp(1),xsamp(2),xsamp(3), vsamp(1),vsamp(2),vsamp(3), &
        vsamp(1),vsamp(2),vsamp(3), time, eoverm,outside_flag,weight,species

        xi = x
        vi = v


        ! print*,"r before keplerian ",sqrt(sum(x**2)), "v ", v
        call keplerian(xi,vi,x,v,time)
        ! print*,"r after keplerian ",sqrt(sum(x**2)), "v ", v
        !call keplerian_rad(xi,vi,x,v,time)
        outside_flag = 1
        vout = v
        xi = x
        vi = v
        ! weight for reentering particles has not been tracked
        !the density of particles at the boundary = #

    else if (outside_flag == 2) then
        call spawn(x,v)
        vout = v
        time = 0.d0
        species = 0
        xi = x
        vi = v
    end if

    write(94,*) x(1),x(2),x(3), v(1),v(2),v(3), vout(1),vout(2),vout(3), time, eoverm, outside_flag,weight,species
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
