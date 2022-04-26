program cosmion
use star
implicit none

double precision :: x,y,z, vx,vy,vz !the DM coordinates
double precision, parameter :: GeV = 1.78266d-24

character*100 :: outfile
! logical antemp, andens, anpot
! double precision mdm
integer Nsteps, i

!set parameters


call random_seed
!set up the star


!masses in g
mdm = 1.*GeV
anTemp = .true.
anDens = .true.
anPot = .true.

Nsteps = 1.d5

outfile = 'positions.dat'

call init_star(anTemp,anDens,anPot,mdm)

open(94,file = outfile)


! call spawn(x,y,z,vx,vy,vz)


! big loop
do i = 1,Nsteps
call spawn(x,y,z,vx,vy,vz)
write(94,*) x,y,z,vx,vy,vz
end do
close(94)
! call propagate()
! call collide()
!
! end do



end program
