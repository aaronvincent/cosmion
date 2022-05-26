program cosmion
use star
implicit none

double precision :: xi(3),vi(3),x(3),v(3),vout(3),r,time !the DM coordinates
double precision, parameter :: GeV = 1.78266d-24
integer outside

character*100 :: outfile
! logical antemp, andens, anpot
! double precision mdm
integer Nsteps, i

!set parameters


call random_seed
!set up the star


!masses in g
mdm = 1*GeV
sigsd = 1.d-34 !cm^2
anTemp = .false.
anDens = .false.
anPot = .true.
outside = 0

Nsteps =1e1

outfile = 'positions.dat'

call init_star(anTemp,anDens,anPot,mdm,sigSD)

open(94,file = outfile)


! call spawn(x,y,z,vx,vy,vz)
call spawn(xi,vi)
print*,xi, vi
vout = vi
time = 0.d0
write(94,*) xi(1),xi(2),xi(3), vi(1),vi(2),vi(3), vout(1),vout(2),vout(3), time,outside

! big loop
call timestamp
do i = 1,Nsteps
    call propagate(xi,vi,x,v,time)

    r = sqrt(sum(x**2))
    if (r<Rsun*0.99949) then  ! Use 0.9995 times the radius because of interp1 subroutine
        call collide(x,v,vout)
        xi = x
        vi = vout
    else
        vout = v
        xi = x
        vi = v
        write(94,*) x(1),x(2),x(3), v(1),v(2),v(3), vout(1),vout(2),vout(3), time,outside
        call keplerian(xi,vi,x,v,time)
        outside = 1
        vout = v
        xi = x
        vi = v
    end if
    write(94,*) x(1),x(2),x(3), v(1),v(2),v(3), vout(1),vout(2),vout(3), time,outside
    outside = 0

end do
close(94)
print*,"Simulation complete"
call timestamp
! call propagate()
! call collide()
!
! end do



end program
