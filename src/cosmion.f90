program cosmion
implicit none
double precision :: x,y,z, vx,vy,vz !the DM coordinates
integer Nsteps

!set parameters

!set up the star
call init_star()

call spawn(x,y,z,vx,vy,vz)

!big loop
do i = 1,Nsteps

call propagate()
call collide()

end do



end program
