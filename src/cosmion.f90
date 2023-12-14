! Cosmion takes 6 arguments and produces a formatted output text file
program cosmion
use star
use init_conds ! used to check energy conservation. Can remove if not needed
implicit none
interface
  subroutine omegaelec(xin,vin,omega_out)
    double precision, intent(in) :: xin(3),vin(3)
    double precision, intent(out) :: omega_out
  end subroutine
  subroutine omeganuc(xin,vin,omega_out,niso)
    double precision, intent(in) :: xin(3),vin(3)
    double precision, intent(out) :: omega_out
    integer, optional :: niso
  end subroutine
end interface
character*100 :: massin, sigmain, Nstepsin, FileNameIn, CrossSectionIn, ParticleIn
double precision :: xi(3),vi(3),x(3),v(3),vout(3),xsamp(3),vsamp(3)
double precision :: r,time,start,finish,weight !the DM coordinates
double precision :: species_precision
character*300 :: outfile, reprofile,spindepin, starfile
! logical antemp, andens, anpot
! double precision mdm

 ! ----- variables for portable seed setting -----
  INTEGER :: i_seed
  INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
  INTEGER, DIMENSION(1:8) :: dt_seed

integer Nsteps, i,ratio


debug_flag = .false. !used for debugging (duh)

IF(COMMAND_ARGUMENT_COUNT().NE.6)THEN
  WRITE(*,*)'ERROR, SIX COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
  STOP
ENDIF

starfile = "data/struct_b16_agss09_modified.dat"

!set parameters

!initialize RNG
CALL RANDOM_SEED(size=i_seed)
ALLOCATE(a_seed(1:i_seed))
CALL RANDOM_SEED(get=a_seed)
CALL DATE_AND_TIME(values=dt_seed)
a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
CALL RANDOM_SEED(put=a_seed)
DEALLOCATE(a_seed)

! call random_seed


CALL GET_COMMAND_ARGUMENT(1,massin)   !first, read in the two values
CALL GET_COMMAND_ARGUMENT(2,sigmain)
CALL GET_COMMAND_ARGUMENT(3,Nstepsin)   !first, read in the two values
CALL GET_COMMAND_ARGUMENT(4,FileNameIn)
CALL GET_COMMAND_ARGUMENT(5,ParticleIn)
CALL GET_COMMAND_ARGUMENT(6,CrossSectionIn)
print*, "m = ", trim(massin), " GeV"
! massin = trim(massin)



read(massin,*)mdm
read(sigmain,* )sigsd
read(Nstepsin,*)Nsteps
! read(FileNameIn,*)outfile
outfile = FileNameIn




print*, "Initializing with :"
print*, "m = ", mdm, " GeV"
print*, "sigma = ", sigsd, " cm^2"

IF(ParticleIn == 'nucleonSD')THEN
  nucleon = .true.
  spinDep = .true.
  print*, "Collision with Nucleons, Spin-Dependent"
ELSEIF(ParticleIn == 'nucleonSI')THEN
  nucleon = .true.
  spinDep = .false.
  print*, "Collision with Nucleons, Spin-Independent"
ELSEIF(ParticleIn == 'electron')THEN
  nucleon = .false.
  spinDep = .false.
  print*, "Collision with Electrons"
ELSE
  WRITE(*,*)'ERROR, particle is not correctly defined (nucleonSD, nucleonSI, electron), STOPPING'
  print*, ParticleIn
  STOP
ENDIF

IF(CrossSectionIn == 'const')THEN
  crossSection = 1
  print*, "Constant cross section"
ELSEIF(CrossSectionIn == 'v2')THEN
  crossSection = 2
  print*, "Cross section proportional to v^2"
ELSEIF(CrossSectionIn == 'q2')THEN
  crossSection = 5
  print*, "Cross section proportional to q^2"
ELSEIF(CrossSectionIn == 'v4')THEN
  crossSection = 3
  print*, "Cross section proportional to v^4"
ELSEIF(CrossSectionIn == 'q4')THEN
  crossSection = 6
  print*, "Cross section proportional to q^4"
ELSEIF(CrossSectionIn == 'vm2')THEN
  crossSection = 4
  print*, "Cross section proportional to v^(-2)"
ELSEIF(CrossSectionIn == 'qm2')THEN
  crossSection = 7
  print*, "Cross section proportional to q^(-2)"
ELSE
  WRITE(*,*)'ERROR, Dependence of the cross section not defined correctly (const, v2, q2, v4, q4, vm2 or qm2), STOPPING'
  STOP
ENDIF

print*, Nsteps, " collisions"
print*, "Printing to ", outfile

mdm = mdm*GeV



!masses in g
! mdm = 5.d0*GeV
! sigsd = 1.d-37 !cm^2
! Nsteps =5e5
species_precision = 1.d-2 ! A precision of 0 uses all species for collisions

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
! spinDep = .false.

!SHO_debug overrides the tabulated phi(r) to provide SHO potential
!for testing of phi(r) and comparison with anPot. Don't set to true for realistic sims
!this flag does nothing if anPot is already true.
SHO_debug = .false.


if (anPot .or. SHO_debug) then
  print*, "Watch out, you are using a SHO potential"
end if

!initialize the star
call init_star(anTemp,anDens,anPot,mdm,sigSD,particle,crossSection,starfile)


if (nucleon) then
  ! Set the elements that the particle will collide with,
  ! based on collision probabilities above the specified precision.
  if (.not. spinDep) then
    call select_species(species_precision)
    print '("Elements:"(29I4.2))', elements
  end if
  end if


!wipe the trajectory history
! open(99,file=reprofile,status='replace')
! close(99)

open(94,file = outfile)

kepflag = 0 !get rid of this noise later

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
weight = 1.d10 !this ensures the initial point is not really counted
write(94,*) xi(1),xi(2),xi(3), vi(1),vi(2),vi(3), vout(1),vout(2),vout(3), time, eoverm,outside_flag,weight


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! MONTE CARLO BEGINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! big loop
call timestamp
call cpu_time(start)
do i = 1,Nsteps
  
    call propagate(xi,vi,x,v,time)
    ! print*, "finished propagate with r/Rsun = ", sqrt(sum(x**2))/Rsun
    ! print*, "outside_flag ", outside_flag
    if (outside_flag == 0) then
    
        call collide(x,v,vout)
        xi = x
        vi = vout
        !get the weight of this position
        if (.not. nucleon) then
          call omegaelec(xi,vi,weight)
        else
          call omeganuc(xi,vi,weight)
        end if
      
    else if (outside_flag == 1) then
    

        ! outside_flag = 3 !this indicates that the weights need to be time/time_total. You need to include this weighting in your analysis script since it can't be done on the fly
!travel to the surface, making friends along the way
        if (anpot) then
        call propagate_to_surface(xi,vi,x,v,time)
        else 
        !these are for testing
        xsamp = xi
        vsamp = vi
        !if in fully numeric mode, propagate to surface needs to be called until the outside flag is triggered again
        !if v < 0, we will propagate the particle to some random optical depths until it turns around and v > 0, then shoot it straight at the surface
        !we'll record positions, with zero momentum transfer (because we are not calling collision) all along the way
        outside_flag = 0
        !this needs to run the full thing below (with the write) while outside-flag is 0, but if it hits -1, we need to restart that step (it means we overshot)
        ! print*, "vr at the beginning of propagato to surf ", sum(xi*vi)/sqrt(sum(xi**2))
        ! print*, "r at the beginning of propagato to surf ", sqrt(sum(xi**2))
        
        do while (outside_flag .ne. 1 ) 
        
          call propagate_to_surface(xi,vi,x,v,time)          
          do while (outside_flag ==-1 ) 
          ! print*,"retrying"
            call propagate_to_surface(xi,vi,x,v,time) !this should work because we're drawing a new random number          
          end do
          
          if (.not. nucleon) then
            call omegaelec(xi,vi,weight)
          else
            call omeganuc(xi,vi,weight)
          end if
          write(94,*) x(1),x(2),x(3), vi(1),vi(2),vi(3), &
          v(1),v(2),v(3), time, eoverm,outside_flag,weight
          vi = v
          xi = x
        end do

        ! call propagate_to_surface(xi,vi,x,v,time)

  ! We've found that we are leaving the star. Current position is now at the surface. We need to record the time that took though
  ! it would also help to take a sample from that path
  ! call omega(xi,vi,weight)
          vout = v
  ! !this counts as in the star, but we'll write some random position sampled from this last trajectory
  !         weight = 1. ! this is wrong
  !         write(94,*) xsamp(1),xsamp(2),xsamp(3), vsamp(1),vsamp(2),vsamp(3), &
  !         vsamp(1),vsamp(2),vsamp(3), time, eoverm,outside_flag,weight

        xi = x
        vi = v

     

        !print*,"r before keplerian ",sqrt(sum(x**2)), "v ", v
        if (anPot) then
          call keplerian(xi,vi,x,v,time)
        else
        ! print*,"calling keplerian at i ", i
        ! print*," before keplerian, xi = ", x, "vi = ", vi
          call keplerian_rad(xi,vi,x,v,time)
          if (time .ne. time) then
            print*, 'time is NaN. params: '
            print*,xi,vi,x,v,time
            print*,'vesc here is ', vescape(sqrt(sum(xi**2)))
            print*, 'positions before the mess'
            print*,xsamp, sqrt(sum(xsamp**2))/Rsun
            print*, 'velocities before the mess'
            print*,vsamp, sqrt(sum(vsamp**2))/vescape(sqrt(sum(xsamp**2)))
            stop
            end if
        end if
        !print*,"r after keplerian ",sqrt(sum(x**2)), "v ", v
        outside_flag = 1
        vout = v
        xi = x
        vi = v
        ! print*," after keplerian, xi = ", x, "vi = ", vi, "i = ", i
        kepflag = 1
        ! weight for reentering particles has not been tracked
        !the density of particles at the boundary = #
  end if !anpot
    else if (outside_flag == 2) then
        call spawn(x,v)
        vout = v
        time = 0.d0
        xi = x
        vi = v
    end if

    write(94,*) x(1),x(2),x(3), v(1),v(2),v(3), vout(1),vout(2),vout(3), time, eoverm, outside_flag,weight
    outside_flag = 0

end do
close(94)
print*,"Simulation complete"
call timestamp
call cpu_time(finish)
Print*, Nsteps, "Collisions perfomed in: ",finish-start, " seconds (", (finish-start)/dble(Nsteps), " sec/coll.)"

! call propagate()
! call collide()

! end do
! Now we reprocess our file
! ratio = 10
! call reprocess(Nsteps,times,ratio,outfile,reprofile)


end program
