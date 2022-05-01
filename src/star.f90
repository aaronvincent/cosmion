
!module with all the star properties
!also contains DM mass and cross section
module star

  implicit none


! Units: length = cm, time = s,energy = erg, mass = g
! Turns on analytic treatment of the temperature, density, potential
  logical :: anTemp, anDens, anPot
  integer nlines
  double precision, allocatable :: tab_mencl(:), tab_starrho(:), tab_mfr(:,:), tab_r(:), tab_vesc(:), tab_dr(:)
  double precision, allocatable :: tab_mfr_oper(:,:), tab_T(:), tab_g(:), tab_atomic(:), vesc_shared_arr(:),tab_phi(:)
  double precision :: rhoSHO,rchi,Rsun,mdm,OmegaSHO,sigSD,mu
  double precision, parameter :: c0=2.99792458d10,GN = 6.672d-8,pi=3.141592653, mnuc = 0.938,mnucg = 1.66054d-24
  double precision, parameter :: hbarc = 1.97d-14,kb = 1.3807d-16
  !this goes with the Serenelli table format
  double precision :: AtomicNumber(29) !29 is is the number from the Serenelli files; if you have fewer it shouldn't matter

!functions
  ! double precision :: ndensity, temperature, potential

  contains

  subroutine init_star(anTempIn,anDensIn,anPotIn,mdm_in,sigSD_in)
    logical, intent(in) :: anTempIn,anDensIn,anPotIn
    double precision, intent(in) :: mdm_in,sigSD_in
    character*300 :: filename
    filename = "data/struct_b16_agss09_nohead.dat" !make this an input later on

    print*,"Initializing star"

    mdm = mdm_in
    sigSD = sigSD_in
    mu = mdm/mnucg !this needs to be fixed for > 1 isotope


!if anything is not analytic, load stellar data
    if ((.not. anTemp) .or. (.not. anDens) .or. (.not. anPot)) then
        call get_solar_params(filename,nlines)
      else
        Rsun = 69.57d9
        AtomicNumber  = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                          18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                          39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                          54.93, 55.845, 58.933, 58.693/)
      end if

    if (anPot) then
      rhoSHO = 148.9d0 !g/cm^3
      rchi = sqrt(3.*kb*temperature(0.d0)/(2.*pi*GN*rhoSHO*mdm))
      OmegaSHO = sqrt(4./3.*pi*GN*rhoSHO)
    end if
    print*,"initialized star"

  end subroutine

  !Get number density of scatterers at radius r
    function ndensity(R,iso)
      double precision ndensity
      double precision, intent(in):: R
      double precision :: nnuc
      integer, intent(in) :: iso
      if (anDens) then
      nnuc = rhoSHO/mnucg
      else
        if (R/Rsun .lt. tab_r(1)) then
          nnuc = tab_starrho(1)*tab_mfr(1,iso)/AtomicNumber(iso)/mnucg
        else
          call interp1(tab_r,tab_starrho*tab_mfr(:,iso)/AtomicNumber(iso)/mnucg,nlines,R/Rsun,nnuc)
        end if
      end if
      ndensity = nnuc
    end function

    function dndr(R,iso)
  ! number density derivative of scatterers at position r. Not used
      double precision dndr
      double precision, intent(in) :: R
      integer, intent(in) :: iso
      if (anDens) then
      dndr = 0.d0
      else
      ! call interp1(tab_r,tab_starrho*tab_mfr(:,iso)/AtomicNumber(iso)/mnucg,nlines,R,nnuc)
      print*,"numerical density not implemented yet"
      end if

  end function

!get temperature at radius r
  function temperature(R)
     double precision, intent(in):: R
    double precision :: T, temperature
    if (anTemp) then
     ! T = 1.5d6 ! Constant case: nuclei equilibrate to this, which is good
     !linear temperature gradient
    T = 1.5d6*(1.-R/Rsun)
    else
      if (R/Rsun < tab_r(1)) then
        T = tab_T(1)
      else

       call interp1(tab_r,tab_T,Nlines,R/Rsun,T)
     end if
    end if
    temperature =  T
  end function

!get potential at radius r
  function potential(R)
    double precision, intent(in):: R
    double precision :: phi, potential
    if (anPot) then
      phi = 2.*pi*GN*rhoSHO*R**2/3.
    else
      call interp1(tab_r,tab_phi,Nlines,R/Rsun,phi)
    end if
    potential =  phi
  end function



!loads solar parameters from file
  subroutine get_solar_params(filename,nlines)
          character*300 :: filename
          double precision :: Pres, Lumi,GMoverR !these aren't used, but dummies are required
          ! double precision, allocatable :: phi(:) !this is used briefly !now up top
          integer :: i,j, nlines,iostatus

          Rsun = 69.57d9 !this is set here, for other stars, this sub is not called
          GMoverR=1.908e15
          !Get number of lines in the file
          open(99,file=filename)
          nlines=0
          do
            read(99,*, iostat=iostatus)
            if(iostatus/=0) then ! to avoid end of file error.
              exit
            else
              nlines=nlines+1
            end if
          end do
          close(99)
          nlines = nlines -1

          !allocate the arrays
          allocate(tab_mencl(nlines))
          allocate(tab_r(nlines))
          allocate(tab_starrho(nlines))
          allocate(tab_mfr(nlines,29)) !we could just allocate niso, but this leads to problems
          allocate(tab_vesc(nlines))
          allocate(tab_phi(nlines))
          allocate(tab_dr(nlines))
          allocate(tab_T(nlines))
          allocate(tab_g(nlines))
          ! allocate(tab_mfr_oper(nlines,16)) ! for the operator method
          allocate(vesc_shared_arr(nlines)) ! for OMP stuff


          !now actually read in the file
          open(99,file=filename)
          do i=1,nlines
            read(99,*) tab_mencl(i),tab_r(i), tab_T(i), tab_starrho(i), Pres, Lumi, tab_mfr(i,:)
          end do
          close(99)

          !we calculate the escape velocity here since all the ingredients are ready
          tab_phi(nlines) = -GMoverR
          tab_vesc(nlines) = sqrt(-2.d0*tab_phi(nlines))
          tab_dr(nlines) = tab_r(nlines)-tab_r(nlines-1)
          do i = 1,nlines-1
            j = nlines-i !trapezoid integral
            tab_phi(j) = tab_phi(j+1) + GMoverR*(tab_r(j)-tab_r(j+1))/2.*(tab_mencl(j)/tab_r(j)**2+tab_mencl(j+1)/tab_r(j+1)**2)
            tab_vesc(j) = sqrt(-2.d0*tab_phi(j)) !escape velocity in cm/s
            tab_dr(j) = -tab_r(j)+tab_r(j+1) !while we're here, populate dr
            ! tab_g(j) = -(-phi(j)+phi(j+1))/tab_dr(j)
            tab_g(i) = -GMoverR*tab_mencl(i)/tab_r(i)**2/Rsun
          end do
          ! tab_g(nlines) = tab_g(nlines-1)
          tab_g(nlines) = -GMoverR*tab_mencl(nlines)/tab_r(nlines)**2/Rsun

            ! Populate the atomic number tables here (because it relies on a specific format)
          AtomicNumber  = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                            18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                            39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                            54.93, 55.845, 58.933, 58.693/)


          return
        end subroutine get_solar_params



        subroutine interp1(xin,yin,lenxin,xout,yout)
          !1d interpolation assuming monotonically increasing vector
          integer, intent(in) :: lenxin
          double precision, intent(in) :: xin(lenxin), yin(lenxin), xout
          double precision :: yout
          integer :: i

          if (xout .lt. xin(1)) then
            print*,"xout ", xout, "xin min ", xin(1)
            stop "Error in interpolation: xout < min(xin)"
          end if
          if (xout .gt. xin(lenxin)) then
          print*,"xout ", xout, "xin max ", xin(lenxin), lenxin
          stop "Error in interpolation: xout > max(xin)"
        end if

          i = 1
          do while (xout .gt. xin(i))
            i = i+1
          end do

          yout = yin(i-1)+(yin(i)-yin(i-1))/(xin(i)-xin(i-1))*(xout-xin(i-1))
          return
        end subroutine interp1

end module
