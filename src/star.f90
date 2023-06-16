
!module with all the star properties
!also contains DM mass and cross section
module star

  implicit none


! Units: length = cm, time = s,energy = erg, mass = g
! Turns on analytic treatment of the temperature, density, potential

  logical :: anTemp, anDens, anPot,fullHistory, SHO_debug,spinDep

  integer nlines
  double precision, allocatable :: tab_mencl(:), tab_starrho(:), tab_mfr(:,:), tab_r(:), tab_vesc(:), tab_dr(:)
  double precision, allocatable :: tab_mfr_oper(:,:), tab_T(:), tab_g(:), tab_atomic(:), vesc_shared_arr(:),tab_phi(:)
  double precision :: rhoSHO,rchi,Rsun,Msun,mdm,OmegaSHO,sigSD,mu
  double precision, parameter :: c0=2.99792458d10,GN = 6.672d-8,pi=3.141592653, mnuc = 0.938,mnucg = 1.66054d-24
  double precision, parameter :: hbarc = 1.97d-14,kb = 1.3807d-16,GeV = 1.78266d-24
  !this goes with the Serenelli table format
  double precision, allocatable :: AtomicNumber(:) !29 is is the number from the Serenelli files; if you have fewer it shouldn't matter
  integer :: outside_flag, species
  integer, allocatable :: elements(:)

!functions
  ! double precision :: ndensity, temperature, potential

  contains

  subroutine init_star(anTempIn,anDensIn,anPotIn,mdm_in,sigSD_in,filename)
    logical, intent(in) :: anTempIn,anDensIn,anPotIn
    double precision, intent(in) :: mdm_in,sigSD_in
    character*300 :: filename

    print*,"Initializing star"

    mdm = mdm_in
    sigSD = sigSD_in
    mu = mdm/mnucg !this needs to be fixed for > 1 isotope
    outside_flag = 0


!if anything is not analytic, load stellar data
    if ((.not. anTemp) .or. (.not. anDens) .or. (.not. anPot)) then
        call get_stellar_params(filename,nlines)
      else
        Rsun = 69.57d9
        AtomicNumber  = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                          18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                          39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                          54.93, 55.845, 58.933, 58.693/)
      end if

    if (anPot .or. SHO_debug) then
      !rhoSHO = 148.9d0 !g/cm^3
      rhoSHO = tab_starrho(1) !set the star's density to that at the centre
      Msun = 4./3.*pi*Rsun**3*rhoSHO
      rchi = sqrt(3.*kb*temperature(0.d0)/(2.*pi*GN*rhoSHO*mdm))
      OmegaSHO = sqrt(4./3.*pi*GN*rhoSHO)
    else
      rchi = sqrt(3.*kb*temperature(0.d0)/(2.*pi*GN*tab_starrho(1)*mdm))
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
      if (R/Rsun <= tab_r(1)) then
        T = tab_T(1)
      else

       call interp1(tab_r,tab_T,Nlines,R/Rsun,T)
     end if
    end if
    temperature =  T
  end function

!get potential at radius r. Only valid inside star
  function potential(R)
    double precision, intent(in):: R
    double precision :: phi, potential
    if (anPot .or. SHO_debug) then
      !making potential negative-definite. This avoid *isssues*
      if (R .ge. Rsun) then
        phi = -4.*pi*GN*rhoSHO*Rsun**3/R
      else
      phi = 2.*pi*GN*rhoSHO*(R**2-2.*Rsun**2)/3. !- 4.*pi/3.*rhoSHO*Rsun**2*GN
    end if
    else
      if (R/Rsun .lt. tab_r(1)) then
        phi = tab_phi(1)
      else if (R .ge. Rsun ) then
        phi = -GN*Msun/R
      else
        call interp1(tab_r,tab_phi,Nlines,R/Rsun,phi)
      end if
    end if
    potential =  phi
  end function

  function vescape(R)
    double precision, intent(in):: R
    double precision ::  vescape
    if (anPot .or. SHO_debug) then
      ! phi = 2.*pi*GN*rhoSHO*R**2/3.
      vescape = sqrt(-2*potential(R))
    else
      if (R/Rsun .le. tab_r(1)) then
        vescape = tab_vesc(1)
      else if (R .ge. Rsun ) then
        vescape = sqrt(2*GN*Msun/R)
      else
        call interp1(tab_r,tab_vesc,Nlines,R/Rsun,vescape)
      end if
    end if

  end function

! get acceleration at R
  function gofr(R)
    double precision, intent(in):: R
    double precision :: gofr
    if (anPot .or. SHO_debug) then
      gofr = -4.*pi*GN*rhoSHO*R/3.
    else
      if (R/Rsun < tab_r(1)) then
        gofr = tab_g(1)
      else if (R .ge. Rsun) then
        gofr = -GN*Msun/R**2
      else
        call interp1(tab_r,tab_g,Nlines,R/Rsun,gofr)
      end if
    end if
  end function


  ! Load stellar parameters from file.
  subroutine get_stellar_params(filename,nlines)
    character*300 :: filename
    integer :: i,nlines
    call find_nlines(filename,nlines)

    !now actually read in the file
    open(99,file=filename)
    read(99,*) Msun,Rsun,AtomicNumber(:)
    !print*,"Msun=",Msun
    !print*,"Rsun=",Rsun
    !print*,"Atomic Number:",AtomicNumber
    do i=1,nlines
      read(99,*) tab_mencl(i),tab_r(i), tab_T(i), tab_starrho(i), tab_mfr(i,:)
      !if (i<5 .or. i>nlines-5) then
      !  print*,tab_r(i)
      !end if
    end do
    close(99)

    call populate_arrays(nlines)
  end subroutine

  subroutine find_nlines(filename,nlines)
    character*300 :: filename
    character*500 :: header
    integer :: nlines,nspecies,iostatus,i
    nspecies = 0
    open(99,file=filename)
    read(99,'(A)') header
    do i=2,len(header)
      ! Check for spaces in the header to find the number of columns.
      if (header(i:i) .eq. " " .and. header(i-1:i-1) .ne. " ") then
        nspecies = nspecies + 1
      end if
    end do
    nspecies = nspecies - 2
  
    nlines=0
    ! Read until the end of the file to find the number of rows.
    do
      read(99,*, iostat=iostatus)
      if(iostatus/=0) then ! to avoid end of file error.
        exit
      else
        nlines=nlines+1
      end if
    end do
    close(99)

    !allocate the arrays
    allocate(AtomicNumber(nspecies))
    allocate(tab_mencl(nlines))
    allocate(tab_r(nlines))
    allocate(tab_starrho(nlines))
    allocate(tab_mfr(nlines,nspecies)) !29we could just allocate niso, but this leads to problems
    allocate(tab_vesc(nlines))
    allocate(tab_phi(nlines))
    allocate(tab_dr(nlines))
    allocate(tab_T(nlines))
    allocate(tab_g(nlines))
    ! allocate(tab_mfr_oper(nlines,16)) ! for the operator method
    allocate(vesc_shared_arr(nlines)) ! for OMP stuff
  end subroutine

  subroutine populate_arrays(nlines)
    integer :: nlines,i,j
    double precision :: GMoverR
    !we calculate the escape velocity here since all the ingredients are ready
    GMoverR = GN*Msun/Rsun
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
    !AtomicNumber  = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
    !                  18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
    !                  39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
    !                  54.93, 55.845, 58.933, 58.693
  end subroutine


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


  subroutine select_species(precision)
    double precision, intent(in) :: precision
    double precision :: numDens,sigma,probs(size(AtomicNumber)),probsTot,pdf
    integer :: i,j,points
    if (precision == 0) then
      ! Consider every element.
      elements = (/(i,i=1,size(AtomicNumber))/)
    else
      ! Compute the average number density times the cross section for each element,
      ! weighted by its estimated radial position.
      points = 1e3 ! Number of radial points used to construct averages.
      pdf = 0.d0
      do i=1,size(AtomicNumber)
        numDens = 0.d0
        do j = 1,points
          pdf = (1.*j/points)**2*exp(-(Rsun*j/points/rchi)**2) ! weight for this radius
          numDens = numDens + ndensity(Rsun*j/points,i) * pdf
        end do
        numDens = numDens/points
        sigma = sigSD * AtomicNumber(i)**4 * ((mdm+mnucg)/(mdm+AtomicNumber(i)*mnucg))**2
        probs(i) = numDens * sigma
      end do
      probs = probs / sum(probs)
      ! Select the elements with collision probabilities above the specified precision.
      allocate(elements(0))
      do i=1,size(AtomicNumber)
        if (probs(i) >= precision) then
          elements = [elements,i]
        end if
      end do
    end if
  end subroutine

end module
