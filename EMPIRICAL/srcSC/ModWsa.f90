!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModWsa
  use ModMagnetogram, ONLY: get_magnetogram_field
  use ModUtilities,      ONLY: CON_stop
#ifdef _OPENACC
  use ModUtilities, ONLY: norm2
#endif
  use ModNumConst, ONLY: cPi,cTwoPi, cDegToRad, cRadToDeg
  implicit none
  PRIVATE ! Except
  SAVE
  type OpenThread
     ! The thread length, \int{ds}, from the photoshere to the given point
     ! Renamed and revised: in meters
     real,pointer :: LengthSi_I(:)
     ! Magnetic field intensity, the inverse of heliocentric distance
     ! Dimensionless.
     real,pointer :: BSi_I(:),RInv_I(:), Coord_DI(:,:)
     ! number of points
     integer :: nPoint
     real    :: SignB
     real    :: ExpansionFactor
     real    :: OpenFlux ! [Gs Rsun**2]
  end type OpenThread

  ! Parameters of the open thread grid

  ! Threads are traced from rMax toward rMin
  real,    public :: rMin = 1.0, rMax = 2.5
  ! Equally spaced cell-centered grid at rMax surface, number of cells
  integer, public :: nLon = 720, nLat = 360

  ! Total number of open threads
  integer :: nOpenThreadAll = 720*360
  ! For parallel computations: last global number
  ! of open thread on the previous PE
  integer, allocatable :: iOpenThread0_P(:)
  ! Number of open threads on the given PE
  integer, allocatable :: nOpenThread_P(:)

  ! Meshes
  real    :: dLon, dLat
  ! Face areas
  real :: FaceAreaAtRmin, FaceAreaAtRmax
  ! Face areas ratio:
  real    :: GeometricExpansionInv = 0.16
  real    :: dS0 = 0.001 ! integration step at R = 1, dS propto R for R > 1
  real    :: BMinSi = 2.0e-7 ! 0.002 G
  ! All threads and points
  type(OpenThread), allocatable :: OpenThreads_II(:,:)
  logical :: DoInit = .true.

  integer, public,  allocatable :: nOpenThread_II(:,:)
  real   , public,  allocatable :: OpenFlux_II(:,:)
  public :: read_wsa_param
  public :: set_open_threads
  ! Shift of the left margin of the map with regards to
  ! the direction of axis
  ! 1. In degrees (integer)
  integer :: iLonShiftDeg = 0
  ! 2. In rads (real)
  real    :: LonShift = 0

  ! Other paramaters of the WSA files
  ! 1. Header of the B0 file used
  character(len=100):: StringDescription = '???'
  ! 2. Number of Carrington rotation
  integer    :: nRot = 0
  ! 3. Time, expressed in terms of the solar (sinodic) rotation period
  real       :: CRFraction = 0

  ! Output
  real,  allocatable :: Buff_VI(:,:), Map_VII(:,:,:)
  integer, parameter :: BrSS_ = 2, Br_ = 1 , OpenFlux_=3
  ! Field unit conversion:
  real, parameter :: Si2Gs = 10000.0
  character(len=*), parameter:: NameMod = 'ModWsa'
contains
  !============================================================================
  subroutine init(iProc,nProc,iComm)

    use ModLookupTable, ONLY: &
         i_lookup_table,  get_lookup_table
    integer, intent(in) :: iProc, nProc, iComm
    integer :: iLookupTable = -1, nParam
    ! Parameters in the B0 file
    real    :: Param_I(4)
    ! Loop variable:
    integer :: iPE
    character(len=*), parameter:: NameSub = 'init'
    !--------------------------------------------------------------------------
    DoInit = .false.
    iLookupTable = i_lookup_table('B0')
    ! Get magnetogram shift and Carrington rotation
    call get_lookup_table(iLookupTable, nParam=nParam, Param_I=Param_I, &
         Time=CRFraction, StringDescription=StringDescription)
    if(nParam > 1)then
       ! Take minimum and maximum radius from the file, if set to negative:
       if(rMin < 0.0)rMin = Param_I(1)
       if(rMax < 0.0)rMax = Param_I(2)
    endif
    if(nParam > 2)then
       if(iLonShiftDeg < 0)iLonShiftDeg = int(Param_I(3))
    else
       if(iLonShiftDeg < 0)call CON_stop(NameMod//':'//NameSub//&
            ': cannot find the offset in longitude from file '//&
            trim(StringDescription))
    end if
    LonShift = iLonShiftDeg*cDegToRad
    if(nParam > 3)nRot      = nint(Param_I(4))

    ! Load balance
    nOpenThreadAll = nLon*nLat
    allocate(iOpenThread0_P(0:nProc-1))
    allocate(nOpenThread_P( 0:nProc-1))

    ! Last global open thread number on iProc-1 PE
    do iPE = 0, nProc - 1
       iOpenThread0_P(iPE) = (iPE*nOpenThreadAll)/nProc
    end do
    do iPE = 0, nProc - 2
       nOpenThread_P(iPE) = iOpenThread0_P(iPE+1) - &
            iOpenThread0_P(iPE)
    end do
    nOpenThread_P(nProc-1) = nOpenThreadAll - &
            iOpenThread0_P(nProc-1)

    ! Meshes
    dLon = cTwoPi/nLon; dLat = cPi/nLat

    ! Factors to convert the field to open flux and back
    FaceAreaAtRmin = dLon*dLat*Rmin**2
    FaceAreaAtRmax = dLon*dLat*Rmax**2

    allocate(OpenThreads_II(nLon,nLat))

    allocate(nOpenThread_II(nLon,nLat))
    allocate(OpenFlux_II(nLon,nLat))
    allocate(Buff_VI(BrSS_,nOpenThreadAll))
    allocate(Map_VII(OpenFlux_,nLon,nLat))

  end subroutine init
  !============================================================================
  subroutine read_wsa_param

    use ModReadParam,   ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('rMin', rMin)
    call read_var('rMax', rMax)
    call read_var('iLonShiftDeg', iLonShiftDeg)
    call read_var('nLon', nLon)
    call read_var('nLat', nLat)
    call read_var('dS0' , dS0 )
    call read_var('BMinSi', BMinSi)
  end subroutine read_wsa_param
  !============================================================================
  subroutine set_open_threads(iProc, nProc, iComm)

    use ModPlotFile, ONLY: save_plot_file
    ! Sets all open threads originating from the source surface
    use ModMpi
    integer, intent(in) :: iProc, nProc, iComm

    ! Loop variables
    integer :: iOpenThreadAll, iLon, iLat, iPe, nPoint ! iLonRMin, iLatRMin, nPoint
    integer :: iError
    ! Generalized and Cartesian coordinates to calculate the open flux
    real    :: Coord_D(3), Xyz_D(3), OpenFlux
    real, parameter :: cHalfPi = cPi/2
    ! Named indexes for open fluxes:
    integer, parameter :: Unsigned_ = 1, Pos_ = 2, Neg_ = 3, Lon_ = 1, Lat_ = 2
    ! Counters for open fluxes
    real :: OpenFluxTotal_I(Unsigned_:Neg_)
    ! Second order nterpolation stencil at R=RMin:
    real :: dCoord_D(Lon_:Lat_), CoordNorm_D(Lon_:Lat_), CosLatDown, CosLatUp
    integer :: iLonDown, iLonUp, iLatDown, iLatUp
    ! Call only once
    !--------------------------------------------------------------------------
    if(DoInit)call init(iProc, nProc, iComm)
    ! Nullify counters
    nOpenThread_II = 0
    OpenFlux_II = 0.0
    OpenFluxTotal_I(Unsigned_:Neg_) = 0.0
    Buff_VI = 0.0
    ! Construct the lines assigned to the given PE
    do iOpenThreadAll = iOpenThread0_P(iProc) + 1, &
         iOpenThread0_P(iProc) + nOpenThread_P(iProc)
       iLat = (iOpenThreadAll - 1)/nLon + 1
       iLon = iOpenThreadAll - nLon*(iLat - 1)
       call set_open_thread(&
            CoordIn_D = [rMax, LonShift + (iLon -0.50)*dLon, &
            -cHalfPi + (iLat - 0.50)*dLat], &
            iLonIn = iLon, iLatIn = iLat, iBuff = iOpenThreadAll)
       nPoint = OpenThreads_II(iLon, iLat)%nPoint
       ! Buff_VI(ExpFactor_,iOpenThreadAll)  = GeometricExpansionInv*&
       !     OpenThreads_II(iLon, iLat)%BSi_I(-nPoint)/&
       !     OpenThreads_II(iLon, iLat)%BSi_I(0)
       Coord_D = OpenThreads_II(iLon, iLat)%Coord_DI(:,-nPoint)
       ! Normalized coords relative to the map margins
       CoordNorm_D(Lon_:Lat_) = [&
            modulo(Coord_D(2) - LonShift,cTwoPi)/dLon, &
            (cHalfPi + Coord_D(3))/dLat]
       iLonDown = int(CoordNorm_D(Lon_) +0.50); iLonUp = iLonDown + 1
       dCoord_D(Lon_) = CoordNorm_D(Lon_) - (iLonDown - 0.50)
       if(iLonDown==0)iLonDown = nLon
       if(iLonUp == nLon + 1)iLonUp = 1
       iLatDown = max(int(CoordNorm_D(Lat_) + 0.50), 1)
       CosLatDown = cos(-cHalfPi + (iLatDown - 0.50)*dLat)
       iLatUp = min(iLatDown + 1, nLat)
       CosLatUp   = cos(-cHalfPi + (iLatUp - 0.50)*dLat)
       dCoord_D(Lat_) = max(CoordNorm_D(Lat_) - (iLatDown - 0.50), 0.0)
       ! iLonRMin = min(int(CoordNorm_D(Lon_)) + 1, nLon)
       ! iLatRMin = min(int(CoordNorm_D(Lat_)) + 1, nLat)
       ! nOpenThread_II(iLonRMin,iLatRMin) = &
       !     nOpenThread_II(iLonRMin,iLatRMin) + 1
       if(any(dCoord_D<0.0).or.any(dCoord_D>1.0))then
          write(*,*)dCoord_D
          call CON_stop('Check weight!')
       end if
       OpenFlux = OpenThreads_II(iLon, iLat)%OpenFlux
       OpenFluxTotal_I(Unsigned_) = OpenFluxTotal_I(Unsigned_) + abs(OpenFlux)
       OpenFluxTotal_I(Pos_) = OpenFluxTotal_I(Pos_) + max(OpenFlux, 0.0)
       OpenFluxTotal_I(Neg_) = OpenFluxTotal_I(Neg_) + min(OpenFlux, 0.0)
       OpenFlux_II(iLonDown,iLatDown) = OpenFlux_II(iLonDown,iLatDown) + &
            OpenFlux/(FaceAreaAtRmin*CosLatDown)* &
            (1 - dCoord_D(Lon_))*(1 - dCoord_D(Lat_)) ! Wieght
       OpenFlux_II(iLonUp,iLatDown) = OpenFlux_II(iLonUp,iLatDown) + &
            OpenFlux/(FaceAreaAtRmin*CosLatDown)* &
            dCoord_D(Lon_)*(1 - dCoord_D(Lat_))       ! Wieght
       OpenFlux_II(iLonDown,iLatUp) = OpenFlux_II(iLonDown,iLatUp) + &
            OpenFlux/(FaceAreaAtRmin*CosLatUp)* &
            (1 - dCoord_D(Lon_))*dCoord_D(Lat_)       ! Wieght
       OpenFlux_II(iLonUp,iLatUp) = OpenFlux_II(iLonUp,iLatUp) + &
            OpenFlux/(FaceAreaAtRmin*CosLatUp)* &
            dCoord_D(Lon_)*dCoord_D(Lat_)             ! Wieght

    end do
    if(nProc > 1) then
       call MPI_reduce_real_array(OpenFluxTotal_I, 3,&
            MPI_SUM, 0, iComm, iError)
       call MPI_allreduce(OpenFlux_II, Map_VII(OpenFlux_,:,:), &
            nOpenThreadAll, MPI_REAL, MPI_SUM, iComm, iError)
       iOpenThreadAll = 0
       do iPE = 0, nProc-1
          if(nOpenThread_P(iPE)==0)CYCLE
          call MPI_BCAST(Buff_VI(:, iOpenThreadAll+1:&
               iOpenThreadAll+nOpenThread_P(iPE)),&
               BrSS_*nOpenThread_P(iPE), MPI_REAL, iPE, iComm, iError)
          iOpenThreadAll = iOpenThreadAll + nOpenThread_P(iPE)
       end do
    else
       Map_VII(OpenFlux_,:,:) = OpenFlux_II
    end if
    if(iProc==0)write(*,*) &
         'Open flux at Rmax: unsigned, positive, negative',&
         OpenFluxTotal_I
    OpenFluxTotal_I = 0
    do iOpenThreadAll = iOpenThread0_P(iProc) + 1, &
         iOpenThread0_P(iProc) + nOpenThread_P(iProc)
       iLat = (iOpenThreadAll - 1)/nLon + 1
       iLon = iOpenThreadAll - nLon*(iLat - 1)
       OpenFlux = Map_VII(OpenFlux_,iLon, iLat)*FaceAreaAtRmin*&
            cos(-cHalfPi + (iLat - 0.50)*dLat)
       OpenFluxTotal_I(Unsigned_) = OpenFluxTotal_I(Unsigned_) + abs(OpenFlux)
       OpenFluxTotal_I(Pos_) = OpenFluxTotal_I(Pos_) + max(OpenFlux, 0.0)
       OpenFluxTotal_I(Neg_) = OpenFluxTotal_I(Neg_) + min(OpenFlux, 0.0)
    end do
    if(nProc > 1)call MPI_reduce_real_array(OpenFluxTotal_I, 3,&
         MPI_SUM, 0, iComm, iError)
    Map_VII(Br_:BrSS_,:,:) = reshape(Buff_VI,[BrSS_,nLon,nLat])
    if(iProc==0)then
       write(*,*)'Open flux at Rmax: unsigned, positive, negative',&
            OpenFluxTotal_I
       call save_plot_file(NameFile = 'wsa_2d.out',  &
            nDimIn  = 2,                       &
            TimeIn = CRFraction,               &
            ParamIn_I= [real(iLonShiftDeg), real(nRot)],&
            VarIn_VII= Map_VII, &
            TypeFileIn    = 'ascii',           &
            CoordMinIn_D  = [180.0/nLon, -90 + 90.0/nLat], &
            CoordMaxIn_D  = [360 - 180.0/nLon, 90 - 90.0/nLat ],&
            StringHeaderIn  = StringDescription, &
            NameUnitsIn  = &
            ' [deg] [deg] [Gs] [Gs] [Gs]  [deg] []', &
            NameVarIn = 'Longitude Latitude Br BrSS OpenFlux'//&
            ' Long0 CRNumber')
    end if
  end subroutine set_open_threads
  !============================================================================
  subroutine set_open_thread(CoordIn_D, iLonIn, iLatIn, iBuff)
    use ModCoordTransform, ONLY: xyz_to_rlonlat, rlonlat_to_xyz
    use ModConst,          ONLY: rSun, MSun, cGravitation
    ! Origin point coordinates (r-lon-lat)
    real, intent(in) :: CoordIn_D(3)
    integer, intent(in) :: iLonIn, iLatIn, iBuff
    ! loop variable
    integer :: iPoint, nPoint
    ! Length interval, ! Heliocentric distance
    real :: Ds, R
    ! coordinates, field vector and modulus
    real :: Xyz_D(3), Coord_D(3), B0_D(3), B0, BMax
    !  for ourward directed field, -1 otherwise
    real ::SignBr
    ! Radial field at the source surface
    real :: BrSS
    ! Coordinates and magnetic field in the midpoint
    ! within the framework of the Runge-Kutta scheme
    real :: XyzAux_D(3), B0Aux_D(3)
    ! Aux
    real :: ROld, Aux, OpenFlux
    real :: XyzOld_D(3), Dir1_D(3), Dir2_D(3), Dir3_D(3), Dir4_D(3)
    integer, parameter :: nPointMax = 10000
    real :: BSi_I(-nPointMax:0),RInv_I(-nPointMax:0), LengthSi_I(-nPointMax:0)
    real :: Coord_DI(3,-nPointMax:0)
    ! Trace the open thread with the given origin point
    character(len=*), parameter:: NameSub = 'set_open_thread'
    !--------------------------------------------------------------------------
    iPoint = 0
    call rlonlat_to_xyz(CoordIn_D,Xyz_D)
    R = RMax
    RInv_I(0) = 1/R
    ! Store Br field at the projection onto Rmin sphere:
    call get_magnetogram_field(Xyz_D*Rmin*RInv_I(0),B0_D)
    Buff_VI(Br_,iBuff) = sum(Xyz_D*B0_D)*RInv_I(0)*Si2Gs
    call get_magnetogram_field(Xyz_D,B0_D)
    BrSS = sum(Xyz_D*B0_D)*RInv_I(0)
    Buff_VI(BrSS_,iBuff) = BrSS*Si2Gs
    OpenFlux = FaceAreaAtRmax*cos(CoordIn_D(3))*Buff_VI(BrSS_,iBuff)
    B0 = norm2(B0_D)
    Coord_DI(:,0) = CoordIn_D
    if(abs(BrSS) < BMinSi)then
       ! Origin point is very close to null point of the potential fieild
       ! Move the origin point down:
       Aux = 1 - Ds0*Rmax*(1 - abs(BrSS)/BMinSi)
       Xyz_D = Xyz_D*Aux
       R = R*Aux
       RInv_I(0) = 1/R
       call xyz_to_rlonlat(Xyz_D,Coord_DI(:,0))
       call get_magnetogram_field(Xyz_D,B0_D)
       B0 = norm2(B0_D )
       Buff_VI(BrSS_,iBuff) = sum(Xyz_D*B0_D)*RInv_I(0)*Si2Gs
       OpenFlux = FaceAreaAtRmax*cos(CoordIn_D(3))*Buff_VI(BrSS_,iBuff)
       Coord_DI(1,0) = R
    end if
    OpenThreads_II(iLonIn, iLatIn)%OpenFlux = OpenFlux
    SignBr = sign(1.0, sum(Xyz_D*B0_D) )
    BMax = max(B0, BMinSi)
    BSi_I(0) = B0
    do
       iPoint = iPoint + 1
       ! If the number of gridpoints in the theads is too
       ! high, coarsen the grid
       if(iPoint > nPointMax)call CON_stop(&
            NameMod//':'//NameSub//': too long open thread')
       ! For the previous point given are Xyz_D, B0_D, B0
       ! R is only used near the photospheric end.
       ! Store R
       ROld = R
       ! Store a point
       XyzOld_D = Xyz_D
       ! Four-stage Runge-Kutta
       BMax = max(BMax, B0)
       Ds = Ds0*R*SignBr/BMax
       Dir1_D = B0_D
       ! 1. Point at the half of length interval:
       XyzAux_D = Xyz_D - 0.50*Ds*Dir1_D
       ! 2. Magnetic field in this point:
       call get_magnetogram_field(XyzAux_D, Dir2_D)
       XyzAux_D = Xyz_D - 0.50*Ds*Dir2_D
       call get_magnetogram_field(XyzAux_D, Dir3_D)
       XyzAux_D = Xyz_D - Ds*Dir3_D
       call get_magnetogram_field(XyzAux_D, Dir4_D)
       ! 3. New grid point:
       Xyz_D = Xyz_D - (Ds/6)*(Dir1_D + 2*Dir2_D + 2*Dir3_D + Dir4_D)
       R = norm2(Xyz_D)
       if(R >  rMax)then
          nPoint = iPoint
          write(*,*)'BrSS=',Buff_VI(BrSS_,iBuff),' Gs, SignBr=', SignBr
          do iPoint = 0, nPoint-1
             write(*,*)iPoint, Coord_DI(:,-iPoint), BSi_I(-iPoint)*Si2Gs
          end do
          call CON_stop(&
               NameMod//':'//NameSub//': thread comes beyond source surface')
       end if
       RInv_I(-iPoint) = 1/R
       call xyz_to_rlonlat(Xyz_D, Coord_D)
       Coord_DI(:,-iPoint) = Coord_D
       call get_magnetogram_field(Xyz_D, B0_D)
       B0 = norm2(B0_D)
       BSi_I(-iPoint) = B0
       if(R <= rMin)EXIT
    end do
    ! Calculate more accurately the intersection point
    ! with the photosphere surface
    Aux = (ROld - RMin) / (ROld - R)
    Xyz_D = (1 - Aux)*XyzOld_D +  Aux*Xyz_D
    ! Store the last point
    RInv_I(-iPoint) = 1/RMin
    call get_magnetogram_field(Xyz_D, B0_D)
    B0 = norm2(B0_D)
    BSi_I(-iPoint) = B0
    ! Allocate thread and store results from tracing
    ! allocate(OpenThreads_II(iLonIn, iLatIn)%LengthSi_I(-iPoint:0))
    ! OpenThreads_II(iLonIn, iLatIn)%LengthSi_I(-iPoint:0) = &
    !     LengthSi_I(-iPoint:0)
    allocate(OpenThreads_II(iLonIn, iLatIn)%BSi_I(-iPoint:0))
    OpenThreads_II(iLonIn, iLatIn)%BSi_I(-iPoint:0) = BSi_I(-iPoint:0)
    allocate(OpenThreads_II(iLonIn, iLatIn)%RInv_I(-iPoint:0))
    OpenThreads_II(iLonIn, iLatIn)%RInv_I(-iPoint:0) = RInv_I(-iPoint:0)
    allocate(OpenThreads_II(iLonIn, iLatIn)%Coord_DI(3,-iPoint:0))
    OpenThreads_II(iLonIn, iLatIn)%Coord_DI(:,-iPoint:0) = &
         Coord_DI(:,-iPoint:0)
    OpenThreads_II(iLonIn, iLatIn)%nPoint = iPoint
  end subroutine set_open_thread
  !============================================================================
end module ModWsa
!==============================================================================
