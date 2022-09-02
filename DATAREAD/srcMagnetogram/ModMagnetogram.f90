!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModMagnetogram

  use ModNumConst,       ONLY: cTwoPi, cRadToDeg, cDegToRad
  use CON_axes,          ONLY: dLongitudeHgrDeg
  use ModUtilities,      ONLY: CON_stop
  use ModCoordTransform, ONLY: rot_xyz_sph, rot_xyz_rlonlat, rot_matrix_z,&
       xyz_to_rlonlat
  use ModLookupTable, ONLY: &
       i_lookup_table, init_lookup_table, make_lookup_table_3d, &
       get_lookup_table, interpolate_lookup_table, Table_I, TableType

  implicit none
  save

  private

  integer, public:: iTableB0    = -1     ! index of B0 lookup table
  integer, public:: iTableB0New = -1     ! index of new B0 lookup table
  public:: read_magnetogram_param        ! read parameters
  public:: init_magnetogram_lookup_table ! initialize lookup table
  public:: get_magnetogram_field         ! get field at a given point

  ! Local variables -------------

  ! True if new harmonics coefficients are to be read from harmonics files
  logical:: DoReadHarmonics = .false., DoReadHarmonicsNew = .false.

  ! Name of the harmonics files
  character(len=100):: NameHarmonicsFile = '???', NameHarmonicsFileNew = '???'

  ! variable names in the lookup table created from harmonics files
  character(len=*), parameter:: &
       NameVarLinear = 'r lon lat bx by bz rMin rMax dLon CR'
  character(len=100):: NameVar = NameVarLinear

  ! Carrington rotation of the magnetogram(s) for temporal interpolation
  real:: CarringtonRot = -1.0, CarringtonRotNew = -1.0

  ! radius of magnetogram and source surface for spatial extrapolation
  real   :: rMagnetogram=1.0, rSourceSurface=2.5
  logical:: IsLogRadius = .false.   ! logarithmic grid in radius
  integer:: nR=30, nLon=72, nLat=30 ! Size of grid created from harmonics

  ! Maximum and actual order of spherical harmonics
  integer:: MaxOrder = 30, nOrder = 30

  ! Weights of the spherical harmonics
  real, allocatable:: g_II(:,:), h_II(:,:)

  ! Powers of radius ratios
  real, allocatable:: rRsPower_I(:), RmRsPower_I(:), RmRPower_I(:)

  ! Legendre polynomials and its derivative
  real, allocatable:: p_II(:,:), Dp_II(:,:)

  ! Azimuthal functions sin(m*phi) and cos(m*phi)
  real, allocatable:: SinPhi_I(:), CosPhi_I(:)

  ! Square roots of integers in spherical harmonics
  real, allocatable:: Sqrt_I(:), SqrtRatio_I(:)

  ! Lookup table related variables
  real:: rMinB0=1.0 ! radial limits of table
  real, public:: rMaxB0=30.0 ! radial limits of table
  real:: LonMinB0 = 0.0    ! starting longitude in the field lookup table (rad)
  real:: dLonB0   = 0.0    ! longitude shift
  real:: RotB0_DD(3,3)     ! rotation matrix due to longitude shift
  real:: dLonB0New   = 0.0 ! longitude shift
  real:: RotB0New_DD(3,3)  ! rotation matrix due to longitude shift

  real:: FactorB0 = 1.0 ! multiplier for the magnetic field

  interface get_magnetogram_field
     module procedure get_magnetogram_field11, get_magnetogram_field31
  end interface get_magnetogram_field

contains
  !============================================================================
  subroutine init_magnetogram_lookup_table(iComm)
    use ModNumConst, ONLY: cDegToRad
    use ModKind,     ONLY: Real4_
    use ModMpi,      ONLY: MPI_comm_rank
    use ModPlotFile, ONLY: save_plot_file
    integer, intent(in), optional :: iComm

    real:: dLon ! longitude shift
    integer:: nParam
    real:: Param_I(4), IndexMin_I(3), IndexMax_I(3)
    logical:: IsLogIndex_I(3)
    ! Two  variables to store magnetogram time.
    ! Carrigton rootation number:
    integer :: nCR = 0
    ! Fraction of Carrington rotation since its beginning till magnetogram time
    real    :: CRFraction  = 0.0
    type(TableType),  pointer :: Ptr
    real(Real4_), allocatable  :: Tmp_VIIIII(:,:,:,:,:,:)
    real, allocatable :: Magnetogram_VII(:,:,:)
    real :: DX_D(3), Lon, Lat ! in degrees
    real :: XyzRLonLat_DD(3,3)
    integer :: iLon, iLat
    ! To get iProc  and determine  the root PE:
    integer :: iProc, iError
    character(len=*), parameter:: NameSub = 'init_magnetogram_lookup_table'
    !--------------------------------------------------------------------------
    ! Make sure these are set
    iTableB0    = i_lookup_table('B0')
    iTableB0New = i_lookup_table('B0New')
    ! Nothing to do if there is no B0 table defined or to be made
    if(iTableB0 < 0 .and. .not.DoReadHarmonics) RETURN
     ! Get processor index and total number of processors
    if(present(iComm))then
       call MPI_comm_rank(iComm,iProc,iError)
    else
       iProc = 0
    end if
    if(DoReadHarmonics)then
       ! Read harmonics coefficient and Carrington rotation info
       call read_harmonics_file(NameHarmonicsFile, CarringtonRot, dLon)
       DoReadHarmonics = .false.
       nCR        = int(CarringtonRot)
       CRFraction = CarringtonRot  - nCR
       ! Initialize lookup table based on #HARMONICSGRID parameters
       call init_lookup_table( &
            NameTable   = 'B0',                    &
            NameCommand = 'save',                  &
            NameVar     = NameVar,                 &
            NameFile    = 'harmonics_bxyz.out',    &
            TypeFile    = 'real4',                 &
            nIndex_I    = [nR+1, nLon+1, nLat+1],  &
            IndexMin_I  = [rMagnetogram,     0.0, -90.0],    &
            IndexMax_I  = [rSourceSurface, 360.0,  90.0],    &
            Param_I = [rMagnetogram, rSourceSurface, dLon, real(nCR)], &
            Time    = CRFraction,                  &
            StringDescription = 'Created from '//trim(NameHarmonicsFile) )
       iTableB0 = i_lookup_table('B0')

       ! If lookup table is not loaded from a file, make it and save it
       call make_lookup_table_3d(iTableB0, calc_b0_table, iComm)

       call deallocate_harmonics_arrays
       if(iProc==0) then
          Ptr => Table_I(iTableB0)
          allocate(Magnetogram_VII(3, nLon, nLat+1))
          allocate(Tmp_VIIIII(3, nR+1, nLon+1,  nLat+1, 1, 1))
          Tmp_VIIIII = Ptr%Value4_VC
          ! Take the radial slice at the photospheric level,
          ! omit the strip at the longitude of 360 degrees:
          Magnetogram_VII  =  real(Tmp_VIIIII(:,1,1:nLon,:,1,1))
          DX_D = Ptr%DIndex_I
          do iLat  = 1, nLat + 1
             Lat = -90 + DX_D(3)*(iLat - 1)
             do iLon = 1, nLon
                Lon = DX_D(2)*(iLon - 1)
                ! Convert to Br, BLon, BLat components
                XyzRLonLat_DD = rot_xyz_rlonlat(lon=cDegToRad*Lon, &
                     lat= cDegToRad*Lat)
                Magnetogram_VII(:, iLon, iLat) =            &
                     matmul(Magnetogram_VII(:, iLon, iLat), &
                     XyzRLonLat_DD)
             end do
          end do
          ! Remesh from latitudinal nodes to pixels:
          Magnetogram_VII(:,:,1:nLat) = 0.5*(&
               Magnetogram_VII(:,:,1:nLat) + Magnetogram_VII(:,:,2:nLat+1))
          call save_plot_file(NameFile = 'field_2d.out',  &
               nDimIn  = 2,                       &
               TimeIn = CRFraction,               &
               ParamIn_I= [real(int(dLon)), real(nCR)], &
               VarIn_VII= Magnetogram_VII(:,1:nLon,1:nLat), &
               TypeFileIn    = 'ascii',           &
               CoordMinIn_D  = [  0.0 + dLon - int(dLon), -90.0 + 90.0/nLat], &
               CoordMaxIn_D  = [360.0 - 360.0/nLon + dLon - int(dLon),  &
               90.0 - 90.0/nLat ],&
               StringHeaderIn  = 'Created from '//trim(NameHarmonicsFile), &
               NameUnitsIn  = ' [deg] [deg] [Gs] [Gs] [Gs] [deg] []', &
               NameVarIn = 'Longitude Latitude Br BTheta BPhi Long0 CRNumber')
          deallocate(Magnetogram_VII)
          deallocate(Tmp_VIIIII)
          nullify(Ptr)
       end if
    end if
       ! Get coordinate limits, Carrington rotation
    call get_lookup_table(iTableB0, nParam=nParam, Param_I=Param_I, &
         IndexMin_I=IndexMin_I, IndexMax_I=IndexMax_I, Time=CRFraction, &
         IsLogIndex_I=IsLogIndex_I)

    if(IsLogIndex_I(1))then
       rMinB0 = 10**IndexMin_I(1); rMaxB0 = 10**IndexMax_I(1)
    else
       rMinB0 = IndexMin_I(1); rMaxB0 = IndexMax_I(1)
    end if
    LonMinB0 = IndexMin_I(2)*cDegToRad

    ! Rotation matrix for longitude shift if needed
    if(nParam > 2) then
       dLonB0 = (Param_I(3) - dLongitudeHgrDeg)*cDegToRad
       RotB0_DD = rot_matrix_z(dLonB0)
    end if
    ! There only 6 meaning digits in the magnetogram metadata fot the
    ! Carrington time:
    CRFraction = real(nint(CRFraction*1000000))/1000000
    ! Get Carrington rotation time of magnetogram
    if(nParam > 3) CarringtonRot = nint(Param_I(4)) + CRFraction
    !
    ! Second lookup table for a different time for temporal interpolation
    !
    if(DoReadHarmonicsNew)then
       ! Read harmonics coefficients and Carrington rotation info
       call read_harmonics_file(NameHarmonicsFileNew, CarringtonRotNew, dLon)
       DoReadHarmonicsNew = .false.
       nCR        = int(CarringtonRotNew)
       CRFraction = CarringtonRotNew  - nCR
       ! Set up lookup table
       call init_lookup_table( &
            NameTable   = 'B0New',                         &
            NameCommand = 'save',                          &
            NameVar     = NameVar,                         &
            NameFile    = 'harmonics_new_bxyz.out',        &
            TypeFile    = 'real4',                         &
            nIndex_I    = [nR+1, nLon+1, nLat+1],          &
            IndexMin_I  = [rMagnetogram,     0.0, -90.0],  &
            IndexMax_I  = [rSourceSurface, 360.0,  90.0],  &
            Param_I = [rMagnetogram, rSourceSurface, dLon, real(nCR)], &
            Time=CRFraction,                               &
            StringDescription = 'Created from '//trim(NameHarmonicsFileNew))

       iTableB0New = i_lookup_table('B0New')

       ! Make second lookup table using the just read harmonics coefficients
       call make_lookup_table_3d(iTableB0New, calc_b0_table, iComm)

       call deallocate_harmonics_arrays
       if(iProc==0) then
          Ptr => Table_I(iTableB0New)
          allocate(Magnetogram_VII(3, nLon, nLat+1))
          allocate(Tmp_VIIIII(3, nR+1, nLon+1,  nLat+1, 1, 1))
          Tmp_VIIIII = Ptr%Value4_VC
          ! Take the radial slice at the photospheric level,
          ! omit the strip at the longitude of 360 degrees:
          Magnetogram_VII  =  real(Tmp_VIIIII(:,1,1:nLon,:,1,1))
          DX_D = Ptr%DIndex_I
          do iLat  = 1, nLat + 1
             Lat = -90 + DX_D(3)*(iLat - 1)
             do iLon = 1, nLon
                Lon = DX_D(2)*(iLon - 1)
                ! Convert to Br, BLon, BLat components
                XyzRLonLat_DD = rot_xyz_rlonlat(lon=cDegToRad*Lon, &
                     lat= cDegToRad*Lat)
                Magnetogram_VII(:, iLon, iLat) =            &
                     matmul(Magnetogram_VII(:, iLon, iLat), &
                     XyzRLonLat_DD)
             end do
          end do
          ! Remesh from latitudinal nodes to pixels:
          Magnetogram_VII(:,:,1:nLat) = 0.5*(&
               Magnetogram_VII(:,:,1:nLat) + Magnetogram_VII(:,:,2:nLat+1))
          call save_plot_file(NameFile = 'field_2d.out',  &
               nDimIn  = 2,                       &
               TimeIn = CRFraction,               &
               ParamIn_I= [real(int(dLon)), real(nCR)], &
               VarIn_VII= Magnetogram_VII(:,1:nLon,1:nLat), &
               TypeFileIn    = 'ascii',           &
               CoordMinIn_D  = [  0.0 + dLon - int(dLon), -90.0 + 90.0/nLat], &
               CoordMaxIn_D  = [360.0 - 360.0/nLon + dLon - int(dLon),  &
               90.0 - 90.0/nLat ],&
               StringHeaderIn  = 'Created from '//trim(NameHarmonicsFileNew), &
               NameUnitsIn  = ' [deg] [deg] [Gs] [Gs] [Gs] [deg] []', &
               NameVarIn = 'Longitude Latitude Br BTheta BPhi Long0 CRNumber')
          deallocate(Magnetogram_VII)
          deallocate(Tmp_VIIIII)
          nullify(Ptr)
       end if
    end if

    if(iTableB0New > 0)then
       ! Get Carrington rotation (time) for new magnetogram
       call get_lookup_table(iTableB0New, nParam=nParam, Param_I=Param_I,&
            Time=CRFraction)
       ! Rotation matrix for longitude shift if needed
       if(nParam > 2) then
          dLonB0New = (Param_I(3) - dLongitudeHgrDeg)*cDegToRad
          RotB0New_DD = rot_matrix_z(dLonB0New)
       end if

       if(nParam > 3)then
          ! Get Carrington rotation time of magnetogram
          CarringtonRotNew = nint(Param_I(4)) + CRFraction
       else
          call CON_stop(NameSub//': missing Carrington rotation info in '// &
               trim(NameHarmonicsFileNew))
       end if
    end if

  end subroutine init_magnetogram_lookup_table
  !============================================================================
  subroutine calc_b0_table(iTable, r, Lon, Lat, b_D)

    ! Calculate B0 at the location

    integer, intent(in):: iTable
    real, intent(in)   :: r, Lon, Lat
    real, intent(out)  :: b_D(:)

    real:: Theta, Phi, Bsph_D(3)
    ! Rotation matrix to convert Br, BTheta, BPhi to  Bx, By, Bz.
    real :: XyzSph_DD(3,3)
    character(len=*), parameter:: NameSub = 'calc_b0_table'
    !--------------------------------------------------------------------------
    Phi   = cDegToRad*Lon
    Theta = cDegToRad*(90 - Lat)

    call get_harmonics_field(r, Theta, Phi, Bsph_D)

    ! Convert to Cartesian components
    XyzSph_DD = rot_xyz_sph(Theta, Phi)

    b_D = matmul(XyzSph_DD, Bsph_D)

  end subroutine calc_b0_table
  !============================================================================
  subroutine get_magnetogram_field31(x, y, z, B0_D, Carrington)

    real, intent(in) ::  x, y, z
    real, intent(out):: B0_D(3)
    real, intent(in), optional :: Carrington

    !--------------------------------------------------------------------------
    call  get_magnetogram_field11([x, y, z], B0_D, Carrington)

  end subroutine get_magnetogram_field31
  !============================================================================
  subroutine get_magnetogram_field11(Xyz_D, B0_D, Carrington)

    ! Return B0_D [Tesla] field at position Xyz_D [Rs]
    ! Interpolat to time Carrington if present and iTableB0New is defined

    real, intent(in) :: Xyz_D(3)
    real, intent(out):: B0_D(3)
    real, intent(in), optional :: Carrington

    real:: rLonLat_D(3), r, B0New_D(3)

    ! Converting to rlonlat (radians)
    character(len=*), parameter:: NameSub = 'get_magnetogram_field11'
    !--------------------------------------------------------------------------
    call xyz_to_rlonlat(Xyz_D, rLonLat_D)

    ! Include the shift in Phi coordinate and make sure that it is
    ! in the range provided by the lookup table
    if(dLonB0 /= 0.0 .or. LonMinB0 /= 0.0) rLonLat_D(2) = &
         modulo(rLonLat_D(2) - dLonB0 - LonMinB0, cTwoPi) + LonMinB0

    ! Lookup table uses degrees
    rLonLat_D(2:3) = cRadToDeg*rLonLat_D(2:3)

    ! Extrapolate for r < rMinB0
    r = rLonLat_D(1)

    call interpolate_lookup_table(iTableB0, rLonLat_D, B0_D, &
         DoExtrapolate=(r<rMinB0) )

    ! Rotate Bx, By based on shifted coordinates
    if(dLonB0 /= 0.0) B0_D = matmul(RotB0_DD, B0_D)

    if(present(Carrington) .and. iTableB0New > 0)then

       if(CarringtonRot < 0.0 .or. CarringtonRotNew < 0.0) call CON_stop( &
            NameSub//': no Carrington time in at least one magnetogram')
       call xyz_to_rlonlat(Xyz_D, rLonLat_D)

       ! Include the shift in Phi coordinate and make sure that it is
       ! in the range provided by the lookup table
       if(dLonB0New /= 0.0 .or. LonMinB0 /= 0.0) rLonLat_D(2) = &
         modulo(rLonLat_D(2) - dLonB0New - LonMinB0, cTwoPi) + LonMinB0

       ! Lookup table uses degrees
       rLonLat_D(2:3) = cRadToDeg*rLonLat_D(2:3)

       if(CarringtonRot < CarringtonRotNew)then
          ! Check if time has passed the original magnetogram time
          if(Carrington > CarringtonRot)then
             call interpolate_lookup_table(iTableB0New, rLonLat_D, B0New_D, &
                  DoExtrapolate=(r<rMinB0) )
             ! Rotate Bx, By based on shifted coordinates
             if(dLonB0New /= 0.0) B0New_D = matmul(RotB0New_DD, B0New_D)
             if(Carrington > CarringtonRotNew)then
                ! Time is beyond new table time, so take that
                B0_D = B0New_D
             else
                ! interpolate in time
                B0_D = ((Carrington       - CarringtonRot)*B0New_D  &
                     +  (CarringtonRotNew - Carrington   )*B0_D   ) &
                     /  (CarringtonRotNew - CarringtonRot)
             end if
          end if
       else
          ! "new" table has earlier time than original table
          if(Carrington < CarringtonRot)then
             call interpolate_lookup_table(iTableB0New, rLonLat_D, B0New_D, &
                  DoExtrapolate=(r<rMinB0) )
             ! Rotate Bx, By based on shifted coordinates
             if(dLonB0New /= 0.0) B0New_D = matmul(RotB0New_DD, B0New_D)
             if(Carrington <= CarringtonRotNew)then
                ! Time is before new table time, so take that
                B0_D = B0New_D
             else
                ! interpolate in time
                B0_D = ((CarringtonRot - Carrington      )*B0New_D  &
                     +  (Carrington    - CarringtonRotNew)*B0_D   ) &
                     /  (CarringtonRot - CarringtonRotNew)
             end if
          end if
       end if
    end if

    ! Scale with r^2 for r > rMaxB0
    if(r > rMaxB0) B0_D = (rMaxB0/r)**2 * B0_D

    ! Multiply with B0 factor and convert from Gauss to Tesla
    B0_D = B0_D*FactorB0*1e-4

  end subroutine get_magnetogram_field11
  !============================================================================
  subroutine read_magnetogram_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    real:: Height ! pointless variable
    character(len=*), parameter:: NameSub = 'read_magnetogram_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#FACTORB0")
       call read_var("FactorB0", FactorB0)

    case("#HARMONICSGRID")
       call read_var('rMagnetogram',   rMagnetogram)
       call read_var('rSourceSurface', rSourceSurface)
       call read_var('IsLogRadius',    IsLogRadius)
       call read_var('MaxOrder',       MaxOrder)
       call read_var('nR',             nR)
       call read_var('nLon',           nLon)
       call read_var('nLat',           nLat)
       if(IsLogRadius)then
          NameVar = 'log'//NameVarLinear
       else
          NameVar = NameVarLinear
       end if

    case("#HARMONICSFILE")
       DoReadHarmonics = .true.
       call read_var('NameHarmonicsFile', NameHarmonicsFile)

    case("#NEWHARMONICSFILE")
       DoReadHarmonicsNew = .true.
       call read_var('NameHarmonicsFileNew', NameHarmonicsFileNew)

    case("#MAGNETOGRAM") ! Kept For backward compatibility only
       call read_var('UseMagnetogram', DoReadHarmonics)
       if(DoReadHarmonics)then
          call read_var('rMagnetogram',        rMagnetogram)
          call read_var('rSourceSurface',      rSourceSurface)
          call read_var('HeightInnerBc',       Height)
          call read_var('NameMagnetogramFile', NameHarmonicsFile)
       end if

    case default
       call CON_stop(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine read_magnetogram_param
  !============================================================================

  subroutine read_harmonics_file(NameFile, Carrington, dLon)

    use ModPlotFile, ONLY: read_plot_file

    character(len=*), intent(in) :: NameFile
    real,             intent(out):: Carrington ! Carrington rotation from file
    real,             intent(out):: dLon       ! longitude offset from file

    ! read properly formatted harmonics coefficient file

    integer:: nHarmonic ! number of harmonics = (nOrder+1)*(nOrder+2)/2
    integer:: nParam, nOrderIn, i, n, m

    real:: Param_I(3), Coef, Coef1, CRFraction
    real, allocatable:: Var_VI(:,:), Coord_DII(:,:,:)

    character(len=*), parameter:: NameSub = 'read_harmonics_file'
    !--------------------------------------------------------------------------

    call read_plot_file(NameFile, TimeOut = CRFraction, &
         nParamOut=nParam, ParamOut_I=Param_I, n1Out=nHarmonic)

    if(nParam < 3) call CON_stop(NameSub// &
         ': not enough parameters in '//trim(NameFile))

    nOrderIn   = nint(Param_I(1))
    nOrder     = min(MaxOrder, nOrderIn)
    Carrington = Param_I(2) + CRFraction
    dLon       = Param_I(3)

    if(nHarmonic /= (nOrderIn+1)*(nOrderIn+2)/2) call CON_stop(NameSub// &
         ': inconsistent n and nOrder in '//trim(NameFile))

    call allocate_harmonics_arrays

    allocate(Var_VI(2,nHarmonic), Coord_DII(2,nHarmonic,1))
    call read_plot_file(NameFile, VarOut_VI=Var_VI, CoordOut_DII=Coord_DII)

    do i = 1, nHarmonic
       n = nint(Coord_DII(1,i,1))
       if (n > nOrder) EXIT  ! the harmonics are listed in increasing order
       m = nint(Coord_DII(2,i,1))
       g_II(n,m) = Var_VI(1,i)
       h_II(n,m) = Var_VI(2,i)
    end do

    deallocate(Var_VI, Coord_DII)

    ! Normalize coefficients
    Coef = rSourceSurface
    do n = 0, nOrder
       Coef  = Coef/rSourceSurface**2
       Coef1 = 1.0/(n+1 + n*Coef)
       g_II(n,0:n) = g_II(n,0:n)*Coef1
       h_II(n,0:n) = h_II(n,0:n)*Coef1
    enddo
    ! Leave out monopole (n=m=0) term::
    g_II(0,0) = 0.0

  end subroutine read_harmonics_file
  !============================================================================

  subroutine allocate_harmonics_arrays

    ! Calculate square roots and their ratios used in the spherical harmonics
    ! functions up to order nOrder for speeding up the calculations

    integer:: m, MaxInt
    !--------------------------------------------------------------------------
    if(allocated(g_II)) RETURN

    ! Spherical harmonics coefficients
    allocate(g_II(0:nOrder,0:nOrder), h_II(0:nOrder,0:nOrder))
    g_II = 0.0; h_II = 0.0

    ! radial functions
    allocate(rRsPower_I(-1:nOrder+2), &
         RmRsPower_I(0:nOrder+2), RmRPower_I(0:nOrder+2))

    ! Legendre polynomlials
    allocate(p_II(0:nOrder,0:nOrder), Dp_II(0:nOrder,0:nOrder))

    ! Azimuthal functions
    allocate(SinPhi_I(0:nOrder), CosPhi_I(0:nOrder))

    MaxInt = max(nOrder**2, 5*nOrder, 10)
    allocate(Sqrt_I(MaxInt), SqrtRatio_I(nOrder+1))

    do m = 1, MaxInt
       Sqrt_I(m) = sqrt(real(m))
    end do

    ! Calculate the ratio sqrt(2m!)/(2^m*m!) recursively
    SqrtRatio_I(1) = 1.0
    do m = 1, nOrder
       SqrtRatio_I(m+1) = SqrtRatio_I(m)*Sqrt_I(2*m-1)/Sqrt_I(2*m)
    enddo

  end subroutine allocate_harmonics_arrays
  !============================================================================
  subroutine deallocate_harmonics_arrays
    !--------------------------------------------------------------------------
    if(.not.allocated(g_II)) RETURN

    deallocate(g_II, h_II, Sqrt_I, SqrtRatio_I, &
         rRsPower_I, RmRsPower_I, RmRPower_I, &
         SinPhi_I, CosPhi_I, p_II, Dp_II)

  end subroutine deallocate_harmonics_arrays
  !============================================================================
  subroutine get_harmonics_field(r, Theta, Phi, Bsph_D)

    ! Calculate the harmonics based potential field Bsph_D at r, Theta, Phi

    real, intent(in):: r, Theta, Phi ! r in Rs, Theta and Phi in radians
    real, intent(out):: Bsph_D(3)    ! Bsph_D with r, Theta, Phi components

    integer:: n, m
    real:: CosTheta, SinTheta

    ! The spherical components of the magnetic field
    real:: Br, Btheta, Bphi

    real:: Coef1, Coef2, Coef3, Coef4
    !--------------------------------------------------------------------------

    ! Calculate the radial part of spherical functions
    call calc_radial_functions(r)

    ! Calculate the set of Legendre polynoms for given Theta
    CosTheta = cos(Theta)
    SinTheta = max(sin(Theta), 1E-10)
    call calc_legendre_polynomial(SinTheta, CosTheta)

    ! Calculate azimuthal functions for given Phi
    call calc_azimuthal_functions(Phi)

    ! Initialize the components
    Br = 0.0; Btheta = 0.0; Bphi = 0.0

    ! Calculate B from spherical harmonics
    do m = 0, nOrder; do n = m, nOrder

       Coef1 = (n+1)*RmRPower_I(n+2) + RmRsPower_I(n+2)*n*rRsPower_I(n-1)
       Coef3 =       RmRPower_I(n+2) - RmRsPower_I(n+2)*rRsPower_I(n-1)
       Coef2 = g_II(n,m)*CosPhi_I(m) + h_II(n,m)*SinPhi_I(m)
       Coef4 = g_II(n,m)*SinPhi_I(m) - h_II(n,m)*CosPhi_I(m)

       ! Br = -d(Psi)/dR
       Br  = Br + p_II(n,m)*Coef1*Coef2

       ! Bt = -(1/r)*d(Psi)/dTheta
       Btheta  = Btheta - Dp_II(n,m)*Coef2*Coef3

       ! Bp = -(1/r)*d(Psi)/dPhi
       Bphi  = Bphi + p_II(n,m)*m/SinTheta*Coef3*Coef4

       ! Potential could be calculated if it was needed:
       ! Potential = Potential + r*p_II(n,m)*Coef2*Coef3

    enddo; enddo

    Bsph_D = [Br, Btheta, Bphi]

  contains
    !==========================================================================

    subroutine calc_legendre_polynomial(SinTheta, CosTheta)

      ! Calculate Legendre polynomials p_II and its derivative Dp_II
      ! with appropriate normalization for Theta
      ! Equation numbers refer to Altschuler et al. 1976

      real, intent(in):: SinTheta, CosTheta

      real:: SinThetaM, SinThetaM1  ! sin(Theta)^m, sin(Theta)^(m-1)
      integer:: m

      ! Cache previous values
      real:: SinThetaOld = -10.0, CosThetaOld = -10.0
      !------------------------------------------------------------------------
      if(SinTheta == SinThetaOld .and. CosTheta == CosThetaOld) RETURN

      SinThetaOld = SinTheta; CosThetaOld = CosTheta

      SinThetaM  = 1.0
      SinThetaM1 = 1.0
      p_II  = 0.0
      Dp_II = 0.0

      do m = 0, nOrder
         if (m == 0) then
            Coef1 = Sqrt_I(2*m+1)
         else
            Coef1 = Sqrt_I(2*(2*m+1))
         endif
         ! Eq.(27)
         p_II(m,m) = SqrtRatio_I(m+1)*Coef1* SinThetaM
         ! Eq.(28)
         if (m < nOrder) p_II(m+1,m) = p_II(m,m)*Sqrt_I(2*m+3)*CosTheta
         ! Eq.(30)
         Dp_II(m,m) = SqrtRatio_I(m+1)*Coef1*m*CosTheta*SinThetaM1
         ! Eq.(31)
         if (m < nOrder) &
              Dp_II(m+1,m) = Sqrt_I(2*m+3)* &
              (CosTheta*Dp_II(m,m) - SinTheta*p_II(m,m))

         ! Increase the powers
         SinThetaM1 = SinThetaM
         SinThetaM  = SinThetaM*SinTheta

      enddo

      ! Recursive rules
      do m = 0, nOrder-2; do n = m+2, nOrder
         ! Eq.(29)
         Coef1 = Sqrt_I(2*n+1)/Sqrt_I(n**2-m**2)
         Coef2 = Sqrt_I(2*n-1)
         Coef3 = Sqrt_I((n-1)**2-m**2)/Sqrt_I(2*n-3)

         p_II(n,m) = Coef1*(Coef2*CosTheta*p_II(n-1,m) - Coef3*p_II(n-2,m))

         ! Eq.(32)
         Dp_II(n,m) = Coef1*(Coef2*(CosTheta*Dp_II(n-1,m) &
              - SinTheta*p_II(n-1,m)) - Coef3*Dp_II(n-2,m))
      enddo; enddo

      ! Apply Schmidt normalization
      do m = 0, nOrder; do n = m, nOrder
         ! Eq.(33)
         Coef1 = 1.0/Sqrt_I(2*n+1)
         ! Eq.(34)
         p_II(n,m)  = p_II(n,m)*Coef1
         Dp_II(n,m) = Dp_II(n,m)*Coef1
      enddo; enddo

    end subroutine calc_legendre_polynomial
    !==========================================================================

    subroutine calc_radial_functions(r)

      real, intent(in):: r

      ! Calculate powers of the ratios of radii up to nOrder

      integer:: m
      real:: RmRs, RmR, rRs

      real:: rOld = -10.0
      !------------------------------------------------------------------------
      if(rOld == r) RETURN
      rOld = r

      RmRs = rMagnetogram/rSourceSurface
      RmR  = rMagnetogram/r
      rRs  = r/rSourceSurface

      ! Zero and -1 powers
      rRsPower_I(-1) = 1.0/rRs
      rRsPower_I(0)  = 1.0
      RmRsPower_I(0) = 1.0
      RmRPower_I(0)  = 1.0

      ! Recursive: x^m = x^m-1 * x
      do m = 1, nOrder+2
         RmRsPower_I(m) = RmRsPower_I(m-1) * RmRs
         RmRPower_I(m)  = RmRPower_I(m-1)  * RmR
         rRsPower_I(m)  = rRsPower_I(m-1)  * rRs
      end do

    end subroutine calc_radial_functions
    !==========================================================================
    subroutine calc_azimuthal_functions(Phi)

      ! Calculate azimuthal harmonics for given Phi

      real,    intent(in):: Phi

      integer:: m
      complex:: z, zM ! Powers of cos(Phi)+i*sin(Phi)

      real   :: PhiOld = -10.0
      !------------------------------------------------------------------------
      if(Phiold == Phi) RETURN
      PhiOld = Phi

      z  = exp( cmplx(0.0, Phi) )
      zM = z

      CosPhi_I(0) = 1.0
      SinPhi_I(0) = 0.0
      CosPhi_I(1) = real(zM)
      SinPhi_I(1) = aimag(zM)

      do m = 2, nOrder
         zM = zM*z
         CosPhi_I(m) = real(zM)
         SinPhi_I(m) = aimag(zM)
      end do

    end subroutine calc_azimuthal_functions
    !==========================================================================
  end subroutine get_harmonics_field
  !============================================================================

end module ModMagnetogram
!==============================================================================
