!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program magnetogram
  use EEE_ModMain, ONLY:  EEE_set_parameters, EEE_get_state_BC, &
       EEE_set_plot_range
  use EEE_ModCommonVariables, ONLY:   Io2Si_V, Si2Io_V, Io2No_V, &
       No2Io_V, Si2No_V, No2Si_V, prefix, x_, y_, z_,            &
       LongitudeCme, LatitudeCme, OrientationCME
  use ModConst, ONLY: cPi, cTwoPi, cDegToRad, cRadToDeg
  use ModPlotFile, ONLY:save_plot_file
  use ModReadParam, ONLY: read_file, read_init, &
       read_line, read_command
  use ModUtilities, ONLY: CON_stop
  use ModCoordTransform, ONLY: rot_xyz_rlonlat, rlonlat_to_xyz, &
       rot_xyz_mercator, rot_matrix_z
  use ModMpi, ONLY: MPI_COMM_SELF
  use ModMagnetogram, ONLY: iTableB0, init_magnetogram_lookup_table, &
       get_magnetogram_field
  use ModLookupTable, ONLY: read_lookup_table_param, get_lookup_table

  implicit none
  ! Number of indexes for magnetogram map
  integer,parameter:: nLong = 360, nLat = 180
  ! Number of indexes for 3D plot
  integer           :: nXY = -1, nZ = -1
  ! Plot variables:
  real, allocatable :: Var_VN(:,:,:,:)
  ! Named indexes, maay be redifined depending on the use of potential and/or
  ! EE field
  integer:: B0x_ = 1, B0y_ = 2, B0z_ = 3, B1x_ = 1, B1y_ = 2, B1z_ = 3
  ! Rotational matrix of coordinate transformation:
  real              :: Rotate_DD(3,3)
  ! Loop indexes    ::
  integer           :: i, j, k
  integer           :: Long0, iLong
  real              :: Longitude, Latitude
  real              :: Longitude_I(nLong), Latitude_I(nLat)
  real              :: BSurface_DC(3,nLong,nLat) = 0
  real              :: XyzRLonLat_DD(3,3)
  real              :: Xyz_D(3), Rho, p, U_D(3), B_D(3)
  character(LEN=3)  :: TypeLatAxis
  character(LEN=100):: StringLine, NameCommand
  !----------------------------------------------------------------------------

  Io2Si_V = 1; Si2Io_V = 1; Io2No_V = 1
  No2Io_V = 1; Si2No_V = 1; No2Si_V = 1

  read(*,*)Long0
  write(*,'(a,i3)')prefix,Long0
  read(*,*)TypeLatAxis
  write(*,'(a)')prefix//TypeLatAxis

  write(*,'(a)')prefix//'Reading CME.in'
  call read_file('CME.in', iCommIn = MPI_COMM_SELF)
  call read_init
  READPARAM: do
     if(.not.read_line(StringLine) ) EXIT READPARAM
     if(.not.read_command(NameCommand)) CYCLE READPARAM
     select case(NameCommand)
     case('#END')
        EXIT READPARAM
     case('#CME')
        write(*,*)'Reading CME Parameters'
        call EEE_set_parameters(NameCommand)
        write(*,'(a,f10.2,a)')'LongitudeCme = ', LongitudeCme, ' [deg]'
        write(*,'(a,f10.2,a)')'LatitudeCme = ', LatitudeCme, ' [deg]'
        write(*,'(a,f10.2,a)')'OrientationCme =', OrientationCme, ' [deg]'
     case('#LOOKUPTABLE')
        call read_lookup_table_param
        call init_magnetogram_lookup_table()
     case default
        call CON_stop('Unknown command:'//NameCommand)
     end select
  end do READPARAM

  do i = 1, nLong
     ! Include the Phishift to correctly place the CME
     iLong = i + Long0
     iLong = 1 + modulo(iLong - 1, nLong)
     Longitude_I(i) = iLong - 0.5
  end do
  do i = 1, nLat
     select case(TypeLatAxis)
     case('uni')
        Latitude_I(i) = i - 90.5
     case('sin')
        Latitude_I(i) = asin((2.0*i - nLat - 1.0)/nLat)*cRadToDeg
     case default
     end select
  end do

  do j = 1,nLat
     Latitude = Latitude_I(j) * cDegToRad
     do i = 1,nLong
        Longitude = Longitude_I(i)*cDegToRad
        call rlonlat_to_xyz(1.0,Longitude,Latitude, Xyz_D)
        call EEE_get_state_BC(Xyz_D, Rho, U_D, B_D, p, &
             Time = 0.0, nStep  = 0, Iteration_Number=0)
        ! Convert to Br, BLon, BLat components
        XyzRLonLat_DD = rot_xyz_rlonlat(lon=Longitude, &
             lat=Latitude)
        BSurface_DC(:,i,j) = matmul(B_D, XyzRLonLat_DD)
     end do
  enddo

  write(*,*)'Saving 2d Flux Rope output file'
  call save_plot_file(&
       NameFile = 'FRMagnetogram.out',&
       TypeFileIn = 'ascii',          &
       StringHeaderIn = 'Flux rope magnetic field at 1 R_s',&
       nDimIn = 2,                    &
       ParamIn_I = [real(Long0)],     &
       CoordIn_I = Longitude_I,       &
       Coord2In_I= Latitude_I,        &
       NameVarIn = 'Longitude Latitude Br BLon BLat Long0',&
       VarIn_VII = BSurface_DC)

  call EEE_set_plot_range(nXY, nZ)
  if(nXY < 1 .and. nZ < 1)STOP
  write(*,*)'nXY, nZ = ', nXY, nZ
  if(iTableB0 > 0)then
     allocate(Var_VN(6, 2*nXY+1, 2*nXY+1, nZ+1))
  else
     allocate(Var_VN(3, 2*nXY+1, 2*nXY+1, nZ+1))
  end if
end program magnetogram
!==============================================================================
