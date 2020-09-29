!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
program magnetogram
  use EEE_ModMain
  use EEE_ModCommonVariables
  use ModConst, ONLY: cPi, cTwoPi, cDegToRad, cRadToDeg
  use ModPlotFile, ONLY:save_plot_file
  use ModReadParam, ONLY: read_file, read_init, &
       read_line, read_command
  use ModUtilities, ONLY: CON_stop
  use ModCoordTransform, ONLY: rlonlat_to_xyz, xyz_to_rlonlat
  use ModMpi

  implicit none


  integer,parameter:: nLong = 360, nLat = 180, nR = 70, nD = 40
  integer:: Long0, i, j, k, iError, iLong
  real:: Radius, Longitude, Latitude, dR
  real:: Radius_I(nR), Longitude_I(nLong), Latitude_I(nLat), Long_I(nD), Lat_I(nD)
  real:: PlotVar_V(3,nR,nD,nD) = 0
  real:: BSurface_DC(3,nLong,nLat) = 0, B_DC(3,nR,nD,nD) = 0
  real::Xyz_D(3),SinLong,CosLong,SinLat,CosLat, Rho, p, U_D(3), B_D(3)
  character(LEN=3)::TypeLatAxis
  character(LEN=100):: StringLine, NameCommand
  logical:: DoDebug = .false.
  !--------------------------------------------------------------------------

  call MPI_init(iError)

  Io2Si_V = 1; Si2Io_V = 1; Io2No_V = 1 
  No2Io_V = 1; Si2No_V = 1; No2Si_V = 1

  read(*,*)Long0
  write(*,'(a,i3)')prefix,Long0
  read(*,*)TypeLatAxis
  write(*,'(a)')prefix//TypeLatAxis

  write(*,*)'Reading CME.in'
  call read_file('CME.in')
  call read_init
  READPARAM: do
     if(.not.read_line(StringLine) ) EXIT READPARAM
     if(.not.read_command(NameCommand)) CYCLE READPARAM
     select case(NameCommand)
     case('#END')
        EXIT READPARAM
     case('#CME')
        write(*,*)'Reading CME Para'
        call EEE_set_parameters(NameCommand)
     case default
        call CON_stop('Unknown command:'//NameCommand)
     end select
  end do READPARAM

  write(*,*)'LongitudeCme, LatitudeCme =',LongitudeCme, LatitudeCme
  dR = 1.0/nR
  do i = 0, nR-1
     Radius_I(i+1) = 1 + i*dR
  enddo
  do i = 1, nLong
     ! Include the Phishift to correctly place the CME
     iLong = i + Long0
     iLong = 1 + modulo(iLong - 1, nLong)
     if (i <= nD) then
        Long_I(i) = iLong + LongitudeCme - nD/2
     endif
     Longitude_I(i) = iLong - 0.5
  end do
  do i = 1, nLat
     select case(TypeLatAxis)
     case('uni')
        if (i <= nD) then
           Lat_I(i) = i + LatitudeCme - nD/2
        endif
        Latitude_I(i) = i - 90.5
     case('sin')
        SinLat = asin((2.0*i - nLat - 1.0)/nLat)
        if (i <= nD) then
           Lat_I(i) = SinLat*cRadToDeg + LatitudeCme - nD/2 + 90.5
        endif
        Latitude_I(i) = SinLat*cRadToDeg
     case default
     end select
  end do
 
  if (DoDebug)then 
     write(*,*)'Radius_I    = ',Radius_I
     write(*,*)'Long_I      =',Long_I
     write(*,*)'Lat_I       =',Lat_I
     write(*,*)'Longitude_I =',Longitude_I
     write(*,*)'Latitude_I  =',Latitude_I     
  endif

  do k = 1,nD
     Latitude = Lat_I(k) * cDegToRad
     SinLat = sin(Latitude)
     CosLat = cos(Latitude)
     do j = 1,nD
        Longitude = Long_I(j) * cDegToRad
        SinLong = sin(Longitude)
        CosLong = cos(Longitude)
        do i = 1, nR
           Radius = Radius_I(i)
	   call rlonlat_to_xyz(Radius,Longitude,Latitude, Xyz_D)
           call EEE_get_state_BC(Xyz_D,Rho,U_D,B_D,p,0.0,0,0)
           B_DC(1,i,j,k) = sum(Xyz_D*B_D)
           B_DC(2,i,j,k) = B_D(y_)*CosLong - B_D(x_)*SinLong
           B_DC(3,i,j,k) = (B_D(x_)*CosLong + B_D(y_)*SinLong)*(-SinLat) &
                + B_D(z_)*CosLat
        end do
     enddo
  end do
  
  PlotVar_V(1:3,:,:,:) = B_DC

  write(*,*)'Saving 3d Flux Rope output file'
  call save_plot_file(&
       NameFile = 'FluxRope.out',&
       TypeFileIn = 'ascii',          &
       StringHeaderIn = '3D Flux rope magnetic field +/- 20 deg of Flux rope location',&
       nDimIn = 3,                    &
       ParamIn_I = (/float(Long0)/),         &
       Coord1In_I = Radius_I, &
       Coord2In_I = Longitude_I,       &
       Coord3In_I = Latitude_I,        &
       NameVarIn  = 'Radius Longitude Latitude Br BLon BLat Long0',&
       VarIn_VIII = PlotVar_V(1:3,:,:,:))

  do j = 1,nLat
     Latitude = Latitude_I(j) * cDegToRad
     SinLat = sin(Latitude)
     CosLat = cos(Latitude)
     do i = 1,nLong
        Longitude = Longitude_I(i) * cDegToRad
        SinLong = sin(Longitude)
        CosLong = cos(Longitude)
        call rlonlat_to_xyz(1.0,Longitude,Latitude, Xyz_D)
        call EEE_get_state_BC(Xyz_D,Rho,U_D,B_D,p,0.0,0,0)
        BSurface_DC(1,i,j) = sum(Xyz_D*B_D)
        BSurface_DC(2,i,j) = B_D(y_)*CosLong - B_D(x_)*SinLong
        BSurface_DC(3,i,j) = (B_D(x_)*CosLong + B_D(y_)*SinLong)*(-SinLat) &
                + B_D(z_)*CosLat
     end do
  enddo

    write(*,*)'Saving 2d Flux Rope output file'
    call save_plot_file(&
         NameFile = 'FRMagnetogram.out',&
         TypeFileIn = 'ascii',          &
         StringHeaderIn = 'Flux rope magnetic field at 1 R_s',&
         nDimIn = 2,                    &
         ParamIn_I = (/float(Long0)/),         &
         CoordIn_I = Longitude_I,       &
         Coord2In_I= Latitude_I,        &
         NameVarIn = 'Longitude Latitude Br BLon BLat Long0',& 
         VarIn_VII = BSurface_DC)

  call MPI_finalize(iError)

end program magnetogram
!============================================================================
