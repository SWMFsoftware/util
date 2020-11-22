!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModMagnetogram

  use ModNumConst,       ONLY: cTwoPi, cRadToDeg, cDegToRad
  use ModIoUnit,         ONLY: io_unit_new
  use CON_axes,          ONLY: dLongitudeHgrDeg
  use ModUtilities,      ONLY: CON_stop
  use ModCoordTransform, ONLY: rot_xyz_sph, rot_matrix_z, xyz_to_rlonlat
  use ModLookupTable, ONLY: &
       i_lookup_table, init_lookup_table, make_lookup_table_3d, &
       get_lookup_table, interpolate_lookup_table

  implicit none
  save

  private

  logical, public :: UseMagnetogram=.false., UseNewMagnetogram = .false.

  integer, public:: iTableB0 = -1 ! index of B0 lookup table

  public :: read_magnetogram_param ! read parameters

  public:: init_magnetogram_lookup_table ! initialize lookup table

  public:: get_magnetogram_field ! field from lookup table at a given point

  public:: interpolate_field ! good old legacy left handed coordinates

  real:: rSourceSurface=2.5, rMagnetogram=1.0, Height=0.0

  ! Units of the magnetic field in the file including corrections
  ! relative to the magnetic units in the SC [Gauss]
  real, public :: UnitB=1.0, UnitBNew=1.0

  public :: read_magnetogram_file, read_new_magnetogram_file

  ! Read H(eliographic) L(ongitude) of the C(entral) M(eridian) of
  ! the M(ap) from the file header. Assign PhiOffset = LonCentral-180
  public :: get_hlcmm

  ! Rotation angle around z-axis, in degrees,
  ! from the coordinate system of the component
  ! towards the coordinate system of the magnetic
  ! map in the positive direction - hence, is positive.
  ! Is set automatically to be equal to the
  ! H(eliographic) L(ongitude) of C(entral) M(eridian)
  ! of the M(ap) minus 180 Deg, if in the input file
  ! PhiOffset is negative.
  real, public :: PhiOffset=-1.0

  ! Local variables -------------

  real:: NewPhiShift = -1.0

  ! spherical harmonics related control parameters

  ! Maximum order of spherical harmonics
  integer, parameter:: nHarmonicsMax=180 ! 90

  ! Weights of the spherical harmonics
  real, dimension(nHarmonicsMax+1,nHarmonicsMax+1):: g_II, h_II

  integer:: nOrder = nHarmonicsMax

  ! Number of header lines in the file
  integer:: nHeadLine=12

  ! Name of the input file
  character (LEN=32):: NameMagnetogram='mf.dat', NameMagnetogram2 = 'newmf.dat'

  integer, parameter:: MaxInt=100000
  real, allocatable:: Sqrt_I(:), SqrtRatio_I(:)

  ! Lookup table related variables
  real:: rMinB0=1.0, rMaxB0=30.0, dLonB0=0.0, FactorB0=1.0, RotB0_DD(3,3)
  real:: LonMinB0 = 0.0

  interface get_magnetogram_field
     module procedure get_magnetogram_field11, get_magnetogram_field31
  end interface get_magnetogram_field

contains
  !============================================================================
  subroutine init_magnetogram_lookup_table(iComm)

    use ModNumConst, ONLY: cDegToRad

    integer, intent(in):: iComm

    integer:: nParam
    real:: Param_I(4), IndexMin_I(3), IndexMax_I(3)
    logical:: IsLogIndex_I(3)
    integer:: nR, nLon, nLat
    !--------------------------------------------------------------------------

    ! Read or create lookup table and get the longitude shift.
    ! Create the rotation matrix based on the shift
    iTableB0 = i_lookup_table('B0')
    if(iTableB0 <= 0)then
       if(.not.UseMagnetogram) RETURN
       ! Set up a default B0 lookup table
       nR = max(30, nOrder); nLon = max(72, nOrder); nLat=max(30, nOrder)
       call init_lookup_table( &
            NameTable   = 'B0',                    &
            NameCommand = 'use',                   &
            NameVar     = 'r lon lat bx by bz',    &
            NameFile    = 'harmonics_bxyz.out',    &
            TypeFile    = 'real4',                 &
            nIndex_I    = [nR+1, nLon+1, nLat+1],  &
            IndexMin_I  = [ 1.0,   0.0, -90.0],    &
            IndexMax_I  = [ 2.5, 360.0,  90.0],    &
            StringDescription = 'Created from '//trim(NameMagnetogram) )
       iTableB0 = i_lookup_table('B0')
    end if

    ! If lookup table is not loaded, make it (and save it)
    call make_lookup_table_3d(iTableB0, calc_b0_table, iComm)

    ! Get coordinate limits
    call get_lookup_table(iTableB0, nParam=nParam, Param_I=Param_I, &
         IndexMin_I=IndexMin_I, IndexMax_I=IndexMax_I, &
         IsLogIndex_I=IsLogIndex_I)

    if(IsLogIndex_I(1))then
       rMinB0 = 10**IndexMin_I(1); rMaxB0 = 10**IndexMax_I(1)
    else
       rMinB0 = IndexMin_I(1); rMaxB0 = IndexMax_I(1)
    end if
    LonMinB0 = IndexMin_I(2)*cDegToRad
    if(nParam > 2) then
       dLonB0 = (Param_I(3) - dLongitudeHgrDeg)*cDegToRad
       RotB0_DD = rot_matrix_z(dLonB0)
    end if

  end subroutine init_magnetogram_lookup_table
  !============================================================================
  subroutine calc_b0_table(iTable, Arg1, Arg2, Arg3, b_D)

    ! Calculate B0 at the location

    integer, intent(in):: iTable
    real, intent(in)   :: Arg1, Arg2, Arg3
    real, intent(out)  :: b_D(:)

    real:: r, Theta, Phi, Bsph_D(3), XyzSph_DD(3,3)

    character(len=*), parameter:: NameSub = 'calc_b0_table'
    !--------------------------------------------------------------------------
    r     = Arg1
    Phi   = Arg2*cDegToRad
    Theta = (90-Arg3)*cDegToRad

    call get_harmonics_field(r, Theta, Phi, Bsph_D)

    ! Convert to Cartesian
    XyzSph_DD = rot_xyz_sph(Theta, Phi)

    b_D = matmul(XyzSph_DD, Bsph_D)

  end subroutine calc_b0_table
  !============================================================================
  subroutine get_magnetogram_field31(x, y, z, B0_D)

    real, intent(in) ::  x, y, z
    real, intent(out):: B0_D(3)

    !--------------------------------------------------------------------------
    call  get_magnetogram_field11([ x, y, z ], B0_D)

  end subroutine get_magnetogram_field31
  !============================================================================
  subroutine get_magnetogram_field11(Xyz_D, B0_D)

    ! Return B0_D [Tesla] field at position Xyz_D [Rs]

    real, intent(in) :: Xyz_D(3)
    real, intent(out):: B0_D(3)

    real:: rLonLat_D(3), r
    !--------------------------------------------------------------------------
    ! Converting to rlonlat (radians)
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

    ! Scale with r^2 for r > rMaxB0
    if(r > rMaxB0) B0_D = (rMaxB0/r)**2 * B0_D

    ! Multiply with B0 factor and convert from Gauss to Tesla
    B0_D = B0_D*FactorB0*1e-4

  end subroutine get_magnetogram_field11
  !============================================================================
  subroutine interpolate_field(rPhiTheta_D, Brphitheta_D)

    ! A great example of left-handed coordinates
    ! Get the field at this location

    real, intent(in) :: rPhiTheta_D(3)
    real, intent(out):: Brphitheta_D(3)

    real:: r, Phi, Theta, Lon, Lat, Bxyz_D(3), Bsph_D(3), XyzSph_DD(3,3)
    !--------------------------------------------------------------------------
    r     = rPhiTheta_D(1)
    Phi   = rPhiTheta_D(2)
    Theta = rPhiTheta_D(3)

    ! Include the shift in Phi coordinate and make sure that it is
    ! in the range provided by the lookup table and convert to degrees
    if(dLonB0 /= 0.0 .or. LonMinB0 /= 0.0) &
         Phi = modulo(Phi - dLonB0 - LonMinB0, cTwoPi) + LonMinB0
    Lon = cRadToDeg*Phi

    ! Convert to latitude in degrees
    Lat = 90 - cRadToDeg*Theta

    ! Get Bxyz
    call interpolate_lookup_table(iTableB0, [r,Lon,Lat], Bxyz_D, &
         DoExtrapolate=(r<rMinB0) )

    ! Scale with r^2 for r > rMaxB0
    if(r > rMaxB0) Bxyz_D = (rMaxB0/r)**2 * Bxyz_D

    ! Multiply with B0 factor and convert from Gauss to Tesla
    Bxyz_D = Bxyz_D*FactorB0*1e-4

    ! Convert to spherical coordinates
    XyzSph_DD = rot_xyz_sph(Theta, Phi)
    Bsph_D = matmul(Bxyz_D, XyzSph_DD)

    ! Convert to left-handed spherical coordinates
    Brphitheta_D = [ Bsph_D(1), Bsph_D(3), Bsph_D(2) ]

  end subroutine interpolate_field
  !============================================================================

  subroutine get_hlcmm(StringHead, Shift)

    ! Read H(eliographic) L(ongitude) of the C(entral) M(eridian) of
    ! the M(ap) from the header line. Assign PhiOffset=HLCMM-180

    character (LEN=80),intent(inout):: StringHead
    real, intent(inout):: Shift

    real:: LonCentral     ! Heliographic Longitude of the central meridian
    integer:: iLonCentral ! LonCentral is integer in WSO magnetograms
    integer:: iErrorRead, iPosition

    !--------------------------------------------------------------------------
    iPosition=index(StringHead,'Centered')
    if (iPosition>0 .and. Shift<0.0)then
       StringHead(1:len(StringHead)-iPosition)=&
            StringHead(iPosition+1:len(StringHead))
       iPosition=index(StringHead,':')
       StringHead(1:len(StringHead)-iPosition)=&
            StringHead(iPosition+1:len(StringHead))
       iPosition=index(StringHead,':')
       StringHead(1:len(StringHead)-iPosition)=&
            StringHead(iPosition+1:len(StringHead))
       read(StringHead,'(i3)',iostat=iErrorRead)iLonCentral
       if(iErrorRead>0)call Con_stop(&
            'Cannot find LonCentral, '//NameMagnetogram//&
            ' is not a true WSO magnetogram')
       ! Rotates based on magnetogram central meridian + HGR system
       Shift = modulo(iLonCentral-180-dLongitudeHgrDeg, 360.0)

       RETURN
    end if

    iPosition=index(StringHead,'Central')
    if(iPosition>0 .and. Shift<0.0)then
       StringHead(1:len(StringHead)-iPosition)=&
            StringHead(iPosition+1:len(StringHead))
       iPosition=index(StringHead,':')
       StringHead(1:len(StringHead)-iPosition)=&
            StringHead(iPosition+1:len(StringHead))
       read(StringHead,*,iostat=iErrorRead)LonCentral
       if(iErrorRead>0)call CON_stop(&
            'Cannot find LonCentral, '//NameMagnetogram//&
            ' is not a true MDI magnetogram')
       ! Rotates based on magnetogram central meridian + HGR system
       Shift=modulo(LonCentral-180-dLongitudeHgrDeg, 360.0)

       RETURN
    end if

    ! Rotation of HGR system is applied
    if (Shift == 0.) &
       Shift = modulo(-dLongitudeHgrDeg, 360.0)

  end subroutine get_hlcmm
  !============================================================================

  subroutine read_magnetogram_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    character(len=*), parameter:: NameSub = 'read_magnetogram_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#FACTORB0")
       call read_var("FactorB0", FactorB0)

    case("#MAGNETOGRAM")
       call read_var('UseMagnetogram', UseMagnetogram)
       if(UseMagnetogram)then
          call read_var('rMagnetogram',        rMagnetogram)
          call read_var('rSourceSurface',      rSourceSurface)
          call read_var('HeightInnerBc',       Height)
          call read_var('NameMagnetogramFile', NameMagnetogram)
          call read_var('nHeaderLine',         nHeadLine)
          call read_var('PhiShift',            PhiOffset)
          call read_var('UnitB',               UnitB)
       end if

    case("#NEWMAGNETOGRAM")
       call read_var('UseNewMagnetogram',   UseNewMagnetogram)
       if(UseNewMagnetogram)then
          call read_var('rMagnetogram',     rMagnetogram)
          call read_var('rSourceSurface',   rSourceSurface)
          call read_var('HeightInnerBc',    Height)
          call read_var('NameMagnetogram2',      NameMagnetogram2)
          call read_var('nHeaderLine',      nHeadLine)
          call read_var('NewPhiShift',      NewPhiShift)
          call read_var('UnitB',            UnitB)
       end if

    case default
       call CON_stop(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine read_magnetogram_param
  !============================================================================

  subroutine read_magnetogram_file
    !--------------------------------------------------------------------------
    ! Initialize once g(n+1,m+1) & h(n+1,m+1) by reading a file
    ! created from Web data::

    call read_harmonics
    call set_sqrt_arrays

  end subroutine read_magnetogram_file
  !============================================================================
  subroutine read_new_magnetogram_file

    character(LEN=32):: NameOldFile
    real:: OldPhiShift

    !--------------------------------------------------------------------------
    NameOldFile = NameMagnetogram
    NameMagnetogram = NameMagnetogram2
    OldPhiShift = PhiOffset
    PhiOffset = NewPhiShift
    call read_magnetogram_file

    PhiOffset = OldPhiShift
    NameMagnetogram = NameOldFile

  end subroutine read_new_magnetogram_file
  !============================================================================

  subroutine read_harmonics

    ! Formats adjusted for wso CR rad coeffs

    integer :: iUnit, nOrderIn, iPosition, iPosition1, n, m, i, iError
    character (LEN=80):: StringHead=''
    real:: gTemp, hTemp, Coef1, Coef
    !--------------------------------------------------------------------------
    iUnit = io_unit_new()
    open(iUnit,file=NameMagnetogram,status='old',iostat=iError)
    if(iError>0)call CON_stop('Cannot open '//NameMagnetogram)
    if (nHeadLine /= 0) then
       do i=1,nHeadLine
          read(iUnit,'(a)') StringHead
          iPosition=index(StringHead,'rder')
          iPosition1=0
          if(iPosition>0)&
               iPosition1=max(&
               index(StringHead(iPosition+4:iPosition+9),'='),&
               index(StringHead(iPosition+4:iPosition+9),':'))
          if(iPosition1>0)then
             read(StringHead(iPosition+4+iPosition1:len(StringHead)),&
                  '(i3)',iostat=iError) nOrderIn
             if(iError>0)call CON_stop('Cannot figure out nOrder')
             if(nOrderIn<nOrder) nOrder=nOrderIn

          end if
          call get_hlcmm(StringHead, PhiOffset)
       enddo
       if(PhiOffset<0.0)call CON_stop(&
            'Did not find central meridian longitude')
    endif

    ! Initialize all coefficient arrays::
    g_II = 0.0; h_II = 0.0
    ! Read file with coefficients, g_II and h_II::
    do
       read(iUnit,*,iostat=iError) n, m, gTemp, hTemp
       if (iError /= 0) EXIT
       if (n > nOrder .or. m > nOrder) CYCLE
       g_II(n+1,m+1) = gTemp
       h_II(n+1,m+1) = hTemp
    enddo
    close(iUnit)
    ! Add correction factor for radial, not LOS, coefficients::
    ! Note old "coefficients" file are LOS, all new coeffs and
    ! files are radial)
    Coef=rSourceSurface
    do n=0,nOrder
       Coef=Coef/rSourceSurface**2
       Coef1 = 1.0/(real(n+1)+real(n)*Coef)
       !       Coef1 = 1.0/real(n+1+(n/(rSourceSurface**(2*n+1))))
       g_II(n+1,1:n+1) = g_II(n+1,1:n+1)*Coef1
       h_II(n+1,1:n+1) = h_II(n+1,1:n+1)*Coef1
    enddo
    ! Leave out monopole (n=0) term::
    g_II(1,1) = 0.0

  end subroutine read_harmonics
  !============================================================================

  subroutine set_sqrt_arrays

    ! Calculate square roots and their ratios used in the spherical harmonics
    ! functions up to order nOrder for speeding up the calculations

    integer:: m
    !--------------------------------------------------------------------------
    if(allocated(Sqrt_I)) RETURN

    allocate(Sqrt_I(MaxInt), SqrtRatio_I(nOrder+1))

    do m = 1, MaxInt
       Sqrt_I(m) = sqrt(real(m))
    end do

    ! Calculate the ratio sqrt(2m!)/(2^m*m!) recursively
    SqrtRatio_I(1) = 1.0
    do m = 1, nOrder
       SqrtRatio_I(m+1) = SqrtRatio_I(m)*Sqrt_I(2*m-1)/Sqrt_I(2*m)
    enddo

  end subroutine set_sqrt_arrays
  !============================================================================

  subroutine get_harmonics_field(r, Theta, Phi, Bsph_D)

    ! Calculate the harmonics based potential field Bsph_D at r, Theta, Phi

    real, intent(in):: r, Theta, Phi ! r in Rs, Theta and Phi in radians
    real, intent(out):: Bsph_D(3)    ! Bsph_D with r, Theta, Phi components

    integer:: n, m
    real:: CosTheta, SinTheta

    ! Azimuthal functions sin(m*phi) and cos(m*phi)
    real, allocatable, save:: SinPhi_I(:), CosPhi_I(:)

    ! The spherical components of the magnetic field
    real:: SumR, SumT, SumP

    ! Legendre polynomials and its derivative
    real, save, allocatable:: p_II(:,:), Dp_II(:,:)

    ! Powers of radius ratios
    real, save, allocatable:: rRsPower_I(:), RoRsPower_I(:), RoRPower_I(:)

    real:: Coef1, Coef2, Coef3, Coef4
    !--------------------------------------------------------------------------

    ! Calculate the radial part of spherical functions
    call calc_radial_functions(r)

    ! Calculate the set of Legendre polynoms for given Theta
    CosTheta = cos(Theta)
    SinTheta = max(sin(Theta), 1E-10)
    call calc_legendre_polynomial(SinTheta, CosTheta)

    call calc_azimuthal_functions(Phi)

    ! Initialize the values for the sums
    SumR = 0.0; SumT = 0.0; SumP = 0.0

    ! Calculate B from spherical harmonics
    do m = 0, nOrder; do n = m, nOrder

       Coef1 = (n+1)*RoRPower_I(n+2) + RoRsPower_I(n+2)*n*rRsPower_I(n-1)
       Coef3 =       RoRPower_I(n+2) - RoRsPower_I(n+2)*rRsPower_I(n-1)
       Coef2 = g_II(n+1,m+1)*CosPhi_I(m) + h_II(n+1,m+1)*SinPhi_I(m)
       Coef4 = g_II(n+1,m+1)*SinPhi_I(m) - h_II(n+1,m+1)*CosPhi_I(m)

       ! Br = -d(Psi)/dR
       SumR  = SumR + p_II(n+1,m+1)*Coef1*Coef2

       ! Bt = -(1/r)*d(Psi)/dTheta
       SumT  = SumT - Dp_II(n+1,m+1)*Coef2*Coef3

       ! Bp = -(1/r)*d(Psi)/dPhi
       SumP  = SumP + p_II(n+1,m+1)*m/SinTheta * Coef3*Coef4

       ! Potential
       ! SumPsi = SumPsi + r*p_II(n+1,m+1)*Coef2*Coef3

    enddo; enddo

    Bsph_D = [ SumR, SumT, SumP ]

  contains
    !==========================================================================

    subroutine calc_legendre_polynomial(SinTheta, CosTheta)

      ! Calculate Legendre polynomials p_II and its derivative Dp_II
      ! with appropriate normalization for Theta

      real, intent(in):: SinTheta, CosTheta

      real:: SinThetaM, SinThetaM1  ! sin(Theta)^m, sin(Theta)^(m-1)
      integer:: m

      integer:: nOrderOld = -10
      real:: SinThetaOld = -10.0, CosThetaOld = -10.0

      ! Cache previous values
      !------------------------------------------------------------------------
      if(nOrder == nOrderOld .and. SinTheta == SinThetaOld &
           .and. CosTheta == CosThetaOld) RETURN

      if(.not.allocated(p_II)) &
           allocate(p_II(nOrder+1,nOrder+1), Dp_II(nOrder+1,nOrder+1))

      nOrderOld = nOrder; SinThetaOld = SinTheta; CosThetaOld = CosTheta

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
         ! Eq.(27) from Altschuler et al. 1976::
         p_II(m+1,m+1) = SqrtRatio_I(m+1)*Coef1* SinThetaM
         ! Eq.(28) from Altschuler et al. 1976::
         if (m < nOrder) p_II(m+2,m+1) = p_II(m+1,m+1)*Sqrt_I(2*m+3)* &
              CosTheta
         ! Eq.(30) from Altschuler et al. 1976::
         Dp_II(m+1,m+1) = SqrtRatio_I(m+1)*Coef1*m*CosTheta*SinThetaM1
         ! Eq.(31) from Altschuler et al. 1976::
         if (m < nOrder) &
              Dp_II(m+2,m+1) = Sqrt_I(2*m+3)* &
              (CosTheta*Dp_II(m+1,m+1) - SinTheta*p_II(m+1,m+1))

         SinThetaM1 = SinThetaM
         SinThetaM  = SinThetaM*SinTheta

      enddo

      ! Recursive rules
      do m = 0, nOrder-2; do n = m+2, nOrder
         ! Eq.(29) from Altschuler et al. 1976::
         Coef1 = Sqrt_I(2*n+1)/Sqrt_I(n**2-m**2)
         Coef2 = Sqrt_I(2*n-1)
         Coef3 = Sqrt_I((n-1)**2-m**2)/Sqrt_I(2*n-3)

         p_II(n+1,m+1) = Coef1*(Coef2*CosTheta*p_II(n,m+1)  &
              - Coef3*p_II(n-1,m+1))

         ! Eq.(32) from Altschuler et al. 1976::
         Dp_II(n+1,m+1) = Coef1*(Coef2*(CosTheta*Dp_II(n,m+1) &
              - SinTheta*p_II(n,m+1)) - Coef3*Dp_II(n-1,m+1))
      enddo; enddo

      ! Apply Schmidt normalization
      do m = 0, nOrder; do n = m, nOrder
         ! Eq.(33) from Altschuler et al. 1976::
         Coef1 = 1.0/Sqrt_I(2*n+1)
         ! Eq.(34) from Altschuler et al. 1976::
         p_II(n+1,m+1)  = p_II(n+1,m+1)*Coef1
         Dp_II(n+1,m+1) = Dp_II(n+1,m+1)*Coef1
      enddo; enddo

    end subroutine calc_legendre_polynomial
    !==========================================================================

    subroutine calc_radial_functions(r)

      real,    intent(in):: r

      ! Calculate powers of the ratios of radii up to nOrder

      integer:: m
      real:: RoRs, RoR, rRs

      integer:: nOrderOld = -10
      real::    rOld = -10.0

      !------------------------------------------------------------------------
      if(nOrderOld == nOrder .and. rOld == r) RETURN

      if(nOrderOld /= nOrder .and. allocated(rRsPower_I)) &
           deallocate(rRsPower_I, RoRsPower_I, RoRPower_I)

      nOrderOld = nOrder; rOld = r

      if(.not.allocated(rRsPower_I)) allocate( &
           rRsPower_I(-1:nOrder+2), &
           RoRsPower_I(0:nOrder+2), &
           RoRPower_I(0:nOrder+2))

      RoRs = rMagnetogram/rSourceSurface
      RoR  = rMagnetogram/r
      rRs  = r/rSourceSurface

      ! Zero and -1 powers
      rRsPower_I(-1) = 1.0/rRs
      rRsPower_I(0)  = 1.0
      RoRsPower_I(0) = 1.0
      RoRPower_I(0)  = 1.0

      ! Recursive: x^m = x^m-1 * x
      do m = 1, nOrder+2
         RoRsPower_I(m) = RoRsPower_I(m-1) * RoRs
         RoRPower_I(m)  = RoRPower_I(m-1)  * RoR
         rRsPower_I(m)  = rRsPower_I(m-1)  * rRs
      end do

    end subroutine calc_radial_functions
    !==========================================================================
    subroutine calc_azimuthal_functions(Phi)

      ! Calculate azimuthal harmonics for given Phi

      real,    intent(in):: Phi

      integer:: m
      real   :: SinPhi, CosPhi

      integer:: nOrderOld = -10
      real   :: PhiOld    = -10.0

      !------------------------------------------------------------------------
      if(nOrderOld == nOrder .and. Phiold == Phi) RETURN

      if(nOrderOld /= nOrder .and. allocated(SinPhi_I)) &
           deallocate(SinPhi_I, CosPhi_I)

      if(.not.allocated(SinPhi_I)) &
           allocate(SinPhi_I(0:nOrder), CosPhi_I(0:nOrder))

      SinPhi_I(0) = 0.0
      CosPhi_I(0) = 1.0
      SinPhi = sin(Phi)
      CosPhi = Cos(Phi)

      ! Recursive: sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
      !            cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
      do m = 1, nOrder
         SinPhi_I(m) = SinPhi*CosPhi_I(m-1) + CosPhi*SinPhi_I(m-1)
         CosPhi_I(m) = CosPhi*CosPhi_I(m-1) - SinPhi*SinPhi_I(m-1)
      end do

    end subroutine calc_azimuthal_functions
    !==========================================================================
  end subroutine get_harmonics_field
  !============================================================================

end module ModMagnetogram
!==============================================================================
