!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModMagnetogram

  use ModNumConst,       ONLY: cTwoPi, cRadToDeg, cDegToRad
  use CON_axes,          ONLY: dLongitudeHgrDeg
  use ModUtilities,      ONLY: CON_stop
  use ModCoordTransform, ONLY: rot_xyz_sph, rot_matrix_z, xyz_to_rlonlat
  use ModLookupTable, ONLY: &
       i_lookup_table, init_lookup_table, make_lookup_table_3d, &
       get_lookup_table, interpolate_lookup_table

  implicit none
  save

  private

  integer, public:: iTableB0    = -1 ! index of B0 lookup table
  integer, public:: iTableB0New = -1 ! index of new B0 lookup table

  public:: read_magnetogram_param ! read parameters

  public:: init_magnetogram_lookup_table ! initialize lookup table

  public:: get_magnetogram_field ! field from lookup table at a given point

  public:: interpolate_field ! good old legacy left handed coordinates

  ! Local variables -------------

  ! Carrington rotation of the magnetogram(s) for temporal interpolation
  real:: CarringtonRot, CarringtonRotNew

  ! radius of magnetogram and source surface for spatial extrapolation
  real:: rMagnetogram=1.0, rSourceSurface=2.5 

  ! True if new harmonics coefficients are to be read
  logical:: DoReadHarmonics = .false., DoReadHarmonicsNew = .false.

  ! Name of the harmonics files
  character(len=100):: NameHarmonicsFile = '???', NameHarmonicsFileNew = '???'

  ! Maximum order of spherical harmonics
  integer:: nOrder = 30

  ! Weights of the spherical harmonics
  real, allocatable:: g_II(:,:), h_II(:,:)

  ! Powers of radius ratios
  real, allocatable:: rRsPower_I(:), RoRsPower_I(:), RoRPower_I(:)

  ! Legendre polynomials and its derivative
  real, allocatable:: p_II(:,:), Dp_II(:,:)

  ! Azimuthal functions sin(m*phi) and cos(m*phi)
  real, allocatable:: SinPhi_I(:), CosPhi_I(:)

  ! Square roots of integers in spherical harmonics
  real, allocatable:: Sqrt_I(:), SqrtRatio_I(:)

  ! Lookup table related variables
  real:: rMinB0=1.0, rMaxB0=30.0 ! radial limits of table
  real:: LonMinB0 = 0.0 ! starting longitude in the table
  real:: dLonB0=0.0     ! longitude shift
  real:: RotB0_DD(3,3)  ! rotation matrix due to longitude shift
  real:: FactorB0=1.0   ! multiplier for the magnetic field

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
    ! Make sure these are set
    iTableB0    = i_lookup_table('B0')
    iTableB0New = i_lookup_table('B0NEW')

    ! Nothing to do if there is no B0 table
    if(iTableB0 < 0 .and. .not.DoReadHarmonics) RETURN
    
    if(DoReadHarmonics)then
       ! Read harmonics coefficient, nOrder and set carrington rotation info
       call read_harmonics_file(NameHarmonicsFile, CarringtonRot)
       DoReadHarmonics = .false.

       ! Create default size lookup table if not specified in PARAM.in
       if(iTableB0 <= 0)then

          ! Set up a default B0 lookup table
          nR = max(30, nOrder); nLon = max(72, nOrder); nLat=max(30, nOrder)
          call init_lookup_table( &
               NameTable   = 'B0',                    &
               NameCommand = 'use',                   &
               NameVar     = 'r lon lat bx by bz',    &
               NameFile    = 'harmonics_bxyz.out',    &
               TypeFile    = 'real4',                 &
               nIndex_I    = [nR+1, nLon+1, nLat+1],  &
               IndexMin_I  = [ rMagnetogram,     0.0, -90.0],    &
               IndexMax_I  = [ rSourceSurface, 360.0,  90.0],    &
               StringDescription = 'Created from '//trim(NameHarmonicsFile) )
          iTableB0 = i_lookup_table('B0')
       end if
       ! If lookup table is not loaded, make it and save it
       call make_lookup_table_3d(iTableB0, calc_b0_table, iComm)

    end if

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

    ! Rotation matrix for longitude shift if needed
    if(nParam > 2) then
       dLonB0 = (Param_I(3) - dLongitudeHgrDeg)*cDegToRad
       RotB0_DD = rot_matrix_z(dLonB0)
    end if

    ! Second lookup table

    if(DoReadHarmonicsNew)then
       ! Read harmonics coefficients and set Carrington rotation info
       call read_harmonics_file(NameHarmonicsFileNew, CarringtonRotNew)
       DoReadHarmonicsNew = .false.

       ! Set up default sized lookup table if not defeined in PARAM.in
       if(iTableB0New < 0)then
          nR = max(30, nOrder); nLon = max(72, nOrder); nLat=max(30, nOrder)
          call init_lookup_table( &
               NameTable   = 'B0New',                 &
               NameCommand = 'use',                   &
               NameVar     = 'r lon lat bx by bz',    &
               NameFile    = 'harmonics_new_bxyz.out',&
               TypeFile    = 'real4',                 &
               nIndex_I    = [nR+1, nLon+1, nLat+1],  &
               IndexMin_I  = [ rMagnetogram,     0.0, -90.0],    &
               IndexMax_I  = [ rSourceSurface, 360.0,  90.0],    &
               StringDescription = 'Created from '//trim(NameHarmonicsFileNew))
          iTableB0New = i_lookup_table('B0')
       end if
       ! Make second lookup table
       call make_lookup_table_3d(iTableB0New, calc_b0_table, iComm)
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

  subroutine read_magnetogram_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    character(len=*), parameter:: NameSub = 'read_magnetogram_param'

    real:: Height ! pointless variable
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#FACTORB0")
       call read_var("FactorB0", FactorB0)

    case("#HARMONICSFILE")
       DoReadHarmonics = .true.
       call read_var('NameHarmonicsFile', NameHarmonicsFile)
       call read_var('rMagnetogram',      rMagnetogram)
       call read_var('rSourceSurface',    rSourceSurface)

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

  subroutine read_harmonics_file(NameFile, Carrington)

    use ModPlotFile, ONLY: read_plot_file

    character(len=*), intent(in) :: NameFile
    real,             intent(out):: Carrington ! carrington rotation from file

    ! read properly formatted harmonics coefficient file

    integer:: nHarmonic ! number of harmonics = (nOrder+1)*(nOrder+2)/2
    integer:: nParam, nOrderIn, i, n, m

    real:: Param_I(3), Coef, Coef1, PhiOffset
    real, allocatable:: Var_VI(:,:), Coord_DII(:,:,:)

    character(len=*), parameter:: NameSub = 'read_harmonics_file'
    !--------------------------------------------------------------------------

    call read_plot_file(NameFile, &
         nParamOut=nParam, ParamOut_I=Param_I, n1Out=nHarmonic)

    if(nParam < 3) call CON_stop(NameSub// &
         ': not enough parameters in '//trim(NameFile))

    nOrderIn   = nint(Param_I(1))
    nOrder     = min(nOrder,nOrderIn)
    Carrington = Param_I(2)
    PhiOffset  = Param_I(3)

    if(nHarmonic /= (nOrderIn+1)*(nOrderIn+2)/2) call CON_stop(NameSub// &
         ': inconsistent n and nOrder in '//trim(NameFile))

    call allocate_harmonics_arrays

    allocate(Var_VI(2,nHarmonic), Coord_DII(2,nHarmonic,1))
    call read_plot_file(NameFile, VarOut_VI=Var_VI, CoordOut_DII=Coord_DII)

    do i = 1, nHarmonic
       n = Coord_DII(1,i,1)
       m = Coord_DII(2,i,1)
       if (n > nOrder .or. m > nOrder) CYCLE
       g_II(n+1,m+1) = Var_VI(1,i)
       h_II(n+1,m+1) = Var_VI(2,i)
    end do

    deallocate(Var_VI, Coord_DII)

    ! Normalize coefficients
    Coef = rSourceSurface
    do n = 0, nOrder
       Coef  = Coef/rSourceSurface**2
       Coef1 = 1.0/(n+1 + n*Coef)
       g_II(n+1,1:n+1) = g_II(n+1,1:n+1)*Coef1
       h_II(n+1,1:n+1) = h_II(n+1,1:n+1)*Coef1
    enddo
    ! Leave out monopole (n=0) term::
    g_II(1,1) = 0.0

  end subroutine read_harmonics_file
  !============================================================================

  subroutine allocate_harmonics_arrays

    ! Calculate square roots and their ratios used in the spherical harmonics
    ! functions up to order nOrder for speeding up the calculations

    integer:: m, MaxInt
    !--------------------------------------------------------------------------
    if(allocated(g_II)) RETURN

    ! Spherical harmonics coefficients
    allocate(g_II(nOrder+1,nOrder+1), h_II(nOrder+1,nOrder+1))
    g_II = 0.0; h_II = 0.0

    ! radial functions
    allocate(rRsPower_I(-1:nOrder+2), &
         RoRsPower_I(0:nOrder+2), RoRPower_I(0:nOrder+2))

    ! Legendre polynomlials
    allocate(p_II(nOrder+1,nOrder+1), Dp_II(nOrder+1,nOrder+1))

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
         rRsPower_I, RoRsPower_I, RoRPower_I, &
         SinPhi_I, CosPhi_I)

  end subroutine deallocate_harmonics_arrays
  !============================================================================
  subroutine get_harmonics_field(r, Theta, Phi, Bsph_D)

    ! Calculate the harmonics based potential field Bsph_D at r, Theta, Phi

    real, intent(in):: r, Theta, Phi ! r in Rs, Theta and Phi in radians
    real, intent(out):: Bsph_D(3)    ! Bsph_D with r, Theta, Phi components

    integer:: n, m
    real:: CosTheta, SinTheta

    ! The spherical components of the magnetic field
    real:: SumR, SumT, SumP

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

      nOrderOld = nOrder; rOld = r

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
