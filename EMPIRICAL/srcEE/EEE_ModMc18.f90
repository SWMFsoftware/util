!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module EEE_ModMc18

  ! Rosenbluth & Bussac 1979 force-free spheromak model
  ! with image dipoles (Lin 2006) and 
  ! uniform field subtraction (Borovikov et al 2018).

#ifdef _OPENACC
  use ModUtilities, ONLY: norm2
#endif
  use EEE_ModCommonVariables

  implicit none

  SAVE

  private ! except

  public :: set_parameters_mc18
  public :: get_mc18_fluxrope
  public :: get_mc18_size
  public :: mc18_init

  ! Local variables

  ! Geometric characteristics of the superimposed configuration:

  ! contraction distance as in   r --> r -a
  real :: Stretch = 0.0
  !$acc declare create(Stretch)

  ! distance from the magnetic configuration center to heliocenter
  real :: rDistance1 = 0.0

  ! Radius of the magnetic configuration (spheromac)
  real :: Radius = 0.0
  !$acc declare create(Radius)

  ! Height of the configuration base above the solar surface [rSun]:
  !   1 + BaseHeight = rDistance1 - Radius
  ! Specifying BaseHeight is more straightforward than specifying rDistance1.
  real :: BaseHeight = 0.0

  ! The derivative of current over flux function
  !(\mu_0)dI/d\psi) has the dimentions of inverse length
  real :: Alpha0
  !$acc declare create(Alpha0)

  ! Characteristic value of magnetic field of the spheromak
  ! configuration: the field in the center of configuration
  ! equals 2(1/3 -\beta_0)B_0 \approx 0.7 B_0:
  real :: B0      ! dimensionless
  real :: B0Dim   ! in Gauss
  !$acc declare create(B0)

  ! Characteristic ratio of plasma pressure to B_0^2/\mu_0
  ! Exact definition in terms of the pressure derivative over
  ! the flux function: \beta_0=\mu_0/(B_0 Alpha_0^2) dp/d\psi
  real :: Beta0
  !$acc declare create(Beta0)

  ! Dimensionless product of R0 by Alpha0. 
  ! First zero of j_1(x): Rosenbluth-Bussac boundary condition k*r0 = Alpha0R0
  ! (GL98 used 5.763854, the first zero of j_2; we use the first zero of j_1)
  real, parameter :: Alpha0R0 = 4.493409457909064

  ! Vector characteristic of the configuration: radius vector of the
  ! configuration center and B0 multiplied by unit vector along
  ! the axis of symmetry
  real :: XyzCenterConf_D(3), bConf_D(3)
  !$acc declare create(XyzCenterConf_D, bConf_D)

  ! Parameter to control self-similar solution
  real :: uCmeSi = 0.0
  !$acc declare create(uCmeSi)
  real, parameter :: Delta = 0.1

  ! Half opening angle of the CME [degrees]:  sin(AlphaAngle) = Radius/rDistance1
  real :: AlphaAngle = 0.0

  ! Lin (2006) image dipole parameters
  logical :: UseImageDipoles = .false.
  !$acc declare create(UseImageDipoles)
  integer :: nDiscDipoles = 200

  ! Pre-computed image dipoles: 1 point image + up to MaxImgDipoles-1 line images
  integer, parameter :: MaxImgDipoles = 1001
  integer :: nImgDipoles = 0
  !$acc declare create(nImgDipoles)
  real :: rImgDipole_DI(3, MaxImgDipoles) = 0.0
  real :: mImgDipole_DI(3, MaxImgDipoles) = 0.0
  !$acc declare create(rImgDipole_DI, mImgDipole_DI)

  ! Plasma beta and ejecta temperature (same convention as TD99)
  ! p = PlasmaBeta * 0.5 * |B|^2 inside; Rho = p / EjectaTemperature
  logical :: UsePlasmaBeta = .false.
  !$acc declare create(UsePlasmaBeta)
  real :: PlasmaBeta = 0.0
  real :: EjectaTemperature = 0.0           ! normalized code units
  real :: EjectaTemperatureDim = 5.0e4      ! [K]
  !$acc declare create(PlasmaBeta, EjectaTemperature, EjectaTemperatureDim)

contains
  !============================================================================
  subroutine mc18_init

    !--------------------------------------------------------------------------
    ! Wave number k = Alpha0R0 / r0  (first zero of j_1 divided by radius)
    Alpha0 = Alpha0R0/Radius

    ! Pure force-free Rosenbluth-Bussac solution: no plasma pressure
    Beta0 = 0.0

    ! Field amplitude and axis from the ambient magnetic field at the CME center.
    ! bAmbientCenterSi_D is filled by the MHD solver before mc18_init is called
    ! (see SC_user_initial_perturbation in ModUserAwsom.f90).
    bConf_D = (3.0/(2.0*cos(Alpha0R0))) * bAmbientCenterSi_D*Si2No_V(UnitB_)
    B0      = norm2(bConf_D)
    B0Dim   = B0*No2Io_V(UnitB_)

    ! Convert self-similar CME speed from km/s to SI
    uCmeSi = uCmeSi*Io2Si_V(UnitU_)

    if(iProc==0)then
       write(*,*) prefix
       write(*,*) prefix, &
            '>>>>>>>>>>>>>>>>>>>                            <<<<<<<<<<<<<<<<<<<<<'
       write(*,*) prefix
       write(*,*) prefix, &
            '     EEGMC Magnetic Cone Model (Rosenbluth-Bussac 1979) is initiated'
       write(*,*) prefix
       write(*,*) prefix, &
            '>>>>>>>>>>>>>>>>>>>                            <<<<<<<<<<<<<<<<<<<<<'
       write(*,*) prefix
       write(*,*) prefix, 'B0Dim          = ', B0Dim,                '[Gauss]'
       write(*,*) prefix, 'Radius         = ', Radius,               '[rSun]'
       write(*,*) prefix, 'AlphaAngle     = ', AlphaAngle,           '[degrees]'
       write(*,*) prefix, 'rDistance1     = ', rDistance1,           '[rSun]'
       write(*,*) prefix, 'BaseHeight     = ', BaseHeight,           '[rSun]'
       write(*,*) prefix, 'LongitudeCme   = ', LongitudeCme,         '[degrees]'
       write(*,*) prefix, 'LatitudeCme    = ', LatitudeCme,          '[degrees]'
       write(*,*) prefix, 'Alpha0         = ', Alpha0,               '[1/rSun]'
       write(*,*) prefix, 'UseImageDipoles= ', UseImageDipoles
       if(UseImageDipoles) &
            write(*,*) prefix, 'nDiscDipoles   = ', nDiscDipoles
       write(*,*) prefix, 'UsePlasmaBeta  = ', UsePlasmaBeta
       if(UsePlasmaBeta)then
          write(*,*) prefix, 'PlasmaBeta     = ', PlasmaBeta
          write(*,*) prefix, 'EjectaTemp     = ', EjectaTemperatureDim, '[K]'
       end if
       write(*,*) prefix, 'Start time     = ', tStartCme,            '[s]'
       write(*,*) prefix, 'CME speed      = ', uCmeSi*Si2Io_V(UnitU_),'[km/s]'
       write(*,*) prefix
    end if

    ! Center position of the configuration in the heliocentric frame
    XyzCenterConf_D = rDistance1*DirCme_D

    EjectaTemperature = EjectaTemperatureDim*Io2No_V(UnitTemperature_)

    !$acc update device(Alpha0, Beta0, XyzCenterConf_D, bConf_D)
    !$acc update device(uCmeSi, B0, Radius, UseImageDipoles)
    !$acc update device(UsePlasmaBeta, PlasmaBeta, EjectaTemperature)

    if(UseImageDipoles) call mc18_compute_image_dipoles

  end subroutine mc18_init
  !============================================================================
  subroutine mc18_compute_image_dipoles

    real :: mDip_D(3), rHat_D(3), dSource
    real :: mR_D(3), mT_D(3)
    real :: dImg, scale3, du, uCtr
    integer :: i
    !--------------------------------------------------------------------------
    ! Equivalent external dipole: m = D*axis, D = -A*r0^3/2
    ! (A = B0, axis = bConf_D/B0)  =>  m = -bConf_D * Radius^3/2
    mDip_D = -bConf_D*(Radius**3/2.0)

    dSource = norm2(XyzCenterConf_D)
    rHat_D  = XyzCenterConf_D/dSource

    mR_D = sum(mDip_D*rHat_D)*rHat_D  ! radial component of moment
    mT_D = mDip_D - mR_D              ! transverse component of moment

    ! Image point at (a^2/d)*rHat = (1/d)*rHat  (solar surface a = 1 R_sun)
    dImg   = 1.0/dSource
    scale3 = dImg**3              ! = (a/d)^3

    nImgDipoles = 1 + nDiscDipoles

    ! Combined point image: moment = scale3*(m_t - m_r)  (Lin 2006 Eq. II.B-C)
    rImgDipole_DI(:,1) = dImg*rHat_D
    mImgDipole_DI(:,1) = scale3*(mT_D - mR_D)

    ! Discretised line images: cell-centred over [0, dImg] along rHat_D
    du = dImg/nDiscDipoles
    do i = 1, nDiscDipoles
       uCtr = (i - 0.5)*du
       rImgDipole_DI(:, 1+i) = uCtr*rHat_D
       mImgDipole_DI(:, 1+i) = -(mT_D/dSource)*uCtr*du
    end do

    !$acc update device(nImgDipoles, rImgDipole_DI, mImgDipole_DI)

  end subroutine mc18_compute_image_dipoles
  !============================================================================
  subroutine set_parameters_mc18(NameCommand)

    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand

    real :: SinAlpha
    character(len=*), parameter:: NameSub = 'set_parameters_mc18'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#CME","#MAGCONE")
       call read_var('AlphaAngle',      AlphaAngle)     ![deg] half opening angle
       call read_var('BaseHeight',      BaseHeight)     ![rSun]
       call read_var('uCmeSi',          uCmeSi)         ![km/s]
       call read_var('UseImageDipoles', UseImageDipoles)
       if(UseImageDipoles)then
          call read_var('nDiscDipoles', nDiscDipoles)
          if(nDiscDipoles > MaxImgDipoles - 1) call CON_stop( &
               NameSub//': nDiscDipoles exceeds MaxImgDipoles-1 = 1000')
       end if
       call read_var('UsePlasmaBeta', UsePlasmaBeta)
       if(UsePlasmaBeta)then
          call read_var('PlasmaBeta', PlasmaBeta)
          call read_var('EjectaTemperature', EjectaTemperatureDim)
       end if

       ! Derive center distance and radius from opening angle and base height.
       ! Geometry: sin(alpha) = Radius/rDistance1
       !           1 + BaseHeight = rDistance1 - Radius = rDistance1*(1 - sin(alpha))
       if(AlphaAngle <= 0.0 .or. AlphaAngle >= 90.0) call CON_stop( &
            NameSub//': AlphaAngle must be in (0, 90) degrees')
       SinAlpha   = sin(AlphaAngle*cDegToRad)
       rDistance1 = (1.0 + BaseHeight)/(1.0 - SinAlpha)
       Radius     = rDistance1*SinAlpha

       ! position of the CME apex
       XyzCmeApexSi_D = DirCme_D*(rDistance1 + Radius)

       ! position of CME center
       XyzCmeCenterSi_D = XyzCmeApexSi_D - DirCme_D*Radius
       DoNormalizeXyz = .true.

    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_mc18
  !============================================================================
  subroutine get_mc18_fluxrope(XyzIn_D, Rho, p, b_D, u_D, TimeNow)
    !$acc routine seq

    ! Magnetic field of the Rosenbluth-Bussac force-free spheromak.
    ! Interior: spherical Bessel (j1) field with Beta0=0.
    ! Exterior: pure dipole of equivalent moment m = -bConf_D * Radius^3/2.
    ! Image dipole corrections (Lin 2006) applied everywhere if UseImageDipoles.

    use ModCoordTransform, ONLY: cross_product

    ! Coordinates of the input point, in rSun
    real, intent(in) :: XyzIn_D(3)

    ! OUTPUTS
    real, intent(out) :: b_D(3)
    real, intent(out), optional :: u_D(3)
    real, optional, intent(in)  :: TimeNow

    ! Density, pressure (non-zero inside if UsePlasmaBeta)
    real, intent(out) :: Rho, p

    real :: XyzConf_D(3), Distance2ConfCenter
    real :: R2CrossB0_D(3), Alpha0R2
    real :: mDip_D(3), MdotR
    real :: PhiInv
    real :: dr_D(3), R2img
    integer :: i
    !--------------------------------------------------------------------------
    if(UseMagCone)then
       PhiInv = 1.0/(1.0 + (TimeNow - tStartCme)*uCmeSi*rCmeApexInvSi)
    else
       PhiInv = 1.0
    end if

    Rho = 0.0
    p   = 0.0
    if(present(u_D)) u_D = 0.0

    ! Position relative to spheromak center (self-similar scaling applied)
    XyzConf_D = XyzIn_D*PhiInv - XyzCenterConf_D

    Distance2ConfCenter = norm2(XyzConf_D)
    if(Distance2ConfCenter <= Delta)then
       XyzConf_D           = XyzConf_D*(Delta/Distance2ConfCenter)
       Distance2ConfCenter = Delta
    end if

    if(Distance2ConfCenter <= Radius)then

       ! INSIDE: force-free spherical Bessel field (Beta0 = 0)
       Alpha0R2    = Alpha0*Distance2ConfCenter
       R2CrossB0_D = cross_product(XyzConf_D, bConf_D)
       b_D = (2*bConf_D - sign(Alpha0, B0)*R2CrossB0_D) &
            *spher_bessel1_over_x(Alpha0R2)              &
            + spher_bessel2(Alpha0R2)/Distance2ConfCenter**2 &
            *cross_product(XyzConf_D, R2CrossB0_D)

       if(UsePlasmaBeta)then
          p   = PlasmaBeta*0.5*sum(b_D**2)
          Rho = p/EjectaTemperature
       end if

       if(present(u_D) .and. UseMagCone)then
          u_D = XyzIn_D*PhiInv  &
               *No2Si_V(UnitX_) &
               *uCmeSi*rCmeApexInvSi
       end if

    else

       ! OUTSIDE: pure dipole  m = -bConf_D * Radius^3/2
       ! (uniform A*axis component is already in the background MHD state)
       mDip_D = -bConf_D*(Radius**3/2.0)
       MdotR  = sum(mDip_D*XyzConf_D)
       b_D = (3.0*MdotR*XyzConf_D/Distance2ConfCenter**2 - mDip_D) &
            /Distance2ConfCenter**3

    end if

    ! Lin (2006) image dipole corrections (heliocentric coordinates, unscaled)
    if(UseImageDipoles)then
       do i = 1, nImgDipoles
          dr_D  = XyzIn_D - rImgDipole_DI(:,i)
          R2img = sum(dr_D**2)
          MdotR = sum(mImgDipole_DI(:,i)*dr_D)
          b_D   = b_D + (3.0*MdotR*dr_D/R2img - mImgDipole_DI(:,i)) &
               /(sqrt(R2img)*R2img)
       end do
    end if

    b_D = b_D*No2Si_V(UnitB_)
    Rho = Rho*No2Si_V(UnitRho_)
    p   = p  *No2Si_V(UnitP_)
    ! u_D is already in SI [m/s] from the interior formula above

  contains
    !==========================================================================
    real function spher_bessel0(x)
      real, intent(in) :: x
      !------------------------------------------------------------------------
      spher_bessel0 = sin(x)/x
    end function spher_bessel0
    !==========================================================================
    real function spher_bessel1(x)
      real, intent(in) :: x
      !------------------------------------------------------------------------
      spher_bessel1 = (sin(x) - x*cos(x))/x**2
    end function spher_bessel1
    !==========================================================================
    real function spher_bessel1_over_x(x)
      real, intent(in) :: x
      !------------------------------------------------------------------------
      if(x == 0)then
         spher_bessel1_over_x = 1.0/3.0
      else
         spher_bessel1_over_x = (sin(x) - x*cos(x))/x**3
      end if
    end function spher_bessel1_over_x
    !==========================================================================
    real function spher_bessel2(x)
      real, intent(in) :: x
      !------------------------------------------------------------------------
      spher_bessel2 = 3*spher_bessel1_over_x(x) - spher_bessel0(x)
    end function spher_bessel2
    !==========================================================================
  end subroutine get_mc18_fluxrope
  !============================================================================
  subroutine get_mc18_size(SizeXY,  SizeZ)
    real,  intent(out) :: SizeXY,  SizeZ
    !--------------------------------------------------------------------------
    SizeXY = Radius                    ! Horizontal size
    SizeZ  = rDistance1 + Radius - 1.0 ! Apex height above solar surface
  end subroutine get_mc18_size
  !============================================================================
end module EEE_ModMc18
!==============================================================================
