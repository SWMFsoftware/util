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
  use EEE_ModCommonVariables, HalfOpeningAngle=>OrientationCme

  implicit none

  SAVE

  private ! except

  public :: set_parameters_mc18
  public :: get_mc18_fluxrope
  public :: get_mc18_size
  public :: mc18_init

  ! Local variables

  ! Geometric characteristics of the superimposed configuration:

  ! distance from the magnetic configuration center to heliocenter
  real :: rDistance1 = 0.0
  !$acc declare create(rDistance1)

  ! Radius of the magnetic configuration (spheromak)
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

  ! Sign of alpha: chosen so the toroidal field at the spheromak bottom
  ! opposes the ambient field component perpendicular to the
  ! plane spanned by DirCme_D and bConf_D.
  real :: iHelicity = 1.0
  !$acc declare create(iHelicity)

  ! Dimensionless product of R0 by Alpha0.
  ! Boundary condition: j1(Alpha0R0)/Alpha0R0 = Beta0
  ! For Beta0=0: Alpha0R0 is the first zero of j1 (~4.4934).
  ! For Beta0>0: solved iteratively in find_alpha0r0.
  ! (GL98 used 5.763854, the first zero of j_2)
  real :: Alpha0R0 = 4.493409457909064
  !$acc declare create(Alpha0R0)

  ! Vector characteristic of the configuration: radius vector of the
  ! configuration center and B0 multiplied by unit vector along
  ! the axis of symmetry
  real :: XyzCenterConf_D(3), bConf_D(3)
  !$acc declare create(XyzCenterConf_D, bConf_D)

  ! Normalized ambient field at the spheromak center: bAmbientCenterSi_D in
  ! code units.  Used for the Borovikov et al. (2018) uniform field subtraction:
  ! interior perturbation = B_Bessel - bAmbientConf_D.
  ! B0 (bConf_D) itself is calibrated against bAmbientConf_D in mc18_get_b0,
  ! optionally including image-dipole feedback if UseImageDipoles is set.
  real :: bAmbientConf_D(3) = 0.0
  !$acc declare create(bAmbientConf_D)

  ! Equivalent external dipole moment of the spheromak (Eq. 3): m = C0*B0,
  ! split into the parts along DirCme_D (mRadDip_D = m_par) and
  ! perpendicular to it (mHorDip_D = m_perp). mDip_D is used for the
  ! exterior field (Sec. 2.1); mRadDip_D/mHorDip_D feed the analytic
  ! image-dipole field (mc18_image_field) once B0 is known.
  real :: C0 = 0.0
  real :: mDip_D(3) = 0.0, mRadDip_D(3) = 0.0, mHorDip_D(3) = 0.0
  !$acc declare create(C0, mDip_D, mRadDip_D, mHorDip_D)

  ! Parameter to control self-similar solution
  real :: uCmeSi = 0.0
  !$acc declare create(uCmeSi)
  real, parameter :: Delta = 0.1

  ! Lin (2006) image dipole parameters
  logical :: UseImageDipoles = .false.
  !$acc declare create(UseImageDipoles)

  ! Spheromak Beta0 and ejecta temperature (same convention as TD99)
  ! Boundary parameter beta0: j1(alpha0*r0)/(alpha0*r0) = beta0.
  ! Pressure from PDF Eq.(5): p = [j1/(a0 r) - b0]*b0*alpha0^2*(r x B0)^2
  logical :: UseBeta0 = .false.
  !$acc declare create(UseBeta0)
  real :: Beta0 = 0.0
  real :: EjectaTemperature = 0.0           ! normalized code units
  real :: EjectaTemperatureDim = 5.0e4      ! [K]
  !$acc declare create(Beta0, EjectaTemperature, EjectaTemperatureDim)

contains
  !============================================================================
  subroutine mc18_init

    use ModCoordTransform, ONLY: cross_product
    !--------------------------------------------------------------------------
    ! Solve boundary condition j1(Alpha0R0)/Alpha0R0 = Beta0.
    ! For force-free (UseBeta0=F), Beta0=0 and Alpha0R0 is
    ! the first zero of j1 (~4.4934).
    Alpha0R0 = find_alpha0r0(merge(Beta0, 0.0, UseBeta0))

    ! Wave number k = Alpha0R0 / r0
    Alpha0 = Alpha0R0/Radius

    ! Center position of the configuration in the heliocentric frame.
    ! Needed already here: mc18_get_b0 below uses it.
    XyzCenterConf_D = rDistance1*DirCme_D

    bAmbientConf_D = bAmbientCenterSi_D*Si2No_V(UnitB_)

    ! Field amplitude and axis: calibrate B0 (bConf_D) from the ambient
    ! field at the spheromak center, optionally including image-dipole
    ! feedback (see mc18_get_b0). bAmbientCenterSi_D is filled by the
    ! MHD solver before mc18_init is called
    ! (see SC_user_initial_perturbation in ModUserAwsom.f90).
    call mc18_get_b0
    B0    = norm2(bConf_D)
    B0Dim = B0*No2Io_V(UnitB_)

    ! Equivalent external dipole moment of the spheromak (Eq. 3): m = C0*B0,
    ! needed for the exterior field regardless of UseImageDipoles; split
    ! into radial/horizontal parts, which is what the analytic image-dipole
    ! field (mc18_image_field) takes as input.
    mDip_D    = C0*bConf_D
    mRadDip_D = sum(mDip_D*DirCme_D)*DirCme_D
    mHorDip_D = mDip_D - mRadDip_D

    ! Helicity: choose sign of alpha so that the toroidal field at the
    ! spheromak bottom opposes the component of the ambient field there
    ! perpendicular to the plane(DirCme_D, bConf_D).
    ! Derivation: b_tor_bot ~ -iHelicity*(DirCme_D x bConf_D), so
    ! iHelicity = sign(bAmbientBottomSi_D · (DirCme_D x bConf_D)).
    iHelicity = sign(1.0, sum(bAmbientBottomSi_D &
         *cross_product(DirCme_D, bConf_D)))

    ! Convert self-similar CME speed from km/s to SI
    uCmeSi = uCmeSi*Io2Si_V(UnitU_)

    if(iProc==0)then
       write(*,*) prefix
       write(*,*) prefix, &
            '>>>>>>>>>>>>>>>>>>>                            '//&
            '<<<<<<<<<<<<<<<<<<<<<'
       write(*,*) prefix
       write(*,*) prefix, &
            '     EEGMC Magnetic Cone Model'//&
            ' (Rosenbluth-Bussac 1979) is initiated'
       write(*,*) prefix
       write(*,*) prefix, &
            '>>>>>>>>>>>>>>>>>>>                            '//&
            '<<<<<<<<<<<<<<<<<<<<<'
       write(*,*) prefix
       write(*,*) prefix, 'B0Dim          = ', B0Dim,               '[Gauss]'
       write(*,*) prefix, 'iHelicity      = ', iHelicity
       write(*,*) prefix, 'Radius         = ', Radius,              '[rSun]'
       write(*,*) prefix, 'HalfOpeningAngle= ', HalfOpeningAngle,   '[degrees]'
       write(*,*) prefix, 'rDistance1     = ', rDistance1,          '[rSun]'
       write(*,*) prefix, 'BaseHeight     = ', BaseHeight,          '[rSun]'
       write(*,*) prefix, 'LongitudeCme   = ', LongitudeCme,        '[degrees]'
       write(*,*) prefix, 'LatitudeCme    = ', LatitudeCme,         '[degrees]'
       write(*,*) prefix, 'Alpha0         = ', Alpha0,              '[1/rSun]'
       write(*,*) prefix, 'UseImageDipoles= ', UseImageDipoles
       write(*,*) prefix, 'UseBeta0  = ', UseBeta0
       if(UseBeta0)then
          write(*,*) prefix, 'Beta0          = ', Beta0
          write(*,*) prefix, 'EjectaTemp     = ', EjectaTemperatureDim, '[K]'
       end if
       write(*,*) prefix, 'Start time     = ', tStartCme,            '[s]'
       write(*,*) prefix, 'CME speed      = ', uCmeSi*Si2Io_V(UnitU_),'[km/s]'
       write(*,*) prefix
    end if

    EjectaTemperature = EjectaTemperatureDim*Io2No_V(UnitTemperature_)

    !$acc update device(Alpha0R0, Alpha0, XyzCenterConf_D, bConf_D, bAmbientConf_D)
    !$acc update device(uCmeSi, B0, Radius, UseImageDipoles)
    !$acc update device(iHelicity, C0, mDip_D, mRadDip_D, mHorDip_D)
    !$acc update device(UseBeta0, Beta0, EjectaTemperature)

  end subroutine mc18_init
  !============================================================================
  subroutine mc18_get_b0

    ! Calibrate B0 (bConf_D). Without image dipoles this reduces to the
    ! original Borovikov et al. (2018) estimate
    !   B0 = -3/(2 j2(a0r0)) Bamb,c.
    ! With image dipoles enabled, the field actually inserted into the MHD
    ! domain is Bamb + (B-Buniform) + Bimg (Eq. 19); requiring this to be
    ! approximately the self-consistent stand-alone configuration B means
    !   Bamb,c + Bimg,c = Buniform(B0),                             (Eq. 21)
    ! i.e. B0 is calibrated so that the uniform term supplies BOTH the
    ! ambient field AND the (separately-inserted) image-dipole field at
    ! the spheromak center - not the other way round. Every image source
    ! (mc18_image_field) is built linearly from either the radial or the
    ! horizontal part of the spheromak's equivalent dipole moment, and all
    ! of them lie on the same solar-center-to-spheromak-center axis as the
    ! spheromak itself, so this calibration decouples exactly into two
    ! independent scalar solves - one along DirCme_D, one in the
    ! perpendicular plane - with no iteration needed.

    real :: bAmbCenter_D(3), Br, bT_D(3), ePerp_D(3)
    real :: GammaUniform, Gr, Gh
    real :: bTestR_D(3), bTestH_D(3)
    real, parameter :: Zero_D(3) = (/0.0, 0.0, 0.0/)
    !--------------------------------------------------------------------------
    bAmbCenter_D = bAmbientCenterSi_D*Si2No_V(UnitB_)

    ! m-per-B0 scalar (Eq. 3): m = C0*B0
    C0 = spher_bessel2(Alpha0R0)*Radius**3/3.0

    ! Response of the uniform-field term alone to B0 (Sec. 2.1)
    GammaUniform = -2.0/3.0*spher_bessel2(Alpha0R0)

    if(.not.UseImageDipoles)then
       bConf_D = bAmbCenter_D/GammaUniform
       RETURN
    end if

    ! Radial/horizontal split of the ambient field at the spheromak center
    Br   = sum(bAmbCenter_D*DirCme_D)
    bT_D = bAmbCenter_D - Br*DirCme_D

    ! Use an arbitrary perpendicular unit vector to obtain the horizontal
    ! image-dipole response.
    ! Any unit vector perpendicular to DirCme_D; the image-dipole response
    ! in the perpendicular plane is isotropic, so the choice doesn't matter.
    ePerp_D = perpendicular_unit_vector(DirCme_D)

    ! Geometric image-dipole response to a unit radial/horizontal moment,
    ! obtained with the same analytic mc18_image_field used for the actual
    ! field (get_mc18_fluxrope), evaluated at the spheromak center.
    bTestR_D = mc18_image_field(XyzCenterConf_D, DirCme_D, Zero_D)
    bTestH_D = mc18_image_field(XyzCenterConf_D, Zero_D, ePerp_D)
    Gr = sum(bTestR_D*DirCme_D)
    Gh = sum(bTestH_D*ePerp_D)

    ! Solve the two decoupled scalar equations for B0 (note the MINUS sign
    ! on the image-dipole term, per Bamb,c + Bimg,c = Buniform, Eq. 21 -
    ! opposite of naively matching Buniform+Bimg to Bamb,c).
    !
    ! Example: r direction:
    ! Bamb,c,r + Bimg,c,r = Buniform(B0),r
    ! => Br*DirCme_D + Gr*DirCme_D*C0*B0,r = GammaUniform*DirCme_D*B0,r
    ! => B0,(Br/(GammaUniform - C0*Gr))
    !
    ! The same for horizontal direction.
    bConf_D = (Br/(GammaUniform - C0*Gr))*DirCme_D &
         + bT_D/(GammaUniform - C0*Gh)

  end subroutine mc18_get_b0
  !============================================================================
  function mc18_image_field(rField_D, mPar_D, mPerp_D) result(b_D)
    !$acc routine seq

    ! Analytic field of the Lin (2006) image-dipole system (Sec. 2.2) at an
    ! arbitrary field point rField_D (heliocentric coordinates, rSun units),
    ! given the radial (mPar_D = m_par, along DirCme_D) and horizontal
    ! (mPerp_D = m_perp, perpendicular to it) parts of the spheromak's
    ! equivalent dipole moment m. Solar radius Rs = 1 in these units.
    ! Combines:
    !  - B1: the radial point image, moment -scale3*mPar_D at
    !    rImg_D = (1/rDistance1)*DirCme_D (Lin 2006 Eq. II.B).
    !  - B2+B3: the horizontal point image plus the continuous line of
    !    horizontal image dipoles between the solar center and rImg_D,
    !    combined and evaluated in closed form via Eq. 24 (Sec. 2.3) rather
    !    than discretizing the line into dipoles or handling B2/B3
    !    separately. Eq. 24, K*gradq*(1/(DS)-1/D^3) +
    !    K*q*(3*gradD/D^4 - gradS/(D*S^2) - gradD/(S*D^2)), is used here
    !    with the gradD coefficients combined algebraically:
    !    3/D^4 - 1/(S*D^2) = (1/D^2)*(3/D^2 - 1/S).

    real, intent(in) :: rField_D(3), mPar_D(3), mPerp_D(3)
    real :: b_D(3)

    real :: dImg, rImg_D(3)
    real :: r, Dvec_D(3), Dist, S, q, gradD_D(3), gradS_D(3), K
    !--------------------------------------------------------------------------
    ! 1.0 is the solar radius in the normalized code units used here.
    dImg   = 1.0/rDistance1
    rImg_D = dImg*DirCme_D

    ! B1: radial point image
    ! -(1.0/rDistance1)**3*mPar_D is the moment of the radial point image (Lin 2006 Eq. II.B)
    b_D = dipole_field(rImg_D, -(1.0/rDistance1)**3*mPar_D, rField_D)

    ! B2+B3: horizontal point image + continuous line, combined
    r       = norm2(rField_D)
    Dvec_D  = rField_D - rImg_D
    Dist    = norm2(Dvec_D)
    S       = r*Dist + r**2 - sum(rField_D*rImg_D)
    q       = sum(mPerp_D*rField_D)
    gradD_D = Dvec_D/Dist
    gradS_D = r*gradD_D + Dist*(rField_D/r) + 2.0*rField_D - rImg_D
    K       = 1.0/rDistance1**3

    b_D = b_D + K*mPerp_D*(1.0/(Dist*S) - 1.0/Dist**3) &
         + K*q*(3.0/Dist**2 - 1.0/S)/Dist**2*gradD_D &
         - K*q/(Dist*S**2)*gradS_D

  end function mc18_image_field
  !============================================================================
  function perpendicular_unit_vector(a_D) result(e_D)
    !$acc routine seq

    ! Returns an arbitrary unit vector perpendicular to a_D, which itself
    ! must be a unit vector.

    real, intent(in) :: a_D(3)
    real :: e_D(3)

    real :: Ref_D(3)
    !--------------------------------------------------------------------------
    if(abs(a_D(1)) < 0.9)then
       Ref_D = (/1.0, 0.0, 0.0/)
    else
       Ref_D = (/0.0, 1.0, 0.0/)
    end if
    e_D = Ref_D - sum(Ref_D*a_D)*a_D
    e_D = e_D/norm2(e_D)

  end function perpendicular_unit_vector
  !============================================================================
  function dipole_field(rSrc_D, m_D, rField_D) result(b_D)
    !$acc routine seq

    ! Field of a point dipole with moment m_D located at rSrc_D, evaluated
    ! at rField_D (all in the normalized code units used throughout this
    ! module, where mu0/(4 pi) is absorbed into m_D, as elsewhere here).

    real, intent(in) :: rSrc_D(3), m_D(3), rField_D(3)
    real :: b_D(3)

    real :: dr_D(3), R2, MdotR
    !--------------------------------------------------------------------------
    dr_D  = rField_D - rSrc_D
    R2    = sum(dr_D**2)
    MdotR = sum(m_D*dr_D)
    b_D   = (3.0*MdotR*dr_D/R2 - m_D)/(sqrt(R2)*R2)

  end function dipole_field
  !============================================================================
  subroutine set_parameters_mc18(NameCommand)

    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand

    real :: SinAlpha
    character(len=*), parameter:: NameSub = 'set_parameters_mc18'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#CME","#MAGCONE")
       call read_var('BaseHeight',      BaseHeight)     ![rSun]
       call read_var('uCmeSi',          uCmeSi)         ![km/s]
       call read_var('UseImageDipoles', UseImageDipoles)
       call read_var('UseBeta0', UseBeta0)
       if(UseBeta0)then
          call read_var('Beta0', Beta0)
          call read_var('EjectaTemperature', EjectaTemperatureDim)
       end if

       ! Derive center distance and radius from opening angle and base height.
       ! Geometry: sin(alpha) = Radius/rDistance1
       !           1 + BaseHeight = rDistance1 - Radius 
       !                          = rDistance1*(1 - sin(alpha))
       if(HalfOpeningAngle <= 0.0 .or. HalfOpeningAngle >= 90.0) call &
            CON_stop(NameSub//': HalfOpeningAngle must be in (0, 90) degrees')
       SinAlpha   = sin(HalfOpeningAngle*cDegToRad)
       rDistance1 = (1.0 + BaseHeight)/(1.0 - SinAlpha)
       Radius     = rDistance1*SinAlpha

       ! position of the CME apex
       XyzCmeApexSi_D = DirCme_D*(rDistance1 + Radius)

       ! position of CME center and bottom
       XyzCmeCenterSi_D = XyzCmeApexSi_D - DirCme_D*Radius
       XyzCmeBottomSi_D = XyzCmeApexSi_D - DirCme_D*2.0*Radius
       DoNormalizeXyz = .true.

    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_mc18
  !============================================================================
  subroutine get_mc18_fluxrope(XyzIn_D, Rho, p, b_D, u_D, TimeNow)
    !$acc routine seq

    ! Magnetic field perturbation of the Rosenbluth-Bussac force-free spheromak.
    ! Interior: spherical Bessel (j1) field minus the uniform ambient field
    !   (Borovikov et al. 2018 uniform-field subtraction for div-B continuity).
    ! Exterior: pure dipole of equivalent moment mDip_D = C0*bConf_D (Eq. 3),
    !   calibrated in mc18_get_b0 (optionally including image-dipole feedback).
    ! Image dipole correction (Lin 2006) applied everywhere if UseImageDipoles.

    use ModCoordTransform, ONLY: cross_product

    ! Coordinates of the input point, in rSun
    real, intent(in) :: XyzIn_D(3)

    ! OUTPUTS
    real, intent(out) :: b_D(3)
    real, intent(out), optional :: u_D(3)
    real, optional, intent(in)  :: TimeNow

    ! Density, pressure (non-zero inside if UseBeta0)
    real, intent(out) :: Rho, p

    real :: XyzConf_D(3), Distance2ConfCenter
    real :: R2CrossB0_D(3), Alpha0R2
    real :: PhiInv
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

       ! INSIDE: force-free spherical Bessel field
       Alpha0R2    = Alpha0*Distance2ConfCenter
       R2CrossB0_D = cross_product(XyzConf_D, bConf_D)
       b_D = (2*bConf_D + iHelicity*Alpha0*R2CrossB0_D) &
            *(spher_bessel1_over_x(Alpha0R2) - Beta0) &
            + spher_bessel2(Alpha0R2)/Distance2ConfCenter**2 &
            *cross_product(XyzConf_D, R2CrossB0_D)          &
            - bAmbientConf_D

       if(UseBeta0)then
          p   = (spher_bessel1_over_x(Alpha0R2) - Beta0) &
               * Beta0 * Alpha0**2 * sum(R2CrossB0_D**2)
          Rho = p/EjectaTemperature
       end if

       if(present(u_D) .and. UseMagCone)then
          u_D = XyzIn_D*PhiInv  &
               *No2Si_V(UnitX_) &
               *uCmeSi*rCmeApexInvSi
       end if

    else

       ! OUTSIDE: pure dipole with equivalent moment mDip_D = C0*bConf_D
       ! (ensures B_r continuity with the subtracted-ambient interior).
       ! XyzConf_D is already relative to the spheromak center.
       b_D = dipole_field((/0.0, 0.0, 0.0/), mDip_D, XyzConf_D)

    end if

    ! Lin (2006) image dipole corrections (heliocentric coordinates, unscaled),
    ! evaluated analytically (Sec. 2.2), not via discretization.
    if(UseImageDipoles) &
         b_D = b_D + mc18_image_field(XyzIn_D, mRadDip_D, mHorDip_D)

    b_D = b_D*No2Si_V(UnitB_)
    Rho = Rho*No2Si_V(UnitRho_)
    p   = p  *No2Si_V(UnitP_)
    ! u_D is already in SI [m/s] from the interior formula above
  end subroutine get_mc18_fluxrope
  !============================================================================
  real function spher_bessel0(x)
    real, intent(in) :: x
    !--------------------------------------------------------------------------
    spher_bessel0 = sin(x)/x
  end function spher_bessel0
  !============================================================================
  real function spher_bessel1(x)
    real, intent(in) :: x
    !--------------------------------------------------------------------------
    spher_bessel1 = (sin(x) - x*cos(x))/x**2
  end function spher_bessel1
  !============================================================================
  real function spher_bessel1_over_x(x)
    real, intent(in) :: x
    !--------------------------------------------------------------------------
    if(x == 0)then
       spher_bessel1_over_x = 1.0/3.0
    else
       spher_bessel1_over_x = (sin(x) - x*cos(x))/x**3
    end if
  end function spher_bessel1_over_x
  !============================================================================
  real function spher_bessel2(x)
    real, intent(in) :: x
    !--------------------------------------------------------------------------
    spher_bessel2 = 3*spher_bessel1_over_x(x) - spher_bessel0(x)
  end function spher_bessel2
  !============================================================================
  real function find_alpha0r0(Beta0In)

    ! Solve j1(x)/x = Beta0In for x = Alpha0R0 by Newton-Raphson.
    ! j1(x)/x is monotone decreasing from 1/3 at x=0 (Taylor limit)
    ! to 0 at x=4.4934 (first zero of j1), so a unique solution
    ! exists for any 0 < Beta0In < 1/3.
    ! Derivative: d/dx[j1(x)/x] = -j2(x)/x (from Bessel recurrence).

    real, intent(in) :: Beta0In

    real    :: x, fx, dfdx
    integer :: i
    integer, parameter :: nIter = 20
    real,    parameter :: Tol   = 1.0e-7
    real,    parameter :: xZero = 4.493409457909064
    character(len=*), parameter :: NameSub = 'find_alpha0r0'
    !--------------------------------------------------------------------------
    if(Beta0In <= 0.0)then
       find_alpha0r0 = xZero
       return
    end if
    if(Beta0In >= 1.0/3.0) call CON_stop( &
         NameSub//': Beta0 >= 1/3, no solution to j1(x)/x = Beta0 exists')

    ! Linear interpolation on (0, xZero) as initial guess
    x = xZero*(1.0 - 3.0*Beta0In)

    do i = 1, nIter
       fx   =  spher_bessel1_over_x(x) - Beta0In
       dfdx = -spher_bessel2(x)/x
       x    =  x - fx/dfdx
       if(abs(fx) < Tol) exit
    end do
    find_alpha0r0 = x

  end function find_alpha0r0
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
