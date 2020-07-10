!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module EEE_ModTDm

  ! The 2014 modified Titov-Demoulin flux rope analytical solution with a 
  ! parabolic current profile. See Titov & Demoulin (1999), and Titov et al. 
  ! (2014 ApJ 790 163) for equations and explanations of the variables (most of
  ! the variable names in this code are consistent with the paper).
  ! Here we also transform the flux rope in Cartesian coordinates into
  ! spherical. See an example of calling this code at the end. Also see Titov
  ! et al. (2018 ApJ 852 L21) for a new method of embeding flux rope into an
  ! already obtained background field.

  use EEE_ModCommonVariables
  implicit none
  save
  private

  public :: set_parameters_tdm14
  public :: calc_tdm14_bfield

  real, public :: AnalyticTdCoordScale=1., AnalyticTdFieldScale=1.
  real, public :: RTube = 1.83 ! torus major radius
  real, public :: ATube = 0.75 ! torus minor radius
  real, public :: Depth = 0.83 ! depth of torus center

  real, public :: XyzApex_D(3), Bstrap_D(3)

  real :: Rotate_DD(3,3)
  real :: rBody = 1.0    ! solar radius when use spherical coordinates

contains
  !============================================================================

  subroutine set_parameters_tdm14(NameCommand)

    use ModUtilities, ONLY: CON_stop
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand
    character(len=*), parameter:: NameSub = 'set_parameters_tdm14'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#CME")
       call read_var('AnalyticTdCoordScale', AnalyticTdCoordScale)
       call read_var('AnalyticTdFieldScale', AnalyticTdFieldScale)
       call read_var('RTube', RTube)
       call read_var('ATube', ATube)
       call read_var('Depth', Depth)
    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

    call init_tdm14

  end subroutine set_parameters_tdm14
  !============================================================================
  subroutine init_tdm14

    use ModCoordTransform, ONLY: rot_matrix_x, rot_matrix_y, rot_matrix_z
    use ModConst, ONLY: cHalfPi, cDegToRad

    real :: CoordRef_D(3)

    character(len=*), parameter:: NameSub = 'init_tdm14'
    !--------------------------------------------------------------------------
    Rotate_DD = rot_matrix_z(-LongitudeCme*cDegToRad)
    Rotate_DD = matmul(rot_matrix_y(LatitudeCme*cDegToRad),Rotate_DD)
    Rotate_DD = matmul(rot_matrix_x(-OrientationCme*cDegToRad),Rotate_DD)
    Rotate_DD = matmul(rot_matrix_y(-cHalfPi),Rotate_DD)

    CoordRef_D = (/ 0., 0., RTube-Depth /)
    XyzApex_D =  &
         matmul(CoordRef_D*AnalyticTdCoordScale+(/0.,0.,rBody/), Rotate_DD)

  end subroutine init_tdm14
  !============================================================================
  subroutine calc_tdm14_bfield(CoordIn_D, bTd_D, bStrap_D)

    real, intent(in) :: CoordIn_D(3)
    real, intent(out) :: bTd_D(3)
    real, intent(in), optional :: bStrap_D(3)

    real :: bStrapRef_D(3), bTdRef_D(3), CoordRef_D(3)
    character(len=*), parameter:: NameSub = 'calc_tdm14_bfield'
    !--------------------------------------------------------------------------

    CoordRef_D = (matmul(Rotate_DD,CoordIn_D)-(/0.,0.,rBody/)) &
         /AnalyticTdCoordScale

    if (present(bStrap_D)) then
       ! Project overlying field to along flux rope, the axial flux term
       ! should be calculated with respect to a straping field that is
       ! perpendicular the flux rope direction 
       bStrapRef_D = matmul(Rotate_DD, bStrap_D)
       call calc_tdm14_cart(CoordRef_D, bTdRef_D, bStrapRef_D(1), .false.)
    else
       call calc_tdm14_cart(CoordRef_D, bTdRef_D, 0., .true.)
    endif

    bTd_D = matmul(bTdRef_D,Rotate_DD) * AnalyticTdFieldScale*No2Si_V(UnitB_)

  end subroutine calc_tdm14_bfield
  !============================================================================
  ! "Modified Titov-Demoulin flux rope" (Titov+ 2014) with a parabolic profile
  ! current distribution inside the rope. The outside straping potential field,
  ! if not given, is bipolar (generated from two monopoles below the surface).
  subroutine calc_tdm14_cart(Coord_D, bTd_D, bVertIn, DoTestWithPointSource)

    use ModConst, ONLY: cTwoPi, cPi
    use ModHyperGeometric, ONLY: calc_elliptic_int_1kind, &
         calc_elliptic_int_2kind

    real, intent(in)    :: Coord_D(3)
    real, intent(out)   :: bTd_D(3)
    real, intent(in)    :: bVertIn
    logical, intent(in) :: DoTestWithPointSource

    real, parameter :: Delta = 0.1
    real, save      :: RTube, aTube, Depth
    real, save      :: TdL, TdQ
    logical, save   :: IsFirst = .true.

    real :: x, y, z
    real :: bTheta_D(3), bQ_D(3), bI_D(3)
    real :: rPlus_D(3), rMinus_D(3), LenRPlus, LenRMinus, bVert
    real :: RPerp, RMinus, RhoStarF, RhoStarG
    real :: RPerpHat_D(3), ThetaHat_D(3), xHat_D(3)
    real :: CurrentI, AxialFlux, Xi, BigTheta
    real :: Kappa, DKappaDx, DKappaDRperp, kStarF, dKSFdX, dKSFdR, kStarG, dKSGdX, dKSGdR
    real :: EllipticKf, EllipticEf, FuncAf, dFuncAf, d2FuncAf, d3FuncAf
    real :: EllipticKg, EllipticEg 
    real :: FuncAg, dFuncAg, d2FuncAg, d3FuncAg
    real :: SewingH, dSewingH, SewingF, dSewingF, SewingG, dSewingG
    real :: Ai, dAIdX, dAIdR
    real :: TmpG, dGdX, dGdR, TmpH, dHdR, Afx, dAFRdX, dAFXdR
    !--------------------------------------------------------------------------

    if (IsFirst) then 
       TdL = RTube/sqrt(2.)/1.03
       TdQ = 1.
       IsFirst = .false.
    endif

    x = Coord_D(1)
    y = Coord_D(2)
    z = Coord_D(3)

    ! ---- point sources field B_q ----
    bQ_D = 0.
    if (DoTestWithPointSource) then
       rPlus_D = (/ x-TdL, y, z+Depth /)
       rMinus_D = (/ x+TdL, y, z+Depth /)
       LenRPlus = sqrt(sum(rPlus_D**2))
       LenRMinus = sqrt(sum(rMinus_D**2))
       bQ_D = TdQ * (rPlus_D/LenRPlus**3 - rMinus_D/LenRMinus**3)
       bVert = - 2*TdQ*TdL/sqrt(TdL**2+RTube**2)**3
    else
       bVert = bVertIn
    endif
    ! In spherical case the straping field magnitude bVert should be provided

    ! ---- variables to calculate toroidal and poloidal vector potentials ----
    RPerp  = sqrt(y**2 + (z + Depth)**2)
    RMinus = sqrt(x**2+(RPerp-RTube)**2)
    !
    !Coordinate unit vectors
    !
    RPerpHat_D = (/ 0., y/RPerp, (z+Depth)/RPerp /)
    ThetaHat_D = (/ 0., -(z+Depth)/RPerp, y/RPerp /)
    xHat_D = (/ 1., 0., 0. /)
    !
    ! Toroidal coordinate, usual argument of special functions
    ! describing solution in the totoidal coordinates
    Kappa = sqrt((RPerp*RTube)/(RPerp*RTube+RMinus**2/4.))
    !
    ! Derivatives over x and over RPerp
    DKappaDx = - (x*Kappa**3) / (4*RPerp*RTube)
    DKappaDRperp = Kappa**3/(8*RPerp**2*RTube) * (RMinus**2-2*RPerp*(RPerp-RTube))

    CurrentI = - (4*cPi*RTube*bVert) / (log(8*RTube/aTube)-25./24.)
    AxialFlux = 3./(5*sqrt(2.)) * CurrentI * aTube

    ! Sewing functions
    Xi = (RMinus-aTube)/(Delta*aTube)
    SewingH = 0.5*(Xi+log(2*cosh(Xi)))
    dSewingH = 0.5*(1+tanh(Xi))
    ! BigTheta = cPi/4*(1+tanh(Xi))
    SewingF = SewingH    ! approximation for parabolic current case
    ! SewingF = SewingH + cF0*exp(cF1*SewingH+cF2*SewingH**2)
    dSewingF = dSewingH
    ! dSewingF = sin(BigTheta)
    ! dSewingF = dSewingH + cF0*exp(cF1*SewingH+cF2*SewingH**2)*(cF1*dSewingH+2*cF2*SewingH*dSewingH)
    SewingG = SewingH
    ! SewingG = SewingH - cF0*exp(cG1*SewingH)
    dSewingG = dSewingH
    ! dSewingG = 1-cos(BigTheta)
    ! dSewingG = dSewingH - cF0*exp(cG1*SewingH)*cG1*dSewingH

    ! curly-A function and its derivatives for k_(six-edged-star)
    RhoStarF = aTube*(1+Delta*SewingF)
    kStarF = sqrt((RPerp*RTube)/(RPerp*RTube+RhoStarF**2/4.))
    dKSFdX = - (x*kStarF**3) / (4*RPerp*RTube) * dSewingF*RhoStarF/RMinus
    dKSFdR = kStarF**3/(8*RPerp**2*RTube) &
         * (RhoStarF**2 - 2*RPerp*(RPerp-RTube)*dSewingF*RhoStarF/RMinus)

    call calc_elliptic_int_1kind(kStarF,EllipticKf)
    call calc_elliptic_int_2kind(kStarF,EllipticEf)

    FuncAf = ((2-kStarF**2)*EllipticKf - 2*EllipticEf) / kStarF
    dFuncAf = (2-kStarF**2)/(kStarF**2*(1-kStarF**2)) * EllipticEf &
         - 2/(kStarF**2) * EllipticKf
    d2FuncAf = -(kStarF**4-7*kStarF**2+4)/(kStarF**3*(1-kStarF**2)**2) &
         *EllipticEf - (5*kStarF**2-4)/(kStarF**3*(1-kStarF**2))*EllipticKf
    d3FuncAf = -(2*kStarF**6-31*kStarF**4+33*kStarF**2-12) &
         /(kStarF**4*(1-kStarF**2)**3) * EllipticEf &
         - (19*kStarF**4-27*kStarF**2+12) &
         /(kStarF**4*(1-kStarF**2)**2) * EllipticKf

    ! curly-A function for k_(five-sided-star)
    RhoStarG = aTube*(1+Delta*SewingG)
    kStarG = sqrt((RPerp*RTube)/(RPerp*RTube+RhoStarG**2/4.))
    dKSGdX = - (x*kStarG**3) / (4*RPerp*RTube) * dSewingG*RhoStarG/RMinus
    dKSGdR = kStarG**3/(8*RPerp**2*RTube) &
         * (RhoStarG**2 - 2*RPerp*(RPerp-RTube)*dSewingG*RhoStarG/RMinus)

    call calc_elliptic_int_1kind(kStarG,EllipticKg)
    call calc_elliptic_int_2kind(kStarG,EllipticEg)

    FuncAg = ((2-kStarG**2)*EllipticKg - 2*EllipticEg) / kStarG
    dFuncAg = (2-kStarG**2)/(kStarG**2*(1-kStarG**2)) * EllipticEg - &
         2/(kStarG**2) * EllipticKg
    d2FuncAg = -(kStarG**4-7*kStarG**2+4)/(kStarG**3*(1-kStarG**2)**2) &
         *EllipticEg - (5*kStarG**2-4)/(kStarG**3*(1-kStarG**2))*EllipticKg
    d3FuncAg = -(2*kStarG**6-31*kStarG**4+33*kStarG**2-12) &
         /(kStarG**4*(1-kStarG**2)**3)*EllipticEg &
         - (19*kStarG**4-27*kStarG**2+12) &
         /(kStarG**4*(1-kStarG**2)**2)*EllipticKg

    ! ---- ring current field B_I ----
    !
    ! A phi-component of vector potential
    Ai = CurrentI/cTwoPi*sqrt(RTube/RPerp)* &
         (FuncAf + dFuncAf*(Kappa - kStarF)+0.5*d2FuncAf*(Kappa - kStarF)**2)
    !
    ! Its partial derivatives
    !
    dAIdX = CurrentI/cTwoPi*sqrt(RTube/RPerp)*&
         (dFuncAf*DKappaDx + d2FuncAf*DKappaDx*(Kappa - kStarF) &
         + 0.5*d3FuncAf*dKSFdX*(Kappa - kStarF)**2)
    dAIdR = CurrentI/cTwoPi*sqrt(RTube/RPerp) &
         *(dFuncAf*DKappaDRperp + d2FuncAf*DKappaDRperp*(Kappa-kStarF) &
         + 0.5*d3FuncAf*dKSFdR*(Kappa-kStarF)**2) - Ai/(2*RPerp)
    !
    ! Poloidal magnetic field
    bI_D = - dAIdX*RPerpHat_D + (dAIdR + Ai/RPerp)*xHat_D

    ! ---- toroidal field B_theta ----

    ! just a temporary variable, same for tmpH below
    TmpG = 3 + 4*dFuncAf*(Kappa-kStarF)   

    dGdX = 4*(d2FuncAf*dKSFdX*(Kappa-kStarF)+dFuncAf*(DKappaDx-dKSFdX))
    dGdR = 4*(d2FuncAf*dKSFdR*(Kappa-kStarF)+dFuncAf*(DKappaDRperp-dKSFdR))

    TmpH = (Kappa**3*(x**2+RTube**2-RPerp**2)-aTube**2*kStarG**3)*dFuncAg + &
         aTube**2*kStarG**3*d2FuncAg*(Kappa-kStarG)
    dHdR = (3*Kappa**2*DKappaDRperp*(x**2+RTube**2-RPerp**2)-2*Kappa**3*RPerp-&
         3*aTube**2*kStarG**2*dKSGdR)*dFuncAg + &
         (Kappa**3*(x**2+RTube**2-RPerp**2)-aTube**2*kStarG**3)*d2FuncAg*dKSGdR + &
         aTube**2*( (3*kStarG**2*dKSGdR*(Kappa - kStarG) &
         + kStarG**3*(DKappaDRperp - dKSGdR))*d2FuncAg &
         + kStarG**3*(Kappa - kStarG)*d3FuncAg*dKSGdR )

    Afx = AxialFlux/(4*cPi*RPerp)*sqrt(RTube/RPerp) * &
         ( FuncAg + (aTube**2*kStarG**3)/(4*RPerp*RTube)*dFuncAg &
         + TmpG**(5./2.)/(30*sqrt(3.)) &
         - 0.3 + TmpG**(3./2.)/(12*sqrt(3.)*RPerp*RTube)*TmpH )

    dAFRdX = AxialFlux/(24*sqrt(3.)*cPi*RPerp)/sqrt(RPerp*RTube) * &
         (1.5*sqrt(TmpG)*dGdX*x*Kappa**3*dFuncAf + &
         sqrt(TmpG)**3*(x*Kappa**3*d2FuncAf*dKSFdX &
         + (Kappa**3+3*x*Kappa**2*DKappaDx)*dFuncAf))
    dAFXdR = (AxialFlux*sqrt(RTube))/(4*cPi)*RPerp**(-3./2.) &
         * ( dFuncAg*dKSGdR + &
         aTube**2/(4*RTube)*((3*kStarG**2*RPerp*dKSGdR-kStarG**3) &
         /(RPerp**2)*dFuncAg + &
         kStarG**3/RPerp*d2FuncAg*dKSGdR) &
         + TmpG**(3./2.)/(12*sqrt(3.))*dGdR + &
         1./(12*sqrt(3.)*RTube)*((1.5*sqrt(TmpG)*dGdR*RPerp-TmpG**(3./2.)) &
         /(RPerp**2)*TmpH + &
         TmpG**(3./2.)/RPerp*dHdR) ) - 3./(2*RPerp)*Afx

    bTheta_D = (dAFRdX-dAFXdR) * ThetaHat_D

    ! ---- combine three parts ----
    bTd_D = bQ_D + bI_D + bTheta_D

  end subroutine calc_tdm14_cart

end module EEE_ModTDm


! ! A minimal example of calling this module is attached

! ! first set parameter for TDm14
! subroutine user_read_inputs
!   case("#TDM14")
!     call set_parameters_tdm14(NameCommand)

! ! This is an example of a tilted dipole as the background field
! subroutine user_get_b0(x, y, z, B0_D)
!   DpM_D = AnalyticDpM*matmul(rot_matrix_z(AnalyticDpLon*cDegToRad), matmul(rot_matrix_y(cHalfPi-AnalyticDpLat*cDegToRad), (/0.,0.,1./)))
!   r = sqrt(x**2+y**2+z**2)
!   rHat_D = (/x,y,z/)/r
!   B0_D = (3*rHat_D*dot_product(DpM_D,rHat_D)-DpM_D)/r**3

! ! Call the initiation function before calling calculation of the flux rope. When setting up the flux rope, we will need the strength of the straping field at the apex of the flux rope torus. This could be calculated directly from user_get_b0, or may be retrieved from a restart file, init file, etc. Ideally the flux rope direction (approximately the direction of the polarity inversion line) should be perpendicular to the direction of the straping field, such that things are symetric. Here I'm not checking that so you can have more flexibility with the flux rope, which is not consistent with the assumptions in the paper and could lead to a crash. Refer to the new paper Titov et al. (2018) "Regularized Biot-Savart Laws" for a better method to generate a flux rope.
! subroutine user_set_ics(iBlock)
!   if (IsFirst) then
!     call init_tdm14
!     Apex_D = get_tdm14_apexlocation()
!     call user_get_b0(Apex_D(1),Apex_D(2),Apex_D(3),bStrap_D)
!     IsFirst = .false.
!   endif
!   do k=MinK, MaxK; do j=MinJ, MaxJ; do i=MinI, MaxI
!     call calc_tdm14_bfield(Xyz_DGB(:,i,j,k,iBlock), B_D, bStrap_D)
!     State_VGB(Bx_:Bz_,i,j,k,iBlock) = B_D
!   enddo; enddo; enddo
