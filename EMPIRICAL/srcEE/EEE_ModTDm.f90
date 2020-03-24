! The 2014 modified Titov-Demoulin flux rope analytical solution with a 
! parabolic current profile. See Titov & Demoulin (1999), and Titov et al. 
! (2014 ApJ 790 163) for equations and explanations of the variables (most of 
! the variable names in this code are consistent with the paper). Here we also 
! transform the flux rope in Cartesian coordinates into spherical. See an 
! example of calling this code at the end. Also see Titov et al. 
! (2018 ApJ 852 L21) for a new method of embeding flux rope into an already 
! obtained background field.

module EEE_ModTDm
  use ModUtilities, ONLY: CON_stop

  implicit none
  save
  private
  public :: set_parameters_tdm14, calc_tdm14_bfield, get_tdm14_apexlocation, init_tdm14

  real, public :: AnalyticTdCoordScale=1., AnalyticTdFieldScale=1.
  real, public :: AnalyticTdCoeffR = 1.83        ! torus major radius
  real, public :: AnalyticTdCoeffA = 0.75        ! torus minor radius
  real, public :: AnalyticTdCoeffD = 0.83        ! depth of torus center
  real, public :: AnalyticTdCoeffL = 1.25        ! magnetic charges distance to center
  ! In spherical: Rotation angle with respect to radial direction (degree), flux rope two foot points align with longitude direction as zero (y axis in Cartesian becomes phi axis), counter clockwise.
  ! In Cartesian: Rotation around z axis, zero is where flux rope lies along y axis
  real, public :: AnalyticTdLocRot = 0.
  ! real, public :: AnalyticTdCoeffQ = 20.         ! strength of magnetic dipole
  real, public :: AnalyticTdLocLon = 50.        ! longitude of flux rope (degree)
  real, public :: AnalyticTdLocLat = 10.        ! Latitude

  logical :: IsInitiated = .false.
  real :: RotateRefGen_DD(3,3)
  logical :: IsCartesian
  real :: rBody = 1.0    ! solar radius when use spherical coordinates

contains

  subroutine set_parameters_tdm14(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand
    character(len=*), parameter:: NameSub = 'set_parameters_tdm14'

    select case(NameCommand)
    case("#TDM14")
      call read_var('AnalyticTdCoordScale', AnalyticTdCoordScale)
      call read_var('AnalyticTdFieldScale', AnalyticTdFieldScale)
      call read_var('AnalyticTdCoeffR', AnalyticTdCoeffR)
      call read_var('AnalyticTdCoeffA', AnalyticTdCoeffA)
      call read_var('AnalyticTdCoeffD', AnalyticTdCoeffD)
      call read_var('Rotation', AnalyticTdLocRot)
      if (.not.IsCartesian) then
        call read_var('Longitude', AnalyticTdLocLon)
        call read_var('Latitude', AnalyticTdLocLat)
      endif
    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select
  end subroutine set_parameters_tdm14


  subroutine init_tdm14(IsCartesianIn)
    use ModCoordTransform, ONLY: rot_matrix_x, rot_matrix_y, rot_matrix_z
    use ModConst, ONLY: cHalfPi, cDegToRad

    logical, intent(in) :: IsCartesianIn
    character(len=*), parameter:: NameSub = 'init_tdm14'

    IsCartesian = IsCartesianIn
    if (IsCartesian) then
      ! Rotation matrix: (TD in reference frame) = RotateRefGen_DD * (TD in requested generic frame)
      RotateRefGen_DD = rot_matrix_z(-AnalyticTdLocRot*cDegToRad)
    else
      ! the same as, but a reversed definition of, the rotation matrix in EEE_ModTd99, easier for me to follow
      RotateRefGen_DD = rot_matrix_z(-AnalyticTdLocLon*cDegToRad)
      RotateRefGen_DD = matmul(rot_matrix_y(AnalyticTdLocLat*cDegToRad),RotateRefGen_DD)
      RotateRefGen_DD = matmul(rot_matrix_x(-AnalyticTdLocRot*cDegToRad),RotateRefGen_DD)
      RotateRefGen_DD = matmul(rot_matrix_y(-cHalfPi),RotateRefGen_DD)
    endif
    IsInitiated = .true.
  end subroutine init_tdm14


  function get_tdm14_apexlocation() result(CoordOut_D)
    real :: CoordOut_D(3)
    real :: CoordRef_D(3)
    character(len=*), parameter:: NameSub = 'get_tdm14_apexlocation'

    if (.not.IsInitiated) then
      call CON_stop(NameSub//' not initiated. Call init_tdm14 first.')
    endif

    CoordRef_D = (/0.,0.,AnalyticTdCoeffR-AnalyticTdCoeffD/)
    CoordOut_D = matmul(CoordRef_D*AnalyticTdCoordScale+(/0.,0.,rBody/),RotateRefGen_DD)
  end function get_tdm14_apexlocation


  subroutine calc_tdm14_bfield(CoordIn_D, bTd_D, bStrap_D)
    real, intent(in) :: CoordIn_D(3)
    real, intent(out) :: bTd_D(3)
    real, intent(in), optional :: bStrap_D(3)

    logical :: IsCartesian
    real :: bStrapRef_D(3), bTdRef_D(3), CoordRef_D(3)
    character(len=*), parameter:: NameSub = 'calc_tdm14_bfield'

    if (.not.IsInitiated) then
      call CON_stop(NameSub//' not initiated. Call init_tdm14 first.')
    endif

    if (IsCartesian) then
      CoordRef_D = matmul(RotateRefGen_DD,CoordIn_D) / AnalyticTdCoordScale
    else
      CoordRef_D = (matmul(RotateRefGen_DD,CoordIn_D)-(/0.,0.,rBody/)) / AnalyticTdCoordScale
    endif
    if (.not.present(bStrap_D)) then
      call calc_tdm14_cart(CoordRef_D, bTdRef_D, 0., .true.)
    else
      ! Project overlying field to along flux rope, the axial flux term should be calculated with respect to a straping field that is perpendicular the flux rope direction 
      bStrapRef_D = matmul(RotateRefGen_DD, bStrap_D)
      call calc_tdm14_cart(CoordRef_D, bTdRef_D, bStrapRef_D(1), .false.)
    endif

    bTd_D = matmul(bTdRef_D,RotateRefGen_DD) * AnalyticTdFieldScale
  end subroutine calc_tdm14_bfield


  ! "Modified Titov-Demoulin flux rope" (Titov+ 2014) with a parabolic profile
  ! current distribution inside the rope. The outside straping potential field,
  ! if not given, is bipolar (generated from two monopoles below the surface).
  subroutine calc_tdm14_cart(Coord_D, bTd_D, bVertIn, DoTestWithPointSource)
    use ModConst, ONLY: cTwoPi, cPi
    use ModHyperGeometric, ONLY: calc_elliptic_int_1kind, calc_elliptic_int_2kind

    real, intent(in) :: Coord_D(3)
    real, intent(out) :: bTd_D(3)
    real, intent(in) :: bVertIn
    logical, intent(in) :: DoTestWithPointSource

    real, parameter :: Delta = 0.1
    ! real, parameter :: cF0 = -0.406982224701535, cF1 = -1.5464309982239, cF2 = -0.249947772314288, cG1 = -2.38261647628
    real, save :: TdR, TdA, TdD
    real, save :: TdL, TdQ
    logical, save :: IsFirst = .true.

    real :: x, y, z
    real :: bTheta_D(3), bQ_D(3), bI_D(3)
    real :: rPlus_D(3), rMinus_D(3), LenRPlus, LenRMinus, bVert
    real :: rVert, Rho, RhoStarF, RhoStarG
    real :: rVertHat_D(3), ThetaHat_D(3), xHat_D(3)
    real :: CurrentI, AxialFlux, Xi, BigTheta
    real :: VarK, dKdX, dKdR, kStarF, dKSFdX, dKSFdR, kStarG, dKSGdX, dKSGdR
    real :: EllipticKf, EllipticEf, FuncAf, dFuncAf, d2FuncAf, d3FuncAf
    real :: EllipticKg, EllipticEg, FuncAg, dFuncAg, d2FuncAg, d3FuncAg
    real :: SewingH, dSewingH, SewingF, dSewingF, SewingG, dSewingG
    real :: Ai, dAIdX, dAIdR
    real :: TmpG, dGdX, dGdR, TmpH, dHdR, Afx, dAFRdX, dAFXdR

    !--------------------------------------------------------------------------

    if (IsFirst) then 
      TdR = AnalyticTdCoeffR
      TdA = AnalyticTdCoeffA
      TdD = AnalyticTdCoeffD
      TdL = TdR/sqrt(2.)/1.03
      TdQ = 1.
      IsFirst = .false.
    endif

    x = Coord_D(1)
    y = Coord_D(2)
    z = Coord_D(3)

    ! ---- point sources field B_q ----
    bQ_D = 0.
    if (DoTestWithPointSource) then
      rPlus_D = (/ x-TdL, y, z+TdD /)
      rMinus_D = (/ x+TdL, y, z+TdD /)
      LenRPlus = sqrt(sum(rPlus_D**2))
      LenRMinus = sqrt(sum(rMinus_D**2))
      bQ_D = TdQ * (rPlus_D/LenRPlus**3 - rMinus_D/LenRMinus**3)
      bVert = - 2*TdQ*TdL/sqrt(TdL**2+TdR**2)**3
    else
      bVert = bVertIn
    endif
    ! In spherical case the straping field magnitude bVert should be provided


    ! ---- variables to calculate toroidal and poloidal vector potentials ----
    rVert = sqrt(y**2+(z+TdD)**2)
    Rho = sqrt(x**2+(rVert-TdR)**2)

    rVertHat_D = (/ 0., y/rVert, (z+TdD)/rVert /)
    ThetaHat_D = (/ 0., -(z+TdD)/rVert, y/rVert /)
    xHat_D = (/ 1., 0., 0. /)

    VarK = sqrt((rVert*TdR)/(rVert*TdR+Rho**2/4.))
    dKdX = - (x*VarK**3) / (4*rVert*TdR)
    dKdR = VarK**3/(8*rVert**2*TdR) * (Rho**2-2*rVert*(rVert-TdR))

    CurrentI = - (4*cPi*TdR*bVert) / (log(8*TdR/TdA)-25./24.)
    AxialFlux = 3./(5*sqrt(2.)) * CurrentI * TdA

    ! Sewing functions
    Xi = (Rho-TdA)/(Delta*TdA)
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
    RhoStarF = TdA*(1+Delta*SewingF)
    kStarF = sqrt((rVert*TdR)/(rVert*TdR+RhoStarF**2/4.))
    dKSFdX = - (x*kStarF**3) / (4*rVert*TdR) * dSewingF
    dKSFdR = kStarF**3/(8*rVert**2*TdR) * (RhoStarF**2-2*rVert*(rVert-TdR)*dSewingF)

    call calc_elliptic_int_1kind(kStarF,EllipticKf)
    call calc_elliptic_int_2kind(kStarF,EllipticEf)

    FuncAf = ((2-kStarF**2)*EllipticKf - 2*EllipticEf) / kStarF
    dFuncAf = (2-kStarF**2)/(kStarF**2*(1-kStarF**2)) * EllipticEf - 2/(kStarF**2) * EllipticKf
    d2FuncAf = -(kStarF**4-7*kStarF**2+4)/(kStarF**3*(1-kStarF**2)**2) * EllipticEf - &
         (5*kStarF**2-4)/(kStarF**3*(1-kStarF**2)) * EllipticKf
    d3FuncAf = -(2*kStarF**6-31*kStarF**4+33*kStarF**2-12)/(kStarF**4*(1-kStarF**2)**3) * EllipticEf - &
         (19*kStarF**4-27*kStarF**2+12)/(kStarF**4*(1-kStarF**2)**2) * EllipticKf

    ! curly-A function for k_(five-sided-star)
    RhoStarG = TdA*(1+Delta*SewingG)
    kStarG = sqrt((rVert*TdR)/(rVert*TdR+RhoStarG**2/4.))
    dKSGdX = - (x*kStarG**3) / (4*rVert*TdR) * dSewingG
    dKSGdR = kStarG**3/(8*rVert**2*TdR) * (RhoStarG**2-2*rVert*(rVert-TdR)*dSewingG)

    call calc_elliptic_int_1kind(kStarG,EllipticKg)
    call calc_elliptic_int_2kind(kStarG,EllipticEg)

    FuncAg = ((2-kStarG**2)*EllipticKg - 2*EllipticEg) / kStarG
    dFuncAg = (2-kStarG**2)/(kStarG**2*(1-kStarG**2)) * EllipticEg - &
         2/(kStarG**2) * EllipticKg
    d2FuncAg = -(kStarG**4-7*kStarG**2+4)/(kStarG**3*(1-kStarG**2)**2) * EllipticEg - &
         (5*kStarG**2-4)/(kStarG**3*(1-kStarG**2)) * EllipticKg
    d3FuncAg = -(2*kStarG**6-31*kStarG**4+33*kStarG**2-12)/(kStarG**4*(1-kStarG**2)**3) * EllipticEg - &
         (19*kStarG**4-27*kStarG**2+12)/(kStarG**4*(1-kStarG**2)**2) * EllipticKg


    ! ---- ring current field B_I ----

    Ai = CurrentI/cTwoPi*sqrt(TdR/rVert) * &
         (FuncAf+dFuncAf*(VarK-kStarF)+0.5*d2FuncAf*(VarK-kStarF)**2)
    dAIdX = CurrentI/cTwoPi*sqrt(TdR/rVert) * &
         (dFuncAf*dKdX+d2FuncAf*dKdX*(VarK-kStarF)+0.5*d3FuncAf*dKSFdX*(VarK-kStarF)**2)
    dAIdR = CurrentI/cTwoPi*sqrt(TdR/rVert) * &
         (dFuncAf*dKdR+d2FuncAf*dKdR*(VarK-kStarF)+0.5*d3FuncAf*dKSFdR*(VarK-kStarF)**2) - Ai/(2*rVert)

    bI_D = - dAIdX * rVertHat_D + (dAIdR+Ai/rVert) * xHat_D

    ! ---- toroidal field B_theta ----

    TmpG = 3+4*dFuncAf*(VarK-kStarF)    ! just a temporary variable, same for tmpH below
    dGdX = 4*(d2FuncAf*dKSFdX*(VarK-kStarF)+dFuncAf*(dKdX-dKSFdX))
    dGdR = 4*(d2FuncAf*dKSFdR*(VarK-kStarF)+dFuncAf*(dKdR-dKSFdR))

    TmpH = (VarK**3*(x**2+TdR**2-rVert**2)-TdA**2*kStarG**3)*dFuncAg + &
         TdA**2*kStarG**3*d2FuncAg*(VarK-kStarG)
    dHdR = (3*VarK**2*dKdR*(x**2+TdR**2-rVert**2)-2*VarK**3*rVert-&
         3*TdA**2*kStarG**2*dKSGdR)*dFuncAg + &
         (VarK**3*(x**2+TdR**2-rVert**2)-TdA**2*kStarG**3)*d2FuncAg*dKSGdR + &
         TdA**2*( (3*kStarG**2*dKSGdR*(VarK-kStarG)+kStarG**3*(dKdR-dKSGdR))*d2FuncAg + &
         kStarG**3*(VarK-kStarG)*d3FuncAg*dKSGdR )

    Afx = AxialFlux/(4*cPi*rVert)*sqrt(TdR/rVert) * &
         ( FuncAg + (TdA**2*kStarG**3)/(4*rVert*TdR)*dFuncAg + TmpG**(5./2.)/(30*sqrt(3.)) &
         - 0.3 + TmpG**(3./2.)/(12*sqrt(3.)*rVert*TdR)*TmpH )

    dAFRdX = AxialFlux/(24*sqrt(3.)*cPi*rVert)/sqrt(rVert*TdR) * &
         (1.5*sqrt(TmpG)*dGdX*x*VarK**3*dFuncAf + &
         sqrt(TmpG)**3*(x*VarK**3*d2FuncAf*dKSFdX+(VarK**3+3*x*VarK**2*dKdX)*dFuncAf))
    dAFXdR = (AxialFlux*sqrt(TdR))/(4*cPi)*rVert**(-3./2.) * ( dFuncAg*dKSGdR + &
         TdA**2/(4*TdR)*((3*kStarG**2*rVert*dKSGdR-kStarG**3)/(rVert**2)*dFuncAg + &
         kStarG**3/rVert*d2FuncAg*dKSGdR) + TmpG**(3./2.)/(12*sqrt(3.))*dGdR + &
         1./(12*sqrt(3.)*TdR)*((1.5*sqrt(TmpG)*dGdR*rVert-TmpG**(3./2.))/(rVert**2)*TmpH + &
         TmpG**(3./2.)/rVert*dHdR) ) - 3./(2*rVert)*Afx

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
!     call init_tdm14(IsCartesian)
!     Apex_D = get_tdm14_apexlocation()
!     call user_get_b0(Apex_D(1),Apex_D(2),Apex_D(3),bStrap_D)
!     IsFirst = .false.
!   endif
!   do k=MinK, MaxK; do j=MinJ, MaxJ; do i=MinI, MaxI
!     call calc_tdm14_bfield(Xyz_DGB(:,i,j,k,iBlock), B_D, bStrap_D)
!     State_VGB(Bx_:Bz_,i,j,k,iBlock) = B_D
!   enddo; enddo; enddo
