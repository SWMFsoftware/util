!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module EEE_ModTD99
!  use ModUtilities, ONLY: norm2
  use EEE_ModCommonVariables
  use ModHyperGeometric
  implicit none
  save
  private

  public :: set_parameters_TD99, get_TD99_fluxrope

  !
  ! Logical variables related to the magnetic field computation:
  logical, public :: DoBqField=.false.

  logical :: DoEquilItube=.false.
  real    :: StartTime

  !
  ! Variables related to the position of the flux rope::
  real, parameter :: Li=0.5

  !
  ! Variables related to the flux rope properties::
  real :: ITube = 0.0, RTube  = 0.0, aTube = 0.0
  real :: Depth = 0.0, aRatio = 0.0
  real :: Mass = 0.0
  real :: InvH0=0.0, Rho0=0.0
  real :: ItubeSaved=0.0

  !
  ! Variables related to the properties of the strapping field, Bq::
  real :: VTransX=1.500E+03             !in [m/s]
  real :: VTransY=-2.900E+04            !in [m/s]
  real :: q=0.0, L=0.0

  ! direction of the magnetic moment of the configuration
  real :: UnitX_D(3)
  ! To be renamed: magnetic field at the center of configuration
  real :: B0
  ! Magnetic configuration center
  real :: XyzCenter_D(3)
contains
  !============================================================================

  subroutine set_parameters_TD99(NameCommand)
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'set_parameters_TD99'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#TD99FLUXROPE")
       call read_var('UseCme',              UseCme)
       call read_var('DoAddFluxRope'      , DoAddFluxRope)
       call read_var('DoEquilItube'       , DoEquilItube)
       call read_var('Itube'              , Itube)
       call read_var('Rtube'              , Rtube)
       call read_var('aTube'              , aTube)
       call read_var('Depth'              , Depth)
       call read_var('Mass'               , Mass)
       call read_var('LongitudeCme'       , LongitudeCme)
       call read_var('LatitudeCme'        , LatitudeCme)
       call read_var('OrientationCme'     , OrientationCme)
       call read_var('DoBqField'          , DoBqField)
       call read_var('q'             , q)
       call read_var('L'             , L)
    case("#CME")
       call read_var('Current',     Itube)
       call read_var('RadiusMajor', Rtube)
       call read_var('RadiusMinor', aTube)
       call read_var('Depth',       Depth)
       call read_var('Mass',        Mass)
    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_TD99
  !=================================
  subroutine get_TD99_fluxrope(Xyz_D,BFRope_D,RhoFRope)
    use ModCoordTransform, ONLY: cross_product
    !\__                                                             __/!
    !    Twisted Magnetic Field Configuration by Titov & Demoulin '99   !
    !                                                                   !
    ! An instability that causes a CME eruption is expected to occur at !
    ! R > L*sqrt(2). For a detailed description of the initial state    !
    ! refer to A&A, 1999, v.351, pp.707-720                             !
    !                                                                   !
    ! ___  This module was written by Ilia Roussev on June 10, 2002 ___ !
    !/                                                                 \!
    real, dimension(3), intent(in) :: Xyz_D(3)
    real, dimension(3), intent(out) :: BFRope_D
    real, intent(out), optional :: RhoFRope
    !\
    ! Coordinates relative to the configuration center:
    !/
    real :: XyzRel_D(3)
    !\
    ! Distance to the configuration center squared:
    !/
    real :: R2 !=sum(XyzRel_D**2)

    real:: xRel, R2Face
    real:: RMinus, RPlus2, Rperp
    real:: ThetaUVy,ThetaUVz
    real:: Kappa,dKappadx,dKappadr
    real:: KappaA,dKappaAdr
    !Values of the hypergeometric functions, 
    !2F1(3/2, 3/2 ; 3; Kappa**2) and 2F1(1/2, 1/2; 1; Kappa**3)
    real:: F32323, F12121

    ! Complete elliptic integrals of related variables::
    real:: KElliptic, EElliptic
    ! Vector potential related variables::
    real:: Ak,dAkdk
    real:: AkA,dAkdkA,d2Akdk2A
    real:: AI,dAIdx,dAIdr
    ! Flux-rope related variables::
    real:: BIPhi_D(3)
    
    logical:: DoFirstCall=.true.
    real   :: B1qField_D(3)
    !--------------------------------------------------------------------------
    ! Initialize the TD99 model parameters once::

    if (DoFirstCall) then
       call init_TD99_parameters
       DoFirstCall=.false.
    endif
   
    !
    ! Calculate the field of thin toroidal current filament
    ! 
    
    !--------------------------------------------------------------------
    ! Assign X,Y,Z coordinates at which to compute the magnetic field::
    !
    ! xRel - coordinate along the axis of symmetry, equal to 0 at the 
    ! plane of symmetry.
    ! zzz - is the heliocetric coordinate along the line passing through 
    ! the center of configuration, which is at the depth d below the 
    ! photosphere level.
    ! yyy coordinate is equal to zero at the axis of symmetry. 
    XyzRel_D = Xyz_D - XyzCenter_D
    R2 = sum(XyzRel_D**2); xRel = sum(XyzRel_D*UnitX_D)
    R2Face = norm2(Xyz_D)

    ! Compute Rperp and TubalDist::

    Rperp = sqrt(R2 - xRel**2)
    RMinus = sqrt(xRel**2 + (Rperp - Rtube)**2)
    RPlus2 = (Rperp + Rtube)**2 + xRel**2


    ! Define the model input, Kappa

    Kappa = 2.0*sqrt(Rperp*Rtube/RPlus2)
    if (RMinus >= aTube) then
       ! Compute the field and density outside the current torus   
       !\
       !Initialize redundant output
       !/
       if (present(RhoFRope))RhoFRope=0.0   

       ! Truncate the value of Kappa::  
       if (abs(1.0-Kappa).lt.cTiny/10.0) &
               Kappa = 1.0-cTiny/10.0

       ! Compute the vector potential, Ak, of the magnetic field 
       ! produced by the ring current Itube and its derivatives
             
       !Calculate F(3/2, 3/2; 3; Kappa**2)
       F32323 = hyper_semi_semi_int(1, 1, 3, Kappa**2)
       !Calculate \tilde{P}^{-1}_{-1/2}(\cosh u)
       Ak     = 0.250*F32323
       F12121 = hyper_semi_semi_int(0, 0, 1, Kappa**2)
       !Calculate d\tilde{P}^{-1}_{-1/2}(\cosh u)/du
       dAkDk  = (F12121 - 0.1250*(2.0 - Kappa**2)*F32323)/&
               (1.0 - Kappa**2)
       
       BFRope_D =B0*(Rtube/sqrt(RPlus2))**3*&
            (dAkDk*(2*xRel*XyzRel_D + (RTube**2 - R2)*UnitX_D)/RPlus2 &
            + Ak*UnitX_D)
       !No toroidal field outside the filament
       BIPhi_D = 0.0
    else
       !\
       ! Compute the field and density inside the current torus
       !/
       ! 1.
       ! Add the prominence material inside the flux rope, assuming that the
       ! total amount mass is 10^13kg, and that the desnity scale-height is
       ! the same as the pressure scale-height, 1/InvH0 (i.e., iso-thermal
       ! atmoshpere)::
       
       if (present(RhoFRope))&
            RhoFRope = Rho0*exp(-10.0*(RMinus/aTube)**6) &
            *exp(-InvH0*abs(R2Face-1.0))
       !2.    
       ! Define the model input, KappaA. A given point is characterized by
       ! two coordinates, say, RPepr and RMinus. In this case,
       ! if we denote Kappa = known_function(RPerp,RMinus), then 
       ! KappaA = known_function(RPepr,ATube)
       KappaA = 2.0*sqrt(Rperp*Rtube &
            /(4.0*Rperp*Rtube + aTube**2))
       call calc_elliptic_int_1kind(KappaA,KElliptic)
       call calc_elliptic_int_2kind(KappaA,EElliptic)
       !\
       ! This is the vector potential and its k derivative at the boundary
       !/
       Ak      = ((2.0-KappaA**2)*KElliptic &
            - 2.0*EElliptic)/KappaA
       dAkdk   = (2.0-KappaA**2)*EElliptic &
            /(KappaA**2*(1.0-KappaA**2)) &
            - 2.0*KElliptic/KappaA**2
       ! Compute the vector potential, AI by smoothly continuing the 
       ! value from the boundary
       AI = Ak + dAkdk*(Kappa - KappaA)
       d2Akdk2A = ((7.0*KappaA**2 - 4.0 &
            - KappaA**4)*EElliptic/(1.0 - KappaA**2) &
            + (4.0 - 5.0*KappaA**2)*KElliptic) &
            /(KappaA**3*(1.0-KappaA**2))
       ! Derive the BI field from the corresponding vector potential,
       ! AI (this involves the comp. of some nasty derivatives)::


       dKappaAdr = KappaA*aTube**2/&
            (4.0*Rperp*Rtube + aTube**2)
       dKappadx  = 2.0*xRel*Kappa/RPlus2
       dKappadr  = Kappa*(Rtube**2 - R2)/RPlus2

       ! Derivative of AI with respect to `x` and `rperp`:: 

       dAIdx   = dAkdk*dKappadx
       dAIdr   = dAkdk*dKappadr
       ! Obtain the BI field in the whole space from the corresponding
       ! vector potential, AI -->
       ! BI = curl(AI*ThetaUV) = BFRope_D(x_:z_)::

       BFRope_D =  B0*(4.0/(Kappa**3*cPi))*(Rtube/sqrt(RPlus2))**3*&
            (dAIdx*XyzRel_D + (dAIdr+d2Akdk2A*dKappaAdr &
            *(Kappa - KappaA)+AI)*UnitX_D)
       ! Compute the toroidal field (BIphix, BIphiy, BIphiz)
       ! produced by the azimuthal current Iphi. This is needed to ensure
       ! that the flux rope configuration is force free. 
       BIPhi_D = abs(Itube)/(2.0*cPi*RPerp*aTube**2) &
            *sqrt(2.0*(aTube**2 - RMinus**2))*&
            cross_product(UnitX_D,XyzRel_D)
    end if
    ! Add the field of the azimuthal current, Iphi::
    ! Compute the field produced by the ring current, Itube, both
    ! inside and outside the torus, BI = BFRope_D(x_:z_)::
    BFRope_D = BFRope_D + BIPhi_D

    !
    ! Add the field of two magnetic charges
    !
    if (DoBqField) then
       call compute_TD99_BqField(Xyz_D, B1qField_D)
       BFRope_D = BFRope_D + B1qField_D
    endif
   
    BFRope_D = BFRope_D*No2Si_V(UnitB_)
    RhoFRope = RhoFRope*No2Si_V(UnitRho_)

  end subroutine get_TD99_fluxrope

  !============================================================================

  subroutine init_TD99_parameters(Time)

    use ModCoordTransform, ONLY: rot_matrix_x,rot_matrix_y,rot_matrix_z
    real, optional, intent(in) :: Time
    real:: AlphaRope,LInduct,WFRope,FootSepar,ItubeDim
    ! Declare the rotational matrix of coordinate transformation::
    real :: Rotate_DD(3,3)
    !--------------------------------------------------------------------------
    !Save start time if used
    if(present(Time))StartTime = Time
    ! Compute the magnetic energy, WFRope, associated with the portion
    ! of the flux rope current that is above the solar surface::

    InvH0 = cGravitation*Msun/Rsun*Si2No_V(UnitU_)**2   ! in [-]
    AlphaRope  = 2.0*acos(Depth/Rtube)                 ! in [rad]
    FootSepar  = Rtube*sin(0.5*AlphaRope)/1.0e6         ! in [Mm]
    LInduct    = cMu*(0.5*AlphaRope/cPi)*Rtube*(log(8.0 &
         *(Rtube - Depth)/aTube) - 1.75)            ! in [H]
    WFRope     = 0.5*LInduct*Itube**2*1.0e7             ! in [ergs]

    ! Compute the average density inside the flux rope assuming that the
    ! total amount of prominence mass is Mass (=10^16g=10^13kg)::

    Rho0  = Mass/(AlphaRope*Rtube*cPi*aTube**2)
    ! in [kg/m^3]

    ! Define the normalized model parameters here::

    ! Flux rope::
    Rtube = Rtube*Si2No_V(UnitX_)
    aTube = aTube*Si2No_V(UnitX_)
    ItubeDim   = Itube
    Itube = ItubeDim*Si2No_V(UnitJ_)*Si2No_V(UnitX_)**2 ! in [A]
    Rho0  = Rho0*Si2No_V(UnitRho_)

    ! Save the maximum value of the current for possible use in
    ! varied_current case::

    ItubeSaved = Itube 

    ! Strapping field::

    Depth     = Depth*Si2No_V(UnitX_)
    L     = L*Si2No_V(UnitX_)
    q     = q*Si2No_V(UnitB_)*Si2No_V(UnitX_)**2

    ! Construct the rotational matrix, Rotate_DD, to position the
    ! flux rope in the desired way on the solar surface::

    Rotate_DD = matmul(rot_matrix_y(-0.5*cPi),&
         rot_matrix_x(-OrientationCme*cDegToRad))
    Rotate_DD = matmul(Rotate_DD,          &
         rot_matrix_y(LatitudeCme*cDegToRad))
    Rotate_DD = matmul(Rotate_DD,          &
         rot_matrix_z(-LongitudeCme*cDegToRad))
    UnitX_D = matmul([1.0, 0.0, 0.0], Rotate_DD)
    B0 = 0.5*Itube/Rtube
    XyzCenter_D = matmul([0.0, 0.0, 1 - Depth], Rotate_DD)
    if (iProc==0) then
       write(*,*) prefix
       write(*,*) prefix,'>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       write(*,*) prefix
       write(*,*) prefix,'    Twisted Flux Rope Model by Titov & Demoulin, 1999.     '
       write(*,*) prefix
       write(*,*) prefix,'>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       write(*,*) prefix
       write(*,*) prefix,'Depth      = ',Depth*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,*) prefix,'Rtube  = ', &
            Rtube*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,*) prefix,'aTube  = ', &
            aTube*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,*) prefix,'atube/Rtube = ',aTube/Rtube,'[-]'
       write(*,*) prefix,'Itube  = ',ItubeDim,'[A]'
       write(*,*) prefix,'Mass   = ',Mass*1.0e3,'[g] '
       write(*,*) prefix,'Rho0   = ',Rho0*No2Io_V(UnitRho_),'[g/cm^3]'
       write(*,*) prefix
       write(*,*) prefix,'q      = ', &
            q*No2Si_V(UnitB_)*No2Si_V(UnitX_)**2,'[T m^2]'
       write(*,*) prefix,'L      = ',L*No2Si_V(UnitX_)/1.0e6,'[Mm]'
       write(*,*) prefix
       write(*,*) prefix,'Free energy of flux rope is ',WFRope,'Ergs.'
       write(*,*) prefix,'Separation of flux rope ends is ',FootSepar,'Mm,'
       write(*,*) prefix,'   or ',cPi*FootSepar*1.0e6/(2.0*Rsun)*cRadToDeg,'deg.'
       write(*,*) prefix
       if (present(Time)) then
          write(*,*) prefix,'>>>>> Time-dependent strapping field is used <<<<<'
          write(*,*) prefix,'StartTime = ',StartTime,'[s]'
          write(*,*) prefix
       endif
    endif
    if (DoEquilItube) then

       ! Compute the equilibrium toroidal current, Itube, based
       ! on the force balance in direction normal to the surface of
       ! the flux tube.

       Itube = 8.0*cPi*q*L*Rtube &
            *(L**2+Rtube**2)**(-1.5) &
            /(alog(8.0*Rtube/aTube) &
            -1.5 + Li/2.0)                           ! in [-]
       WFRope    = 0.5*LInduct*(ItubeDim)**2*1.0e7      ! in [ergs]
    endif
    if (DoEquilItube.and.iProc==0) then
       write(*,*) prefix,'The strapping field, Bq, is added and the EQUILIBRIUM value'
       write(*,*) prefix,'of Itube is computed!!!'
       write(*,*) prefix
       write(*,*) prefix,'The value of Itube is reset to :: ',Itube
       write(*,*) prefix,'The free energy of the flux rope is :: ',WFRope,'Ergs.'
       write(*,*) prefix
    endif

  end subroutine init_TD99_parameters

  !=====================================================================!

  subroutine compute_TD99_BqField(Xyz_D, BqField_D, TimeNow)

    real, intent(in)  :: Xyz_D(3)
    real, intent(out) :: BqField_D(3)
    real, intent(in), optional :: TimeNow

    ! Variables related to coordinates::
    real:: R2Plus,R2Mins
    real, dimension(3):: RPlus_D,RMins_D
    !--------------------------------------------------------------------

    ! Compute the locations, RMins_D and RPlus_D, of the two magnetic
    ! charges, -/+q::

    if (present(TimeNow)) then
       RPlus_D(x_) = Xyz_D(x_)-L &
            - VTransX*(TimeNow - StartTime)*Si2No_V(UnitX_)
       RMins_D(x_) = Xyz_D(x_)+L &
            + VTransX*(TimeNow - StartTime)*Si2No_V(UnitX_)
       RPlus_D(y_) = Xyz_D(y_) &
            - VTransY*(TimeNow - StartTime)*Si2No_V(UnitX_)
       RMins_D(y_) = Xyz_D(y_) &
            + VTransY*(TimeNow - StartTime)*Si2No_V(UnitX_)
    else
       RPlus_D(x_) = Xyz_D(x_)-L
       RMins_D(x_) = Xyz_D(x_)+L
       RPlus_D(y_) = Xyz_D(y_)
       RMins_D(y_) = RPlus_D(y_)
    endif
    RPlus_D(z_) = Xyz_D(z_) + Depth -1.0
    RMins_D(z_) = RPlus_D(z_)
    R2Plus = norm2(RPlus_D)
    R2Mins = norm2(RMins_D)

    ! Compute the field of the strapping magnetic field, BqField_D::

    BqField_D = q*(RPlus_D/R2Plus**3 - RMins_D/R2Mins**3)
  end subroutine compute_TD99_BqField

  !=====================================================================!
end module EEE_ModTD99
