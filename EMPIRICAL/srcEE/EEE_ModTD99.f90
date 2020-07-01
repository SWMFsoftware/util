!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module EEE_ModTD99
!  use ModUtilities, ONLY: norm2
  use EEE_ModCommonVariables
  use ModHyperGeometric
  implicit none
  save
  private

  public :: set_parameters_TD99, get_TD99_fluxrope, compute_TD99_BqField

  !
  ! Logical variables related to the magnetic field computation:
  logical, public :: UseBqField=.false.

  logical :: DoEquilItube=.false.
  real    :: StartTime
  
  logical :: UsePlasmaBeta = .false. 
  ! If ON, the parameter Beta, having a meaning of a constant 
  ! gaskinetic-to-magnetic pressure
  ! ratio, is set for the ejecta, as well as the ejecta temperature. In this 
  ! case both density and pressure are set inside the flux rope. Otherwise, the
  ! density only is calculated in agreeement with the input parameter of mass.
  real    :: PlasmaBeta = 0.0
  real    :: EjectaTemperature, EjectaTemperatureDim = 5.0e4 ! 50,000 K

  ! Variables related to the position of the flux rope:
  !
  ! internal inductance releted to \mu_0 R_\infty
  real, parameter :: Li=0.5

  ! Variables related to the flux rope properties::
  real :: ITubeSi = 0.0, ITube = 0.0
  !Major radius (in fact - R_\infty, the radius of the circumference
  !at which the toroidal coordinate, u, tends to infinity
  real :: RTube  = 0.0
  !Minor radius:
  real :: aTube = 0.0
  !Negative height (depth) of the current tube center with respect to
  !the photosphere level.
  real :: Depth = 0.0

  real :: MassSi, MassDim
  real :: Rho0=0.0

  ! Magnetic field at the center of configuration is alway positive,
  ! Starpping field is negative. Phi-conponent of the toroidal current
  ! is positive. However, the sign of the field toroidal component may be
  ! both positive and negative. To account for this, we introduce the
  ! variable to store the helicity sign.
  !
  real :: HelicitySign = 1.0
  !
  ! To set the negative helicity, the input parameter, BcTube, should negative
  ! This choice affects only the sign of helicity, but not the direction
  ! of the nagnetic field at the center of configuration.

  ! MAGNETIC CHARGES
  !
  character(LEN=10) :: TypeCharge = 'none'

  ! magnetic charge and its distance from the current ring
  real :: q = 0.0, BqFieldDim, qDistance = 0.0
  !
  ! coordinate vectors of the charge location
  !
  real :: qPlusLoc_D(3), qMinusLoc_D(3)
  real :: VTransXSi = 0.0   
  real :: VTransYSi = 0.0
  !
  ! The velocity vector, for a negative charge
  !
  real :: qVelMinusSi_D(3)
 

  ! direction of the magnetic moment of the configuration
  real :: UnitX_D(3)
  ! magnetic field of the current tube at the center of configuration
  real :: BcTube, BcTubeDim
  ! Magnetic configuration center
  real :: XyzCenter_D(3)

contains
  !============================================================================
  subroutine set_parameters_TD99(NameCommand)

    use ModReadParam, ONLY: read_var

    integer :: iError
    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'set_parameters_TD99'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#TD99FLUXROPE")
       call read_var('UseCme',              UseCme)
       call read_var('DoAddFluxRope'      , DoAddFluxRope)
       call read_var('DoEquilItube'       , DoEquilItube)
       call read_var('Itube'              , ItubeSi)
       call read_var('Rtube'              , Rtube)
       call read_var('aTube'              , aTube)
       call read_var('Depth'              , Depth)
       call read_var('MassDim'            , MassDim)
       MassSi = MassDim*1.0e-3
       call read_var('LongitudeCme'       , LongitudeCme)
       call read_var('LatitudeCme'        , LatitudeCme)
       call read_var('OrientationCme'     , OrientationCme)
       call read_var('DoBqField'          , UseBqField)
       call read_var('q'             , q)
       call read_var('L'             , qDistance)
    case("#CME")
       call read_var('BcTubeDim', BcTubeDim)
       !
       ! Set negative helicity, if the input is negative:
       !
       HelicitySign = sign(1.0, BcTubeDim)
       BcTubeDim = abs(BcTubeDim)
       call read_var('RadiusMajor', Rtube)
       call read_var('RadiusMinor', aTube)
       call read_var('Depth',       Depth)
       call read_var('UsePlasmaBeta', UsePlasmaBeta)
       if(UsePlasmaBeta)then
          call read_var('PlasmaBeta', PlasmaBeta)
          call read_var('EjectaTemperature', EjectaTemperatureDim)
       else
          PlasmaBeta = 0.0
          call read_var('MassDim',     MassDim)
          MassSi = MassDim*1.0e-3  !g to kg
       end if
       call read_var('TypeCharge', TypeCharge, iError)
       if(iError/=0)then
          TypeCharge = 'none'; UseBqField = .false.
          RETURN
       elseif(len_trim(TypeCharge)==0)then
          TypeCharge = 'none'; UseBqField = .false.
          RETURN
       elseif(trim(TypeCharge)=='none')then
          UseBqField = .false.
          RETURN
       end if
       call read_var('BqFieldDim', BqFieldDim)
       call read_var('qDistance',  qDistance)
       if(trim(TypeCharge)=='steady')then
          UseBqField = .true.
          RETURN
       end if
       call CON_stop(NameSub//&
            ':moving charges and flux cancelation is work in progress') 
    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_TD99
  !===========================================================================
  subroutine get_TD99_fluxrope(Xyz_D, BFRope_D, RhoFRope, pFluxRope)

    use ModCoordTransform, ONLY: cross_product

    !    Twisted Magnetic Field Configuration by Titov & Demoulin '99   !
    !                                                                   !
    ! An instability that causes a CME eruption is expected to occur at !
    ! R > L*sqrt(2). For a detailed description of the initial state    !
    ! refer to A&A, 1999, v.351, pp.707-720                             !
    !                                                                   !
    ! ___  This module was written by Ilia Roussev on June 10, 2002 ___ !

    real, intent(in)  :: Xyz_D(3)
    real, intent(out) :: BFRope_D(3)
    real, intent(out) :: RhoFRope
    real, intent(out) :: pFluxRope

    ! Coordinates relative to the configuration center:
    real :: XyzRel_D(3)

    ! Distance to the configuration center squared:
    real :: R2Rel !=sum(XyzRel_D**2)

    real:: xRel
    real:: RMinus, RPlus2, Rperp
    real:: Kappa, Kappa2, DKappaDx, DKappaDr
    real:: KappaA, KappaA2, DKappaAdr
 
    ! Vector potential related variables::
    real:: Ak, dAkdk, D2Akdk2A
    real:: AI, dAIdx, dAIdr
    ! Flux-rope related variables::
    real:: BIPhi_D(3), B1qField_D(3)
    
    logical:: DoFirstCall=.true.
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
    R2Rel = sum(XyzRel_D**2); xRel = sum(XyzRel_D*UnitX_D)

    ! Compute Rperp and TubalDist::

    Rperp = sqrt(R2Rel - xRel**2)
    RMinus = sqrt(xRel**2 + (Rperp - Rtube)**2)
    RPlus2 = (Rperp + Rtube)**2 + xRel**2


    ! Define the model input, Kappa

    Kappa2 = 4.0*Rperp*Rtube/RPlus2; Kappa = sqrt(Kappa2)
    if (RMinus >= aTube) then
       !
       ! No pressure and density outside flux rope
       RhoFRope = 0.0; pFluxRope = 0.0
       ! Compute the field outside the current torus   
       !


       ! Compute the vector potential, Ak, of the magnetic field 
       ! produced by the ring current Itube and its derivatives

       !Calculate \tilde{P}^{-1}_{-1/2}(\cosh u)
       Ak     = toroid_P(0, Kappa2In=Kappa2)
       !Calculate d\tilde{P}^{-1}_{-1/2}(\cosh u)/du
       dAkDk  = toroid_dpdu(0, Kappa2In=Kappa2)/(1.0 - Kappa2)
       
       BFRope_D = BcTube*(Rtube/sqrt(RPlus2))**3*&
            (dAkDk*(2*xRel*XyzRel_D + (RTube**2 - R2Rel)*UnitX_D)/RPlus2 &
            + Ak*UnitX_D)
       !No toroidal field outside the filament
       BIPhi_D = 0.0
    else
       !
       ! Compute the field and density inside the current torus
       !
       ! 1.   
       ! Define the model input, KappaA. A given point is characterized by
       ! two coordinates, say, RPepr and RMinus. In this case,
       ! if we denote 
       ! Kappa = function(RPerp,RMinus): i.e. Kappa**2 = 4*RPerp*RTube *R_+**2
       ! then 
       ! KappaA = function(RPepr,ATube), or
       KappaA2 = 4.0*Rperp*Rtube/(4.0*Rperp*Rtube + aTube**2)
       KappaA  = sqrt(KappaA2)
       ! 
       ! 2.
       ! The vector potential of the externap field, Ak and its derivative,
       ! dAk/dk, are both prolonged to the tube interior using the following 
       ! sewing function (below we use A_k by a factor of Kappa**3 different
       ! from what we used above:
       ! A_i = A_k(KappaA) + dA_k/dKappaA*(Kappa - Kappa_A) (*)
       !
       Ak      = toroid_P(0, Kappa2In=KappaA2)*KappaA**3
       dAkdk   = toroid_dpdu(0, Kappa2In=KappaA2)*KappaA2/(1 - KappaA2) 
       AI = Ak + dAkdk*(Kappa - KappaA)
       !
       ! 3. 
       ! Now, we derive the field components in terms of x and rPerp 
       ! derivatives. Function A_i depends on x only via Kappa, so that:
       ! dA_i/dx = dA_k/dKappaA*dKappa/dx:
       !
       dKappadx  = 2.0*xRel*Kappa/RPlus2
       dAIdx   = dAkdk*dKappadx      
       !
       ! 4.
       ! Analogously we account for the dependence of Kappa on the radial
       ! coordinate:
       ! dA_i/dr = dA_k/dKappaA*dKappaA/dr:
       !
       dKappadr  = Kappa*(Rtube**2 - R2Rel)/RPlus2
       dAIdr   = dAkdk*dKappadr
       !
       ! 5.
       ! Now, we account for the dependence of KappaA on rPerp. From (*), the 
       ! contributions from the first derivative dA_k(KappaA)/dKappaA cancel 
       ! each other, so that only the second derivative matters, which is 
       ! equal to  d^2A_k(KappaA)/dKappaA^2:
       !
       d2Akdk2A = KappaA/(1 - KappaA**2)*(3*toroid_P(0, Kappa2In=KappaA**2) +&
           toroid_dpdu(0, Kappa2In=KappaA**2)*(1 + KappaA**2)/(1 - KappaA**2))
       dKappaAdr = KappaA*aTube**2/&
            (4.0*Rperp*Rtube + aTube**2)
       dAIdr = dAIdr + d2Akdk2A*dKappaAdr*(Kappa - KappaA)

       BFRope_D =  BcTube/(Kappa**3)*(Rtube/sqrt(RPlus2))**3*&
            (dAIdx*XyzRel_D + (dAIdr + AI)*UnitX_D)
       ! Compute the toroidal field (BIphix, BIphiy, BIphiz)
       ! produced by the azimuthal current Iphi. This is needed to ensure
       ! that the flux rope configuration is force free. 
       BIPhi_D = HelicitySign*abs(Itube)/(2.0*cPi*RPerp*aTube**2) &
            *sqrt(2.0*(aTube**2 - RMinus**2))*&
            cross_product(UnitX_D,XyzRel_D)
       ! Add the prominence material inside the flux rope, assuming that the
       ! given total amount of mass
       
       if (.not.UsePlasmaBeta)then
          !Cold plasma density is applied with density estimated from
          !the total mass of eject
          RhoFRope = Rho0*exp(-10.0*(RMinus/aTube)**6)
          pFluxRope = 0
       else
          !
          !Rescale BIPhi, which is not only a part of total pressure:
          !
          BIPhi_D = BIPhi_D/sqrt(1 + PlasmaBeta)
          pFluxRope = 0.50*sum(BIPhi_D**2)*PlasmaBeta
          RhoFRope = pFluxRope/EjectaTemperature
       end if
    end if
    ! Add the field of the azimuthal current, Iphi::
    ! Compute the field produced by the ring current, Itube, both
    ! inside and outside the torus, BI = BFRope_D(x_:z_)::
    BFRope_D = BFRope_D + BIPhi_D
    
    !
    ! Add the field of two magnetic charges
    !
    if (UseBqField) then
       call compute_TD99_BqField(Xyz_D, B1qField_D)
       BFRope_D = BFRope_D + B1qField_D
    endif

   
    BFRope_D  = BFRope_D *No2Si_V(UnitB_)
    RhoFRope  = RhoFRope *No2Si_V(UnitRho_)
    pFluxRope = pFluxRope*No2Si_V(UnitP_)
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

    AlphaRope  = 2.0*acos(Depth/Rtube)                 ! in [rad]
    FootSepar  = Rtube*No2Si_V(UnitX_)*sin(0.5*AlphaRope)/1.0e6  ! in [Mm]
    LInduct    = cMu*(0.5*AlphaRope/cPi)*Rtube*No2Si_V(UnitX_)*(log(8.0 &
         *(Rtube - Depth)/aTube) - 1.75)            ! in [H]
    ! Compute the average density inside the flux rope assuming that the
    ! total amount of prominence mass is Mass (=10^16g=10^13kg)::

    Rho0  = MassSi/(AlphaRope*Rtube*cPi*aTube**2*No2Si_V(UnitX_)**3)
    ! in [kg/m^3]
    Rho0  = Rho0*Si2No_V(UnitRho_)
    
    !Alternatively, if non-zero PlasmaBeta is used, the pressure
    !is profiled proportionally to the total pressure, and the density 
    !is found from the plasma pressure and constant ejecta temperature:
    EjectaTemperature =  EjectaTemperatureDim*Io2No_V(UnitTemperature_)


    ! Define the normalized model parameters here::

    ! Flux rope::
    Rtube = Rtube*Io2No_V(UnitX_)
    aTube = aTube*Io2No_V(UnitX_)
    BcTube= BcTubeDim*Io2No_V(UnitB_)
    
    !Invert the equation, BcTubeSi = 0.5 \Mu_0 * ITubeSi/ RTubeSi
    ITubeSi = 2*(BcTube*No2Si_V(UnitB_)) * (RTube*No2Si_V(UnitX_))/cMu ![A]
    ITube   = ITubeSi*Si2No_V(UnitJ_)*Si2No_V(UnitX_)**2
    WFRope     = 0.5*LInduct*ItubeSi**2*1.0e7             ! in [ergs]

    ! Strapping field::

    Depth     = Depth*Io2No_V(UnitX_)
    qDistance = qDistance*Io2No_V(UnitX_)
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
    XyzCenter_D = matmul([0.0, 0.0, 1 - Depth], Rotate_DD)
    if (iProc==0) then
       write(*,'(a)') prefix
       write(*,'(a)') prefix//'>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       write(*,'(a)') prefix
       write(*,'(a)') prefix//'    Twisted Flux Rope Model by Titov & Demoulin, 1999.     '
       write(*,'(a)') prefix
       write(*,'(a)') prefix//'>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       write(*,'(a)') prefix
       write(*,'(a,es13.5,a)') prefix//'Depth  = ',Depth*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,'(a,es13.5,a)') prefix//'Rtube  = ', &
            Rtube*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,'(a,es13.5,a)') prefix//'aTube  = ', &
            aTube*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,'(a,es13.5,a)') prefix//'atube/Rtube = ',aTube/Rtube,'[-]'
       write(*,'(a,es13.5,a)') prefix//'Itube  = ',ITubeSi,'[A]'
       if(UsePlasmaBeta)then
          write(*,'(a,es13.5,a)') prefix//'Beta   =',PlasmaBeta
          write(*,'(a,es13.5,a)') prefix//'Tejecta=',EjectaTemperatureDim,'[K]'
       else
          write(*,'(a,es13.5,a)') prefix//'Mass   = ',MassDim,'[g] '
          write(*,'(a,es13.5,a)') prefix//'Rho0   = ',Rho0*No2Io_V(UnitRho_),'[g/cm^3]'
       end if
       write(*,'(a)') prefix
       write(*,'(a,es13.5,a)') prefix//'q      = ', &
            q*No2Si_V(UnitB_)*No2Si_V(UnitX_)**2,'[T m^2]'
       write(*,'(a,es13.5,a)') prefix//'L      = ',qDistance*No2Si_V(UnitX_)/1.0e6,'[Mm]'
       write(*,'(a)') prefix
       write(*,'(a,es13.5,a)') prefix//'Free energy of flux rope is ',WFRope,'[erg].'
       write(*,'(a,es13.5,a)') prefix//'Separation of flux rope ends is ',FootSepar,'[Mm],'
       write(*,'(a,es13.5,a)') prefix//'   or ',cPi*FootSepar*1.0e6/(2.0*Rsun)*cRadToDeg,'[deg].'
       write(*,'(a)') prefix
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

       ItubeSi = 8.0*cPi*q*qDistance*Rtube &
            *(qDistance**2+Rtube**2)**(-1.5) &
            /(alog(8.0*Rtube/aTube) &
            -1.5 + Li/2.0)                           ! in [-]
       WFRope    = 0.5*LInduct*(ItubeDim)**2*1.0e7      ! in [ergs]
    endif
    if (DoEquilItube.and.iProc==0) then
       write(*,'(a)') prefix//'The strapping field, Bq, is added and the EQUILIBRIUM value'
       write(*,'(a)') prefix//'of Itube is computed!!!'
       write(*,'(a)') prefix
       write(*,'(a,es13.5)') prefix//'The value of Itube is reset to :: ',ItubeSi
       write(*,'(a,es13.5,a)') prefix,'The free energy of the flux rope is :: ',WFRope,'Ergs.'
       write(*,'(a)') prefix
    endif

  end subroutine init_TD99_parameters
  !===========================================================================
  subroutine compute_TD99_BqField(Xyz_D, BqField_D, TimeNow)

    real, intent(in)  :: Xyz_D(3)
    real, intent(out) :: BqField_D(3)
    real, intent(in), optional :: TimeNow

    ! Variables related to coordinates::
    real:: R2Plus, R2Mins
    real:: RPlus_D(3), RMins_D(3)
    logical :: DoSaveStartTime = .true.
    !-------------------------------------------------------------------------
    if(DoSaveStartTime.and.present(TimeNow))then
       StartTime = TimeNow
       DoSaveStartTime = .false.
    end if

    ! Compute the locations, RMins_D and RPlus_D, of the two magnetic
    ! charges, -/+q::

    if (present(TimeNow)) then
       RPlus_D(x_) = Xyz_D(x_) - qDistance &
            - VTransXSi*(TimeNow - StartTime)*Si2No_V(UnitX_)
       RMins_D(x_) = Xyz_D(x_) + qDistance &
            + VTransXSi*(TimeNow - StartTime)*Si2No_V(UnitX_)
       RPlus_D(y_) = Xyz_D(y_) &
            - VTransYSi*(TimeNow - StartTime)*Si2No_V(UnitX_)
       RMins_D(y_) = Xyz_D(y_) &
            + VTransYSi*(TimeNow - StartTime)*Si2No_V(UnitX_)
    else
       RPlus_D(x_) = Xyz_D(x_) - qDistance
       RMins_D(x_) = Xyz_D(x_) + qDistance
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

  !===========================================================================
end module EEE_ModTD99
