!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module EEE_ModTD99

#ifdef OPENACC
  use ModUtilities, ONLY: norm2 
#endif
  use EEE_ModCommonVariables
  use ModHyperGeometric
  implicit none
  save
  private

  public :: set_parameters_TD99, get_TD99_fluxrope, compute_TD99_BqField
  !
  !
  ! Variables related to the position of the flux rope:
  !
  ! Major radius (in fact - R_\infty, the radius of the circumference
  ! at which the toroidal coordinate, u, tends to infinity
  real :: RTube  = 0.0
  ! Minor radius:
  real :: aTube = 0.0
  ! Negative height (depth) of the current tube center with respect to
  ! the photosphere level.
  real :: Depth = 0.0
  ! direction of the magnetic moment of the configuration
  real :: UnitX_D(3)
  ! Magnetic configuration center
  real :: XyzCenter_D(3)

  ! magnetic field of the current tube at the center of configuration
  real :: BcTube, BcTubeDim, ITube
  ! Magnetic field at the center of configuration is always parallel,
  ! while starpping field is antiparallel, to the x-axis. Phi-conponent
  ! of the toroidal current is positive. However, the sign of the toroidal
  ! field component may be both positive and negative. To account for
  ! this, we introduce the variable to store the helicity sign.
  !
  real :: HelicitySign = 1.0
  !
  ! To set the negative helicity, the input parameter, BcTube, should be
  ! negative. This choice affects only the sign of helicity, but not the
  ! direction of the nagnetic field at the center of configuration.

  !-------------------Parameters of plasma inside the rope------------
  logical :: UsePlasmaBeta = .false.
  ! If ON, the parameter Beta, having a meaning of a constant
  ! gaskinetic-to-magnetic pressure
  ! ratio, is set for the ejecta, as well as the ejecta temperature. In this
  ! case both density and pressure are set inside the flux rope. Otherwise, the
  ! density only is calculated in agreeement with the input parameter of mass.
  real    :: PlasmaBeta = 0.0
  real    :: EjectaTemperature, EjectaTemperatureDim = 5.0e4 ! 50,000 K
  ! If the ejecta is cold ( PlasmaBeta = 0.0), the total mass of ejecta is
  ! set as the initial condition.
  real :: MassSi
  real :: Rho0=0.0
  !
  ! There is an alternative way to parameterize the current in the loop:
  ! namely, via the strength of the overarching (strapping) field. With
  ! this choice, the current is determined by the condition of equilibrium
  !`with the given strapping field intensity. The possible options are:
  ! readbstrap (the starapping field intensity is provided in the parameter
  ! file, getbstrap (in this case the field is taken from the MHD solution).
  character(LEN=10) :: TypeBStrap = 'none'
  !
  ! The following logical is set to .true. if TypeBStrap = 'getbstrap'
  logical :: DoGetBStrap = .false.
  !
  ! If .true. the current of the configuration is found from the
  ! equilibrium condition at the given strapping field
  logical :: UseEquilibriumCurrent = .false.
  !
  !
  real    :: bStrapping = 0.0, bStrappingDim = 0.0
  ! MAGNETIC CHARGES
  !
  character(LEN=10) :: TypeCharge = 'none'

  ! magnetic charge and its distance from the current ring
  real :: q = 0.0,  qDistance = 0.0, bQStrappingDim, bQStrapFraction
  !
  ! coordinate vectors of the charge location
  !
  real    :: RPlus_D(3), RMinus_D(3)
  real    :: UChargeX = 0.0
  !
  ! If .true. in the initial distribution of the magnetic field the
  ! strapping field from two magnetic charges
  logical :: UseStaticCharge = .false.
  !
  ! If .true., the magnetic field from two magnetic charges is varied
  ! dynamically, via varied B0Field
  logical :: UseDynamicStrapping = .false.
contains
  !============================================================================

  subroutine init_TD99_parameters

    use ModCoordTransform, ONLY: rot_matrix_x,rot_matrix_y,rot_matrix_z

    real :: AlphaRope, LInduct, WFRope, FootSepar, ITubeSi, bQStrapping
    ! Declare the rotational matrix of coordinate transformation:
    real :: Rotate_DD(3,3)
    ! internal inductance releted to \mu_0 R_\infty
    real, parameter :: Li = 0.5

    ! Define the normalized model parameters here::

    ! Flux rope::
    character(len=*), parameter:: NameSub = 'init_TD99_parameters'
    !--------------------------------------------------------------------------
    Rtube = Rtube*Io2No_V(UnitX_)
    aTube = aTube*Io2No_V(UnitX_)
    Depth = Depth*Io2No_V(UnitX_)

    if (iProc==0) then
       write(*,'(a)') prefix
       write(*,'(a)') prefix//&
            '>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       write(*,'(a)') prefix
       if(UseTD14)then
          write(*,'(a)') prefix//&
               '    Twisted Flux Rope Model by Titov, 2014.     '
       else
          write(*,'(a)') prefix//&
               '    Twisted Flux Rope Model by Titov & Demoulin, 1999.     '
       end if
       write(*,'(a)') prefix
       write(*,'(a)') prefix//&
            '>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       write(*,'(a)') prefix
       write(*,'(a,es13.5,a)') prefix//'Depth  = ', &
            Depth*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,'(a,es13.5,a)') prefix//'Rtube  = ', &
            Rtube*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,'(a,es13.5,a)') prefix//'aTube  = ', &
            aTube*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,'(a,es13.5,a)') prefix//'atube/Rtube = ',aTube/Rtube,'[-]'
       write(*,'(a)') prefix
    end if

    ! Computation of the portion of the flux rope current that is above
    ! the solar surface:
    AlphaRope  = 2.0*acos(Depth/Rtube)                 ! in [rad]
    FootSepar  = Rtube*No2Si_V(UnitX_)*sin(0.5*AlphaRope)/1.0e6  ! in [Mm]
    if(iProc==0)then
       write(*,'(a,es13.5,a)') prefix//'Separation of flux rope ends is ',&
            FootSepar,'[Mm],'
       write(*,'(a,es13.5,a)') prefix//'   or ',AlphaRope*cRadToDeg,'[deg].'
       write(*,'(a)') prefix
    end if
    !
    !
    ! Construct the rotational matrix, Rotate_DD, to position the
    ! flux rope in the desired way on the solar surface::

    Rotate_DD = matmul(rot_matrix_y(-0.5*cPi),&
         rot_matrix_x(-OrientationCme*cDegToRad))
    Rotate_DD = matmul(Rotate_DD,          &
         rot_matrix_y(LatitudeCme*cDegToRad))
    Rotate_DD = matmul(Rotate_DD,          &
         rot_matrix_z(-LongitudeCme*cDegToRad))
    UnitX_D = matmul([1.0, 0.0, 0.0], Rotate_DD)
    XyzCenter_D = (1 - Depth)*DirCme_D
    if(.not.UseTD14)then
       if(.not.UsePlasmaBeta)then
          !
          ! Compute the average density inside the flux rope assuming that the
          ! total amount of prominence mass is Mass (=10^16g=10^13kg):
          Rho0  = MassSi/(AlphaRope*Rtube*cPi*aTube**2*No2Si_V(UnitX_)**3)
          !
          ! Convert density from [kg/m^3] to dimensionless units
          Rho0  = Rho0*Si2No_V(UnitRho_)
          if(iProc==0)then
             write(*,'(a,es13.5,a)') prefix//'Mass   = ',MassSi*1.0e3,'[g] '
             write(*,'(a,es13.5,a)') prefix//'Rho0   = ',&
                  Rho0*No2Io_V(UnitRho_),'[g/cm^3]'
             write(*,'(a)') prefix
          end if
       else
          ! Alternatively, if non-zero PlasmaBeta is used, the pressure
          ! is profiled proportionally to the total pressure, and the density
          ! is found from the plasma pressure and constant ejecta temperature:
          EjectaTemperature =  EjectaTemperatureDim*Io2No_V(UnitTemperature_)
          if(iProc==0)then
             write(*,'(a,es13.5,a)') prefix//'Beta   =',PlasmaBeta
             write(*,'(a,es13.5,a)') prefix//'Tejecta=',&
                  EjectaTemperatureDim,'[K]'
             write(*,'(a)') prefix
          end if
       end if
    end if

    if(UseEquilibriumCurrent)then
       if(DoGetBStrap)then
          !
          ! Calculate projection of the MHD field at the apex point
          ! anti-parallel to the X-axis
          bStrapping = - sum(BAmbientApexSi_D*UnitX_D)*Si2No_V(UnitB_)
          if(bStrapping < 0.0)call CON_stop(NameSub//&
               ': Strapping MHD field is negative')
       else
          !
          ! Normalize field read from PARAM.in file
          bStrapping = bStrappingDim*Io2No_V(UnitB_)
       end if
       if(.not.UseTD14)then
          ! Compute the equilibrium toroidal current, Itube, based
          ! on the force balance in direction normal to the surface of
          ! the flux tube.
          ! Itube = 8*cPi*q*qDistance*Rtube &
          !     *(qDistance**2+Rtube**2)**(-1.5) &
          !     /(alog(8.0*Rtube/aTube) -1.5 + Li/2.0)    =
          !     4*cPi*(2*q*qDistance/(qDistance**2 + RTube**2))*RTube/&
          !     /(alog(8.0*Rtube/aTube) -1.5 + Li/2.0)    =
          !     4*cPi*bStrapping*RTube/(alog(8.0*Rtube/aTube) -1.5 + Li/2.0)
          ! From here, we can calculate BcTube = ITube/(2*RTube)
          ! in [-]
          BcTube  = 2*cPi*bStrapping/ &
               (log(8.0*Rtube/aTube) - 1.5 + Li/2.0)            ! in [-]
          ! Invert the equation, BcTubeSi = 0.5 \Mu_0 * ITubeSi/ RTubeSi
          ITubeSi = 2*(BcTube*No2Si_V(UnitB_))*&
               (RTube*No2Si_V(UnitX_))/cMu ![A]
       else
          ! Use equation from TD14:
          ITube =  (4*cPi*RTube*bStrapping) / (log(8*RTube/aTube)-25./24.)
          ITubeSi = ITube*No2Si_V(UnitB_)*No2Si_V(UnitX_)/cMu ![A]
       end if
       if (iProc==0) then
          write(*,'(a,es13.5,a)') prefix//' The strapping field, bStrapping=',&
               bStrapping*No2Io_V(UnitB_),' [Gs] is introduced and'
          write(*,'(a)') prefix//' EQUILIBRIUM value of Itube is computed!!!'
          write(*,'(a)') prefix
          write(*,'(a,es13.5,a)') prefix//' The value of Itube is reset to: ',&
               ItubeSi,' [A]'
          write(*,'(a)') prefix
       endif
       if(UseStaticCharge.or.UseDynamicStrapping)&
            bQStrapping = bStrapping*bQStrapFraction
    else
       ! Normalize field read from parameter file
       BcTube= BcTubeDim*Io2No_V(UnitB_)
       ! Invert the equation, BcTube = 0.5  * ITube/ RTube
       ITube   = 2*BcTube*RTube
       ITubeSi = ITube*No2Si_V(UnitB_)*No2Si_V(UnitX_)/cMu ![A]
       if(iProc==0)then
          write(*,'(a,es13.5,a)') prefix//'Itube  = ',ITubeSi,'[A]'
          write(*,'(a)') prefix
       end if
       if(UseStaticCharge.or.UseDynamicStrapping)&
            bQStrapping = bQStrappingDim*Io2No_V(UnitB_)
    end if
    if(iProc==0)then
       write(*,'(a,f4.0)') prefix//'Sign of helicity is ', HelicitySign
       write(*,'(a)') prefix
    end if

    LInduct    = cMu*(0.5*AlphaRope/cPi)*Rtube*No2Si_V(UnitX_)*(log(8.0 &
         *(Rtube - Depth)/aTube) - 1.75)               ! in [H]
    WFRope  = 0.5*LInduct*ItubeSi**2*1.0e7             ! in [ergs]

    if(iProc==0)then
       write(*,'(a,es13.5,a)') prefix//'Free energy of flux rope is ',&
            WFRope,'[erg].'
       write(*,'(a)') prefix
    end if

    if(UseStaticCharge.or.UseDynamicStrapping)then
       qDistance = qDistance*Io2No_V(UnitX_)
       ! Strapping field of charges::
       ! Invert formula:
       ! bQStrapping = 2*q*qDistance &
       !     *(qDistance**2+Rtube**2)**(-1.5)
       q = bQStrapping**(qDistance**2 + Rtube**2)**1.50/(2*qDistance)
       if(iProc==0)then
          write(*,'(a)') prefix//&
               ' The field of point magnetic charges is introduced'
          write(*,'(a,es13.5,a)') prefix//'q      = ', &
               q*No2Si_V(UnitB_)*No2Si_V(UnitX_)**2,'[T m^2]'
          write(*,'(a,es13.5,a)') prefix//'L      = ',&
               qDistance*No2Si_V(UnitX_)/1.0e6,'[Mm]'
          write(*,'(a)') prefix
       end if
       RPlus_D  = XyzCenter_D + qDistance*UnitX_D
       RMinus_D = XyzCenter_D - qDistance*UnitX_D
       if(UseDynamicStrapping)then
          UChargeX = UChargeX*Io2Si_V(UnitU_) &   ! From km/s to m/s
               *Si2No_V(UnitX_)                        ! To R_sun per s
          if(iProc==0)then
             write(*,'(a)') prefix,&
                  '>>>>> Time-dependent strapping field is used <<<<<'
             write(*,'(a,es13.5,a)') prefix//&
                  'StartTime = ',StartTime,'[s]'
             write(*,'(a,es13.5,a)') prefix//&
                  'Velocity of converging charges = ',&
                  UChargeX,'[R_s/s]'
             write(*,*) prefix
          end if ! iProc==0
       endif  ! UseDynamicStrapping
    end if  ! UseCharge
  end subroutine init_TD99_parameters
  !============================================================================
  subroutine set_parameters_TD99(NameCommand)

    use ModReadParam, ONLY: read_var

    integer :: iError
    real    :: MassDim
    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'set_parameters_TD99'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#CME")
       ! Set defaults:
       UseEquilibriumCurrent = .false.
       DoGetBStrap = .false.
       UseStaticCharge = .false.
       UseDynamicStrapping = .false.
       TypeBStrap = 'none'
       TypeCharge = 'none'
       call read_var('BcTubeDim', BcTubeDim)
       !
       ! Set negative helicity, if the input is negative:
       !
       HelicitySign = sign(1.0, BcTubeDim)
       BcTubeDim = abs(BcTubeDim)
       call read_var('RadiusMajor', Rtube)
       call read_var('RadiusMinor', aTube)
       call read_var('Depth',       Depth)
       !
       ! Save coords of the configuration center...
       XyzCmeCenterSi_D = DirCme_D*(1 - Depth)
       !
       ! ...and those of apex:
       XyzCmeApexSi_D = XyzCmeCenterSi_D + DirCme_D*RTube
       DoNormalizeXyz = .true.
       if(.not.UseTD14)then
          call read_var('UsePlasmaBeta', UsePlasmaBeta)
          if(UsePlasmaBeta)then
             call read_var('PlasmaBeta', PlasmaBeta)
             call read_var('EjectaTemperature', EjectaTemperatureDim)
          else
             PlasmaBeta = 0.0
             call read_var('MassDim',     MassDim)
             MassSi = MassDim*1.0e-3  ! g to kg
          end if
          ! If UseTD14 all these parameters are not effective
       end if
       call read_var('TypeBStrap', TypeBStrap, iError)
       if(iError/=0)then
          ! The line is empty, return
          TypeBStrap = 'none'
          RETURN
       elseif(len_trim(TypeBStrap)==0)then
          ! The line is empty, return
          TypeBStrap = 'none'
          RETURN
       end if
       select case(trim(TypeBStrap))
       case('readbstrap')
          UseEquilibriumCurrent  = .true.
          call read_var('bStrappingDim', bStrappingDim)
          call read_var('TypeCharge', TypeCharge, iError)
          if(iError/=0)then
             ! The line is empty, return
             TypeCharge = 'none'
             RETURN
          elseif(len_trim(TypeCharge)==0)then
             ! The line is empty, return
             TypeCharge = 'none'
             RETURN
          end if
       case('getbstrap')
          UseEquilibriumCurrent  = .true.
          DoGetBStrap = .true.
          call read_var('TypeCharge', TypeCharge, iError)
          if(iError/=0)then
             ! The line is empty, return
             TypeCharge = 'none'
             RETURN
          elseif(len_trim(TypeCharge)==0)then
             ! The line is empty, return
             TypeCharge = 'none'
             RETURN
          end if
       case('none')
          ! No TypeCharge, no TypeBStrap. Do nothing
          RETURN
       case('steady','moving','cancelflux')
          ! In this line, there is TypeCharge, not TypeBStrap
          UseEquilibriumCurrent  = .false.
          TypeCharge = trim(TypeBStrap)
          if(iProc==0)write(*,'(a)')TypeCharge//'                TypeCharge'
          TypeBStrap = 'none'
       case default
          if(iProc==0)call CON_stop(NameSub// ': '//&
              trim(TypeBStrap)//' is unknown as TypeCharge or TypeBStrap')
       end select
       ! Now, TypeCharge is read, either after TypeBStrap or instead of it'
       if(UseEquilibriumCurrent)then
          ! In this case we read, which fraction of the above
          ! defined strapping field is due to magnetic charges
          call read_var('bQStrapFraction', bQStrapFraction)
       else
          ! The magnetude of magnetic charges is characterized in terms
          ! of the strapping field they produce at the apex of flux rope
          call read_var('bQStrappingDim', bQStrappingDim)
       end if
       call read_var('qDistance',  qDistance)

       select case(trim(TypeCharge))
       case('steady')
          UseStaticCharge = .true.
       case('cancelflux')
          UseDynamicStrapping = .true.
          StartTime = -1.0
          call read_var('UChargeX', UChargeX)
       case('moving')
          UseStaticCharge = .true.
          UseDynamicStrapping = .true.
          StartTime = -1.0
          call read_var('UChargeX', UChargeX)
       case default
          if(iProc==0)call CON_stop(&
               NameSub//': TypeCharge='//trim(TypeCharge)//' is unknown')
       end select

    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_TD99
  !============================================================================
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
    real:: Kappa, Kappa2
    real:: B1qField_D(3)

    ! Initialize the TD99 model parameters once::

    !--------------------------------------------------------------------------
    if (DoInit) then
       call init_TD99_parameters
       DoInit=.false.
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
    if(.not.UseTD14)then
       call td99
    else
       call tdm
    end if
    !
    ! Add the field of two magnetic charges
    !
    if (UseStaticCharge) then
       call compute_TD99_BqField(Xyz_D, B1qField_D)
       BFRope_D = BFRope_D + B1qField_D
    endif

    BFRope_D  = BFRope_D *No2Si_V(UnitB_)
    RhoFRope  = RhoFRope *No2Si_V(UnitRho_)
    pFluxRope = pFluxRope*No2Si_V(UnitP_)
  contains
    !==========================================================================
    subroutine td99
      real:: DKappaDx, DKappaDr
      real:: KappaA, KappaA2, DKappaAdr

      ! Vector potential related variables::
      real:: Ak, dAkdk, D2Akdk2A
      real:: AI, dAIdx, dAIdr
      ! Flux-rope related variables::
      real:: BIPhi_D(3)

      !------------------------------------------------------------------------
      if (RMinus >= aTube) then
         !
         ! No pressure and density outside flux rope
         RhoFRope = 0.0; pFluxRope = 0.0
         ! Compute the field outside the current torus
         !

         ! Compute the vector potential, Ak, of the magnetic field
         ! produced by the ring current Itube and its derivatives

         ! Calculate \tilde{P}^{-1}_{-1/2}(\cosh u)
         Ak     = toroid_P(0, Kappa2In=Kappa2)
         ! Calculate d\tilde{P}^{-1}_{-1/2}(\cosh u)/du
         dAkDk  = toroid_dpdu(0, Kappa2In=Kappa2)/(1.0 - Kappa2)

         BFRope_D = BcTube*(Rtube/sqrt(RPlus2))**3*&
              (dAkDk*(2*xRel*XyzRel_D + (RTube**2 - R2Rel)*UnitX_D)/RPlus2 &
              + Ak*UnitX_D)
         ! No toroidal field outside the filament
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
         BIPhi_D = HelicitySign*BcTube*RTube/(cPi*RPerp*aTube**2) &
              *sqrt(2.0*(aTube**2 - RMinus**2))*&
              cross_product(UnitX_D,XyzRel_D)
         ! Add the prominence material inside the flux rope, assuming that the
         ! given total amount of mass

         if (.not.UsePlasmaBeta)then
            ! Cold plasma density is applied with density estimated from
            ! the total mass of eject
            RhoFRope = Rho0*exp(-10.0*(RMinus/aTube)**6)
            pFluxRope = 0
         else
            !
            ! Rescale BIPhi, which is not only a part of total pressure:
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
    end subroutine td99
    !==========================================================================
    subroutine tdm
       use ModConst, ONLY: cTwoPi, cPi
       use ModHyperGeometric, ONLY: calc_elliptic_int_1kind, &
         calc_elliptic_int_2kind

       real, parameter :: Delta = 0.1
       real :: bTheta_D(3),  bI_D(3)
       real :: RhoStarF, RhoStarG
       real :: RPerpHat_D(3), ThetaHat_D(3)
       real :: AxialFlux, Xi, BigTheta
       real :: DKappaDx, DKappaDRperp, kStarF, dKSFdX, dKSFdR, kStarG, dKSGdX, dKSGdR
       real :: EllipticKf, EllipticEf, FuncAf, dFuncAf, d2FuncAf, d3FuncAf
       real :: EllipticKg, EllipticEg
       real :: FuncAg, dFuncAg, d2FuncAg, d3FuncAg
       real :: SewingH, dSewingH, SewingF, dSewingF, SewingG, dSewingG
       real :: Ai, dAIdX, dAIdR
       real :: TmpG, dGdX, dGdR, TmpH, dHdR, Afx, dAFRdX, dAFXdR
      !------------------------------------------------------------------------

       ! In spherical case the straping field magnitude bVert should be provided

       !
       ! Coordinate unit vectors
       !
       RPerpHat_D = (XyzRel_D - xRel*UnitX_D)/RPerp
       ThetaHat_D = cross_product(UnitX_D, RPerpHat_D)
       !
       ! Toroidal coordinate, usual argument of special functions
       ! describing solution in the totoidal coordinates
       !
       ! Derivatives over x and over RPerp
       DKappaDx = - (xRel*Kappa**3) / (4*RPerp*RTube)
       DKappaDRperp = Kappa**3/(8*RPerp**2*RTube) * (RMinus**2-2*RPerp*(RPerp-RTube))

       AxialFlux = 3./(5*sqrt(2.)) * ITube * aTube

       ! Sewing functions
       Xi = (RMinus - aTube)/(Delta*aTube)
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
       RhoStarF = aTube*(1 + Delta*SewingF)
       kStarF = sqrt((RPerp*RTube)/(RPerp*RTube + RhoStarF**2/4.))
       dKSFdX = - (xRel*kStarF**3) / (4*RPerp*RTube)*dSewingF*RhoStarF/RMinus
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
       RhoStarG = aTube*(1 + Delta*SewingG)
       kStarG = sqrt((RPerp*RTube)/(RPerp*RTube+RhoStarG**2/4.))
       dKSGdX = - (xRel*kStarG**3) / (4*RPerp*RTube) * dSewingG*RhoStarG/RMinus
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
       Ai = ITube/cTwoPi*sqrt(RTube/RPerp)* &
            (FuncAf + dFuncAf*(Kappa - kStarF)+0.5*d2FuncAf*(Kappa - kStarF)**2)
       !
       ! Its partial derivatives
       !
       dAIdX = ITube/cTwoPi*sqrt(RTube/RPerp)*&
            (dFuncAf*DKappaDx + d2FuncAf*DKappaDx*(Kappa - kStarF) &
            + 0.5*d3FuncAf*dKSFdX*(Kappa - kStarF)**2)
       dAIdR = ITube/cTwoPi*sqrt(RTube/RPerp) &
            *(dFuncAf*DKappaDRperp + d2FuncAf*DKappaDRperp*(Kappa-kStarF) &
            + 0.5*d3FuncAf*dKSFdR*(Kappa-kStarF)**2) - Ai/(2*RPerp)
       !
       ! Poloidal magnetic field
       bI_D = - dAIdX*RPerpHat_D + (dAIdR + Ai/RPerp)*UnitX_D

       ! ---- toroidal field B_theta ----

       ! just a temporary variable, same for tmpH below
       TmpG = 3 + 4*dFuncAf*(Kappa-kStarF)

       dGdX = 4*(d2FuncAf*dKSFdX*(Kappa-kStarF)+dFuncAf*(DKappaDx-dKSFdX))
       dGdR = 4*(d2FuncAf*dKSFdR*(Kappa-kStarF)+dFuncAf*(DKappaDRperp-dKSFdR))

       TmpH = (Kappa**3*(xRel**2 + RTube**2 - RPerp**2) - &
            aTube**2*kStarG**3)*dFuncAg + &
            aTube**2*kStarG**3*d2FuncAg*(Kappa-kStarG)
       dHdR = (3*Kappa**2*DKappaDRperp*(xRel**2 + RTube**2 -RPerp**2) - &
            2*Kappa**3*RPerp-&
            3*aTube**2*kStarG**2*dKSGdR)*dFuncAg + &
            (Kappa**3*(xRel**2 + RTube**2 - RPerp**2) - &
            aTube**2*kStarG**3)*d2FuncAg*dKSGdR + &
            aTube**2*( (3*kStarG**2*dKSGdR*(Kappa - kStarG) &
            + kStarG**3*(DKappaDRperp - dKSGdR))*d2FuncAg &
            + kStarG**3*(Kappa - kStarG)*d3FuncAg*dKSGdR )

       Afx = AxialFlux/(4*cPi*RPerp)*sqrt(RTube/RPerp) * &
            ( FuncAg + (aTube**2*kStarG**3)/(4*RPerp*RTube)*dFuncAg &
            + TmpG**(5./2.)/(30*sqrt(3.)) &
            - 0.3 + TmpG**(3./2.)/(12*sqrt(3.)*RPerp*RTube)*TmpH )

       dAFRdX = AxialFlux/(24*sqrt(3.)*cPi*RPerp)/sqrt(RPerp*RTube) * &
            (1.5*sqrt(TmpG)*dGdX*xRel*Kappa**3*dFuncAf + &
            sqrt(TmpG)**3*(xRel*Kappa**3*d2FuncAf*dKSFdX &
            + (Kappa**3+3*xRel*Kappa**2*DKappaDx)*dFuncAf))
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
       BFRope_D = bI_D + bTheta_D
       RhoFRope  = 0.0
       pFluxRope = 0.0
    end subroutine tdm
    !==========================================================================
  end subroutine get_TD99_fluxrope
  !============================================================================
  subroutine compute_TD99_BqField(Xyz_D, BqField_D, TimeNow)

    real, intent(in)  :: Xyz_D(3)
    real, intent(out) :: BqField_D(3)
    real, intent(in), optional :: TimeNow

    ! Variables related to coordinates::
    real:: RPlus, RMinus
    !
    ! coordinate vectors of the charge location
    !
    real    :: RPlusSteady_D(3), RMinusSteady_D(3)
    real    :: RPlusMoving_D(3), RMinusMoving_D(3)
    !--------------------------------------------------------------------------
    ! if UseDynamicStrapping is .false. the cal may be only accidental
    if(present(TimeNow).and.(.not.UseDynamicStrapping))RETURN
    ! Compute the locations, RMins_D and RPlus_D, of the two magnetic
    ! charges, -/+q::
    RPlusSteady_D  = Xyz_D - RPlus_D
    RMinusSteady_D = Xyz_D - RMinus_D
    RPlus  = norm2(RPlusSteady_D)
    RMinus = norm2(RMinusSteady_D)
    !
    ! Compute the field of the strapping magnetic field, BqField_D::
    !
    BqField_D = q*(RPlusSteady_D/RPlus**3 - RMinusSteady_D/RMinus**3)
    if (.not.present(TimeNow)) RETURN
    ! In the steady-state location there are charges of the opposite
    ! sign, therefore the calculated field should be flipped
    BqField_D = - BqField_D
    ! When the time is long enough, the moving charges annihilate
    if((TimeNow - StartTime)*UChargeX >= qDistance)then
       BqField_D = BqField_D*No2Si_V(UnitB_)
       RETURN
    end if

    ! Compute the locations, RMins_D and RPlus_D, of the two magnetic
    ! charges, -/+q::
    RPlusMoving_D  = Xyz_D - RPlus_D  + (TimeNow - StartTime)*UChargeX*UnitX_D
    RMinusMoving_D = Xyz_D - RMinus_D - (TimeNow - StartTime)*UChargeX*UnitX_D
    RPlus  = norm2( RPlusMoving_D)
    RMinus = norm2(RMinusMoving_D)
    BqField_D = (BqField_D + &
         q*(RPlusMoving_D/RPlus**3 - RMinusMoving_D/RMinus**3))*No2Si_V(UnitB_)
  end subroutine compute_TD99_BqField
  !============================================================================
end module EEE_ModTD99
!==============================================================================
