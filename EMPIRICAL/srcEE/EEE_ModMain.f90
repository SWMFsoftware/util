!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module EEE_ModMain

  use EEE_ModGetB0, ONLY: EEE_get_B0
  use ModUtilities, ONLY: CON_stop

  implicit none
  save

  private

  public :: EEE_initialize
  public :: EEE_set_parameters
  public :: EEE_get_state_init
  public :: EEE_get_state_bc
  public :: EEE_get_b0
  public :: EEE_do_not_add_cme_again
  public :: EEE_set_plot_range

contains
  !============================================================================
  subroutine EEE_initialize(BodyNDim, BodyTDim, Gamma, iCommIn, TimeNow)
    use EEE_ModCommonVariables
    use EEE_ModGL98, ONLY: gl98_init
#ifdef _OPENACC
    use ModUtilities, ONLY: norm2
#endif

    real, intent(in)             :: BodyNDim, BodyTDim, Gamma
    integer, optional, intent(in):: iCommIn
    real, optional, intent(in)   :: TimeNow

    integer :: iComm, iError
    !--------------------------------------------------------------------------
    if(present(iCommIn))then
       iComm = iCommIn
    else
       iComm = MPI_COMM_WORLD
    end if
    call MPI_COMM_RANK(iComm,iProc,iError)
    if(present(TimeNow) .and. tStartCme < 0.0)then
       tStartCme = TimeNow
       if(iProc==0)write(*,*) prefix,' tStartCme=',tStartCme,'[s]'
    end if

    g = Gamma
    inv_g = 1.0/g
    gm1 = g - 1.0
    inv_gm1 = 1.0/(g - 1.0)

    ! assume MassIon_I(1) = 1.0
    No2Si_V(UnitX_)   = rSun
    No2Si_V(UnitU_)   = sqrt(g*cBoltzmann*BodyTDim/cProtonMass)
    No2Si_V(UnitRho_) = 1000000*cProtonMass*BodyNDim

    ! Set other normalizing SI variables from the independent ones.
    !
    ! For sake of convenience
    !  units of B are chosen to satisfy v_A = B/sqrt(rho)       (mu = 1)
    !  units of n are chosen to satisfy  n  = rho/(ionmass/amu) (mp = 1)
    !  units of T are chosen to satisfy  T  = p/n               (kBoltzmann=1)
    !
    ! Note that No2Si_V(UnitN_) is NOT EQUAL TO 1/No2Si_V(UnitX_)^3 !!!
    No2Si_V(UnitT_)          = No2Si_V(UnitX_)/No2Si_V(UnitU_)         ! s
    No2Si_V(UnitN_)          = No2Si_V(UnitRho_)/cProtonMass           ! #/m^3
    No2Si_V(UnitP_)          = No2Si_V(UnitRho_)*No2Si_V(UnitU_)**2    ! Pa
    No2Si_V(UnitB_)          = No2Si_V(UnitU_) &
         *sqrt(cMu*No2Si_V(UnitRho_))                                  ! T
    No2Si_V(UnitRhoU_)       = No2Si_V(UnitRho_)*No2Si_V(UnitU_)     ! kg/m^2/s
    No2Si_V(UnitEnergyDens_) = No2Si_V(UnitP_)                         ! J/m^3
    No2Si_V(UnitPoynting_)   = No2Si_V(UnitEnergyDens_)*No2Si_V(UnitU_)! J/m^2/s
    No2Si_V(UnitJ_)          = No2Si_V(UnitB_)/( No2Si_V(UnitX_)*cMu ) ! A/m^2
    No2Si_V(UnitElectric_)   = No2Si_V(UnitU_)*No2Si_V(UnitB_)         ! V/m
    No2Si_V(UnitTemperature_)= No2Si_V(UnitP_) &
         /( No2Si_V(UnitN_)*cBoltzmann )                               ! K
    No2Si_V(UnitDivB_)       = No2Si_V(UnitB_)/No2Si_V(UnitX_)         ! T/m
    No2Si_V(UnitAngle_)      = 1.0                                     ! radian

    ! Set inverse conversion SI -> normalized
    Si2No_V = 1.0/No2Si_V

    ! As a default use SI units, so below only the differences need to be set
    Io2Si_V = 1.0
    No2Io_V = No2Si_V

    Io2Si_V(UnitX_)           = rSun                      ! R
    Io2Si_V(UnitRho_)         = 1.0E+3                    ! g/cm^3
    Io2Si_V(UnitN_)           = 1.0E+6                    ! #/cm^3
    Io2Si_V(UnitU_)           = 1.0E+3                    ! km/s
    Io2Si_V(UnitP_)           = 1.0E-1                    ! dyne/cm^2
    Io2Si_V(UnitB_)           = 1.0E-4                    ! Gauss
    Io2Si_V(UnitRhoU_)        = 1.0E+1                    ! g/cm^2/s
    Io2Si_V(UnitEnergydens_)  = 1.0E-1                    ! erg/cm^3
    Io2Si_V(UnitJ_)           = 1.0E-6                    ! uA/m^2
    Io2Si_V(UnitDivB_)        = 1.0E-2                    ! Gauss/cm
    Io2Si_V(UnitAngle_)       = cRadToDeg                 ! degrees

    ! Calculate the remaining unit conversions
    Si2Io_V = 1/Io2Si_V
    No2Io_V = No2Si_V*Si2Io_V
    Io2No_V = 1/No2Io_V

    Gbody  = -cGravitation*mSun*(Si2No_V(UnitU_)**2 * Si2No_V(UnitX_))
    if(DoNormalizeXyz)then
       XyzCmeCenterSi_D = XyzCmeCenterSi_D*Io2Si_V(UnitX_)
       XyzCmeApexSi_D   = XyzCmeApexSi_D  *Io2Si_V(UnitX_)
       rCmeApexInvSi    = 1/norm2(XyzCmeApexSi_D)
       ! To avoid repeated scaling in subsequent sessions
       DoNormalizeXyz = .false.
    end if

    DoInit = .true.

    if(UseGL .or. UseSpheromak) then 
      call gl98_init
    end if

   !$acc update device(UseCme, UseGL, UseTD, UseTD14, UseTD22)
   !$acc update device(UseSpheromak, DoAddTD, DoAddGL, DoAddSpheromak)

   !$acc update device(Io2Si_V, Si2Io_V, Io2No_V, No2Io_V, Si2No_V, No2Si_V)
   !$acc update device(rCmeApexInvSi, tStartCme, tDecayCmeDim, tDecayCme)
   !$acc update device(Gbody)
  end subroutine EEE_initialize
  !============================================================================
  subroutine EEE_set_parameters(NameCommand)

    use ModReadParam,     ONLY: read_var
    use EEE_ModGL98,      ONLY: set_parameters_GL98
    use EEE_ModTD99,      ONLY: set_parameters_TD99
    use EEE_ModArch,      ONLY: set_parameters_arch
    use EEE_ModShearFlow, ONLY: set_parameters_shearflow
    use EEE_ModCms,       ONLY: set_parameters_cms
    use EEE_ModCommonVariables, ONLY: &
         UseCme, DoAddFluxRope, UseTD, UseTD14, UseGL, UseShearFLow, UseArch, &
         UseTD22, LongitudeCme, LatitudeCme, OrientationCme, DirCme_D, &
         UseCms, UseSpheromak, DoAddTD, DoAddGL, DoAddSpheromak, &
         tStartCme, tDecayCmeDim
    use ModNumConst,      ONLY: cDegToRad
    use ModCoordTransform, ONLY: lonlat_to_xyz
    use CON_axes,          ONLY: dLongitudeHgrDeg

    character(len=*), intent(in) :: NameCommand

    character(len=20):: TypeCme

    character(len=*), parameter:: NameSub = 'EEE_set_parameters'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#CMETIME")
       call read_var("tStartCme", tStartCme)
    case("#CME")
       call read_var("UseCme", UseCme)
       if(UseCme)then
          call read_var("DoAddFluxRope",  DoAddFluxRope)
          ! Reset tStartCme to -1 (say for a second CME)
          if(DoAddFluxRope) tStartCme = -1.0
          call read_var("tDecayCme",      tDecayCmeDim)
          call read_var("LongitudeCme",   LongitudeCme)
          call read_var("LatitudeCme",    LatitudeCme)
          ! Correct for Rotation
          LongitudeCme = LongitudeCme - dLongitudeHgrDeg
          ! Get direction vector for the CME
          call lonlat_to_xyz(LongitudeCme*cDegToRad, &
               LatitudeCme*cDegToRad, DirCme_D)
          call read_var("OrientationCme", OrientationCme)
          call read_var("TypeCme", TypeCme, IsUpperCase=.true.)
          select case(TypeCme)
          case("TD")
             UseTD = .true.
             DoAddTD = DoAddFluxRope
             call set_parameters_TD99(NameCommand)
          case("TD14")
             UseTD = .true.
             UseTD14 = .true.
             DoAddTD = DoAddFluxRope
             call set_parameters_TD99(NameCommand)
          case("TD22")
             UseTD = .true.
             UseTD22 = .true.
             DoAddTD = DoAddFluxRope
             call set_parameters_TD99(NameCommand)
          case("GL")
             UseGL = .true.
             DoAddGL = DoAddFluxRope
             call set_parameters_GL98(NameCommand)
          case("SPHEROMAK")
             UseSpheromak = .true.
             DoAddSpheromak = DoAddFluxRope
             call set_parameters_GL98("#SPHEROMAK")
          case("BREAKOUT")
             UseShearFlow = .true.
             call set_parameters_shearflow(NameCommand)
             UseArch = .true.
             call set_parameters_arch(NameCommand)
          case default
             call CON_stop(NameSub//': invalid value for TypeCme='//TypeCme)
          end select
       else
          UseTD        = .false.
          UseTD14      = .false.
          UseGL        = .false.
          UseSpheromak = .false.
          UseShearFlow = .false.
          UseArch      = .false.
       end if

       ! The remaining commands are preserved for backwards compatibility
       ! and possibly for expert use (more options than #CME command)
    case("#ARCH")
       call set_parameters_arch(NameCommand)
    case("#SHEARFLOW")
       UseCme       = .true.
       UseShearFlow = .true.
       call set_parameters_shearflow(NameCommand)
    case("#CMS")
       UseCme = .true.
       UseCms = .true.
       call set_parameters_cms(NameCommand)
    end select

  end subroutine EEE_set_parameters
  !============================================================================
  subroutine EEE_get_state_bc(Xyz_D, Rho, U_D, B_D, p, Time, nStep, nIter)
    !$acc routine seq
    use EEE_ModCommonVariables, ONLY: UseCme, UseTD, UseShearFlow, UseGL, &
         UseCms, UseSpheromak, tStartCme, tDecayCmeDim
    use EEE_ModTD99, ONLY: get_TD99_fluxrope
    use EEE_ModShearFlow, ONLY: get_shearflow
    use EEE_ModGL98, ONLY: get_GL98_fluxrope
    use EEE_ModCms, ONLY: get_cms

    real, intent(in) :: Xyz_D(3), Time
    real, intent(out) :: Rho, U_D(3), B_D(3), p
    integer, intent(in):: nStep, nIter

    ! Perturbations due to CME
    real :: Rho1, U1_D(3), B1_D(3), p1

    ! Coefficient for perturbations (less than 1 for decay)
    real :: Coeff
    ! initialize perturbed state variables
    !--------------------------------------------------------------------------
    Rho = 0.0; U_D = 0.0; B_D = 0.0; p = 0.0

    if(.not.UseCme) RETURN

    ! linearly decay the perturbed magnetic field to 0 during tDecay time
    Coeff = 1.0
    if(tDecayCmeDim > 0.0) Coeff = max(0.0, &
         1 - (Time - tStartCme)/tDecayCmeDim)
#ifndef _OPENACC
    if (UseTD) then
       call get_TD99_fluxrope(Xyz_D, B1_D, Rho1, p1)
       Rho = Rho + Coeff*Rho1; B_D = B_D + Coeff*B1_D; p = p + Coeff*p1
    end if

    if(UseGL)then
       ! Add Gibson & Low (GL98) flux rope
       call get_GL98_fluxrope(Xyz_D, Rho1, p1, B1_D, U1_D, Time) !! send Time
       B_D = B_D + Coeff*B1_D
    end if
#endif

    if(UseSpheromak)then
       call get_GL98_fluxrope(Xyz_D, Rho1, p1, B1_D, U1_D, Time)
       B_D = B_D + Coeff*B1_D; U_D = U_D + Coeff*U1_D
    endif

#ifndef _OPENACC    
    if(UseShearFlow)then
       call get_shearflow(Xyz_D, Time, U1_D, nIter)
       U_D = U_D + Coeff*U1_D
    end if

    if(UseCms) call get_cms(Xyz_D, B_D)
#endif
  end subroutine EEE_get_state_BC
  !============================================================================
  subroutine EEE_get_state_init(Xyz_D, Rho, B_D, p, nStep, nIter)

    use EEE_ModCommonVariables, ONLY: UseCme, DoAddFluxRope, DoAddTD, &
         DoAddGL, UseCms, DoAddSpheromak, tStartCme
    use EEE_ModGL98, ONLY: get_GL98_fluxrope
    use EEE_ModTD99, ONLY: get_TD99_fluxrope
    use EEE_ModCms,  ONLY: get_cms

    real, intent(in) :: Xyz_D(3)
    real, intent(out) :: Rho, B_D(3), p
    integer, intent(in) :: nStep, nIter
    real :: U_D(3)
    real :: Rho1, U1_D(3), B1_D(3), p1
    !--------------------------------------------------------------------------
    ! initialize perturbed state variables
    Rho = 0.0; U_D = 0.0; B_D = 0.0; p = 0.0

    if(.not. (UseCme .and. (DoAddFluxRope .or. UseCms))) RETURN

    if(DoAddTD)then
       ! Add Titov & Demoulin (TD99) flux rope
       call get_TD99_fluxrope(Xyz_D, B1_D, Rho1, p1)
       Rho = Rho + Rho1; B_D = B_D + B1_D; p = p + p1
    endif

    if(DoAddGL)then
       ! Add Gibson & Low (GL98) flux rope
       call get_GL98_fluxrope(Xyz_D, Rho1, p1, B1_D, U1_D)
       Rho = Rho + Rho1; B_D = B_D + B1_D; p = p + p1
    end if

    if(DoAddSpheromak)then
       call get_GL98_fluxrope(Xyz_D, Rho1, p1, B1_D, U1_D, tStartCme)
       Rho = Rho + Rho1; B_D = B_D + B1_D
       p = p + p1; U_D = U_D + U1_D
    end if
    if(UseCms) call get_cms(Xyz_D, B_D)

  end subroutine EEE_get_state_init
  !============================================================================
  subroutine EEE_do_not_add_cme_again

    ! Designed to avoid repeatedly adding CME in subsequent sessions

    use EEE_ModCommonVariables, ONLY: DoAddFluxRope, &
         DoAddTD, DoAddGL, DoAddSpheromak
    !--------------------------------------------------------------------------
    DoAddFluxRope = .false.; DoAddGL = .false.; DoAddTD = .false.
    DoAddSpheromak = .false.

  end subroutine EEE_do_not_add_cme_again
  !============================================================================
  subroutine EEE_set_plot_range(nXY, nZ)
    use EEE_ModCommonVariables, ONLY: UseCme, UseTD, UseGL, UseSpheromak, &
         DXyzPlot, ExtensionFactor
    use EEE_ModTD99, ONLY: get_TD99_size
    use EEE_ModGL98, ONLY: get_GL98_size
    integer, intent(out) :: nXY, nZ
    real :: SizeXY, SizeZ
    !--------------------------------------------------------------------------
    nXY = -1; nZ = -1
    if(.not.UseCme)RETURN
    if(UseGL.or.UseSpheromak)then
       call get_GL98_size(SizeXY,  SizeZ)
    elseif(UseTD) then
       call get_TD99_size(SizeXY,  SizeZ)
    else
       RETURN
    end if
    nXY = nint(ExtensionFactor*SizeXY/DXyzPlot)
    nZ  = nint(ExtensionFactor*SizeZ /DXyzPlot)
  end subroutine EEE_set_plot_range
  !============================================================================
end module EEE_ModMain
!==============================================================================
