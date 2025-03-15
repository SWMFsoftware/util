!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModTransitionRegion

  use ModExactRS
  use ModUtilities, ONLY: CON_stop
  use ModConst,       ONLY: cBoltzmann, ceV, cProtonMass, cGravitation, mSun,&
       rSun, cMu
  implicit none
  SAVE
  PRIVATE ! except for a random set of 100 things below

  ! Public members
  public :: init_tr      ! Initialize module, check analytical TR model table
  public :: get_trtable_value     ! Access the TR model table
  public :: integrate_emission    ! Emission from the analytic TR model
  public :: plot_tr               ! Plot emission from the analytic TR model
  public :: read_tr_param         ! Read TR model parameters
  public :: check_tr_param        ! Synchtonize TR and SC model parameters
  public :: allocate_thread_arr   ! Threads originating from face centers
  public :: deallocate_thread_arr ! Deallocate threads
  public :: set_thread            ! Trace a thread from the face center
  public :: advance_thread_expl   ! Advance hd motion on thread
  public :: advance_thread_semi_impl   ! Advance heat conduction on thread
  public :: save_plot_thread      ! Plot state vector distribution on thread

  type OpenThread
     ! The Length along the thread, from the face i to the face i+1
     ! Exceptions: Ds_G(-nCell-1) is the length from
     ! the chromosphere to the face with the index -nCell
     ! Ds_G(0) is the distance from the center of interface
     ! between the threaded gap and the computational domain to
     ! simulate SC, and the point used to interpolate the state variables
     ! at the said interface. In meters
     real, pointer :: Ds_G(:)
     ! Magnetic field intensity (SI) in the points further associated to
     ! being the faces of the control volume scheme
     real, pointer :: B_F(:)
     ! (Dimensionless) spherical coordinates (of the faces)
     real, pointer :: Coord_DF(:,:)
     ! Heliocentric distance at the face center.
     real, pointer :: R_F(:)
     ! CELL (nor face!) centered array of the physical quantities
     ! to be solved for the given thread
     real, pointer :: State_VG(:,:)
     ! Cell-centered array of B1 (induction) field
     real, pointer :: B1_DG(:,:)
     ! Unit direction vector of the cell-centered magnetic field
     real, pointer :: DirB_DG(:,:)
     real, pointer :: ConservativeOld_VC(:,:)
     real, pointer :: Te_G(:)
     ! Since the arrays are further updated with different routines,
     ! we store the intermediate results of simulation:
     ! 1. Local time step (set be the hydro update procedure)
     real, pointer :: Dt_C(:)
     ! 2. The transition region parameters, solved in
     ! the heat conduction routine, and the global time step
     real :: TeTr = -1.0, uTr = -1.0, PeTr = -1.0, Dt = -1.0
     ! number of cells, the _C arrays has the index range
     ! from -nCell to -1
     integer :: nCell = -1
     real    :: OpenFlux = 0.0 ! [T Rsun**2]
  end type OpenThread
  public :: OpenThread

  ! Named indexes for state variables
  !
  ! State variables:
  ! in the state VG array:+++++++++++++++
  ! 1 - density
  ! 2 - velocity
  ! 3 - ion par pressure
  ! 4 - ion pressure = 2/3 PerpPressure + 1/3 ParPressure
  ! 5 - electron par pressure
  ! 6 - electron pressure = 2/3 Perp pressure + 1/3 ParPressure
  ! 7 - "Major wave" - propagating from the Sun outward/representative
  ! 8 - "Minor wave" - propagating toward the Sun/representative
  ! Primitive variables:+++++++++++++++++
  ! 1 - density
  ! 2 - velocity
  ! 3 - PparIon
  ! 4 - PperpIon
  ! 5 - PparElectron
  ! 6 - PperpElectron
  ! 7 - "Major wave" - propagating from the Sun outward/representative
  ! 8 - "Minor wave" - propagating toward the Sun/representative
  ! Conservative variables:++++++++++++++
  ! 1 - mass density
  ! 2 - momentum density
  ! 3 - ion energy density (kinetic+parallel internal)
  ! 4 - PperpIon
  ! 5 - Half of PparElectron (internal energy in parallel component)
  ! 6 - PperpElectron
  ! 7 - "Major wave" - propagating from the Sun outward/representative
  ! 8 - "Minor wave" - propagating toward the Sun/representative
  integer, public, parameter :: Rho_ = 1, U_ = 2, RhoU_ = 2, Ppar_ = 3, &
       Energy_ = 3, P_ = 4, Pperp_= 4, PePar_ = 5, PePerp_ = 6, Pe_ = 6,&
       Wmajor_ = 7, Wminor_ = 8 !, Wdiff_ = 9

  ! The Poynting flux to magnetic field ratio (one of the input parameters
  ! in SI unins)
  real, public :: PoyntingFluxPerBSi = 1.0e6 ! W/(m^2 T)
  real, public :: LperpTimesSqrtBSi = 7.5e4  ! m T^(1/2)
  real :: rMinReflectionTr = 0.0

  ! Normalization as used in the radcool table
  real, parameter :: RadNorm = 1.0E+22

  ! Correction coefficient to transform the units as used in the radcool
  ! table to SI system: erg = 1e-7 J / cm3 = 1e-6 m3 = 0.1
  real, parameter :: Cgs2SiEnergyDens = 0.1

  ! To find the volumetric radiative cooling rate in J/(m^3 s)
  ! the value found from radcool table should be multiplied by
  ! RadcoolSi and by n_e[m^{-3}]*n_i[m^{-3}]
  real, parameter :: Radcool2Si = 1.0e-12 & ! (cm-3=>m-3)**2
                          *Cgs2SiEnergyDens/RadNorm

  ! A constant factor to calculate the electron heat conduction
  real, public :: HeatCondParSi

  real :: cExchangeRateSi

  ! Gravity potential at the solar surface (m^2 s^-2), negative
  real, parameter :: cGravityPotAt1Rs = -cGravitation*mSun/rSun

  ! Coulomb logarithm
  real, public  :: CoulombLog = 20.0

  real, parameter :: cTwoSevenths = 2.0/7.0
  real,   public  :: cTolerance   = 1.0e-6

  ! Correspondent named indexes: meaning of the columns in the table
  integer, parameter, public :: LengthPavrSi_ = 1, uHeat_ = 2, &
       HeatFluxLength_ = 3, dHeatFluxXOverU_ = 4, LambdaSi_=5, &
       DlogLambdaOverDlogT_ = 6
  ! Tabulated analytical solution:
  real, public :: TrTable_V(LengthPAvrSi_:DlogLambdaOverDlogT_)
  ! integer, parameter :: R_=1 ! to mark the radial coordinate, in R_sun

  ! Control parameter: below TeSiMin the observables are not calculated
  real, public :: TeSiMin = 5.0e4
  ! Average ion charge number and its square root
  real, public :: SqrtZ   = 1.0
  real :: Z = 1.0

  ! Table numbers needed to use lookup table
  integer         :: iTableRadCool = -1
  integer, public :: iTableTr = -1

  ! Chromosphere top boundary, in Rsun
  real, public, parameter :: rChromo = 1.0

  ! By default, this logical is .false. If logical set is true, the energy
  ! flux to/from the first control volume on the thread is accounted for.
  logical, public :: UseChromoEvap  = .false.
  ! To apply CromoEvaporation, the factor below is non-zero.
  real,    public :: ChromoEvapCoef = 0.0
  ! Parameters of the TR table:
  real, parameter :: TeTrMin = 1.0e4
  real, parameter :: TeTrMax = real(2.9988442312309672D+06)
  integer, parameter :: nPointTe = 310, nPointU = 89
  ! Global array used in calculating the table
  ! To be allocated as Value_VII (DlogLambdaOverDlogT_,nPointTe,nPointU)
  real, allocatable :: Value_VII(:,:,:)
  real            :: DeltaLogTe
  ! Max speed
  real, parameter :: uMax = 2000.0, uMin = -200.0
  real, parameter :: DeltaU = (uMax - uMin)/(nPointU - 1)
  ! Parameter of the chromosphere energy budget:
  ! How many hydrogen atoms should be ionized by the heat conduction flux
  ! to lift a single pair of proton+electron to the solar wind
  real, parameter :: cEnergyEfficiency = 1.0
  ! Needed for initialization:
  logical :: DoInit = .true.
  interface advance_heat_conduction
     module procedure advance_heat_conduction_ta
     module procedure advance_heat_conduction_ss
  end interface advance_heat_conduction
  ! Inner boundary
  ! To be read from the parameter file
  real :: TempInner = 2.0e6, PressInner = 4.0e-3
  ! Dimensionless parameters for stochastic heating
  logical :: UseStochasticHeating = .true.
  real :: StochasticExponent   = 0.21
  real :: StochasticAmplitude  = 0.18
  real :: StochasticExponent2  = 0.21
  real :: StochasticAmplitude2 = 0.0 ! 1.17
  ! Threads are traced from rMax toward rMin
  real,    public :: rMin = rChromo, rMax = 2.5
  integer, public :: nPointMax = 10000
  ! integration step at R = 1, dS propto R for R > 1
  real, public :: dS0 = 0.001
  ! Field unit conversion:
  real, public, parameter :: Si2Gs = 10000.0, Gs2Si = 1.0e-4
  ! Minimum field at the source surface and at the inner boundary
  real, public :: BssMinSi = 0.002*Gs2Si
  real, public :: BMinSi  = 0.0125*Gs2Si

  ! Minimum pressure and density
  real, public :: MinPress = 0.0
  real, public :: MinRho = 0.0
  ! Constant values read from the parameter file,
  ! if the uniform partitioning is assumed for dissipated turbulent energy
  real :: QparPerQtotal, QperpPerQtotal, QePerQtotal
  ! Time stepping
  real, public :: CflLocal = 0.5
  logical :: DoLimitLogVar = .false.
  integer, parameter :: iLogVar_V(5) = &
       [Rho_, Ppar_, Pperp_, PePar_, PePerp_]
  character(len=*), parameter:: NameMod = 'ModTransitionRegion'
contains
  !============================================================================
  subroutine init_tr(zIn, TeChromoSi, iComm)
    use ModConst,       ONLY: kappa_0_e, te_ti_exchange_rate
    use ModLookupTable, ONLY: i_lookup_table
    !    use BATL_lib,       ONLY: test_start, test_stop, iProc
    ! INPUTS
    ! Average ion charge
    real, intent(in) :: zIn
    ! Temperature on top of chromosphere
    real, intent(in) :: TeChromoSi
    integer, optional, intent(in) :: iComm
    logical:: DoTest
    !   call test_start(NameSub, DoTest)
    character(len=*), parameter:: NameSub = 'init_tr'
    !--------------------------------------------------------------------------
    if(.not.DoInit)RETURN
    DoInit = .false.
    ! Set model parameters
    Z = zIn
    SqrtZ = sqrt(Z)
    TeSiMin = TeChromoSi
    DeltaLogTe = log(TeTrMax/TeTrMin)/(nPointTe - 1)
    ! electron heat conduct coefficient for single charged ions
    ! = 9.2e-12 W/(m*K^(7/2))
    HeatCondParSi   = kappa_0_e(CoulombLog)
    cExchangeRateSi = te_ti_exchange_rate(CoulombLog)
    ! Init Riemann solver:
    call exact_rs_set_gamma(3.0)
    iTableTr = i_lookup_table('TR')
    if(iTableTr<=0)then
      iTableRadCool = i_lookup_table('radcool')
       if(iTableRadCool <=0 )&
            call CON_stop('To create TR table, the radcool table is needed')
       call check_tr_table(iComm=iComm)
    end if
    ! call test
    ! call CON_stop('Check test')
    !  To make a unit test:
    ! 1. Uncomment the lines below as well as the declaration and the use
    ! of DoTest, test_start and test_stop
    ! 2. Add init_tr identifier to the #TEST command for SC component in
    !    Param/PARAM.in.test.SCIH_threadbc
    !    Param/PARAM.in.test.restart.SCIH_threadbc
    ! 3. run test8
    ! In run_test there will be files test_tr.out and test_tr_AiaXrt.out
    ! if(DoTest.and.iProc==0)then
    !   ! Plot test solution
    !   iTableAiaXrt = i_lookup_table('AiaXrt')
    !   if(iTableAiaXrt>0)then
    !      call plot_tr(&
    !           NamePlotFile='test_tr_AiaXrt.out',&
    !           nGrid       = 100,         &
    !           TeSi        = 5.0e5,       &
    !           PeSi        = 5.0e-3,      &
    !           iTable      = iTableAiaXrt)
    !   else
    !      call plot_tr(&
    !           NamePlotFile='test_tr.out',&
    !           nGrid       = 100,         &
    !           TeSi        = 5.0e5,       &
    !           PeSi        = 5.0e-3       )
    !   end if
    ! end if
  end subroutine init_tr
  !============================================================================
  subroutine read_tr_param(NameCommand)

    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand
    character(len=*), parameter:: NameSub = 'read_tr_param'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case('#CHROMOEVAPORATION')
       DoLimitLogVar = .true.
       MinPress = 1.e-14     ! Pa
       MinRho   = 1.e-23     ! kg/m3
       call read_var('UseChromoEvap',UseChromoEvap)
       if(UseChromoEvap)then
          ChromoEvapCoef = TeTrMin
       else
          ChromoEvapCoef = 0.0
       end if
    case("#STOCHASTICHEATING")
       UseStochasticHeating = .true.
       ! Stochastic heating when Beta_proton is below 1
       call read_var('StochasticExponent', StochasticExponent)
       call read_var('StochasticAmplitude', StochasticAmplitude)
    case("#HIGHBETASTOCHASTIC")
       ! Correction for stochastic heating when Beta_proton is between 1 and 30
       ! KAWs are non-propagating for Beta_proton > 30.
       call read_var('StochasticExponent2', StochasticExponent2)
       call read_var('StochasticAmplitude2', StochasticAmplitude2)
    case("#TURBULENCE")
       call read_var("PoyntingFluxPerBSi",PoyntingFluxPerBSi)
       call read_var("LperpTimesSqrtBSi",LperpTimesSqrtBSi)
       call read_var("rMinReflectionTr", rMinReflectionTr)
    case("#INNERBOUNDARY")
       call read_var('TempInner', TempInner)
       call read_var('PressInner', PressInner)
    case("#MINPRESSTR")
       call read_var('MinPress', MinPress)
    case("#MINRHOTR")
       call read_var('MinRho', MinRho)
    case("#LIMITLOGVAR")
       call read_var('DoLimitLogVar', DoLimitLogVar)
    case('#UNIFORMPARTITION')
       UseStochasticHeating = .false.
       call read_var('QparPerQtotal', QparPerQtotal)
       call read_var('QperpPerQtotal', QperpPerQtotal)
       QePerQtotal = 1 - (QparPerQtotal + QperpPerQtotal)
    case default
       call CON_stop(NameSub//": unknown command="//trim(NameCommand))
    end select
  end subroutine read_tr_param
  !============================================================================
  subroutine check_tr_param(rMaxIn, rMinReflectionTrIn, &
       TempInnerIn, PressInnerIn,                       &
       UseStochasticHeatingIn,                          &
       StochasticExponentIn, StochasticAmplitudeIn,     &
       StochasticExponent2In, StochasticAmplitude2In,   &
       QparPerQtotalIn, QperpPerQtotalIn)
    real, intent(in) :: rMaxIn, rMinReflectionTrIn
    real, intent(in) :: TempInnerIn, PressInnerIn ! Coronal base parameters SI
    ! Dimensionless parameters for turbulent heating
    logical, intent(in) :: UseStochasticHeatingIn
    ! If UseStochastic heating, assign the following:
    real, optional, intent(in) :: StochasticExponentIn
    real, optional, intent(in) :: StochasticAmplitudeIn
    real, optional, intent(in) :: StochasticExponent2In
    real, optional, intent(in) :: StochasticAmplitude2In
    ! Else
    real, optional, intent(in) :: QparPerQtotalIn
    real, optional, intent(in) :: QperpPerQtotalIn
    !--------------------------------------------------------------------------
    rMax = rMaxIn
    rMinReflectionTr = rMinReflectionTrIn
    TempInner = TempInnerIn
    PressInner = PressInnerIn
    UseStochasticHeating = UseStochasticHeatingIn
    if(UseStochasticHeating)then
       if(present(StochasticExponentIn))&
            StochasticExponent = StochasticExponentIn
       if(present(StochasticAmplitudeIn))&
            StochasticAmplitude = StochasticAmplitudeIn
       if(present(StochasticExponent2In))&
            StochasticExponent2 = StochasticExponent2In
       if(present(StochasticAmplitude2In))&
            StochasticAmplitude2 = StochasticAmplitude2In
    else
       if(present(QparPerQtotalIn))QparPerQtotal = QparPerQtotalIn
       if(present(QperpPerQtotalIn))QperpPerQtotal = QperpPerQtotalIn
       QePerQtotal = 1 - (QparPerQtotal + QperpPerQtotal)
    end if
  end subroutine check_tr_param
  !============================================================================
  subroutine check_tr_table(TypeFileIn,iComm)
    use ModLookupTable, ONLY: i_lookup_table, &
         init_lookup_table, make_lookup_table, interpolate_lookup_table

    character(LEN=*),optional,intent(in)::TypeFileIn
    integer, optional, intent(in):: iComm

    character(len=5)::TypeFile
    ! Ionization potential for hydrogen
    ! cPotential_II(1,1) from util/CRASH/src/ModIonizPotential.f90:
    real, parameter :: cIonizPotentialH =  13.59844  ! eV
    ! Ratio of enthalpy flux, divided by (u(Tmin)*Pe):
    real, parameter :: FluxTouPe = cEnergyEfficiency*&
         cIonizPotentialH*ceV/(cBoltzmann*TeTrMin)   & ! ionization loss,
         + 5                                           ! enthalpy, per ion
    ! Misc:
    ! Loop variables
    integer            :: iTe, iU
    ! Arguments corresponding to these indexes:
    real   :: TeSi_I(nPointTe), uTr

    real, dimension(nPointTe) :: LambdaSi_I, DLogLambdaOverDLogT_I, &
         LengthPe_I, UHeat_I, dFluxXLengthOverDU_I, dHeatFluxOverVel_I
    real            :: SemiIntUheat_I(1:nPointTe-1)
    real            :: DeltaLogTeCoef
    real   :: LambdaCgs_V(1)

    ! at the moment, radcool not a function of Ne, but need a dummy 2nd
    ! index, and might want to include Ne dependence in table later.
    ! Table variable should be normalized to radloss_cgs * 10E+22
    ! since we don't want to deal with such tiny numbers
    ! real, parameter :: RadNorm = 1.0E+22
    ! real, parameter :: Cgs2SiEnergyDens = &
    !     1.0e-7&   ! erg = 1e-7 J
    !     /1.0e-6    ! cm3 = 1e-6 m3
    ! real, parameter :: Radcool2Si = 1.0e-12 & ! (cm-3=>m-3)**2
    !     /RadNorm*Cgs2SiEnergyDens
    ! Misc:
    real :: FactorStep, uOverTmin5, SqrtOfU2, KinEnergyFlux

    character(len=*), parameter:: NameSub = 'check_tr_table'
    !--------------------------------------------------------------------------
    if(present(TypeFileIn))then
       TypeFile = TypeFileIn
    else
       TypeFile = 'real8'
    end if
    ! initialize the TR table
    call init_lookup_table(                                      &
         NameTable = 'TR',                                       &
         NameCommand = 'save',                                   &
         Param_I = [cEnergyEfficiency, uMin, uMax],              &
         NameVar =                                               &
         'logTe u LPe UHeat FluxXLength '//                      &
         'dFluxXLegthOverDU Lambda dLogLambdaOverDLogT '//&
         'EnergyEff uMin uMax',                                  &
         nIndex_I = [nPointTe,nPointU],                          &
         IndexMin_I = [TeTrMin, uMin],                           &
         IndexMax_I = [TeTrMax, uMax],                           &
         NameFile = 'TR8.dat',                                   &
         TypeFile = TypeFile,                                    &
         StringDescription =                                     &
         'Model for transition region: '//                       &
         '[K] [m/s] [N/m] [m/s] [W/m] [1] [W*m3/(k_B2)] [1] [W/m]'//&
         ' [1] [m/s] [m/s]')

    ! The table is now initialized.
    iTableTr = i_lookup_table('TR')
    ! Fill in the array of tabulated values:
    allocate(Value_VII(DlogLambdaOverDlogT_,nPointTe,nPointU))

    FactorStep = exp(DeltaLogTe)
    ! Fill in TeSi array
    TeSi_I(1) = TeTrMin
    do iTe = 2, nPointTe
       TeSi_I(iTe) = TeSi_I(iTe-1)*FactorStep
    end do
    ! Fill in LambdaSi
    do iTe = 1, nPointTe
       call interpolate_lookup_table(iTableRadCool,&
            TeSi_I(iTe), LambdaCgs_V)
       LambdaSi_I(iTe) = LambdaCgs_V(1)*Radcool2Si
    end do
    ! Calculate dLogLambda/DLogTe
    DLogLambdaOverDLogT_I(1) = &
         log(LambdaSi_I(2)/LambdaSi_I(1))/DeltaLogTe
    do iTe = 2,nPointTe-1
       DLogLambdaOverDLogT_I(iTe) = &
            log(LambdaSi_I(iTe+1)/LambdaSi_I(iTe-1))/(2*DeltaLogTe)
    end do
    DLogLambdaOverDLogT_I(nPointTe) = &
         log(LambdaSi_I(npointTe)/LambdaSi_I(nPointTe-1))/DeltaLogTe

    DeltaLogTeCoef = DeltaLogTe*HeatCondParSi/cBoltzmann**2
    do iU = 1, nPointU
       uTr = (iU - 1)*DeltaU + uMin
       uOverTmin5 = 5*(uTr/TeTrMin)*DeltaLogTe
       KinEnergyFlux = (cProtonMass/cBoltzmann)*(uTr/TeTrMin)**3*DeltaLogTe
       LengthPe_I = 0.0; uHeat_I = 0.0
       UHeat_I(1) = (max(uTr, 0.0)*(FluxTouPe + &
            0.50*cProtonMass*uTr**2/(cBoltzmann*TeTrMin)))**2
       do iTe = 2, nPointTe
          ! Integrate \sqrt{2\int{\kappa_0\Lambda Te**1.5 d(log T)}}/k_B
          ! Predictor at half step
          SqrtOfU2 = sqrt(UHeat_I(iTe-1) )
          SemiIntUheat_I(iTe-1) = sqrt( UHeat_I(iTe-1) + &
               SqrtOfU2*(uOverTmin5*TeSi_I(iTe-1)      + &
               KinEnergyFlux*TeSi_I(iTe-1)**2)         + &
               LambdaSi_I(iTe-1)*TeSi_I(iTe-1)**1.50*DeltaLogTeCoef)
          UHeat_I(iTe) = UHeat_I(iTe-1) + SemiIntUheat_I(iTe-1)*&
               (uOverTmin5*(TeSi_I(iTe-1) + TeSi_I(iTe)) +&
               KinEnergyFlux*(TeSi_I(iTe-1)**2 + TeSi_I(iTe)**2) ) &
               + (LambdaSi_I(iTe-1)*TeSi_I(iTe-1)**1.50 + &
               LambdaSi_I(iTe)*TeSi_I(iTe)**1.50)*DeltaLogTeCoef
          UHeat_I(iTe-1) = SqrtOfU2
       end do
       UHeat_I(nPointTe) = sqrt(UHeat_I(nPointTe))

       do iTe = 2, nPointTe
          ! Integrate \int{\kappa_0\Lambda Te**3.5 d(log T)/UHeat}
          LengthPe_I(iTe) = LengthPe_I(iTe-1) + 0.5*DeltaLogTe* &
               ( TeSi_I(iTe-1)**3.5 + TeSi_I(iTe)**3.5 )&
               /SemiIntUheat_I(iTe-1)
          ! Not multiplied by \kappa_0
       end do
       dFluxXLengthOverDU_I(1) = &
            (LengthPe_I(2)*UHeat_I(2) - &
            LengthPe_I(1)*UHeat_I(1))/&
            (DeltaLogTe*TeSi_I(1)**3.5)
       do iTe = 2, nPointTe - 1
          dFluxXLengthOverDU_I(iTe) = &
               ( LengthPe_I(iTe+1)*UHeat_I(iTe+1)   &
               - LengthPe_I(iTe-1)*UHeat_I(iTe-1) ) &
               /(2*DeltaLogTe*TeSi_I(iTe)**3.5)
       end do
       dFluxXLengthOverDU_I(nPointTe) = &
            (LengthPe_I(nPointTe)*UHeat_I(nPointTe) - &
            LengthPe_I(nPointTe-1)*UHeat_I(nPointTe-1))/&
            (DeltaLogTe*TeSi_I(nPointTe)**3.5)

       LengthPe_I(:) = LengthPe_I(:)*HeatCondParSi
       do iTe = 1, nPointTe
          Value_VII(LengthPAvrSi_:DlogLambdaOverDlogT_, iTe, iU) = &
               [ LengthPe_I(iTe), UHeat_I(iTe), &
               LengthPe_I(iTe)*UHeat_I(iTe), dFluxXLengthOverDU_I(iTe),   &
               LambdaSi_I(iTe)/cBoltzmann**2, DLogLambdaOverDLogT_I(iTe)]
       end do
    end do
    ! Fill in the velocity derivative
    do iU = 2, nPointU - 1
       uTr = (iU - 1)*DeltaU + uMin
       dHeatFluxOverVel_I = &
            (Value_VII(Uheat_,:,iU+1)*Value_VII (LengthPavrSi_,:,iU+1) - &
            Value_VII (Uheat_,:,iU-1)*Value_VII (LengthPavrSi_,:,iU-1))/ &
            (2*DeltaU)
       Value_VII(dHeatFluxXOverU_,:,iU) = Value_VII(dHeatFluxXOverU_,:,iU) - &
            uTr*dHeatFluxOverVel_I/(HeatCondParSi*TeSi_I(:)**3.5)
    end do
    ! Shape the table using the filled in array
    call make_lookup_table(iTableTr, calc_tr_table, iComm)
    deallocate(Value_VII)
  end subroutine check_tr_table
  !============================================================================
  subroutine calc_tr_table(iTableIn, Arg1, Arg2, Value_V)

    integer, intent(in):: iTableIn
    real, intent(in)   :: Arg1, Arg2
    real, intent(out)  :: Value_V(:)
    integer            :: iTe, iU
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_tr_table'
    !--------------------------------------------------------------------------
    iTe = 1 + nint(log(Arg1/TeTrMin)/DeltaLogTe)
    iU = nint( (Arg2 - uMin)/ DeltaU) + 1
    Value_V(:) = Value_VII(:, iTe, iU)
  end subroutine calc_tr_table
  !============================================================================
  subroutine get_trtable_value(Te, uFace)

    ! For two entries, Te and u on top of the transition region,
    ! the routine calculates TrTable_V array of the tabulated values
    ! If uIn is not provided, it is taken to be 0.
    ! If optional uOut is requested, the limited value of u is provided
    ! calculated with the value uMin < uTr < uMax at the bottom of
    ! transition region.
    ! All inputs and outputs are in SI units
    use ModLookupTable,  ONLY: interpolate_lookup_table
    real, intent(in) :: Te  ! Temperature on top of the transition region
    ! Speed on top of the transition region
    real, OPTIONAL, intent(inout) :: uFace
    ! Misc:
    ! Actual table entry, speed at the bottom of TR:
    real :: uTr
    !--------------------------------------------------------------------------
    if(present(uFace))then
       uTr = min(uMax, max(ChromoEvapCoef*uFace/Te, uMin))
       uFace = uTr*Te/TeTrMin
    else
       uTr = 0.0
    end if
    call interpolate_lookup_table(iTableTr, Te, uTr, TrTable_V, &
         DoExtrapolate=.false.)
  end subroutine get_trtable_value
  !============================================================================
  subroutine set_thread(XyzIn_D, FaceArea, OpenThread1, &
       xyz_to_coord, get_field)
    ! Origin point coordinates (Xyz)
    real, intent(in) :: XyzIn_D(3)
    ! Face area, to calculate flux
    real, intent(in) :: FaceArea
    ! Thread to set
    type(OpenThread), intent(inout) :: OpenThread1
    interface
       subroutine xyz_to_coord(Xyz_D, Coord_D)
         implicit none
         real, intent(in) :: Xyz_D(3)
         real, intent(out):: Coord_D(3)
       end subroutine xyz_to_coord
       subroutine get_field(Xyz_D, B_D, B1Out_D)
         implicit none
         real, intent(in) :: Xyz_D(3)
         real, intent(out):: B_D(3)
         real, OPTIONAL, intent(out) :: B1Out_D(3)
       end subroutine get_field
    end interface
    ! loop variable
    integer :: iPoint, nCell
    ! Length interval, ! Heliocentric distance
    real :: Ds, R
    ! coordinates, field vector and modulus
    real :: XyzStart_D(3), Xyz_D(3), B_D(3), B, BMax
    !  for ourward directed field, -1 otherwise
    real :: SignBr
    ! Radial field at the source surface
    real :: BrSS
    ! Coordinates and magnetic field in the midpoint
    ! within the framework of the Runge-Kutta scheme
    real :: XyzAux_D(3)!, Baux_D(3)
    ! Aux
    real :: ROld, Aux
    real :: XyzOld_D(3), Dir1_D(3), Dir2_D(3), Dir3_D(3), Dir4_D(3)
    real :: B_F(-nPointMax:0), Length_I(-nPointMax:0), &
         B1_DF(3,-nPointMax:0), DirB_DG(3,-nPointMax:0)
    real :: Coord_DI(3,-nPointMax:0), R_F(-nPointMax:0)
    real :: CosBrMin = -1.0
    integer :: iRefine ! Factor of step refinement near null point
    ! Trace the open thread with the given origin point
    character(len=*), parameter:: NameSub = 'set_thread'
    !--------------------------------------------------------------------------
    XyzStart_D = XyzIn_D
    R = norm2(XyzStart_D)
    CosBrMin = -1.0   ! Arbitrary angle between magnetic field and radial dir
    iRefine = 10      ! Step refinement by a factor of 1/10 is allowed
    call get_field(XyzStart_D, B_D, B1_DF(:,0))
    BrSS = sum(XyzStart_D*B_D)/R
    OpenThread1%OpenFlux = FaceArea*BrSS
    ! Unit direction vector with positive radial component
    DirB_DG(:,0) = B_D/sign(norm2(B_D),OpenThread1%OpenFlux)
    SignBr = sign(1.0, BrSS)
    IBEGINLOOP: do
       B = norm2(B_D)
       call xyz_to_coord(XyzStart_D,Coord_DI(:,0))
       iPoint = 0
       B_F(0) = B
       BMax = max(B, BssMinSi)
       Xyz_D = XyzStart_D
       R_F(0) = R
       POINTS: do
          iPoint = iPoint + 1
          ! For the previous point given are Xyz_D, B_D, B
          ! Store R
          ROld = R
          ! Store a point
          XyzOld_D = Xyz_D
          ! Four-stage Runge-Kutta
          BMax = max(BMax, B)
          Ds = Ds0*R*SignBr/min(BMax, iRefine*B)
          Dir1_D = B_D
          if(CosBrMin /= -1.0)call limit_cosbr(Xyz_D, Dir1_D)
          ! 1. Point at the half of length interval:
          XyzAux_D = Xyz_D - 0.50*Ds*Dir1_D
          ! 2. Magnetic field in this point:
          call get_field(XyzAux_D, Dir2_D)
          if(CosBrMin /= -1.0)call limit_cosbr(XyzAux_D, Dir2_D)
          XyzAux_D = Xyz_D - 0.50*Ds*Dir2_D
          call get_field(XyzAux_D, Dir3_D)
          if(CosBrMin /= -1.0)call limit_cosbr(XyzAux_D, Dir3_D)
          XyzAux_D = Xyz_D - Ds*Dir3_D
          call get_field(XyzAux_D, Dir4_D)
          if(CosBrMin /= -1.0)call limit_cosbr(XyzAux_D, Dir4_D)
          ! 3. New grid point:
          Xyz_D = Xyz_D - (Ds/6)*(Dir1_D + 2*Dir2_D + 2*Dir3_D + Dir4_D)
          R = norm2(Xyz_D)
          if(R >  rMax.or.iPoint==nPointMax)then
             R = norm2(XyzStart_D)
             call get_field(XyzStart_D, B_D)
             ! Don't allow magnetic field too much decline from radial dir
             CosBRMin = ( (R**2-rChromo**2)/nPointMax +Ds0**2)/&
                  (2*rChromo*Ds0)
             if(CosBRMin>0.9)call CON_stop('Increase nPointThreadMax')
             ! Prohibit step refinement
             iRefine = 1
             CYCLE IBEGINLOOP
          end if
          R_F(-iPoint) = R
          call xyz_to_coord(Xyz_D, Coord_DI(:,-iPoint))
          call get_field(Xyz_D, B_D, B1Out_D=B1_DF(:,-iPoint))
          B = norm2(B_D)
          B_F(-iPoint) = B
          Length_I(-iPoint) = norm2(XyzOld_D - Xyz_D)
          DirB_DG(:,-iPoint) = (XyzOld_D - Xyz_D)/Length_I(-iPoint)
          if(R <= rMin + Ds0*max(1 - abs(B)/BminSi,0.0))EXIT IBEGINLOOP
       end do POINTS
    end do IBEGINLOOP
    ! Calculate more accurately the intersection point
    ! with the photosphere surface
    Aux = (ROld - RMin) / (ROld - R)
    Xyz_D = (1 - Aux)*XyzOld_D +  Aux*Xyz_D
    ! Store the last point
    call xyz_to_coord(Xyz_D, Coord_DI(:,-iPoint))
    call get_field(Xyz_D, B_D)
    B = norm2(B_D)
    B_F(-iPoint) = B
    Length_I(-iPoint) = norm2(Xyz_D - XyzOld_D)
    ! R_F(-iPoint) = norm2(Xyz_D)
    R_F(-iPoint) = norm2(Xyz_D)
    nCell = iPoint - 2
    ! Allocate threads if needed
    if(nCell /= OpenThread1%nCell)then
       if(.not.associated(OpenThread1%Ds_G))then
          call allocate_pointer
       else
          call deallocate_pointer
          call allocate_pointer
       end if
    end if
    ! Store results from tracing
    OpenThread1%Ds_G(-nCell:-1) = &
         Length_I(-nCell:-1)*Rsun
    ! Length of the analytical transition region,
    ! combined from the last two intervals
    OpenThread1%Ds_G(-nCell-1) = &
         (Length_I(-iPoint) + Length_I(-iPoint+1))*Rsun
    ! Distance from the "ghost point" to the interface (0 index)
    OpenThread1%Ds_G(0) = 0.0
    OpenThread1%B_F(-iPoint:0) = B_F(-iPoint:0)
    OpenThread1%Coord_DF(:,-iPoint:0) = Coord_DI(:,-iPoint:0)
    allocate(OpenThread1%R_F(-nCell:0))
    OpenThread1%R_F(-nCell:0) = R_F(-nCell:0)
    OpenThread1%B1_DG(:,-nCell:-1) = 0.50*(&
         B1_DF(:,-nCell:-1) + B1_DF(:,1-nCell:0))
    OpenThread1%B1_DG(:,0) = B1_DF(:,0)
    OpenThread1%DirB_DG(:,-nCell:0) = DirB_DG(:,-nCell:0)
    ! Re-initialize state vector if needed
    if(nCell /= OpenThread1%nCell)then
       OpenThread1%nCell = nCell
       call init_thread_variables(OpenThread1)
    end if
  contains
    !==========================================================================
    subroutine allocate_pointer
      !------------------------------------------------------------------------
      allocate(OpenThread1%B_F(-iPoint:0))
      allocate(OpenThread1%Ds_G(-nCell-1:0))
      allocate(OpenThread1%Coord_DF(3,-iPoint:0))
      allocate(OpenThread1%Dt_C(-nCell:-1))
      allocate(OpenThread1%R_F(-nCell:0))
      allocate(OpenThread1%Te_G(-nCell:0))
      allocate(OpenThread1%B1_DG(3,-nCell:0))
      allocate(OpenThread1%DirB_DG(3,-nCell:0))
      allocate(OpenThread1%ConservativeOld_VC(Rho_:Wminor_,-nCell:-1))
      allocate(OpenThread1%State_VG(Rho_:Wminor_,-nCell:0))
    end subroutine allocate_pointer
    !==========================================================================
    subroutine deallocate_pointer
      !------------------------------------------------------------------------
      deallocate(OpenThread1%B_F)
      deallocate(OpenThread1%Ds_G)
      deallocate(OpenThread1%Coord_DF)
      deallocate(OpenThread1%Dt_C)
      deallocate(OpenThread1%R_F)
      deallocate(OpenThread1%Te_G)
      deallocate(OpenThread1%B1_DG)
      deallocate(OpenThread1%DirB_DG)
      deallocate(OpenThread1%ConservativeOld_VC)
      deallocate(OpenThread1%State_VG)
    end subroutine deallocate_pointer
    !==========================================================================
    subroutine init_thread_variables(OpenThread1)

      use ModConst,only : cBoltzmann, cProtonMass

      type(OpenThread),intent(inout) :: OpenThread1

      real :: RInvCenter_I(-OpenThread1%nCell:-1)
      integer :: nCell

      !------------------------------------------------------------------------
      nCell = OpenThread1%nCell
      RInvCenter_I(-nCell:-1) = &
           (1/OpenThread1%R_F(-nCell:-1) + &
           1/OpenThread1%R_F(-nCell+1: 0))*0.50

      ! initial velocity is zero

      OpenThread1%State_VG(U_,-nCell:-1) = 0.0

      ! exponential pressure

      OpenThread1%State_VG(P_,-nCell:-1) = PressInner * 0.5 * exp( &
           cProtonMass * mSun * cGravitation / cBoltzmann / &
           TempInner * (RInvCenter_I - 1.0/rMin) / rSun)
      OpenThread1%State_VG(Ppar_,-nCell:-1) = &
           OpenThread1%State_VG(P_,-nCell:-1)
      OpenThread1%State_VG(Pe_,-nCell:-1) = &
           OpenThread1%State_VG(P_,-nCell:-1)
      OpenThread1%State_VG(PePar_,-nCell:-1) = &
           OpenThread1%State_VG(Pe_,-nCell:-1)

      ! Pi = n * kB * T, n = Rho / Mp
      ! -> Rho = Mp * Pi / (kB * T)
      OpenThread1%State_VG(Rho_,-nCell:-1) = cProtonMass * &
           OpenThread1%State_VG(P_,-nCell:-1) / (cBoltzmann*TempInner)

      OpenThread1%State_VG(Wmajor_,-nCell:-1) = 1.0
      OpenThread1%State_VG(Wminor_,-nCell:-1) = 1.0e-8
      OpenThread1%Te_G(-nCell:0) = TempInner
      OpenThread1%TeTr = TempInner
      OpenThread1%uTr  = 0.0
      OpenThread1%PeTr = 0.50*PressInner

    end subroutine init_thread_variables
    !==========================================================================
    subroutine limit_cosbr(Xyz_D, B_D)

      real, intent(in)    :: Xyz_D(3) ! Location
      real, intent(inout) :: B_D(3)   ! Magnetic field to be corrected
      real                :: CosBR, B, DirR_D(3), DirB_D(3)
      !------------------------------------------------------------------------
      DirR_D = Xyz_D/norm2(Xyz_D)
      B = norm2(B_D)
      DirB_D = SignBr*B_D/B
      CosBR = sum(DirB_D*DirR_D)
      if(CosBR>=CosBRMin)RETURN
      DirB_D = (DirB_D - CosBR*DirR_D)*&      ! Tangential componenets
           sqrt((1 - CosBRMin**2)/(1 - CosBR**2))+& ! Reduced in magnitude
           DirR_D*CosBRMin          ! Plus increased radial comp
      B_D = SignBr*B*DirB_D
    end subroutine limit_cosbr
    !==========================================================================
  end subroutine set_thread
  !============================================================================
  subroutine allocate_thread_arr(Threads_II, nI, nJ)
    type(OpenThread), allocatable, intent(inout) :: Threads_II(:,:)
    integer, intent(in) :: nI, nJ
    integer :: i, j
    character(len=*), parameter:: NameSub = 'allocate_thread_arr'
    !--------------------------------------------------------------------------
    allocate(Threads_II(nI,nJ))
    do j = 1, nJ; do i = 1, nI
       call init_thread(Threads_II(i,j))
    end do; end do
  contains
    !==========================================================================
    subroutine init_thread(OpenThread1)

      ! Thread to set
      type(OpenThread), intent(inout) :: OpenThread1
      !------------------------------------------------------------------------

      nullify(OpenThread1%ConservativeOld_VC)
      nullify(OpenThread1%Dt_C)
      nullify(OpenThread1%State_VG)
      nullify(OpenThread1%B1_DG)
      nullify(OpenThread1%DirB_DG)
      nullify(OpenThread1%Te_G)
      nullify(OpenThread1%R_F)
      nullify(OpenThread1%Coord_DF)
      nullify(OpenThread1%B_F)
      nullify(OpenThread1%Ds_G)
      OpenThread1%TeTr = -1.0
      OpenThread1%uTr  = -1.0
      OpenThread1%PeTr = -1.0
      OpenThread1%Dt   = -1.0
      OpenThread1%OpenFlux = 0.0
      OpenThread1%nCell = -1
    end subroutine init_thread
    !==========================================================================
  end subroutine allocate_thread_arr
  !============================================================================
  subroutine deallocate_thread_arr(Threads_II, nI, nJ)
    type(OpenThread), allocatable, intent(inout) :: Threads_II(:,:)
    integer, intent(in) :: nI, nJ
    integer :: i, j
    character(len=*), parameter:: NameSub = 'deallocate_thread_arr'
    !--------------------------------------------------------------------------
    do j = 1, nJ; do i = 1, nI
       call deallocate_thread(Threads_II(i,j))
    end do; end do
    deallocate(Threads_II)
  contains
    !==========================================================================
    subroutine deallocate_thread(OpenThread1)

      ! Thread to set
      type(OpenThread), intent(inout) :: OpenThread1
      !------------------------------------------------------------------------
      OpenThread1%TeTr = -1.0
      OpenThread1%uTr  = -1.0
      OpenThread1%PeTr = -1.0
      OpenThread1%Dt   = -1.0
      OpenThread1%nCell = -1
      OpenThread1%OpenFlux = 0.0
      if(.not.associated(OpenThread1%Coord_DF))RETURN
      deallocate(OpenThread1%ConservativeOld_VC)
      deallocate(OpenThread1%Dt_C)
      deallocate(OpenThread1%State_VG)
      deallocate(OpenThread1%B1_DG)
      deallocate(OpenThread1%DirB_DG)
      deallocate(OpenThread1%Te_G)
      deallocate(OpenThread1%R_F)
      deallocate(OpenThread1%Coord_DF)
      deallocate(OpenThread1%B_F)
      deallocate(OpenThread1%Ds_G)
    end subroutine deallocate_thread
    !==========================================================================
  end subroutine deallocate_thread_arr
  !============================================================================
  subroutine integrate_emission(TeSi, PeSi, iTable, nVar, Integral_V)
    use ModLookupTable,  ONLY: interpolate_lookup_table
    ! INPUTS:
    ! The plasma parameters on top of the transition region:
    real,    intent(in)  :: TeSi, PeSi
    integer, intent(in)  :: iTable, nVar
    real,    intent(out) :: Integral_V(nVar)
    ! The model is parameterized in terms of PAvr=sqrt(Pi*Pe) at Ti=Te
    ! In terms of SqrtZ: PAvr = Pi*SqrtZ = Pe/SqrtZ
    real    :: PAvrSi
    ! 1D Grid across the TR
    ! Number of uniform meshes in the range from TeSiMin to Te on the TR top
    integer, parameter :: nTrGrid = 20
    ! mesh-centered temperature and the spatial length of intervals, in cm!!
    real    :: TeAvrSi_I(nTrGrid + 1), DeltaLCgs_I(nTrGrid + 1)
    ! Gen table values:
    real    :: Value_VI(nVar, nTrGrid +1)
    real    :: DeltaTe      ! Mesh of a temperature
    real    :: LengthPavrSi_I(nTrGrid + 1), Te_I(nTrGrid + 1)
    real    :: DeltaLengthPavrSi_I(nTrGrid + 1)
    integer ::  i, iVar ! Loop variables
    ! Electron density in particles per cm3:
    real    :: NeCgs, NiCgs
    !--------------------------------------------------------------------------
    PAvrSi = PeSi/SqrtZ
    DeltaTe = (TeSi - TeSiMin)/nTrGrid
    Te_I(1) = TeSiMin
    call get_trtable_value(TeSiMin)
    ! First value is now the product of the thread length in meters times
    ! a geometric mean pressure, so that
    LengthPavrSi_I(1) = TrTable_V(LengthPAvrSi_)
    do i = 1, nTrGrid
       Te_I(i +1) = Te_I(i) + DeltaTe
       call get_trtable_value(Te_I(i + 1))
       LengthPavrSi_I(i + 1) = TrTable_V(LengthPAvrSi_)
       DeltaLengthPavrSi_I(i) = LengthPavrSi_I(i + 1) - LengthPavrSi_I(i)
       TeAvrSi_I(i) = (Te_I(i + 1) + Te_I(i))*0.50
    end do

    TeAvrSi_I(nTrGrid + 1) = TeSi
    DeltaLengthPavrSi_I(nTrGrid + 1) = &
         LengthPavrSi_I(1) - LengthPavrSi_I(nTrGrid + 1)

    ! Now, DeltaL_I is the length interval multiplied by PAvrSi.
    ! Get rid of the latter factor and convert from meters to cm:
    DeltaLCgs_I = DeltaLengthPavrSi_I*100.0/PAvrSi

    do i = 1, nTrGrid + 1
       ! Calculate mesh-centered electron density:
       NeCgs = 1.0e-6*PAvrSi*SqrtZ/(cBoltzmann*TeAvrSi_I(i))
       NiCgs = 1.0e-6*PAvrSi/(cBoltzmann*TeAvrSi_I(i)*SqrtZ)
       call interpolate_lookup_table(iTable, TeAvrSi_I(i), NeCgs, &
            Value_VI(:,i), DoExtrapolate=.true.)
       ! Multiply by a scale factor 1e-26 of the table times Ne*Ni*DeltaL
       Value_VI(:,i) = Value_VI(:,i)*NeCgs*NiCgs*1.0e-26*&
            DeltaLCgs_I(i)
    end do

    ! Sum up the contributions from all temperature intervals:
    do iVar = 1, nVar
       Integral_V(iVar) = sum(Value_VI(iVar,:))
    end do

  end subroutine integrate_emission
  !============================================================================
  subroutine advance_thread_expl(iStage, OpenThread1, IsTimeAccurate, &
       RightFace0_V, LeftFace0_V, DtIn)

    integer, intent(in) :: iStage
    type(OpenThread),intent(inout) :: OpenThread1
    logical, intent(in) :: IsTimeAccurate
    real, OPTIONAL, intent(in) :: RightFace0_V(Rho_:Wminor_)
    real, OPTIONAL, intent(out) :: LeftFace0_V(Rho_:Wminor_)
    real, OPTIONAL, intent(in)  :: DtIn

    integer :: nCell ! # of cells = OpenThread1%nCell
    ! _C(ell) - are cell-centered values
    ! _G(host) - same, but with ghost cells
    ! _F(ace) - are face centered
    ! for each cell the left face has the same index as the cell,
    ! the right face has the index by 1 greater
    ! _C arrays have indexes from -nCell to -1, with nPoint-2 being called
    ! nCell, which is the total number of cells
    ! _G arrays have indexes from -nCell-1 to 0
    ! _F arrays have indexes from -nCell to 0
    !
    !
    ! State variables:
    ! in the state VG array:
    ! 1 - density
    ! 2 - velocity
    ! 3 - ion par pressure
    ! 4 - ion pressure = 2/3 PerpPressure + 1/3 ParPressure
    ! 5 - electron par pressure
    ! 6 - electron pressure = 2/3 Perp pressure + 1/3 ParPressure
    ! 7 - "Major wave" - propagating from the Sun outward/representative
    ! 8 - "Minor wave" - propagating toward the Sun/representative
    ! Primitive variables:
    ! 1 - density
    ! 2 - velocity
    ! 3 - PparIon
    ! 4 - PperpIon
    ! 5 - PparElectron
    ! 6 - PperpElectron
    ! 7 - "Major wave" - propagating from the Sun outward/representative
    ! 8 - "Minor wave" - propagating toward the Sun/representative
    ! Conservative variables:
    ! 1 - mass density
    ! 2 - momentum density
    ! 3 - ion energy density (kinetic+parallel internal)
    ! 4 - PperpIon
    ! 5 - Half of PparElectron (internal energy in parallel component)
    ! 6 - PperpElectron
    ! 7 - "Major wave" - propagating from the Sun outward/representative
    ! 8 - "Minor wave" - propagating toward the Sun/representative

    ! Primitive variables (inputs for the limited interpolation procedure)
    ! Include one layer of ghost cells
    real    :: Primitive_VG(Rho_:Wminor_,-OpenThread1%nCell-1:0)
    ! Conservative variables, in the physical cells only (old)
    real, pointer    :: ConservativeOld_VC(:,:)
    ! State vector:
    real, pointer :: State_VG(:,:)
    ! Electron temperature for semi-implicit scheme
    real, pointer :: Te_G(:)
    ! Local time step for a given cell
    real, pointer :: Dt_C(:)
    ! Cell size along the magnetic field line. In the leftmost ghost cell
    ! there is the length of analytical transition
    real, pointer :: Ds_G(:)
    !
    ! Elements of the control volume scheme:
    !

    ! Fluxes for all variables, thorough a given face
    real :: Flux_VF(Rho_:Wminor_,-OpenThread1%nCell:0)
    ! Face area (divided by the flux)
    real    :: FaceArea_F(-OpenThread1%nCell:0)
    ! Control volume (the inverse of)
    real    :: vInv_C(-OpenThread1%nCell:-1)
    !
    ! What we need for the source calculation:
    !
    ! Face area difference divided per control volume, needed to calculate
    ! the contributions from the perp pressures (Xianyu, why is it called like
    ! this?)
    real :: NablaHatB_C(-OpenThread1%nCell:-1)
    ! Projection of the gravity acceleration, onto the field direction
    real :: GravityAcc_C(-OpenThread1%nCell:-1)
    ! Pseudo-divergence of velocity, needed to calculate contributions from
    ! parallel pressures
    real :: DuDs_C(-OpenThread1%nCell:-1)
    ! Interpolated left and right states at the face
    real    :: pLeft_VF(Rho_:Wminor_,-OpenThread1%nCell:0)
    real    :: pRight_VF(Rho_:Wminor_,-OpenThread1%nCell:0)
    ! Differences across the face, to be limited
    real    :: dVarUp_V(Rho_:Wminor_), dVarDown_V(Rho_:Wminor_)
    ! To get the radial component of spherical coordinates
    real, parameter :: cThird = 1.0/3.0, cTwoThird = 2.0/3.0, Beta = 1.50
    ! Staging:
    integer, parameter :: nStage = 2
    real    :: StageCoef
    ! Misc:
    real :: Source_V(Rho_:Wminor_)
    real :: GravitySource = 0.0
    real :: WavePress, WavePperpSource, WavePpar
    real :: PperpSource,PePerpSource,PePar
    real :: Un_F(-OpenThread1%nCell:0), Rho_F(-OpenThread1%nCell:0),&
         B_F(-OpenThread1%nCell:0)
    real :: Cleft_F(-OpenThread1%nCell:0)
    real :: Cright_F(-OpenThread1%nCell:0)
    real :: VaFace_F(-OpenThread1%nCell:0)
    real :: VaCell_C(-OpenThread1%nCell:-1)
    real :: VaDtOverDs_C(-OpenThread1%nCell:-1)
    real :: DvaDs_C(-OpenThread1%nCell:-1)
    real :: DissMajor_C(-OpenThread1%nCell:-1)
    real :: DissMinor_C(-OpenThread1%nCell:-1)
    real :: Heating_C(-OpenThread1%nCell:-1)
    real :: Reflection_C(-OpenThread1%nCell:-1)
    ! Calculated using the constants above or by applying apportion_heating
    real :: PparHeating, PperpHeating, PeHeating
    ! Loop variables:
    integer :: iCell, iFace
    character(len=*), parameter:: NameSub = 'advance_thread_expl'
    !--------------------------------------------------------------------------
    if(present(DtIn).and.(.not.IsTimeAccurate))call CON_stop(&
         NameSub//':DtIn input is only allowed in time accurate mode')
    nCell = OpenThread1%nCell
    Dt_C=>OpenThread1%Dt_C
    State_VG => OpenThread1%State_VG
    ConservativeOld_VC => OpenThread1%ConservativeOld_VC
    Te_G => OpenThread1%Te_G
    ! Mesh size along the line
    Ds_G => OpenThread1%Ds_G
    ! Face centered field, in Si
    B_F = OpenThread1%B_F(-nCell:0)
    ! Face area
    FaceArea_F(-nCell:0) = 1/B_F
    do iCell = -nCell, -1
       ! Control volume (the inverse of): 1/(dS*0.5*(FaceLeft + FaceRight))
       vInv_C(iCell) = 2/(Ds_G(iCell)*&
            (FaceArea_F(iCell) + FaceArea_F(iCell+1)) )
       ! B d(1/B)/ds, in m^{-1}
       NablaHatB_C(iCell) = vInv_C(iCell)* &
            (FaceArea_F(iCell+1) - FaceArea_F(iCell) )
       ! Gravity force projection onto the field direction
       ! Is calculated as the negative of increment in gravity potential
       ! divided by the mesh length.
       GravityAcc_C(iCell) =  cGravityPotAt1Rs* & ! Calc -\delta Grav.Pot.
            (1/OpenThread1%R_F(iCell) -    &
            1/OpenThread1%R_F(iCell+1)) /  &
            OpenThread1%Ds_G(iCell)  ! Divide by \delta s
    end do
    Primitive_VG(:,-nCell:0) = State_VG(:,-nCell:0)
    Primitive_VG(Pperp_,-nCell:0) = Primitive_VG(P_,-nCell:0) + 0.5* &
         (Primitive_VG(P_,-nCell:0) - Primitive_VG(Ppar_,-nCell:0))
    Primitive_VG(PePerp_,-nCell:0) = Primitive_VG(Pe_,-nCell:0) + 0.5* &
         (Primitive_VG(Pe_,-nCell:0) - Primitive_VG(PePar_,-nCell:0))
    if(iStage==1)then
       do iCell = -nCell, -1
          !
          ! Initial values for the conserved variables:
          !
          ConservativeOld_VC(:,iCell) = Primitive_VG(:,iCell)
          ! Except:
          ! 1. Momentum density
          ConservativeOld_VC(RhoU_,iCell) = &
               Primitive_VG(U_,iCell)*Primitive_VG(Rho_,iCell)
          ! 2. Energy density
          ConservativeOld_VC(Energy_,iCell) = 0.50*Primitive_VG(Ppar_,iCell)&
               + 0.50*Primitive_VG(Rho_,iCell)*Primitive_VG(U_,iCell)**2
          ! 3. "Conserved" electron parallel energy density is a half of PePar)
          ConservativeOld_VC(PePar_,iCell) = 0.50*Primitive_VG(PePar_,iCell)
       end do
    end if
    !
    !
    ! Initial distributions are all initialized, start the 2-stage advance
    ! procedure
    !
    !
    StageCoef = 0.50*iStage ! Half time step is applied at first stage
    !
    ! Boundary conditions:
    !
    ! Apply left BC in the ghostcell #=-nCell-1
    call set_thread_inner_bc(&
         TeIn = OpenThread1%TeTr, &
         uIn  = OpenThread1%uTr,  &
         pIn  = OpenThread1%PeTr, &
         Primitive_V= Primitive_VG(:,-nCell-1))
    !
    ! limited interpolation procedure:
    ! 1. logarithm of density/pressure is better to be limited
    if(DoLimitLogVar.and.any(Primitive_VG(iLogVar_V,-nCell-1:0)<=0.0))then
       write(*,*)'-ncell-1=',-nCell - 1
       write(*,*)'OpenThread TeTr=', OpenThread1%TeTr
       write(*,*)'OpenThread uTr=', OpenThread1%uTr
       write(*,*)'OpenThread PeTr=', OpenThread1%PeTr
       if(OpenThread1%PeTr <= 0.0)write(*,*)'TrTable_V=', TrTable_V
       do iCell = -nCell-1, 0
          if(any(Primitive_VG(iLogVar_V,iCell)<=0.0))&
               write(*,*)'iCell=',iCell,' state:',&
               Primitive_VG(:,iCell)
       end do
       call save_plot_thread(OpenThread1,'failed_thread.out')
       call CON_stop('Negative pressure/density')
    end if
    if(DoLimitLogVar)Primitive_VG(iLogVar_V,-nCell-1:0) = log(&
         Primitive_VG(iLogVar_V,-nCell-1:0))
    ! 2. Left boundary:
    pLeft_VF(:,-nCell) = Primitive_VG(:,-nCell-1)
    ! 3. Calculate leftmost unlimited slope. The half of this (as assumed in
    ! the limiter function) approximates the face value
    dVarDown_V = 2*(Primitive_VG(:,-nCell) - Primitive_VG(:,-nCell-1))
    ! Calculate the up slope
    dVarUp_V = Primitive_VG(:,-nCell+1) - Primitive_VG(:,-nCell)
    ! Calculate and apply the limited slopes, to get the limited
    ! reconstructed values:
    pRight_VF(:,-nCell) = Primitive_VG(:,-nCell) - &
         (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
         min(Beta*abs(dVarUp_V), abs(dVarDown_V), &
         cThird*abs(2*dVarDown_V+dVarUp_V))
    pLeft_VF(:,-nCell+1) = Primitive_VG(:,-nCell) + &
         (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
         min(Beta*abs(dVarUp_V), abs(dVarDown_V), &
         cThird*abs(dVarDown_V+2*dVarUp_V))
    ! 4. Pass through all cells, limit slopes
    do iCell = -nCell+1, -2
       ! Propagate the old up slope to become the down slope
       dVarDown_V = dVarUp_V
       ! Calculate the up slope
       dVarUp_V = Primitive_VG(:,iCell+1) - Primitive_VG(:,iCell)
       ! Calculate and apply the limited slopes, to get the limited
       ! reconstructed values:
       pRight_VF(:,iCell) = Primitive_VG(:,iCell) - &
            (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
            min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
            cThird*abs(2*dVarDown_V+dVarUp_V))
       pLeft_VF(:,iCell+1) = Primitive_VG(:,iCell) + &
            (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
            min(Beta*abs(dVarUp_V), Beta*abs(dVarDown_V), &
            cThird*abs(dVarDown_V+2*dVarUp_V))
    end do
    ! Propagate the old up slope to become the down slope
    dVarDown_V = dVarUp_V
    ! Calculate the up slope. The half of this (as assumed in the
    ! limiter function) approximates the face value
    dVarUp_V = Ds_G(-1)/(0.50*Ds_G(-1) + Ds_G(0))*&
         (Primitive_VG(:,0) - Primitive_VG(:,-1))
    ! Calculate and apply the limited slopes, to get the limited
    ! reconstructed values:
    pRight_VF(:,-1) = Primitive_VG(:,-1) - &
         (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
         min(abs(dVarUp_V), Beta*abs(dVarDown_V), &
         cThird*abs(2*dVarDown_V+dVarUp_V))
    pLeft_VF(:,0) = Primitive_VG(:,-1) + &
         (sign(0.25, dVarUp_V) + sign(0.25, dVarDown_V))* &
         min(abs(dVarUp_V),Beta*abs(dVarDown_V), &
         cThird*abs(dVarDown_V+2*dVarUp_V))
    if(DoLimitLogVar)then
       Primitive_VG(iLogVar_V,-nCell-1:0) = exp(&
            Primitive_VG(iLogVar_V,-nCell-1:0))
       pRight_VF(iLogVar_V,-nCell:-1) = exp(&
            pRight_VF(iLogVar_V,-nCell:-1))
       pLeft_VF(iLogVar_V,-nCell:0) = exp(&
            pLeft_VF(iLogVar_V,-nCell:0))
    end if
    ! 2. Apply right BC on the external boundary (# = 0)
    if(present(RightFace0_V))then
       pRight_VF(:,0) = RightFace0_V
    else
       ! First order BC
       pRight_VF(:,0) = Primitive_VG(:,0)
    end if
    ! 2. Save left BC on the external boundary (# = 0)
    if(present(LeftFace0_V))LeftFace0_V =  pLeft_VF(:,0)
    ! Get fluxes
    ! Loop over faces:
    do iFace = -nCell, 0
       call get_thread_flux(pLeft_VF(:,iFace), pRight_VF(:,iFace), &
            Flux_VF(:,iFace), Cleft_F(iFace), Cright_F(iFace), &
            Un_F(iFace), Rho_F(iFace))
       VaFace_F(iFace) = B_F(iFace)/sqrt(cMu*Rho_F(iFace))
    end do
    ! Set the speed on top of the transition region from the RS:
    ! Fix 2025-3-13
    OpenThread1%uTr = 0.0 ! Un_F(-nCell)*Rho_F(-nCell)/(&
         ! OpenThread1%PeTr*cProtonMass/(OpenThread1%TeTr*cBoltzmann))
    call get_trtable_value(OpenThread1%TeTr, OpenThread1%uTr)
    ! Correct pressure for updated plasma speed
    OpenThread1%PeTr = TrTable_V(LengthPavrSi_)/Ds_G(-nCell-1)
    ! du/ds, dV_A/ds, V_A
    do iCell = -nCell,-1
       DuDs_C(iCell)  = (Un_F(iCell+1) - Un_F(iCell)) / &
            Ds_G(iCell)
       DvaDs_C(iCell) = (VaFace_F(iCell+1) - VaFace_F(iCell)) / &
            Ds_G(iCell)
       VaCell_C(iCell) = 0.50*(VaFace_F(iCell+1) + VaFace_F(iCell))
       DissMajor_C(iCell) = 2.0*sqrt(PoyntingFluxPerBsi*cMu*&
            Primitive_VG(Wminor_,iCell)*VaCell_C(iCell))/LperpTimesSqrtBsi
       DissMinor_C(iCell) = 2.0*sqrt(PoyntingFluxPerBsi*cMu*&
            Primitive_VG(Wmajor_,iCell)*VaCell_C(iCell))/LperpTimesSqrtBsi
       Reflection_C(iCell) = sign(min(abs(DvaDs_C(iCell)),0.50*abs(&
            DissMajor_C(iCell) - DissMinor_C(iCell))),&
            Primitive_VG(Wmajor_,iCell) - Primitive_VG(Wminor_,iCell))
       Heating_C(iCell) = &
            (DissMajor_C(iCell)*Primitive_VG(Wmajor_,iCell)&
            + DissMinor_C(iCell)*Primitive_VG(Wminor_,iCell))*&
            PoyntingFluxPerBsi*sqrt(cMu*Primitive_VG(Rho_,iCell))
    end do
    where(OpenThread1%R_F(-nCell:-1)<rMinReflectionTr)&
         Reflection_C = 0.0
    if(iStage==1)then
       do iCell = -nCell, -1
          ! Get dt
          ! Dt_C(iCell) = CflLocal * OpenThread1%Ds_G / Cmax_I
          ! Cmax is the maximum of the right wave
          ! speed from the left face and the left wave speed
          ! from the right face
          Dt_C(iCell) = CflLocal/max(&
               vInv_C(iCell)* &
               abs(max(Cright_F(iCell)*FaceArea_F(iCell),&
               -Cleft_F(iCell+1)*FaceArea_F(iCell+1)))&
               + &  ! CFL
               max(DissMajor_C(iCell),DissMinor_C(iCell))&
               ,& ! + diss rate
               cooling_rate(RhoSi=Primitive_VG(Rho_,iCell),&  ! or cooling
               PeSi=Primitive_VG(PePar_,iCell)))              ! rate
       end do
       if(present(DtIn))then
          Dt_C(-nCell:-1) = min(Dt_C(-nCell:-1), DtIn)
          OpenThread1%Dt = DtIn
       else
          OpenThread1%Dt = minval(Dt_C(-nCell:-1))
          if(IsTimeAccurate)Dt_C(-nCell:-1) = OpenThread1%Dt
       end if
    end if
    !
    ! Advance the conserved variables to the half or full time step,
    ! depending on the stage
    !
    do iCell = -nCell, -1
       ! Apply fluxes

       ! dU = dt / Volume * (F_{i-1/2}*FaceArea_{i-1/2} &
       !                        -F_{i+1/2}*FaceArea_{i+1/2} )
       State_VG(:,iCell) = ConservativeOld_VC(:,iCell) + &
            StageCoef*Dt_C(iCell)*vInv_C(iCell)* &
            (Flux_VF(:,iCell)*FaceArea_F(iCell) - &
            Flux_VF(:,iCell+1)*FaceArea_F(iCell+1) )
       ! Calculate source terms
       ! 1. Pperp source term
       PperpSource = &
            Primitive_VG(Pperp_,iCell)*NablaHatB_C(iCell)
       ! 2. Wave pressure, convert from the representative functions
       WavePress = 0.50*&
            sum(Primitive_VG(Wmajor_:Wminor_,iCell))*PoyntingFluxPerBsi*&
            sqrt(cMu*Primitive_VG(Rho_,iCell))
       ! 3. Wave pressure source terms
       WavePperpSource = WavePress*NablaHatB_C(iCell)
       WavePpar = WavePress
       ! 4. PePar source term
       PePar = Primitive_VG(PePar_,iCell)
       ! 5. PePerp source term
       PePerpSource = &
            Primitive_VG(PePerp_,iCell)*NablaHatB_C(iCell)
       ! 6. Heating and partitioning
       if(UseStochasticHeating)then
          call apportion_heating(&
               ! Inputs, all in SI:
               PparIn = Primitive_VG(Ppar_,iCell)      ,&
               PperpIn= Primitive_VG(Pperp_,iCell)     ,&
               PeIn   = cThird*Primitive_VG(PePar_,iCell) + &
               cTwoThird*Primitive_VG(PePerp_,iCell)   ,&
               RhoIn  = Primitive_VG(Rho_,iCell)       ,&
               BIn    = 0.50*(B_F(iCell) + B_F(iCell+1)),&
               WmajorIn = Primitive_VG(Wmajor_,iCell)  ,&
               WminorIn = Primitive_VG(Wminor_,iCell)  ,&
               DissRateMajorIn = DissMajor_C(iCell)     ,&
               DissRateMinorIn = DissMinor_C(iCell)     ,&
               ! Outputs:
               QparPerQtotal     = QparPerQtotal        ,&
               QperpPerQtotal    = QperpPerQtotal       ,&
               QePerQtotal       = QePerQtotal)
       end if
       PparHeating  = QparPerQtotal *Heating_C(iCell)
       PperpHeating = QperpPerQtotal*Heating_C(iCell)
       PeHeating    = QePerQtotal   *Heating_C(iCell)
       ! 7. Gravity force density, rho*g:
       GravitySource = Primitive_VG(Rho_,iCell)*GravityAcc_C(iCell)
       !
       ! Get full source vector
       !
       Source_V     = [0.0, &  ! No density source
            PperpSource + PePerpSource + WavePperpSource +&
            GravitySource,                                & ! Momentum
            (PperpSource + PePerpSource + WavePperpSource + GravitySource)&
            *Primitive_VG(U_,iCell) + &
            (PePar + WavePpar)*DuDs_C(iCell) + &
            PparHeating,                  & ! Energy
            - PperpSource*Primitive_VG(U_,iCell) + PperpHeating, & ! Pperp
            - PePar*DuDs_C(iCell) + cThird*PeHeating, & ! <-PePar / PePerp
            - PePerpSource*Primitive_VG(U_,iCell)+ cTwoThird*PeHeating, &
            -DissMajor_C(iCell)*Primitive_VG(Wmajor_,iCell) - &
            Reflection_C(iCell)*sqrt(Primitive_VG(Wmajor_,iCell)*&
            Primitive_VG(Wminor_,iCell)),               & ! Wmajor
            -DissMinor_C(iCell)*Primitive_VG(Wminor_,iCell) + &
            Reflection_C(iCell)*sqrt(Primitive_VG(Wmajor_,iCell)*&
            Primitive_VG(Wminor_,iCell))                 ] ! Wminor
       ! Apply source
       State_VG(:,iCell) = State_VG(:,iCell) +&
            StageCoef*Dt_C(iCell)*Source_V
    end do
    !
    ! Finalize the stage
    !
    ! Transform to primitives, set minimum pressure and density
    do iCell = -nCell, -1
       State_VG(Rho_,iCell) = max(MinRho, State_VG(Rho_,iCell) )
       State_VG(U_,iCell)   = State_VG(RhoU_,iCell)/State_VG(Rho_,iCell)

       ! Ppar = 2*Energy - Rho U**2
       State_VG(Ppar_,iCell)   = max(MinPress, 2*State_VG(Energy_,iCell) - &
            State_VG(Rho_,iCell)*State_VG(U_,iCell)**2 )

       ! Pressure = 2/3 Pperp +1/3 Ppar
       State_VG(P_,iCell)   = max(MinPress, cThird*State_VG(Ppar_,iCell)&
            + cTwoThird*State_VG(Pperp_,iCell))

       State_VG(PePar_,iCell)   = max(MinPress, 2*State_VG(PePar_,iCell) )

       State_VG(Pe_,iCell)  = max(MinPress, cThird*State_VG(PePar_,iCell)&
            + cTwoThird*State_VG(PePerp_,iCell))
       Te_G(iCell) = State_VG(Pe_,iCell)*cProtonMass/&
         (cBoltzmann*State_VG(Rho_,iCell))
    end do
    if(iStage==1)RETURN
    ! Semi-implicit stage
    ! First, solve implicit equation
    ! Wmajor_i^{n+1} -Wmajor_i^n +(Va*Dt/Ds)*(-W^{n+1}_{i-1} + W^{n+1)_i)
    Primitive_VG(Wmajor_:Wminor_,-nCell:-1) = &
         State_VG(Wmajor_:Wminor_,-nCell:-1)
    do iCell = -nCell, -1
       VaDtOverDs_C(iCell) = VaCell_C(iCell)*Dt_C(iCell)/Ds_G(iCell)
       Primitive_VG(Wmajor_,iCell) = ( Primitive_VG(Wmajor_,iCell) +&
            Primitive_VG(Wmajor_,iCell-1)*VaDtOverDs_C(iCell)) / &
            (1 + VaDtOverDs_C(iCell))
    end do
    ! Second, solve implicit equation
    ! Wminor_i^{n+1} -Wminor_i^n -(Va*Dt/Ds)*(+W^{n+1}_{i+1} - W^{n+1)_i) = 0
    do iCell = -1, -nCell, -1
       Primitive_VG(Wminor_,iCell) = ( Primitive_VG(Wminor_,iCell) +&
            Primitive_VG(Wminor_,iCell+1)*VaDtOverDs_C(iCell)) / &
            (1 + VaDtOverDs_C(iCell))
    end do
    State_VG(Wmajor_:Wminor_,-nCell:-1) = &
         Primitive_VG(Wmajor_:Wminor_,-nCell:-1)
  contains
    !==========================================================================
    subroutine set_thread_inner_bc(TeIn, uIn, pIn, Primitive_V)

      real, intent(in) :: TeIn, uIn, pIn
      real, intent(out):: Primitive_V(Rho_:Wminor_)
      !------------------------------------------------------------------------
      Primitive_V = [&
           cProtonMass*pIn/(cBoltzmann*TeIn), &
           uIn, &
           pIn, & ! Ppar
           pIn, & ! Pperp
           pIn, & ! PePar
           pIn, & ! PePerp
           1.0, & ! Wmajor
           1e-8]  ! Wminor
    end subroutine set_thread_inner_bc
    !==========================================================================
    subroutine get_thread_flux(pLeft_V, pRight_V, &
         Flux_V, Cleft, Cright, UnFace, RhoFace)

      real,   intent(in) :: pLeft_V(Rho_:Wminor_), pRight_V(Rho_:Wminor_)
      real,  intent(out) :: Flux_V(Rho_:Wminor_)
      real,  intent(out) :: Cleft, Cright, UnFace, RhoFace

      real    :: UpwindState_V(Rho_:Wminor_), &
           TotalP, DensityRatio, PparFace, ExtraP, ExtraPright, ExtraPleft

      !------------------------------------------------------------------------
      RhoL = pLeft_V(Rho_ )
      ExtraPleft = pLeft_V(PePar_) + sum(pLeft_V(Wmajor_:Wminor_))*0.5*&
           sqrt(cMu*RhoL)*PoyntingFluxPerBsi
      pL   = pLeft_V(Ppar_) + ExtraPleft
      UnL  = pLeft_V(U_)

      RhoR = pRight_V(Rho_ )
      ExtraPright = pRight_V(PePar_) + sum(pRight_V(Wmajor_:Wminor_))*0.5*&
           sqrt(cMu*RhoR)*PoyntingFluxPerBsi
      pR   = pRight_V(Ppar_) + ExtraPright
      UnR  = pRight_V(U_)

      ! exact solver

      call exact_rs_pu_star
      call exact_rs_sample(0.0, RhoFace, UnFace, TotalP)

      ! Store WL WR

      Cleft  = WL
      Cright = WR

      ! Upwind scheme to get pPerpFace

      if (UnStar  >  0) then
         ! Factor for variables linearly scaling with density
         DensityRatio = RhoFace / pLeft_V(Rho_)
         UpwindState_V = pLeft_V
         ExtraP = ExtraPleft
      else
         DensityRatio = RhoFace / pRight_V(Rho_)
         UpwindState_V = pRight_V
         ExtraP = ExtraPright
      end if

      ! get the flux

      Flux_V(Rho_) = RhoFace * UnFace
      Flux_V(RhoU_) = RhoFace*UnFace**2 + TotalP
      ! All contributions to pressure except for ion one are assumed
      ! to scale adiabatically with density:
      PparFace = TotalP - ExtraP*DensityRatio**3
      Flux_V(Energy_) = (0.50*(RhoFace*UnFace**2 + PparFace) + TotalP)*UnFace
      ! All passive scalars should be upwinded, scaled proportional to density
      ! and multiplied by Uface, to get their flux:
      Flux_V(Pperp_:Wminor_) = DensityRatio*UpwindState_V(Pperp_:Wminor_)&
           *UnFace
      ! Exception: both the conservative variable is 1/2 PePar and its flux
      ! while the left and right state vertors include just PePar, so that
      ! this flux should be multiplied by 1/2
      Flux_V(PePar_) = 0.50*Flux_V(PePar_)
    end subroutine get_thread_flux
    !==========================================================================
  end subroutine advance_thread_expl
  !============================================================================
  subroutine advance_thread_semi_impl(OpenThread1)

    type(OpenThread),intent(inout) :: OpenThread1
    ! Solver for electron heat conduction
    real :: Ti_C(-OpenThread1%nCell:-1)
    real :: N_C(-OpenThread1%nCell:-1)
        ! State vector:
    real, pointer :: State_VG(:,:)
    ! Electron temperature for semi-implicit scheme
    real, pointer :: Te_G(:)
    ! Local time step for a given cell
    real, pointer :: Dt_C(:)
    ! Cell size along the magnetic field line. In the leftmost ghost cell
    ! there is the length of analytical transition
    real, pointer :: Ds_G(:)
    integer :: nCell
    character(len=*), parameter:: NameSub = 'advance_thread_semi_impl'
    !--------------------------------------------------------------------------
    nCell = OpenThread1%nCell
    Dt_C=>OpenThread1%Dt_C
    State_VG => OpenThread1%State_VG
    Te_G => OpenThread1%Te_G
    ! Mesh size along the line
    Ds_G => OpenThread1%Ds_G
    ! Electron heat conduction and losses
    N_C(-nCell:-1) = State_VG(Rho_,-nCell:-1)/cProtonMass
    Ti_C(-nCell:-1) = State_VG(P_,-nCell:-1)/(cBoltzmann*N_C)
    call advance_heat_conduction(&
         nPoint = nCell, &
         Dt_I = Dt_C,    &
         Te_I = Te_G,    &
         Ti_I = Ti_C,    &
         Ni_I = N_C,     &
         Ds_I = Ds_G,    &
         uFace   = OpenThread1%uTr, &
         B_I  = OpenThread1%B_F(-nCell:0), &
         PeFaceOut = OpenThread1%PeTr, & ! Store state for the top of TR
         TeFaceOut = OpenThread1%TeTr, &
         DoLimitTimestep=.true.)
    State_VG(P_ ,-nCell:-1) = cBoltzmann*N_C*Ti_C(-nCell:-1)
    State_VG(Pe_,-nCell:-1) = cBoltzmann*N_C*Te_G(-nCell:-1)
    ! If no anysotropic pressure is used
    State_VG(Ppar_ ,-nCell:-1) = State_VG(P_ ,-nCell:-1)
    State_VG(PePar_ ,-nCell:-1) = State_VG(Pe_ ,-nCell:-1)
  end subroutine advance_thread_semi_impl
  !============================================================================
  subroutine solve_tr_face(TeCell, & ! cell-centered, in SI
       uFace,                      & ! face-centered
       LengthTr, DeltaS, & ! Distances photosphere-to-face, face-to-center
       TeFaceOut, HeatFluxOut, DFluxOverDConsOut, &
       PeFaceOut) ! SI, face centered

    ! For given values of Te and u*Pavr/Te at the cell center, all inputs
    ! are in SI and expressed in K and W/(m.K). The distance from the cell
    ! to the face is DeltaS, from the face to photosphere is LengthTr.
    ! All distances are measured along the magnetic field line in meters.
    ! Outputs provided are TeFace, PeFace, uFace, HeatFluxFace, all in SI:
    ! [K], [N/m^2], [m/s], [W/m^2], at the face.

    ! The plasma temperature at the cell center:
    real,    intent(in)    :: TeCell
    ! Velosity at the face
    real,    intent(inout) :: uFace
    ! Distances photosphere-to-face, face-to-center
    real,    intent(in)  :: LengthTr, DeltaS
    ! Parameters at the face, connected by the analytical TR solution to the
    ! photosphere:
    real, optional, intent(out) :: TeFaceOut
    real, optional, intent(out) :: HeatFluxOut, dFluxOverdConsOut
    real, optional, intent(out) :: PeFaceOut

    real    :: TeFace, HeatFluxFace
    real    :: ConsFace, dConsFace, ConsCell, Tolerance
    ! Newton-Rapson coefficients: dConsFace = -Res/dResdCons
    real    :: Res, dResdCons
    integer :: iCounter
    integer, parameter :: iCounterMax=20
    character(len=*), parameter:: NameSub = 'solve_tr_face'
    !--------------------------------------------------------------------------
    ConsCell = cTwoSevenths*HeatCondParSi*TeCell**3.5
    Tolerance = cTolerance*ConsCell
    ! In the first approximation, take TeFace = TeCell
    ConsFace = ConsCell
    TeFace = TeCell
    iCounter = 0
    do
       call get_trtable_value(TeFace, uFace)
       HeatFluxFace = TrTable_V(HeatFluxLength_)/LengthTr
       ! Solve equation HeatFlux = (ConsCell - ConsFace)/DeltaS
       Res = DeltaS*HeatFluxFace + ConsFace - ConsCell
       if(abs(Res) < Tolerance) EXIT
       if(iCounter == iCounterMax)then
          write(*,*)'ChromoEvapCoef=',ChromoEvapCoef
          write(*,*)'Lengthtr',LengthTr
          write(*,*)'TeCell=', TeCell
          write(*,*)'uFace=',uFace
          write(*,*)'DeltaS=',DeltaS
          call CON_stop(NameSub//': No convergence')
       end if
       ! New iteration:
       iCounter = iCounter + 1
       dResdCons = DeltaS*TrTable_V(dHeatFluxXOverU_)/LengthTr + 1
       dConsFace = -Res/dResdCons
       ConsFace = ConsFace + dConsFace
       TeFace = (ConsFace*3.50/HeatCondParSi)**cTwoSevenths
    end do
    if(present(TeFaceOut))TeFaceOut = TeFace
    if(present(HeatFluxOut))HeatFluxOut = HeatFluxFace
    if(present(dFluxOverdConsOut))then
       dFluxOverdConsOut = TrTable_V(dHeatFluxXOverU_)/LengthTr
       dFluxOverdConsOut = dFluxOverdConsOut/(1 + DeltaS*dFluxOverdConsOut)
    end if
    if(present(PeFaceOut))PeFaceOut = TrTable_V(LengthPavrSi_)*SqrtZ/LengthTr
  end subroutine solve_tr_face
  !============================================================================
  real function cooling_rate(RhoSi, PeSi)
    ! Volumetric radiation energy loss divided by the electron internal energy

    real, intent(in) :: RhoSi, PeSi
    real :: TeSi, NiSi, Cooling

    !--------------------------------------------------------------------------
    NiSi = RhoSi/cProtonMass
    TeSi = PeSi/(cBoltzmann*Z*NiSi)
    call get_trtable_value(TeSi)
    Cooling =TrTable_V(LambdaSI_)*Z*(cBoltzmann*NiSi)**2
    cooling_rate =  (3.50 + max(0.0, -TrTable_V(DLogLambdaOverDLogT_)))*&
               Cooling/(1.5*PeSi)
  end function cooling_rate
  !============================================================================
  subroutine advance_heat_conduction_ta(nPoint, Dt, Te_I, Ti_I, Ni_I, Ds_I, &
       uFace, B_I, BCellIn_I, PeFaceOut, TeFaceOut, DoLimitTimestep)
    integer, intent(in)    :: nPoint
    real,    intent(in)    :: Dt
    real,    intent(inout) :: Te_I(nPoint+1), Ti_I(nPoint)
    real,    intent(in)    :: Ni_I(nPoint), Ds_I(0:nPoint+1), B_I(nPoint+1)
    real,    intent(inout) :: uFace
    real, optional, intent(in) :: BCellIn_I(nPoint)
    real, optional, intent(out) :: PeFaceOut, TeFaceOut
    logical, optional, intent(in) :: DoLimitTimestep
    ! MISC:
    real :: Dt_I(nPoint)
    !--------------------------------------------------------------------------
    Dt_I = Dt
    call advance_heat_conduction_ss(nPoint, Dt_I, Te_I, Ti_I, Ni_I, Ds_I, &
       uFace, B_I, BCellIn_I, PeFaceOut, TeFaceOut, DoLimitTimestep)
  end subroutine advance_heat_conduction_ta
  !============================================================================
  subroutine advance_heat_conduction_ss(nPoint, Dt_I, Te_I, Ti_I, Ni_I, Ds_I,&
       uFace, B_I, BCellIn_I, PeFaceOut, TeFaceOut, DoLimitTimestep)
    use ModConst, ONLY: cBoltzmann
    integer, intent(in)    :: nPoint
    real,    intent(in)    :: Dt_I(nPoint)
    real,    intent(inout) :: Te_I(nPoint+1), Ti_I(nPoint)
    real,    intent(in)    :: Ni_I(nPoint), Ds_I(0:nPoint+1), B_I(nPoint+1)
    real,    intent(inout) :: uFace
    real, optional, intent(in) :: BCellIn_I(nPoint)
    real, optional, intent(out) :: PeFaceOut, TeFaceOut
    logical, optional, intent(in) :: DoLimitTimestep

    real :: Res_VI(2,nPoint), Weight_VI(2,nPoint), bFaceInv_I(1:nPoint+1),  &
         Lower_VVI(2,2,nPoint), Main_VVI(2,2,nPoint), Upper_VVI(2,2,nPoint),&
         SpecHeat_I(nPoint), ExchangeRate_I(nPoint), Cons_I(nPoint+1),      &
         TeStart_I(nPoint), TiStart_I(nPoint), DeltaEnergy_I(nPoint),       &
         DeltaIonEnergy_I(nPoint), dCons_VI(2,nPoint), DtInv_I(nPoint),     &
         HeatFlux2Tr, dFluxOverdCons, Cooling, BcellInv_I(nPoint)
    real, parameter :: QuasiCfl = 0.85
    integer, parameter :: Cons_ = 1, Ti_=2
    ! Misc:
    ! Loop variable
    integer :: iPoint, iIter
    integer :: nIterMax = 40
    !--------------------------------------------------------------------------
    BfaceInv_I = 1/B_I
    if(present(BcellIn_I))then
       BcellInv_I = 1/BcellIn_I
    else
       BcellInv_I = 0.50*(BFaceInv_I(1:nPoint) + BFaceInv_I(2:nPoint+1))
    end if
    ! Initialization
    TeStart_I(1:nPoint) = Te_I(1:nPoint)
    TiStart_I(1:nPoint) = Ti_I(1:nPoint)
    ! dCons = kappa(Te)dTe=> Cons = 2/7 kappa*Te
    Cons_I(1:nPoint+1) = cTwoSevenths*HeatCondParSi*Te_I(1:nPoint+1)**3.5

    SpecHeat_I     = 1.50*cBoltzmann*Ni_I*Ds_I(1:nPoint)*BcellInv_I
    do iIter = 1,nIterMax
       Main_VVI = 0.0; Upper_VVI = 0.0; Lower_VVI = 0.0
       ! Contribution from heat conduction fluxes
       ! Flux linearizations over small dCons
       Upper_VVI(Cons_,Cons_,nPoint) = &
            -BFaceInv_I(nPoint+1)/(0.50*Ds_I(nPoint) + Ds_I(nPoint+1))
       Upper_VVI(Cons_,Cons_,1:nPoint-1) = &
            -BFaceInv_I(2:nPoint)*2/(Ds_I(1:nPoint-1) + Ds_I(2:nPoint))
       Lower_VVI(Cons_,Cons_,2:nPoint) = Upper_VVI(Cons_,Cons_,1:nPoint-1)
       Main_VVI(Cons_,Cons_,2:nPoint) = &
            - Upper_VVI(Cons_,Cons_,2:nPoint) &
            - Lower_VVI(Cons_,Cons_,2:nPoint)
       ! Right heat fluxes
       Res_VI(Cons_,1:nPoint) = &
            (Cons_I(1:nPoint) - Cons_I(2:nPoint+1)) &
            *Upper_VVI(Cons_,Cons_,1:nPoint)
       ! Add other left heat fluxes
       Res_VI(Cons_,2:nPoint) = Res_VI(Cons_,2:nPoint) + &
           (Cons_I(2:nPoint) - Cons_I(1:nPoint-1))*&
           Lower_VVI(Cons_,Cons_,2:nPoint)
       call solve_tr_face(TeCell =        Te_I(1), &
            uFace = uFace,                         &
            LengthTr      = Ds_I(0),               &
            DeltaS        = 0.50*Ds_I(1),          &
            HeatFluxOut   = HeatFlux2Tr,           &
            DFluxOverDConsOut = dFluxOverdCons,    &
            PeFaceOut = PeFaceOut,                 &
            TeFaceOut = TeFaceOut)
       ! Add left heat flux to the TR
       Res_VI(Cons_,1) = Res_VI(Cons_,1) - HeatFlux2Tr*BFaceInv_I(1)
       ! Linearize left heat flux to the TR
       Main_VVI(Cons_,Cons_,1) = &
            - Upper_VVI(Cons_,Cons_,1) + dFluxOverdCons*BFaceInv_I(1)
       ! Radiative cooling, limit timestep
       do iPoint = 1, nPoint
          call get_trtable_value(Te_I(iPoint))
          Cooling =Ds_I(iPoint)*BcellInv_I(iPoint)*TrTable_V(LambdaSI_)*Z*&
               (cBoltzmann*Ni_I(iPoint))**2
          Res_VI(Cons_,iPoint) = Res_VI(Cons_,iPoint) - Cooling
          ! linearized -dCooling/dCons
          Main_VVI(Cons_,Cons_,iPoint) = Main_VVI(Cons_,Cons_,iPoint) +&
               TrTable_V(DLogLambdaOverDLogT_)*Cooling/(3.5*Cons_I(iPoint))
          ! Limit time step:
          if(iIter==1)then
             if(present(DoLimitTimestep))then
                DtInv_I(iPoint) = max(1/Dt_I(iPoint),   &
                     (3.50 + max(0.0, -TrTable_V(DLogLambdaOverDLogT_)))&
                     *Cooling/(QuasiCfl*Z*Te_I(iPoint)*SpecHeat_I(iPoint)))
             else
                DtInv_I(iPoint) = 1/Dt_I(iPoint)
             end if
          end if
       end do
       ! Change in the internal energy (to correct the energy source
       ! for the time-accurate mode):
       DeltaEnergy_I(1:nPoint) = Z*SpecHeat_I(1:nPoint)*DtInv_I* &
            (Te_I(1:nPoint) - TeStart_I(1:nPoint))
       DeltaIonEnergy_I(1:nPoint) = SpecHeat_I(1:nPoint)*DtInv_I* &
            (Ti_I(1:nPoint) - TiStart_I(1:nPoint))
       ! Energy evolution:
       Main_VVI(Cons_,Cons_,1:nPoint) = Main_VVI(Cons_,Cons_,1:nPoint) + &
              DtInv_I*Z*SpecHeat_I(1:nPoint)*Te_I(1:nPoint)/   &
              (3.5*Cons_I(1:nPoint))
       Main_VVI(Ti_,Ti_,1:nPoint) = Main_VVI(Ti_,Ti_,1:nPoint) + &
            DtInv_I*SpecHeat_I(1:nPoint)
       Res_VI(Cons_,1:nPoint) = Res_VI(Cons_,1:nPoint) - &
            DeltaEnergy_I(1:nPoint)
       Res_VI(Ti_,1:nPoint) =  - DeltaIonEnergy_I(1:nPoint)
       ExchangeRate_I = cExchangeRateSi*Z**2*Ni_I/Te_I(1:nPoint)**1.5
       ! Point implicit limiter:
       ExchangeRate_I = Z*SpecHeat_I*ExchangeRate_I/&
            (1 + 2*ExchangeRate_I/DtInv_I)
       Res_VI(Cons_,1:nPoint) = Res_VI(Cons_,1:nPoint) + &
            ExchangeRate_I(1:nPoint)*&
            (Ti_I(1:nPoint) - Te_I(1:nPoint))
       Main_VVI(Cons_,Ti_,1:nPoint) = Main_VVI(Cons_,Ti_,1:nPoint) -&
            ExchangeRate_I(1:nPoint)
       Main_VVI(Cons_,Cons_,1:nPoint) = Main_VVI(Cons_,Cons_,1:nPoint) +&
            ExchangeRate_I(1:nPoint)*Te_I(1:nPoint)/&
            (3.5*Cons_I(1:nPoint))
       Res_VI(Ti_,1:nPoint) = Res_VI(Ti_,1:nPoint) &
            + ExchangeRate_I(1:nPoint)*&
            (Te_I(1:nPoint) - Ti_I(1:nPoint))
       Main_VVI(Ti_,Ti_,1:nPoint) = Main_VVI(Ti_,Ti_,1:nPoint) +&
            ExchangeRate_I(1:nPoint)
       Main_VVI(Ti_,Cons_,1:nPoint) = Main_VVI(Ti_,Cons_,1:nPoint) -&
            ExchangeRate_I(1:nPoint)*Te_I(1:nPoint)/&
            (3.5*Cons_I(1:nPoint))
       call tridiag_block22(n=nPoint,  &
            Lower_VVI=Lower_VVI(:,:,1:nPoint),&
            Main_VVI=Main_VVI(:,:,1:nPoint),&
            Upper_VVI=Upper_VVI(:,:,1:nPoint),&
            Res_VI=Res_VI(:,1:nPoint),  &
            Weight_VI=DCons_VI(:,1:nPoint))
       Cons_I(1:nPoint) = Cons_I(1:nPoint) + DCons_VI(Cons_,1:nPoint)
       if(any(Cons_I(1:nPoint)<=0))then
          do iPoint = 1, nPoint
             write(*,*)'iPoint Cons_VI=',iPoint,Cons_I(iPoint), &
                  Ti_I(iPoint) + DCons_VI(Ti_,iPoint)
          end do
          call CON_stop('Reduce time step')
       end if
       ! Recover temperature
       Te_I(1:nPoint) = max((3.5*Cons_I(1:nPoint)/HeatCondParSi)**cTwoSevenths&
            , TeSiMin)
       Ti_I(1:nPoint) = Ti_I(1:nPoint) + DCons_VI(Ti_,1:nPoint)
       if(all(abs(DCons_VI(Cons_,1:nPoint))<cTolerance*Cons_I(1:nPoint)))EXIT
    end do
  end subroutine advance_heat_conduction_ss
  !============================================================================
  subroutine tridiag_block22(n,Lower_VVI,Main_VVI,Upper_VVI,Res_VI,Weight_VI)
    ! This routine solves three-diagonal system of equations:
    !
    !  ||m_1 u_1  0....        || ||w_1|| ||r_1||
    !  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||
    !  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||
    !  ||...                   || ||...|| ||...||
    !  ||.............0 l_n m_n|| ||w_n|| ||r_n||
    !
    ! Prototype: Numerical Recipes, Chapter 2.6, p.40.
    ! Here each of the components w_i and r_i are 3-component vectors and
    ! m_i, l_i, u_i are 2*2 matrices                                       !

    integer, intent(in):: n
    real, intent(in):: &
         Lower_VVI(2,2,n), Main_VVI(2,2,n), Upper_VVI(2,2,n),Res_VI(2,n)
    real, intent(out):: &
         Weight_VI(2,n)

    integer:: j
    real   :: TildeM_VV(2,2), TildeMInv_VV(2,2), TildeMInvDotU_VVI(2,2,2:n)

    ! If tilde(M)+L.Inverted(\tilde(M))\dot.U = M, then the equation
    !      (M+L+U)W = R
    ! may be equivalently written as
    ! (tilde(M) +L).(I + Inverted(\tilde(M)).U).W=R

    character(len=*), parameter:: NameSub = 'tridiag_block22'
    !--------------------------------------------------------------------------
    TildeM_VV = Main_VVI(:,:,1)
    TildeMInv_VV = inverse_matrix(TildeM_VV)
    ! First 2-vector element of the vector, Inverted(tilde(M) + L).R
    Weight_VI(:,1) = matmul(TildeMInv_VV,Res_VI(:,1))
    do j=2,n
       ! Next 3*3 blok element of the matrix, Inverted(Tilde(M)).U
       TildeMInvDotU_VVI(:,:,j) = matmul(TildeMInv_VV,Upper_VVI(:,:,j-1))
       ! Next 3*3 block element of matrix tilde(M), obeying the eq.
       ! tilde(M)+L.Inverted(\tilde(M))\dot.U = M
       TildeM_VV = Main_VVI(:,:,j) - &
            matmul(Lower_VVI(:,:,j), TildeMInvDotU_VVI(:,:,j))
       ! Next element of inverted(Tilde(M))
       TildeMInv_VV = inverse_matrix(TildeM_VV)
       ! Next 2-vector element of the vector, Inverted(tilde(M) + L).R
       ! satisfying the eq. (tilde(M) + L).W = R
       Weight_VI(:,j) = matmul(TildeMInv_VV,Res_VI(:,j) - &
            matmul(Lower_VVI(:,:,j),Weight_VI(:,j-1)))
    end do
    do j=n-1,1,-1
       ! Finally we solve equation
       ! (I + Inverted(Tilde(M)).U).W =  Inverted(tilde(M) + L).R
       Weight_VI(:,j) = Weight_VI(:,j) &
            - matmul(TildeMInvDotU_VVI(:,:,j+1),Weight_VI(:,j+1))
    end do
  contains
    !==========================================================================
    function inverse_matrix(A_II)
      real :: inverse_matrix(2,2)
      real, intent(in) :: A_II(2,2)
      real :: Determinant
      !------------------------------------------------------------------------
      Determinant = A_II(1,1)*A_II(2,2) - A_II(1,2)*A_II(2,1)
      if (Determinant == 0.0) call CON_stop(NameSub//' failed')
      inverse_matrix(1,1) =  A_II(2,2)
      inverse_matrix(2,2) =  A_II(1,1)
      inverse_matrix(1,2) = -A_II(1,2)
      inverse_matrix(2,1) = -A_II(2,1)
      inverse_matrix = inverse_matrix/Determinant
    end function inverse_matrix
    !==========================================================================
  end subroutine tridiag_block22
  !============================================================================
  subroutine apportion_heating(&
       ! Inputs, all in SI:
    !--------------------------------------------------------------------------
       PparIn, PperpIn, PeIn, RhoIn, BIn, &
       WmajorIn, WminorIn,                &
       DissRateMajorIn, DissRateMinorIn,  &
       ! Outputs:
       QparPerQtotal, QperpPerQtotal, QePerQtotal)

    use ModConst, ONLY: cMu, cProtonMass, cElectronCharge
    ! Input pressures
    real, intent(in) :: PparIn, PperpIn, PeIn
    ! Input density and magnetic field
    real, intent(in) :: RhoIn, Bin
    ! Input wave energy densities
    real, intent(in) :: WmajorIn, WminorIn
    ! Dissipation rates:
    real, intent(in) :: DissRateMajorIn, DissRateMinorIn
    ! Outputs
    real, intent(out) :: QparPerQtotal, QperpPerQtotal, QePerQtotal

    real :: SqrtRho, Wmajor, Wminor, Qmajor, Qminor, P
    real :: BetaProton, BetaElectron, TeByTp, Vperp, Qtotal
    real :: DampingElectron, DampingPar, DampingPerp, DampingProton
    real :: GyroRadiusTimesB, InvGyroRadius, LperpInvGyroRad
    real :: WmajorGyro, WminorGyro, Wgyro
    real :: CascadeTimeMajor, CascadeTimeMinor, DeltaU, Epsilon, DeltaB, Delta
    real :: Qproton, QminorFraction, QmajorFraction
    real, parameter :: cTwoThird = 2.0/3.0, cThird = 1.0/3.0

#ifndef SCALAR
    character(len=*), parameter:: NameSub = 'apportion_heating'
    !--------------------------------------------------------------------------

    SqrtRho = PoyntingFluxPerBsi*sqrt(cMu*RhoIn)
    ! Recover wave energy densities from representatives
    Wmajor = WmajorIn*SqrtRho
    Wminor = WminorIn*SqrtRho
    ! Energy dissipation rates
    Qmajor = Wmajor*DissRateMajorIn
    Qminor = Wminor*DissRateMinorIn
    Qtotal = Qmajor + Qminor

    P = cTwoThird*PperpIn + cThird*PparIn
    BetaProton = 2.0*cMu*P/(Bin*Bin)
    TeByTp = PeIn/P

    BetaElectron = 2.0*cMu*PeIn/(Bin*Bin)

    DampingElectron = 0.01*sqrt(TeByTp/BetaProton) &
         *(1.0 + 0.17*BetaProton**1.3) &
         /(1.0 +(2800.0*BetaElectron)**(-1.25))
    DampingPar = 0.08*sqrt(sqrt(TeByTp))*BetaProton**0.7 &
               *exp(-1.3/BetaProton)

    ! Stochasting heating contributes to perpendicular ion heating.

    ! Perpendicular ion thermal speed
    Vperp = sqrt(2.0*PperpIn/RhoIn)

    GyroRadiusTimesB = Vperp*cProtonMass/cElectronCharge

    InvGyroRadius = Bin/GyroRadiusTimesB

    LperpInvGyroRad = InvGyroRadius*LperpTimesSqrtBsi/sqrt(Bin)

    WmajorGyro = Wmajor/sqrt(LperpInvGyroRad)
    WminorGyro = Wminor/sqrt(LperpInvGyroRad)

    Wgyro = WmajorGyro + WminorGyro

    ! Cascade timescale at the gyroscale
    CascadeTimeMajor = WmajorGyro/max(Qmajor,1e-30)
    CascadeTimeMinor = WminorGyro/max(Qminor,1e-30)

    ! For protons the following would be DeltaU and DeltaB at ion gyro
    ! radius, except that we assumed that the Alfven ratio is one.
    DeltaU = sqrt(Wgyro/RhoIn)
    DeltaB = sqrt(cMu*Wgyro)

    Epsilon = DeltaU/Vperp
    Delta = DeltaB/Bin

    ! Damping rate for stochastic heating.
    ! It interpolates between the beta<1 and 1<beta<30 version.
    ! This formula is at the moment only suitable for protons.
    DampingPerp = (StochasticAmplitude &
         *exp(-StochasticExponent/max(Epsilon,1e-15)) &
         + StochasticAmplitude2*sqrt(BetaProton) &
         *exp(-StochasticExponent2/max(Delta,1e-15))) &
         *RhoIn*DeltaU**3 &
         *InvGyroRadius/max(Wgyro,1e-15)

    ! Set k_parallel*V_Alfven = 1/t_minor (critical balance)
    DampingElectron = DampingElectron/max(CascadeTimeMinor,1e-30)
    DampingPar = DampingPar/max(CascadeTimeMinor, 1e-30)

    ! Total damping rate around proton gyroscale
    DampingProton = DampingElectron + DampingPar + DampingPerp

    QmajorFraction = DampingProton*CascadeTimeMajor &
         /(1.0 + DampingProton*CascadeTimeMajor)
    QminorFraction = DampingProton*CascadeTimeMinor &
         /(1.0 + DampingProton*CascadeTimeMinor)

    Qproton = (QmajorFraction*Qmajor &
         + QminorFraction*Qminor)/Qtotal

    QparPerQtotal = DampingPar/DampingProton*Qproton

    QperpPerQtotal = DampingPerp/DampingProton*Qproton

    QePerQtotal = 1 - QparPerQtotal - QperpPerQtotal
#endif
  end subroutine apportion_heating
  !============================================================================
  subroutine plot_tr(NamePlotFile, nGrid, TeSi, PeSi, iTable)

    use ModPlotFile,    ONLY: save_plot_file
    use ModLookupTable, ONLY: interpolate_lookup_table, get_lookup_table
    use ModConst,       ONLY: cBoltzmann
    use ModUtilities,   ONLY: split_string, join_string

    character(LEN=*), intent(in) :: NamePlotFile

    ! Number of grid points. The grid is uniform if Te, but not in X
    integer, intent(in)  :: nGrid

    ! The plasma parameters on top of the transition region:
    real,    intent(in)  :: TeSi, PeSi

    ! The TR is mostly used to account for the integral of the spectral
    ! intensity across the transition region, which is tabulated in the
    ! lookup table. In this routine we can visualize the integrand
    integer, optional, intent(in)  :: iTable

    ! The model is parameterized in terms of PAvr=sqrt(Pi*Pe) at Ti=Te
    ! In terms of SqrtZ: PAvr = Pi*SqrtZ = Pe/SqrtZ
    real    :: PAvrSi

    ! 1D Grid across the TR
    real :: LengthSi_I(nGrid)

    ! Plot variables: Electron temperature and density in particles per cm3
    integer, parameter :: TeSi_ = 1, NeCgs_ = 2
    real, allocatable  :: Value_VI(:,:)
    real               :: DeltaTe      ! Mesh of a temperature
    real               :: LengthPavrSi_I(nGrid)
    integer            :: i             ! Loop variable

    ! Ion density in particles per cm3:
    real    :: NiCgs

    ! To work with the emissivity table, if any
    logical :: DoPlotEmissivity !=present(iTable)
    integer :: nValue, nVar

    ! Array for variable names
    integer, parameter:: MaxVar = 200
    character(len=20) :: NameVar_I(MaxVar)
    character(len=500):: NameVarTable, NameVarPlot, NameUnitPlot
    !--------------------------------------------------------------------------
    DoPlotEmissivity = present(iTable)
    if(DoPlotEmissivity)then
       call  get_lookup_table(&
            iTable=iTable,    &
            nValue=nValue,    &
            NameVar=NameVarTable)
       call split_string(NameVarTable, MaxVar, NameVar_I, nVar, &
         UseArraySyntaxIn=.true.)
       call join_string(NameVar_I(3:nValue+2),NameVarTable) ! only var names!
       NameVarPlot = 'Length Te Ne '//trim(NameVarTable)//' TeTop PeTop'
       NameVar_I(1:nValue) = 'Response/m'
       call join_string(NameVar_I(:nValue),NameUnitPlot)
       NameUnitPlot='m K 1/cm3 '//trim(NameUnitPlot)//' K N/m2'
       allocate(Value_VI(TeSi_:NeCgs_+nValue,1:nGrid))
    else
       NameVarPlot = 'Length Te Ne TeTop PeTop'
       NameUnitPlot= 'm K 1/cm3 K N/m2'
       allocate(Value_VI(TeSi_:NeCgs_,1:nGrid))
    end if

    PAvrSi = PeSi/SqrtZ
    DeltaTe = (TeSi - TeSiMin)/(nGrid - 1)
    Value_VI(TeSi_,1) = TeSiMin
    Value_VI(NeCgs_,1) = 1.0e-6*PeSi/(cBoltzmann*TeSiMin)

    call interpolate_lookup_table(iTableTr, TeSiMin, 0.0, TrTable_V, &
         DoExtrapolate=.false.)

    ! First value is now the product of the thread length in meters times
    ! a geometric mean pressure, so that
    LengthPavrSi_I(1) = TrTable_V(LengthPAvrSi_)

    do i = 2, nGrid
       Value_VI(TeSi_,i) = Value_VI(TeSi_,i-1) + DeltaTe
       call interpolate_lookup_table(iTableTr, Value_VI(TeSi_,i), 0.0, &
            TrTable_V, DoExtrapolate=.false.)
       LengthPavrSi_I(i) = TrTable_V(LengthPAvrSi_)
       Value_VI(NeCgs_,i) = 1.0e-6*PeSi/(cBoltzmann*Value_VI(TeSi_,i))
    end do

    ! Now, LPAvrSi is the length interval multiplied by PAvrSi.
    ! Get rid of the latter factor and offset coordinate so that]
    ! LengthSi_I(1) = 0
    LengthSi_I = (LengthPavrSi_I - LengthPavrSi_I(1))/PAvrSi
    if(DoPlotEmissivity)then
       do i = 1, nGrid
          call interpolate_lookup_table(iTable, Value_VI(TeSi_,i), 0.0, &
            Value_VI(NeCgs_,i), Value_VI(NeCgs_+1:NeCgs_+nValue,i),    &
            DoExtrapolate=.false.)
          ! Multiply response function by: (1) NeCgs*NiCgs;
          ! (2) factor 1e-26 for response functions in CHIANTI tables
          ! (3) factor 100 to convert response/(cm of LoS length) to that per m
          NiCgs = 1.0e-6*PAvrSi/(SqrtZ*cBoltzmann*Value_VI(TeSi_,i))
          Value_VI(NeCgs_+1:NeCgs_+nValue,i) =    &
               Value_VI(NeCgs_+1:NeCgs_+nValue,i)*&
               Value_VI(NeCgs_,i)*NiCgs*          &
               1.0e-26*                           & ! response function per cm
               1.0e2                                ! response function per m
       end do
    end if
    call save_plot_file(&
         NameFile=NamePlotFile,   &
         ParamIn_I = [TeSi, PeSi],&
         NameVarIn = NameVarPlot, &
         nDimIn=1,                &
         Coord1In_I = LengthSi_I, &
         VarIn_VI = Value_VI,     &
         StringHeaderIn = 'Analytical model for transition region: ['//&
         trim(NameUnitPlot)//']')
    deallocate(Value_VI)

  end subroutine plot_tr
  !============================================================================
  subroutine save_plot_thread(OpenThread1, NameFile)

    use ModPlotFile, ONLY: save_plot_file

    type(OpenThread),intent(inout) :: OpenThread1
    character(LEN=*), intent(in) :: NameFile
    integer, parameter :: sWPlus_ = 5, sWminus_ = 6,sTi_=7,  sTe_=8, &
         iPlotVar_V(Rho_:sWminus_) = [Rho_, U_, P_, Pe_, Wmajor_, Wminor_]
    real :: Value_VI(Rho_:sTe_+1,-nPointMax:-1)
    real, pointer :: State_VG(:,:)
    integer :: nCell
    !--------------------------------------------------------------------------
    nCell = OpenThread1%nCell
    State_VG => OpenThread1%State_VG
    Value_VI(Rho_:sWminus_,-nCell:-1) = State_VG(iPlotVar_V,-nCell:-1)
    if(sign(1.0, OpenThread1%OpenFlux) > 0)then
       Value_VI(sWplus_,-nCell:-1) = PoyntingFluxPerBsi*sqrt(cMu*&
            State_VG(Rho_,-nCell:-1))*State_VG(Wmajor_,-nCell:-1)
       Value_VI(sWminus_,-nCell:-1) = PoyntingFluxPerBsi*sqrt(cMu*&
            State_VG(Rho_,-nCell:-1))*State_VG(Wminor_,-nCell:-1)

    else
       Value_VI(sWplus_,-nCell:-1) = PoyntingFluxPerBsi*sqrt(cMu*&
            State_VG(Rho_,-nCell:-1))*State_VG(Wminor_,-nCell:-1)
       Value_VI(sWminus_,-nCell:-1) = PoyntingFluxPerBsi*sqrt(cMu*&
            State_VG(Rho_,-nCell:-1))*State_VG(Wmajor_,-nCell:-1)
    end if
    Value_VI(sTi_,-nCell:-1) = State_VG(P_,-nCell:-1)*cProtonMass/&
         (cBoltzmann*State_VG(Rho_,-nCell:-1))
    Value_VI(sTe_,-nCell:-1) = OpenThread1%Te_G(-nCell:-1)
    Value_VI(sTe_+1,-nCell:-1) = &
         OpenThread1%B_F(-nCell:-1)*Si2Gs
    call save_plot_file(NameFile = NameFile,  &
         nDimIn  = 1,                       &
         ParamIn_I= [OpenThread1%TeTr,&
         OpenThread1%uTr,&
         OpenThread1%PeTr,&
         OpenThread1%Dt],&
         VarIn_VI= Value_VI(:,-nCell:-1), &
         TypeFileIn    = 'ascii',           &
         CoordIn_I  = 0.5*(OpenThread1%R_F(-nCell:-1) +&
         OpenThread1%R_F(-nCell+1:0)), &
         StringHeaderIn  = 'Thread file ', &
         NameUnitsIn  = &
         '[Rsun] '//&
         '[kg/m3] [m/s] [J/m3] [J/m3] [J/m3] [J/m3] [K] [K] [Gs] '//&
         '[K] [m/s] [J/m3] [s]',&
         NameVarIn = 'R '//&
         'Rho U P Pe Wp Wm Ti Te B'//&
         ' TeTR UTR PeTR Dt' )
  end subroutine save_plot_thread
  !============================================================================
  subroutine test
    use ModConst, ONLY: Rsun, cBoltzmann
    real :: PeFace, TeFace, HeatFluxFace
    real :: Te_I(100) = 1.0e5  ! K
    real :: Ni_I(99) = 3.0e14 ! m-3
    real :: Ti_I(99) = 1.0e5  ! K
    real :: uFace = 0.0
    real :: Ds_I(0:100) = 1.0e-3*Rsun
    real :: B_I(100) = 5.0e-4      ! T
    integer :: iPoint, iTime
    !--------------------------------------------------------------------------
    call solve_tr_face(TeCell = 1.0e6, &
         uFace = uFace,                   & ! cell-centered, in SI
         LengthTr = 0.05*Rsun,            &
         DeltaS   = 0.0025*Rsun,          & ! Distance face-to-center
         TeFaceOut   = TeFace,            &
         PeFaceOut   = PeFace,            &
         HeatFluxOut = HeatFluxFace)
    write(*,'(a,es13.6)')'TeFace = ',TeFace
    write(*,'(a,es13.6)')'PeFace = ',PeFace
    write(*,'(a,es13.6)')'uFace  = ',uFace
    write(*,'(a,es13.6)')'Heat flux = ', HeatFluxFace
    Te_I(100) = 2.0e6
    do iTime = 1,36
       call advance_heat_conduction_ta(99, 100.0, Te_I, Ti_I, Ni_I, Ds_I, &
            uFace, B_I, PeFaceOut = PeFace, DoLimitTimestep=.true.)
       Ni_I = PeFace/(cBoltzmann*Te_I(1:99))
    end do
    write(*,'(a,es13.6)')'PeFace = ',PeFace,' PeFace*Length=', PeFace*0.1*Rsun
    write(*,'(a)')'iPoint   Te     Ti    Ni'
    do iPoint = 1,99
       write(*,'(i2,3es13.6)')iPoint, Te_I(iPoint), Ti_I(iPoint), Ni_I(iPoint)
    end do
  end subroutine test
  !============================================================================
end module ModTransitionRegion
!==============================================================================
