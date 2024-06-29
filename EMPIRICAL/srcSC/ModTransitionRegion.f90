!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModTransitionRegion

  use ModUtilities, ONLY: CON_stop
  use ModConst,       ONLY: cBoltzmann, ceV
  implicit none
  SAVE
  PRIVATE
  ! The Poynting flux to magnetic field ratio (one of the input parameters
  ! in SI unins)
  real, public :: PoyntingFluxPerBSi = 1.0e6 ! W/(m^2 T)
  real, public :: LperpTimesSqrtBSi = 7.5e4  ! m T^(1/2)
  
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

  real, public :: cExchangeRateSi
  
  ! Coulomb logarithm
  real, public  :: CoulombLog = 20.0
  
  real, parameter :: cTwoSevenths = 2.0/7.0
  real, parameter :: cTolerance   = 1.0e-6

  ! Correspondent named indexes: meaning of the columns in the table
  integer, parameter, public :: LengthPavrSi_ = 1, uHeat_ = 2, &
       HeatFluxLength_ = 3, dHeatFluxXOverU_ = 4, LambdaSi_=5, &
       DlogLambdaOverDlogT_ = 6, uLengthPavrSi_ = 7

  ! Control parameter: minimum temerature in the TR table
  real, public :: TeSiMin = 5.0e4
  real, public :: SqrtZ   = 1.0
  real :: Z = 1.0
  public :: init_tr, check_tr_table, integrate_emission, plot_tr, &
       read_tr_param, solve_tr_face, test

  ! Table numbers needed to use lookup table
  integer         :: iTableRadCool = -1
  integer, public :: iTableTr = -1

  ! Chromosphere top boundary
  real, public, parameter :: rChromo = 1.0

  ! By default, this logical is .false. the entholpy increase needed to
  ! heat the plasma flow across the transition region to the top temperature
  ! If logical set is true, the energy frlux to/from the first control
  ! volume is accounted for
  logical, public :: UseChromoEvap  = .false.
  ! To apply CromoEvaporation, the factor below is non-zero.
  real,    public :: ChromoEvapCoef = 0.0
  real, parameter :: TeTrMin = 1.0e4
  real, parameter :: TeTrMax = 1.0e8
  integer, parameter :: nPointTe = 500, nPointU = 31
   ! Global arrays used in calculating the tables
  real, dimension(nPointTe) :: TeSi_I, LambdaSi_I, DLogLambdaOverDLogT_I, &
       LengthPe_I, UHeat_I, dFluxXLengthOverDU_I
  real            :: SemiIntUheat_I(1:nPointTe-1)
  real            :: DeltaLogTe, DeltaLogTeCoef
  ! Max speed
  real, parameter :: uMax = 200.0
  real, parameter :: DeltaU = 1.5*uMax/(nPointU - 1)
  ! Ionization potential for hydrogen
  ! cPotential_II(1,1) from util/CRASH/src/ModIonizPotential.f90:
  real, parameter :: cIonizPotentialH =  13.59844  ! eV
  ! Ratio of enthalpy flux, divided by (u(Tmin)*Pe):
  real, parameter :: FluxTouPe = cIonizPotentialH*ceV/(cBoltzmann*TeTrMin) + 5
  ! Needed for initialization:
  logical, private :: DoInit = .true.
  interface advance_heat_conduction
     module procedure advance_heat_conduction_ta
     module procedure advance_heat_conduction_ss
  end interface advance_heat_conduction
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
       call read_var('UseChromoEvap',UseChromoEvap)
       if(UseChromoEvap)then
          ChromoEvapCoef = TeTrMin
       else
          ChromoEvapCoef = 0.0
       end if
    case("#TURBULENCE")
       call read_var("PoyntingFluxPerBSi",PoyntingFluxPerBSi)
       call read_var("LperpTimesSqrtBSi",LperpTimesSqrtBSi)
    case default
       call CON_stop(NameSub//": unknown command="//trim(NameCommand))
    end select
  end subroutine read_tr_param
  !============================================================================
  subroutine check_tr_table(TypeFileIn,iComm)
    use ModLookupTable, ONLY: i_lookup_table, &
         init_lookup_table, make_lookup_table

    character(LEN=*),optional,intent(in)::TypeFileIn
    integer, optional, intent(in):: iComm

    character(len=5)::TypeFile

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
         NameVar =                                               &
         'logTe u LPe UHeat FluxXLength '//                      &
         'dFluxXLegthOverDU Lambda dLogLambdaOverDLogT uLPe',    &
         nIndex_I = [nPointTe,nPointU],                          &
         IndexMin_I = [TeTrMin,-0.5*uMax],                       &
         IndexMax_I = [TeTrMax, uMax],                           &
         NameFile = 'TR8.dat',                                   &
         TypeFile = TypeFile,                                    &
         StringDescription =                                     &
         'Model for transition region: '//                       &
         '[K] [m/s] [N/m] [m/s] [W/m] [1] [W*m3/(k_B2)] [1] [W/m]')

    ! The table is now initialized.
    iTableTr = i_lookup_table('TR')
    ! Fill in the table
    call make_lookup_table(iTableTr, calc_tr_table, iComm)
  end subroutine check_tr_table
  !============================================================================
  subroutine calc_tr_table(iTableIn, Arg1, Arg2, Value_V)

    use ModLookupTable, ONLY: interpolate_lookup_table
    integer, intent(in):: iTableIn
    real, intent(in)   :: Arg1, Arg2
    real, intent(out)  :: Value_V(:)
    integer            :: iTe, iU
    integer, save      :: iUlast = -1
    real   :: LambdaCgs_V(1)

    logical, save:: IsFirstCall = .true.
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
    real,save :: FactorStep, uOverTmin5, SqrtOfU2
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_tr_table'
    !--------------------------------------------------------------------------
    if(IsFirstCall)then
       ! write(*,*)'Start init'
       IsFirstCall=.false.
       FactorStep = exp(DeltaLogTe)
       do iTe = 1, nPointTe
          if(iTe==1)then
             TeSi_I(1) = TeTrMin
          else
             TeSi_I(iTe) = TeSi_I(iTe-1)*FactorStep
          end if
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
    end if
    iU = nint( (Arg2 + 0.50*uMax)/ DeltaU) + 1
    if(iU/=iUlast)then
       uOverTmin5 = 5*(Arg2/TeTrMin)*DeltaLogTe
       LengthPe_I = 0.0; uHeat_I = 0.0
       UHeat_I(1) = (max(Arg2, 0.0)*FluxTouPe)**2
       do iTe = 2, nPointTe
          ! Integrate \sqrt{2\int{\kappa_0\Lambda Te**1.5 d(log T)}}/k_B
          ! Predictor at half step
          SqrtOfU2 = sqrt(UHeat_I(iTe-1) )
          SemiIntUheat_I(iTe-1) = sqrt( UHeat_I(iTe-1) + &
               SqrtOfU2*uOverTmin5*TeSi_I(iTe-1)       + &
               LambdaSi_I(iTe-1)*TeSi_I(iTe-1)**1.50*DeltaLogTeCoef)
          UHeat_I(iTe) = UHeat_I(iTe-1) + &
               SemiIntUheat_I(iTe-1)*uOverTmin5*(TeSi_I(iTe-1) + TeSi_I(iTe))&
               + (LambdaSi_I(iTe-1)*TeSi_I(iTe-1)**1.50 + &
               LambdaSi_I(iTe)*TeSi_I(iTe)**1.50)*DeltaLogTeCoef
          UHeat_I(iTe-1) = SqrtOfU2
       end do
       UHeat_I(nPointTe) = sqrt(UHeat_I(nPointTe))

       ! Temporary fix to avoid a singularity in the first point
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
       iUlast = iU
    end if
    iTe = 1 + nint(log(Arg1/TeTrMin)/DeltaLogTe)
    Value_V(LengthPAvrSi_:uLengthPAvrSi_) = &
         [ LengthPe_I(iTe), UHeat_I(iTe), &
         LengthPe_I(iTe)*UHeat_I(iTe), dFluxXLengthOverDU_I(iTe),   &
         LambdaSi_I(iTe)/cBoltzmann**2, DLogLambdaOverDLogT_I(iTe), &
         LengthPe_I(iTe)*Arg2]
  end subroutine calc_tr_table
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
    ! Tabulated analytical solution:
    real    :: TrTable_V(LengthPAvrSi_:uLengthPAvrSi_)
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
    call interpolate_lookup_table(iTableTr, TeSiMin, 0.0, TrTable_V, &
         DoExtrapolate=.false.)
    ! First value is now the product of the thread length in meters times
    ! a geometric mean pressure, so that
    LengthPavrSi_I(1) = TrTable_V(LengthPAvrSi_)
    do i = 1, nTrGrid
       Te_I(i +1) = Te_I(i) + DeltaTe
       call interpolate_lookup_table(iTableTr, Te_I(i + 1), 0.0, &
            TrTable_V, DoExtrapolate=.false.)
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
  subroutine solve_tr_face(TeCell, uPeOverTeCell, & ! cell-centered, in SI
       LengthTr, DeltaS, & ! Distances photosphere-to-face, face-to-center
       TeFaceOut, HeatFluxOut, DFluxOverDConsOut, &
       PeFaceOut, uFaceOut) ! SI, face centered

    ! For given values of Te and u*Pavr/Te at the cell center, all inputs
    ! are in SI and expressed in K and W/(m.K). The distance from the cell
    ! to the face is DeltaS, from the face to photosphere is LengthTr.
    ! All distances are measured along the magnetic field line in meters.
    ! Outputs provided are TeFace, PeFace, uFace, HeatFluxFace, all in SI:
    ! [K], [N/m^2], [m/s], [W/m^2], at the face.
    use ModLookupTable, ONLY: interpolate_lookup_table

    ! The plasma parameters at the cell center region:
    real,    intent(in)  :: TeCell, uPeOverTeCell
    ! Distances photosphere-to-face, face-to-center
    real,    intent(in)  :: LengthTr, DeltaS
    ! Parameters at the face, connected by the analytical TR solution to the
    ! photosphere:
    real, optional, intent(out) :: TeFaceOut
    real, optional, intent(out) :: HeatFluxOut, dFluxOverdConsOut
    real, optional, intent(out) :: PeFaceOut, uFaceOut

    ! Tabulated analytical solution:
    real    :: Value_V(LengthPAvrSi_:uLengthPAvrSi_), TeFace, HeatFluxFace
    real    :: ConsFace, dConsFace, ConsCell, uMinTr, Tolerance
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
       call interpolate_lookup_table(iTableTr, Arg1In=TeFace,  &
            iVal=uLengthPAvrSi_,                               &
            ValIn=ChromoEvapCoef*LengthTr*uPeOverTeCell/SqrtZ, &
            Value_V=Value_V,                                   &
            Arg2Out=uMinTr,                                    &
            DoExtrapolate=.false.)
       HeatFluxFace = Value_V(HeatFluxLength_)/LengthTr
       ! Solve equation HeatFlux = (ConsCell - ConsFace)/DeltaS
       Res = DeltaS*HeatFluxFace + ConsFace - ConsCell
       if(abs(Res) < Tolerance) EXIT
       if(iCounter == iCounterMax)call CON_stop(NameSub//': No convergence')
       ! New iteration:
       iCounter = iCounter + 1
       dResdCons = DeltaS*Value_V(dHeatFluxXOverU_)/LengthTr + 1
       dConsFace = -Res/dResdCons
       ConsFace = ConsFace + dConsFace
       TeFace = (ConsFace*3.50/HeatCondParSi)**cTwoSevenths
    end do
    if(present(TeFaceOut))TeFaceOut = TeFace
    if(present(HeatFluxOut))HeatFluxOut = HeatFluxFace
    if(present(dFluxOverdConsOut))then
       dFluxOverdConsOut = Value_V(dHeatFluxXOverU_)/LengthTr
       dFluxOverdConsOut = dFluxOverdConsOut/(1 + DeltaS*dFluxOverdConsOut)
    end if
    if(present(PeFaceOut))PeFaceOut = Value_V(LengthPavrSi_)*SqrtZ/LengthTr
    if(present(uFaceOut)) uFaceOut  = uMinTr*TeFace/TeTrMin
  end subroutine solve_tr_face
  !============================================================================
  subroutine advance_heat_conduction_ta(nPoint, Dt, Te_I, Ti_I, Ni_I, Ds_I, &
       U1, B_I, BCellIn_I, PeFaceOut, TeFaceOut, UfaceOut)
    integer, intent(in)    :: nPoint
    real,    intent(in)    :: Dt
    real,    intent(inout) :: Te_I(nPoint+1), Ti_I(nPoint)
    real,    intent(in)    :: Ni_I(nPoint), Ds_I(0:nPoint+1), U1, B_I(nPoint+1)
    real, optional, intent(in) :: BCellIn_I(nPoint)
    real, optional, intent(out) :: PeFaceOut, TeFaceOut, UfaceOut
    ! MISC:
    real :: Dt_I(nPoint)
    !--------------------------------------------------------------------------
    Dt_I = Dt
    call advance_heat_conduction_ss(nPoint, Dt_I, Te_I, Ti_I, Ni_I, Ds_I, &
       U1, B_I, BCellIn_I, PeFaceOut, TeFaceOut, UfaceOut)
  end subroutine advance_heat_conduction_ta
  !============================================================================
  subroutine advance_heat_conduction_ss(nPoint, Dt_I, Te_I, Ti_I, Ni_I, Ds_I,&
       U1, B_I, BCellIn_I, PeFaceOut, TeFaceOut, UfaceOut)
    use ModConst, ONLY: cBoltzmann
    use ModLookupTable, ONLY: interpolate_lookup_table
    integer, intent(in)    :: nPoint
    real,    intent(in)    :: Dt_I(nPoint)
    real,    intent(inout) :: Te_I(nPoint+1), Ti_I(nPoint)
    real,    intent(in)    :: Ni_I(nPoint), Ds_I(0:nPoint+1), U1, B_I(nPoint+1)
    real, optional, intent(in) :: BCellIn_I(nPoint)
    real, optional, intent(out) :: PeFaceOut, TeFaceOut, UfaceOut

    real :: Res_VI(2,nPoint), Weight_VI(2,nPoint), bFaceInv_I(1:nPoint+1),  &
         Lower_VVI(2,2,nPoint), Main_VVI(2,2,nPoint), Upper_VVI(2,2,nPoint),&
         SpecHeat_I(nPoint), ExchangeRate_I(nPoint), Cons_I(nPoint+1),      &
         TeStart_I(nPoint), TiStart_I(nPoint), DeltaEnergy_I(nPoint),       &
         DeltaIonEnergy_I(nPoint), dCons_VI(2,nPoint), DtInv_I(nPoint),     &
         HeatFlux2Tr, dFluxOverdCons, Value_V(LengthPAvrSi_:uLengthPAvrSi_),&
         Cooling, PeFace, TeFace, uFace, BcellInv_I(nPoint)
    real, parameter :: QuasiCfl = 0.85
    integer, parameter :: Cons_ = 1, Ti_=2, nIterMax = 40
    ! Misc:
    ! Loop variable
    integer :: iPoint, iIter
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
            uPeOverTeCell = u1*cBoltzmann*Ni_I(1), & 
            LengthTr      = Ds_I(0),               &
            DeltaS        = 0.50*Ds_I(1),          &
            HeatFluxOut   = HeatFlux2Tr,           &
            DFluxOverDConsOut = dFluxOverdCons,    &
            PeFaceOut = PeFace,                    &
            TeFaceOut = TeFace,                    &
            uFaceOut  = uFace)
       ! Add left heat flux to the TR
       Res_VI(Cons_,1) = Res_VI(Cons_,1) - HeatFlux2Tr*BFaceInv_I(1)
       ! Linearize left heat flux to the TR
       Main_VVI(Cons_,Cons_,1) = &
            - Upper_VVI(Cons_,Cons_,1) + dFluxOverdCons*BFaceInv_I(1)
       ! Radiative cooling, limit timestep
       do iPoint = 1, nPoint
          call interpolate_lookup_table(iTableTR, Te_I(iPoint), 0.0, &
               Value_V, &
               DoExtrapolate=.false.)
          Cooling =Ds_I(iPoint)*BcellInv_I(iPoint)*Value_V(LambdaSI_)*Z*&
               (cBoltzmann*Ni_I(iPoint))**2
          Res_VI(Cons_,iPoint) = Res_VI(Cons_,iPoint) - Cooling
          ! linearized -dCooling/dCons
          Main_VVI(Cons_,Cons_,iPoint) = Main_VVI(Cons_,Cons_,iPoint) +&
               Value_V(DLogLambdaOverDLogT_)*Cooling/(3.5*Cons_I(iPoint))
          ! Limit time step:
          if(iIter==1)DtInv_I(iPoint) = max(1/Dt_I(iPoint),   &
               (3.50 + max(0.0, -Value_V(DLogLambdaOverDLogT_)))*Cooling/&
               (QuasiCfl*Z*Te_I(iPoint)*SpecHeat_I(iPoint)))
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
             write(*,*)'iPoint Cons_VI(:,iPoint)=',iPoint,Cons_I(iPoint),&
                  Ti_I(iPoint) + DCons_VI(Ti_,iPoint)
          end do
       end if
       ! Recover temperature
       Te_I(1:nPoint) = (3.5*Cons_I(1:nPoint)/HeatCondParSi)**cTwoSevenths
       Ti_I(1:nPoint) = Ti_I(1:nPoint) + DCons_VI(Ti_,1:nPoint)
       if(all(abs(DCons_VI(Cons_,1:nPoint))<cTolerance*Cons_I(1:nPoint)))EXIT
    end do
    if(present(PeFaceOut))PeFaceOut = PeFace
    if(present(TeFaceOut))TeFaceOut = TeFace
    if(present(UfaceOut))UfaceOut = uFace
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

    ! Tabulated analytical solution:
    real    :: TrTable_V(LengthPAvrSi_:uLengthPAvrSi_)

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
         Coord1In_I = LengthSi_i, &
         VarIn_VI = Value_VI,     &
         StringHeaderIn = 'Analytical model for transition region: ['//&
         trim(NameUnitPlot)//']')
    deallocate(Value_VI)

  end subroutine plot_tr
  !============================================================================
  subroutine test
    use ModConst, ONLY: Rsun, cBoltzmann
    real :: PeFace, uFace, TeFace, HeatFluxFace
    real :: Te_I(100) = 1.0e5  ! K
    real :: Ni_I(99) = 3.0e14 ! m-3
    real :: Ti_I(99) = 1.0e5  ! K
    real :: U1 = 0.0
    real :: Ds_I(0:100) = 1.0e-3*Rsun
    real :: B_I(100) = 5.0e-4      ! T
    integer :: iPoint, iTime
    !--------------------------------------------------------------------------
    call solve_tr_face(TeCell = 1.0e6, &
         uPeOverTeCell = 0.0,          & ! cell-centered, in SI
         LengthTr = 0.05*Rsun,         &
         DeltaS   = 0.0025*Rsun,        & ! Distance face-to-center
         TeFaceOut   = TeFace,            &
         PeFaceOut   = PeFace,            &
         uFaceOut    = uFace,             &
         HeatFluxOut = HeatFluxFace)
    write(*,'(a,es13.6)')'TeFace = ',TeFace
    write(*,'(a,es13.6)')'PeFace = ',PeFace
    write(*,'(a,es13.6)')'uFace  = ',uFace
    write(*,'(a,es13.6)')'Heat flux = ', HeatFluxFace
    Te_I(100) = 2.0e6
    do iTime = 1,36
       call advance_heat_conduction_ta(99, 100.0, Te_I, Ti_I, Ni_I, Ds_I, &
            U1, B_I, PeFaceOut = PeFace)
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
