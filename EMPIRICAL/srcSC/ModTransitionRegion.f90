!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModTransitionRegion

  use ModUtilities, ONLY: CON_stop
  use ModConst,       ONLY: cBoltzmann, ceV
  implicit none
  SAVE
  PRIVATE
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

  ! Coulomb logarithm
  real, public  :: CoulombLog = 20.0

  ! Correspondent named indexes: meaning of the columns in the table
  integer, parameter, public :: LengthPavrSi_ = 1, uHeat_ = 2, &
       HeatFluxLength_ = 3, dHeatFluxXOverU_ = 4, LambdaSi_=5, &
       DlogLambdaOverDlogT_ = 6, uLengthPavrSi_ = 7

  ! Control parameter: minimum temerature in the TR table
  real, public :: TeSiMin = 5.0e4
  real, public :: SqrtZ   = 1.0
  public :: init_tr, check_tr_table, integrate_emission, & ! plot_tr, &
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

contains
  !============================================================================
  subroutine init_tr(Z, TeChromoSi, iComm)
    use ModConst,       ONLY: kappa_0_e
    use ModLookupTable, ONLY: i_lookup_table
    !    use BATL_lib,       ONLY: test_start, test_stop, iProc
    ! INPUTS
    ! Average ion charge
    real, intent(in) :: Z
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
    SqrtZ = sqrt(Z)
    TeSiMin = TeChromoSi
    DeltaLogTe = log(TeTrMax/TeTrMin)/(nPointTe - 1)
    ! electron heat conduct coefficient for single charged ions
    ! = 9.2e-12 W/(m*K^(7/2))
    HeatCondParSi = kappa_0_e(CoulombLog)
    iTableTr = i_lookup_table('TR')
    if(iTableTr<=0)then
      iTableRadCool = i_lookup_table('radcool')
       if(iTableRadCool <=0 )&
            call CON_stop('To create TR table, the radcool table is needed')
       call check_tr_table(iComm=iComm)
    end if
    call test
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
    real    :: LengthPavrSi_I(nTrGrid + 1), TeSi_I(nTrGrid + 1)
    real    :: DeltaLengthPavrSi_I(nTrGrid + 1)
    integer ::  i, iVar ! Loop variables
    ! Electron density in particles per cm3:
    real    :: NeCgs, NiCgs
    !--------------------------------------------------------------------------
    PAvrSi = PeSi/SqrtZ
    DeltaTe = (TeSi - TeSiMin)/nTrGrid
    TeSi_I(1) = TeSiMin
    call interpolate_lookup_table(iTableTr, TeSiMin, 0.0, TrTable_V, &
         DoExtrapolate=.false.)
    ! First value is now the product of the thread length in meters times
    ! a geometric mean pressure, so that
    LengthPavrSi_I(1) = TrTable_V(LengthPAvrSi_)
    do i = 1, nTrGrid
       TeSi_I(i +1) = TeSi_I(i) + DeltaTe
       call interpolate_lookup_table(iTableTr, TeSi_I(i + 1), 0.0, &
            TrTable_V, DoExtrapolate=.false.)
       LengthPavrSi_I(i + 1) = TrTable_V(LengthPAvrSi_)
       DeltaLengthPavrSi_I(i) = LengthPavrSi_I(i + 1) - LengthPavrSi_I(i)
       TeAvrSi_I(i) = (TeSi_I(i + 1) + TeSi_I(i))*0.50
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
       TeFace, PeFace, uFace, HeatFluxFace) ! SI, face centered

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
    real,    intent(out) :: TeFace, PeFace, uFace, HeatFluxFace

    ! Tabulated analytical solution:
    real    :: Value_V(LengthPAvrSi_:uLengthPAvrSi_)
    real, parameter :: cTwoSevenths = 2.0/7.0
    real    :: ConsFace, dConsFace, ConsCell, uMinTr, Tolerance
    ! Newton-Rapson coefficients: dConsFace = -Res/dResdCons
    real    :: Res, dResdCons
    integer :: iCounter
    integer, parameter :: iCounterMax=20
    real, parameter    :: cTolerance = 1.0e-6
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
    PeFace = Value_V(LengthPavrSi_)*SqrtZ/LengthTr
    uFace  = uMinTr*TeFace/TeTrMin
    write(*,*)log10(TeFace),Value_V
  end subroutine solve_tr_face
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
    use ModConst, ONLY: Rsun
    real :: PeFace, uFace, TeFace, HeatFluxFace
    !--------------------------------------------------------------------------
    call solve_tr_face(TeCell = 1.0e6, &
         uPeOverTeCell = 0.0,          & ! cell-centered, in SI
         LengthTr = 0.05*Rsun,         &
         DeltaS   = 0.0025*Rsun,        & ! Distance face-to-center
         TeFace   = TeFace,            &
         PeFace   = PeFace,            &
         uFace    = uFace,             &
         HeatFluxFace = HeatFluxFace)
    write(*,'(a,e13.6)')'TeFace = ',TeFace
    write(*,'(a,e13.6)')'PeFace = ',PeFace
    write(*,'(a,e13.6)')'uFace  = ',uFace
    write(*,'(a,e13.6)')'Heat flux = ', HeatFluxFace
  end subroutine test
  !============================================================================
end module ModTransitionRegion
!==============================================================================
