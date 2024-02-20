!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModCosmicRay
  implicit none
  interface local_interstellar_spectrum
     module procedure local_interstellar_spectrum_s
     module procedure local_interstellar_spectrum_a
  end interface local_interstellar_spectrum
contains
  !============================================================================
  subroutine local_interstellar_spectrum_s(MomentumSi, XyzSi_D,    & ! Input
       DistTimesP2Si,                                              & ! Output
       TypeLisBcIn, A, Z)                                          ! Optional

    ! INPUTS
    real, intent(in) :: MomentumSi    ! Particle momentum in Si units
    real, intent(in) :: XyzSi_D(3)    ! Where to set BC
    ! OUTPUT
    real, intent(out):: DistTimesP2Si ! Distribution function*p**2, Si
    ! OPTIONAL
    real, OPTIONAL, intent(in) :: A   ! To convert to energy per nuc
    real, OPTIONAL, intent(in) :: Z   ! To convert to rigidity
    character(LEN=18),OPTIONAL,intent(in):: TypeLisBcIn ! Decide which LIS
    ! LOCAL VARS
    real :: FakeMomentumArr_I(1), FakeDistTimesP2_I(1)
    !--------------------------------------------------------------------------
    FakeMomentumArr_I(1) = MomentumSi
    call local_interstellar_spectrum_a(1,FakeMomentumArr_I, XyzSi_D, & ! Input
         FakeDistTimesP2_I,                                          & ! Output
         TypeLisBcIn, A, Z)                                          ! Optional
    DistTimesP2Si = FakeDistTimesP2_I(1)
  end subroutine local_interstellar_spectrum_s
  !============================================================================
  subroutine local_interstellar_spectrum_a(nP,MomentumSi_I, XyzSi_D, & ! Input
       DistTimesP2Si_I,                                              & ! Output
       TypeLisBcIn, A, Z)                                            ! Optional
    ! Formulae for LISM flux by Burger et al. 2000 (doi: 10.1029/2000JA000153)
    ! and modified by Usoskin et al. 2005 (doi: 10.1029/2005JA011250)
    ! Here we refer to Eq.(2) in Usoskin et al. 2005 as default
    ! We also include Eq.(13) in Corti et al. 2019 for potential use

    use ModConst, ONLY: cGeV, cLightSpeed, cRmeProton,  &
         cElectronCharge
    use ModUtilities, ONLY: CON_stop
    ! INPUTS
    integer, intent(in) :: nP  ! The number of meshes in the logarithmic grid
    real, intent(in) :: MomentumSi_I(1:nP)   ! Particle momenta in Si units
    real, intent(in) :: XyzSi_D(3)           ! Where to set BC
    ! OUTPUT
    real, intent(out):: DistTimesP2Si_I(1:nP) ! Distribution function*p**2, Si
    ! OPTIONAL
    real, OPTIONAL, intent(in) :: A          ! To convert to energy per nuc
    real, OPTIONAL, intent(in) :: Z          ! To convert to rigidity
    character(LEN=18),OPTIONAL,intent(in):: TypeLisBcIn ! Decide which LIS
    ! LOCAL VARS
    character(LEN=18) :: TypeLisBc               ! Decide LIS to use finally
    real :: EnergySi_I(1:nP), RigiditySi_I(1:nP) ! E_K and R(E_K) in SI units
    real :: RigidityGv_I(1:nP)                   ! R(E_K) in the unit of GV
    real :: DistTimesP2Gn_I(1:nP) ! Distribution function*p**2, [.../(GeV/nuc)]
    real :: Ai, Zi                ! Default values for A and Z if not present
    ! For Eq.(13) in Corti et al., 2019 (doi: 10.3847/1538-4357/aafac4)
    real, parameter :: nFit = 5658.0, Gamma0Fit = 1.669
    real, parameter :: RigidityFit(3) = [0.572, 6.2, 540.0]
    real, parameter :: sFit(3) = [1.78, 3.89, 1.53]
    real, parameter :: DeltaFit(3) = [-4.117, -0.423, -0.26]
    !--------------------------------------------------------------------------
    Ai = 1.0; Zi = 1.0
    if(present(A)) Ai = A
    if(present(Z)) Zi = Z

    ! p[Si] -> E_K[Si]
    EnergySi_I = sqrt((MomentumSi_I*cLightSpeed)**2 +  &
         (Ai*cRmeProton)**2) - Ai*cRmeProton
    ! E_K[Si] -> R(E_K)[Si]
    RigiditySi_I = Ai/abs(Zi*cElectronCharge)*         &
         sqrt(EnergySi_I*(EnergySi_I+2*cRmeProton))
    ! R(E_K)[Si] -> R(E_K)[GV]
    RigidityGv_I = RigiditySi_I/cGeV

    if(present(TypeLisBcIn)) then
       TypeLisBc = TypeLisBcIn
    else
       TypeLisBc = 'default'
    end if

    select case(TypeLisBc)
    case('default', 'usoskin2005')
       ! For Eq.(2) in Usoskin et al. 2005 (doi: 10.1029/2005JA011250)
       ! R(E_K)[GV] -> j_LIS[.../(GeV/nuc)]
       DistTimesP2Gn_I = 1.9E+4*RigidityGv_I**(-2.78)/    &
            (1.0+0.4866*RigidityGv_I**(-2.51))
       ! j_LIS[.../(GeV/nuc)] -> p^2*f|_{infty}[Si]
       DistTimesP2Si_I = DistTimesP2Gn_I/(Ai*cGeV)
    case('corti2019')
       ! For Eq.(13) in Corti et al. 2019 (doi: 10.3847/1538-4357/aafac4)
       ! R(E_K)[GV] -> j_LIS[.../(GV)] -> j_LIS[.../(GeV/nuc)]
       DistTimesP2Gn_I = nFit*(RigidityGv_I**Gamma0Fit) *                  &
            ((1.0+(RigidityGv_I/RigidityFit(1))**sFit(1))/                 &
            (1.0+RigidityFit(1)**(-sFit(1))))**(DeltaFit(1)/sFit(1)) *     &
            ((1.0+(RigidityGv_I/RigidityFit(2))**sFit(2))/                 &
            (1.0+RigidityFit(2)**(-sFit(2))))**(DeltaFit(2)/sFit(2)) *     &
            ((1.0+(RigidityGv_I/RigidityFit(3))**sFit(3))/                 &
            (1.0+RigidityFit(3)**(-sFit(3))))**(DeltaFit(3)/sFit(3)) *     &
            Ai/abs(Zi*cElectronCharge)*(EnergySi_I+cRmeProton)/            &
            sqrt(EnergySi_I*(EnergySi_I+2*cRmeProton))
       ! j_LIS[.../(GeV/nuc)] -> p^2*f|_{infty}[Si]
       DistTimesP2Si_I = DistTimesP2Gn_I/(Ai*cGeV)
    case default
       call CON_stop('Unknown type of LIS of: '//TypeLisBc)
    end select
  end subroutine local_interstellar_spectrum_a
  !============================================================================
end module ModCosmicRay
!==============================================================================
