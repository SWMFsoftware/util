!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModTransitionRegion
  implicit none
  interface local_interstellar_spectrum
     module procedure local_interstellar_spectrum_s
     module procedure local_interstellar_spectrum_a
  end interface local_interstellar_spectrum
contains
  !============================================================================
  subroutine local_interstellar_spectrum_s(MomentumSi, XyzSi_D,    & ! Input
       DistTimesP2Si,                                              & ! Output
       A, Z)                                                         ! Optional

    ! INPUTS
    real, intent(in) :: MomentumSi    ! Particle momentum in Si units
    real, intent(in) :: XyzSi_D(3)    ! Where to set BC
    ! OUTPUT
    real, intent(out):: DistTimesP2Si ! Distribution function*p**2, Si
    ! OPTIONAL
    real, OPTIONAL, intent(in) :: A   ! To convert to energy per nuc
    real, OPTIONAL, intent(in) :: Z   ! To convert to rigidity
    ! LOCAL VARS
    real :: FakeMomentumArr_I(1), FakeDistTimesP2_I(1)
    !--------------------------------------------------------------------------
    FakeMomentumArr_I(1) = MomentumSi
    call local_interstellar_spectrum_a(1,FakeMomentumArr_I, XyzSi_D, & ! Input
         FakeDistTimesP2_I,                                          & ! Output
         A, Z)                                                       ! Optional
    DistTimesP2Si = FakeDistTimesP2_I(1)
  end subroutine local_interstellar_spectrum_s
  !============================================================================
  subroutine local_interstellar_spectrum_a(nP,MomentumSi_I, XyzSi_D, & ! Input
       DistTimesP2Si_I,                                              & ! Output
       A, Z)                                                         ! Optional
    ! Formulae for LISM flux by Burger et al. 2000 (doi: 10.1029/2000JA000153)
    ! and modified by Usoskin et al. 2005 (doi: 10.1029/2005JA011250)
    ! Here we refer to Eq. (2) in Usoskin et al. 2005

    use ModConst, ONLY: cGEV, cLightSpeed, cProtonMass, &
         cElectronCharge, cRMEProton
    ! INPUTS
    integer, intent(in) :: nP  ! The number of meshes in the logarithmic grid
    real, intent(in) :: MomentumSi_I(1:nP)   ! Particle momenta in Si units
    real, intent(in) :: XyzSi_D(3)           ! Where to set BC
    ! OUTPUT
    real, intent(out):: DistTimesP2Si_I(1:nP) ! Distribution function*p**2, Si
    ! OPTIONAL
    real, OPTIONAL, intent(in) :: A          ! To convert to energy per nuc
    real, OPTIONAL, intent(in) :: Z          ! To convert to rigidity
    ! LOCAL VARS
    real :: EnergySi_I(1:nP), RigiditySi_I(1:nP) ! E_K and R(E_K) in SI units
    real :: RigidityGV_I(1:nP)                   ! R(E_K) in the unit of GV
    real :: DistTimesP2Gn_I(1:nP) ! Distribution function*p**2, [.../(GeV/nuc)]
    real :: Ai, Zi                ! Default values for A and Z if not present
    !--------------------------------------------------------------------------
    Ai = 1.0; Zi = 1.0
    if(present(A)) Ai = A
    if(present(Z)) Zi = Z

    ! p[Si] -> E_K[Si]
    EnergySi_I = sqrt((MomentumSi_I * cLightSpeed)**2 + &
         (Ai * cProtonMass * cLightSpeed**2)**2)
    ! E_K[Si] -> R(E_K)[Si]
    RigiditySi_I = Ai / abs(Zi*cElectronCharge) *       &
         sqrt(EnergySi_I*(EnergySi_I+2*cRMEProton))
    ! R(E_K)[Si] -> R(E_K)[GV]
    RigidityGV_I = RigiditySi_I / cGEV
    ! R(E_K)[GV] -> j_LIS[.../(GeV/nuc)]
    DistTimesP2Gn_I = 1.9E+4 * RigidityGV_I**(-2.78) /  &
         (1.0 + 0.4866 * RigidityGV_I**(-2.51))
    ! j_LIS[.../(GeV/nuc)] -> p^2*f|_{infty}[Si]
    DistTimesP2Si_I = DistTimesP2Gn_I/(Ai * cGEV)
  end subroutine local_interstellar_spectrum_a
  !============================================================================
end module ModTransitionRegion
!==============================================================================
