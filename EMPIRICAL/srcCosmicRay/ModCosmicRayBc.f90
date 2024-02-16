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
  subroutine local_interstellar_spectrum_a(nP,MomentumSi_I, XyzSi_D, & ! Input
       DistTimesP2Si_I,                                              & ! Output
       A, Z)                                                         ! Optional

    ! INPUTS
    integer, intent(in) :: nP  ! The number of meshes in the logarithmic grid
    real, intent(in) :: MomentumSi_I(1:nP)   ! Particle momenta in Si units
    real, intent(in) :: XyzSi_D(3)           ! Where to set BC
    ! OUTPUT
    real, intent(out):: DistTimesP2Si_I(1:nP) ! Distribution function*p**2, Si
    ! OPTIONAL
    real, OPTIONAL, intent(in) :: A          ! To convert to energy per nuc
    real, OPTIONAL, intent(in) :: Z          ! To convert to rigidity
    !----------------------------
  end subroutine local_interstellar_spectrum_a
  !============================================================================
  subroutine local_interstellar_spectrum_s(MomentumSi, XyzSi_D,    & ! Input
       DistTimesP2Si,                                              & ! Output
       A, Z)                                                         ! Optional

    use ModConst, ONLY: cGEV
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
    !----------------------------
    FakeMomentumArr_I(1) = MomentumSi
    call local_interstellar_spectrum_a(1,FakeMomentumArr_I, XyzSi_D, & ! Input
         FakeDistTimesP2_I,                                          & ! Output
         A, Z)                                                       ! Optional
    DistTimesP2Si = FakeDistTimesP2_I(1)
  end subroutine local_interstellar_spectrum_s
  !============================================================
end module ModTransitionRegion
