! Merging FDIPS_module.f90 & ModMagHarmonics.f90 to read the raw
! magnetogram & modifying the field.

module ModReadMagnetogram

  use ModNumConst
  use ModUtilities, ONLY: CON_stop
  use ModConst, ONLY: cPi, cTwoPi, cHalfPi

  implicit none

  ! Input parameters
  character (len=100):: NameFileIn = 'fitsfile.out'
  real:: BrMax = 3500.0

  ! Optional enhancement of the polar magnetic field with a factor
  !  1 + (PolarFactor-1)*abs(sin(Latitude))^PolarExponent
  logical           :: DoChangePolarField = .false.
  real              :: PolarFactor = 1.0
  real              :: PolarExponent = 2.0

  ! Parameters used to apply modified scaling to the raw magnetogram
  ! Apply a scaling factor to small magnetic fields, to compensate
  ! the measurement error for the coronal hole (=very low) field
  logical :: DoChangeWeakField = .false.
  real    :: BrFactor = 1.0, BrMin = 0.0

  ! To check the latitude grid
  logical:: UseCosTheta = .true.
  ! Parameters for uniformTheta_Transform
  logical:: UseChebyshevNode = .false.
  real, allocatable:: ChebyshevWeightE_I(:), ChebyshevWeightW_I(:)

  ! For Reading the original Magnetogram & Modifying
  real, allocatable   :: Br0_II(:,:)
  real, allocatable   :: Phi0_I(:), Lat0_I(:)
  integer :: nTheta0, nPhi0, nThetaorig, nPhiorig, nThetaAll, nPhiAll

  ! To indicate if fdips is done or harmonics
  logical :: IsFdips = .false.

  ! Carrington rotation #. If the map covers parts of two rotations,
  ! the carrington rotation for the central meridian is provided
  integer:: iCarringtonRotation = 0
  ! Garrington longitude of the left margin of the map
  ! ("leading longitude")
  integer:: iLong0 = 0

  ! For the associated functions
  real:: dR=1.0, dPhi=1.0, dTheta=1.0, dSinTheta=1.0

contains

  subroutine read_magnetogram_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand
    character(len=*), parameter:: NameSub = 'read_magnetogram_param'
    !-------------------------------------------------------------------------
    select case(NameCommand)
    case("#MAGNETOGRAMFILE")
       call read_var('NameFileIn', NameFileIn)
       call read_var('BrMax',      BrMax)
    case("#CHANGEPOLARFIELD")
       DoChangePolarField = .true.
       call read_var('PolarFactor',   PolarFactor)
       call read_var('PolarExponent', PolarExponent) 
    case("#CHANGEWEAKFIELD")
       DoChangeWeakfield = .true.
       call read_var('BrFactor', BrFactor)
       call read_var('BrMin', BrMin)
    case("#CHEBYSHEV")
       call read_var('UseChebyshevNode', UseChebyshevNode)
    case default
       call CON_stop(NameSub//': unknown command='//trim(NameCommand))
    end select

  end subroutine read_magnetogram_param

  !============================================================================
  subroutine read_orig_magnetogram(IsPhiThetaOrder, UseWedge, DoRemoveMonopole, &
       nThetaCoarse, nPhiCoarse)

    use ModPlotFile, ONLY: read_plot_file, save_plot_file
    use ModNumConst, ONLY: cDegToRad

    logical, optional,  intent(in):: IsPhiThetaOrder 
    logical, optional,  intent(in):: UseWedge
    logical, optional,  intent(in):: DoRemoveMonopole
    integer, optional,  intent(in):: nThetaCoarse, nPhiCoarse

    real:: Param_I(2)
    real:: BrAverage, weight
    real, allocatable :: BrTmp_II(:,:), BrTrans_II(:,:)

    integer:: iError, nParam, iTheta, iPhi, nThetaRatio, nPhiRatio
    integer:: iTheta0, iTheta1, jPhi0, jPhi1, jPhi, kPhi

    logical :: IsInputLatReverse = .false.
    logical :: RemoveMonopole = .true.

    character(len=200) :: NameVarOut

    character(len=*), parameter:: NameSub = 'read_orig_magnetogram'
    !------------------------------------------------------------------------
    call read_plot_file(NameFileIn, n1Out = nPhi0, n2Out = nTheta0, &
         ParamOut_I=Param_I, iErrorOut=iError, nParamOut=nParam)

    if(iError /= 0) call CON_stop(NameSub// &
         ': could not read header from file '//trim(NameFileIn))

    if(nParam>0)iLong0 = nint(Param_I(1)-0.50*360.0/nPhi0)
    if(nParam>1)iCarringtonRotation = nint(Param_I(2))

    write(*,*)'nTheta0, nPhi0, LongitudeShift: ', nTheta0, nPhi0, iLong0

    allocate(Phi0_I(nPhi0), Lat0_I(nTheta0)) 
    allocate(Br0_II(nPhi0,nTheta0))

    call read_plot_file(NameFileIn, &
         Coord1Out_I=Phi0_I, Coord2Out_I=Lat0_I, VarOut_II = Br0_II, &
         iErrorOut=iError, NameVarOut=NameVarOut)

    if(iError /= 0) call CON_stop(NameSub// &
         ': could not read data from file'//trim(NameFileIn))

    ! Check if the latitude coordinate is uniform or not
    ! There is no point using Chebyshev transform if the original grid                   
    ! is already uniform in theta                                                        
    UseCosTheta = abs(Lat0_I(3) - 2*Lat0_I(2) + Lat0_I(1)) > 1e-6
    if(.not.UseCosTheta) UseChebyshevNode = .false.

    IsInputLatReverse = Lat0_I(1) > Lat0_I(nTheta0)

    ! For #CHANGEPOLARFIELD
    if(DoChangePolarField)then
       do iTheta = 1, nTheta0
          Br0_II(:,iTheta) = Br0_II(:,iTheta) * (1 + &
               (PolarFactor-1)*abs(sin(cDegToRad*Lat0_I(iTheta)))**PolarExponent)
       end do
    end if

    ! options not supported with UseWedge
    if(present(UseWedge))then
       if(DoChangePolarField .and. UseWedge)then
          call CON_stop('UseWedge currently does not support DoChangePolarField.')
       endif
       if(UseCosTheta .and. UseWedge)then
          call CON_stop(NameSub//&
               ': Currently UseWedge only works with uniform latitude grid')
       endif
       if(present(DoRemoveMonopole))then
          if(DoRemoveMonopole .and. UseWedge)then
             call CON_stop(NameSub//&
                  ': Currently UseWedge does not work with DoRemoveMonopole')
          elseif (.not. DoRemoveMonopole .and. .not. UseWedge)then
             RemoveMonopole = .false.
          endif
       endif
    endif

    if(present(DoRemoveMonopole))then
       if(.not.DoRemoveMonopole)&
            RemoveMonopole = .false.
    endif

    ! For #CHANGEWEAKFIELD
    if(DoChangeWeakField)then
       if(BrMin > 0.0 .or. BrFactor > 1.0) &
            Br0_II = sign(min(abs(Br0_II) + BrMin, BrFactor*abs(Br0_II)), Br0_II)
    endif

    ! Fix too large values of Br
    where (abs(Br0_II) > BrMax) Br0_II = sign(BrMax, Br0_II)

    ! Done for both Harmonics & Fdips unless specified otherwise
    ! for FDIPS when wedge is used
    if(RemoveMonopole)then
       if (UseCosTheta) then
          BrAverage = sum(Br0_II)/(nTheta0*nPhi0)
       else
          BrAverage = 0.0
          do iTheta = 1, nTheta0
             BrAverage = BrAverage + &
                  sum(Br0_II(:,nTheta0+1-iTheta)) * cos(cDegToRad*Lat0_I(nTheta0+1-iTheta))
          enddo
          BrAverage = BrAverage / (nTheta0*nPhi0)
       endif
       Br0_II = Br0_II - BrAverage
    endif

    ! Save 2D file
    call save_plot_file('field_2d.out', TypeFileIn='real4',&
         StringHeaderIn='Longitude, Latitude [Deg], Br[G]', &
         NameVarIn='Longitude Latitude Br',&
         VarIn_II=Br0_II, &
         Coord1In_I=Phi0_I, &
         Coord2In_I=Lat0_I)

    ! to be passed to harmonics before remeshing
    nThetaorig = nTheta0
    nPhiorig   = nPhi0

    dPhi      = cTwoPi/nPhi0
    dTheta    = cPi/nTheta0
    dSinTheta = 2.0/nTheta0

    if(present(nThetaCoarse) .and. present(nPhiCoarse))then
       IsFdips = .true.
    else if(UseChebyshevNode) then
       write(*,*)'Transforming to Uniform Theta for Harmonics'
       call uniformTheta_transform
       write(*,*)"Transformed to Uniform Theta grid for Harmonics"
    endif

    ! To be done for FDIPS only
    ! Step 1 : Coarsening for FDIPS based on Input Theta & Phi
    if(IsFdips)then
       if(nThetaorig > nThetaCoarse .and. nThetaCoarse > 1)then
          if(modulo(nThetaorig, nThetaCoarse) /= 0)then
             write(*,*) NameSub,' nTheta in file    =', nThetaorig
             write(*,*) NameSub,' nTheta in FDIPS.in=', nThetaCoarse
             call CON_stop(NameSub//': not an integer coarsening ratio')
          end if
          ! Set integer coarsening ratio
          nThetaRatio = nThetaorig / nThetaCoarse
          nThetaAll   = nThetaorig / nThetaRatio
       else
          nThetaRatio = 1
          nThetaAll   = nThetaorig
       end if

       if(nPhiorig > nPhiCoarse .and. nPhiCoarse > 1)then
          if(modulo(nPhiorig, nPhiCoarse) /= 0)then
             write(*,*) NameSub,' nPhi in file    =', nPhiorig
             write(*,*) NameSub,' nPhi in FDIPS.in=', nPhiCoarse
             call CON_stop(NameSub//': not an integer coarsening ratio')
          end if
          nPhiRatio = nPhiorig / nPhiCoarse
          nPhiAll   = nPhiorig / nPhiRatio
       else
          nPhiRatio = 1
          nPhiAll   = nPhiorig
       end if

       allocate(BrTmp_II(nPhiAll,nThetaAll))
       BrTmp_II=0.0

       do iPhi = 1, nPhiAll
          if (.not. UseWedge) then
             jPhi0 = nPhiRatio*(iPhi-1) - nPhiRatio/2 + 1
             jPhi1 = nPhiRatio*(iPhi-1) + nPhiRatio/2 + 1
          else
             ! phi binning is different for UseWedge because no periodicity
             jPhi0 = nPhiRatio*(iPhi-1) + 1
             jPhi1 = jPhi0 + nPhiRatio -1
          endif

          do jPhi = jPhi0, jPhi1
             if( .not. UseWedge .and. modulo(nPhiRatio,2) == 0 .and. &
                  (jPhi == jPhi0 .or. jPhi == jPhi1) )then
                ! For even coarsening ratio use 0.5 weight at the two ends
                Weight = 0.5
             else
                Weight = 1.0
             end if

             if (.not. UseWedge) then
                ! Apply periodicity
                kPhi = modulo(jPhi-1,nPhi0) + 1
             else
                kPhi = jPhi
             endif

             do iTheta = 1, nThetaAll
                iTheta0 = nThetaRatio*(iTheta-1) + 1
                iTheta1 = iTheta0 + nThetaRatio - 1

                BrTmp_II(iPhi,iTheta) = BrTmp_II(iPhi,iTheta) &
                     + Weight * sum( Br0_II(kPhi,iTheta0:iTheta1))
             end do
          end do
       end do
       BrTmp_II = BrTmp_II / (nThetaRatio*nPhiRatio)
       if(allocated(Br0_II)) deallocate(Br0_II)
       allocate(Br0_II(nPhiAll,nThetaAll))
       Br0_II = BrTmp_II
       deallocate(BrTmp_II)

       ! Br0_II is phi,theta till here and not reverse
       if(UseChebyshevNode)then
          call uniformTheta_transform
          nThetaAll = nTheta0
          nPhiAll = nPhi0
       endif
       ! For FDIPS ONLY, transpose
       allocate(BrTrans_II(nThetaAll,nPhiAll))
       if(present(IsPhiThetaOrder)) then
          if(.not. IsPhiThetaOrder)then
             if (IsInputLatReverse) then
                BrTrans_II = transpose(Br0_II)
             else
                BrTrans_II = transpose(Br0_II(:,nThetaAll:1:-1))
             endif
          endif
       endif
       if(allocated(Br0_II)) deallocate(Br0_II)
       allocate(Br0_II(nThetaAll,nPhiAll))
       Br0_II = BrTrans_II
       ! Br0_II is now (Theta, Phi) not lon,Lat
       deallocate(BrTrans_II)
    endif

  end subroutine read_orig_magnetogram
  !============================================================================== 
  real function sin_latitude(iTheta)
    integer,intent(in)::iTheta

    sin_latitude = (real(iTheta)+0.5)*dSinTheta-1.0

  end function sin_latitude
  !==============================================================================
  real function r_latitude(iTheta)
    integer,intent(in)::iTheta
    !----------------------------------------------------------------------------
    if(UseCosTheta)then
       r_latitude = asin(sin_latitude(iTheta))
    else
       r_latitude = (iTheta + 0.50)*dTheta - cPi*0.50
    end if
  end function r_latitude
  !==============================================================================
  real function colatitude(iTheta)
    integer,intent(in)::iTheta
    colatitude = cPi*0.5 - r_latitude(iTheta)
  end function colatitude
  !==============================================================================
  subroutine uniformTheta_transform

    ! The magnetogram is remeshed onto a uniform in theta (co-latitude) grid.
    ! This corresponds to the Chebyshev points with respect to cos theta,
    ! which is the argument of the associated Legendre polynomials.
    ! The equal spacing in the theta direction makes the calculation much
    ! more accurate, especially in the polar regions.
    ! The number of points in the theta direction of the remeshed grid
    ! is always an odd number.

    integer :: nThetaIn, nThetaOut, iw, iu, iTheta, iPhi
    integer :: iLower, iUpper
    real    :: dThetaChebyshev, dThetaInterpolate, BrSouthPole, BrNorthPole
    real, allocatable :: ThetaIn_I(:), ThetaOut_I(:), NewBr_II(:,:), &
         ChebyshevWeightEu_I(:), Br_II(:,:)
    !--------------------------------------------------------------------------
    write(*,*) ' Remeshing to Uniform Theta'

    if(IsFdips)then
       nTheta0 = nThetaAll
       nPhi0 = nPhiAll
       write(*,*)'ISFdips, nTheta0, nPhi0',IsFdips, nTheta0, nPhi0
    endif

    nThetaIn = nTheta0
    nThetaOut = ceiling(nThetaIn * cPi/2)

    write(*,*)'nTHeta0,nTHetaIn,nTHetaOut',nTheta0,nThetaIn,nThetaOut

    !Notice the number of point in theta direction is an odd number
    if (mod(nThetaOut,2) == 0) then
       nThetaOut=nThetaOut+1
    else
       nThetaOut=nThetaOut
    end if

    write(*,*) 'Original nTheta=', nThetaIn
    write(*,*) 'New nTheta=     ', nThetaOut

    ! Note that the indices start from 0
    ! for both FDIPS & Harmonics
    ! Lat grid to Theta grid 
    allocate(NewBr_II(0:nPhi0-1,0:nThetaOut-1))
    allocate(ThetaIn_I(0:nThetaIn-1))
    allocate(ThetaOut_I(0:nThetaOut-1))

    do iTheta=0,nThetaIn-1
       ThetaIn_I(iTheta) = colatitude(iTheta)
    end do

    dThetaChebyshev = cPi/(nThetaOut-1)

    ! +90 is 0 deg Theta
    do iTheta = 0, nThetaOut-1
       ThetaOut_I(iTheta) = cPi - iTheta*dThetaChebyshev
    end do

    iLower = -1
    iUpper = -1
    NewBr_II = 0

    ! Allocate a new Magnetic field variable which has 
    ! first index as 0
    ! common for both FDIPS & Harmonics
    allocate(Br_II(0:nPhi0-1,0:nTheta0-1))
    Br_II = Br0_II

    ! Still in PhiTheta order,
    ! Assumes that the Lat(0) is -90 (south pole)
    BrNorthPole = sum(Br_II(:,nThetaIn-1))/nPhi0
    BrSouthPole = sum(Br_II(:,0))/nPhi0

    write(*,*)'BrNorth,BrSouth',BrNorthPole,BrSouthPole
    ! Use linear interpolation to do the data remesh
    do iPhi=0,nPhi0-1
       do iTheta=0, nThetaOut-1
          ! A search is needed in case the sin(latitude) grid is not uniform
          !do iTheta_search=0,nThetaIn-2
          !   if(ThetaOut_I(iTheta) <= ThetaIn_I(iTheta_search) .and. &
          !      ThetaOut_I(iTheta) >= ThetaIn_I(iTheta_search+1)) then
          !      iLower=iTheta_search+1
          !      iUpper=iTheta_search
          !      exit
          !   else
          !      iLower=-1
          !      iUpper=-1
          !   end if
          !end do
          iUpper = &
               floor((cos(ThetaOut_I(iTheta)) - cos(ThetaIn_I(0)))/dSinTheta)
          iLower = iUpper+1

          if (iUpper /= -1 .and. iUpper /= nThetaIn-1 ) then
             dThetaInterpolate = ThetaIn_I(iUpper) - ThetaIn_I(iLower)
             NewBr_II(iPhi,iTheta) = Br_II(iPhi,iLower)* &
                  (ThetaIn_I(iUpper)-ThetaOut_I(iTheta))/dThetaInterpolate &
                  +Br_II(iPhi,iUpper)* &
                  (ThetaOut_I(iTheta)-ThetaIn_I(iLower))/dThetaInterpolate
          else
             if (iUpper == nThetaIn-1) then
                dThetaInterpolate=ThetaIn_I(nThetaIn-1)
                NewBr_II(iPhi,iTheta) = BrNorthPole & 
                     *(ThetaIn_I(nThetaIn-1)-ThetaOut_I(iTheta)) &
                     /dThetaInterpolate &
                     + Br_II(iPhi,nThetaIn-1) &
                     *(ThetaOut_I(iTheta))/dThetaInterpolate
             end if
             if (iUpper == -1) then
                dThetaInterpolate = cPi-ThetaIn_I(0)
                NewBr_II(iPhi,iTheta) = Br_II(iPhi,0)* &
                     (cPi - ThetaOut_I(iTheta))/dThetaInterpolate &
                     +BrSouthPole* &
                     (ThetaOut_I(iTheta) - ThetaIn_I(0))/dThetaInterpolate
             end if
          end if
       end do
    end do

    ! Copying the remeshed grid after interpolation into the previous theta
    ! Copy the remeshed grid size into nTheta
    nTheta0 = nThetaOut

    ! Copy the remeshed NewBr_II into Br0_II
    ! Br0_II indices begin from 1
    if(allocated(Br0_II)) deallocate(Br0_II)
    allocate(Br0_II(1:nPhi0,1:nTheta0))
    Br0_II = NewBr_II

    ! Calculate the weights for the transformation in theta direction.
    allocate(ChebyshevWeightE_I(0: nThetaOut-1))
    allocate(ChebyshevWeightW_I(0: nThetaOut-1))
    allocate(ChebyshevWeightEu_I(0:(nThetaOut-1)/2))

    ChebyshevWeightW_I    = 0.0
    ChebyshevWeightE_I    = 1.0
    ChebyshevWeightE_I(0) = 0.5
    ChebyshevWeightE_I(nThetaOut-1) = 0.5

    ChebyshevWeightEu_I    = 1.0
    ChebyshevWeightEu_I(0) = 0.5
    ChebyshevWeightEu_I((nThetaOut-1)/2) = 0.5

    do iw=0,nThetaOut-1
       do iu=0,(nThetaOut-1)/2
          ChebyshevWeightW_I(iw) = ChebyshevWeightW_I(iw) + &
               ChebyshevWeightEu_I(iu)*(-2.0)/(4*(iu)**2-1)* & 
               cos(iw*iu*cPi/((nThetaOut-1)/2))
       end do
    end do
    ChebyshevWeightW_I = ChebyshevWeightW_I/(nThetaOut-1)

    deallocate(NewBr_II, ThetaIn_I, ThetaOut_I, ChebyshevWeightEu_I)

  end subroutine uniformTheta_transform

end module ModReadMagnetogram
