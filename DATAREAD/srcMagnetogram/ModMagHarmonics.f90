!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModMagHarmonics

  use ModNumConst
  use ModUtilities, ONLY: CON_stop
  use ModReadMagnetogram, ONLY: dPhi, dTheta, dSinTheta

  implicit none  

  logical:: UseSinLatitudeGrid = .true.
  logical:: UseChebyshevNode   = .true.
  ! **********************Choice of this parameter**********************
  ! *Sin(Latitude): WSO : http://wso.stanford.edu                      *
  ! *   MDI: see http://soi.stanford.edu/magnetic/index6.html          * 
  ! *   SOLIS:http://solis.nso.edu/vsm/map_info.html                   *
  ! *   GONG: http://gong.nso.edu/data/magmap/index.html               *
  ! ********************************************************************

  ! Name of output file
  character (len=100):: NameFileOut='harmonics.dat'

  ! This Module reads a raw (RADIAL, not LOS!!!) magnetogram data file and
  ! generates a magnetogram file in the form of spherical
  ! harmonics to be used by the SWMF. 

  ! ************************ Data Links ********************************
  ! * MDI:   http://soi.stanford.edu/magnetic/index6.html              *
  ! * WSO:   http://wso.stanford.edu/forms/prgs.html                   *  
  ! * GONG:  http://gong.nso.edu/data/magmap/QR/mqs/                   *
  ! * SOLIS: ftp://solis.nso.edu/synoptic/level3/vsm/merged/carr-rot   *
  ! * MWO:   ftp://howard.astro.ucla.edu/pub/obs/synoptic_charts       *
  ! ********************************************************************
  ! * Field in Gauss: MDI,GONG,SOLIS                                   *
  ! * Field in microTesla(0.01Gs): WSO, MWO                            *
  ! ********************************************************************


  integer:: i,n,m,iTheta,iPhi,iR,mm,nn
  real, allocatable, dimension(:,:):: g_nm, h_nm, Br_II
  real:: dR=1.0 !, dPhi=1.0, dTheta, dSinTheta=1.0
  integer:: nPhi=72, nTheta=29

  integer, parameter:: MaxHarmonics = 180
  integer:: nHarmonics = MaxHarmonics, nHarmonicsIn = MaxHarmonics

  real:: CosTheta,SinTheta
  real:: stuff1,stuff2,stuff3
  real:: Theta,Phi

  real :: SumArea,da
  real :: NormalizationFactor
  integer :: iUnit2, iNM
  real,allocatable,dimension(:) :: gArray, hArray
  integer,allocatable,dimension(:) :: nArray, mArray
  integer :: SizeOfnm, ArrPerProc,EndProc,SizeLastProc
  real, allocatable, dimension(:,:)   :: p_nm, CosMPhi_II,SinMPhi_II
  real, allocatable, dimension(:,:,:) :: PNMTheta_III

  real:: SinThetaM, SinThetaM1 
  integer:: delta_m0
  real, allocatable:: FactRatio1(:)
  integer, parameter:: MaxInt=100000
  real:: Sqrt_I(MaxInt)

  real, allocatable, dimension(:) :: ChebyshevWeightE_I, ChebyshevWeightW_I

contains

  subroutine read_harmonics_param

    use ModReadParam
    use ModReadMagnetogram, ONLY: read_magnetogram_param

    character(len=lStringLine) :: NameCommand

    character(len=*), parameter:: NameSub = 'read_harmonics_param'
    !-----------------------------------------------------------------------
    call read_file('HARMONICS.in')
    call read_init
    call read_echo_set(.true.)

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#HARMONICS")
          call read_var('nHarmonics', nHarmonicsIn)
       case("#OUTPUT")
          call read_var('NameFileOut', NameFileOut)
       case("#MAGNETOGRAMFILE", "#CHEBYSHEV", '#CHANGEWEAKFIELD', '#CHANGEPOLARFIELD')
          call read_magnetogram_param(NameCommand)
       case default
          call CON_stop(NameSub//': unknown command='//trim(NameCommand))
       end select
    end do

  end subroutine read_harmonics_param
  !=================================================================
  subroutine read_modified_magnetogram

    use ModReadMagnetogram, ONLY: read_orig_magnetogram, nTheta0, nPhi0, &
         nThetaorig, nPhiorig

    character(len=*), parameter:: NameSub = 'read_modified_magnetogram'
    !----------------------------------------------------------------
    call read_orig_magnetogram

    nTheta = nTheta0
    nPhi   = nPhi0
    write(*,*)' Original nTheta, nPhi =', nThetaorig, nPhiorig
    write(*,*)' Remeshed nTheta, nPhi =', nTheta0, nPhi0

    ! Setting the order on harmonics to be equal to the 
    ! latitudinal resolution.
    if(nHarmonicsIn > 0 .and. nHarmonicsIn < MaxHarmonics)then
       nHarmonics = nHarmonicsIn
    else
       nHarmonics = min(nThetaorig, MaxHarmonics)
    endif
    write(*,*)'Order of harmonics: ',nHarmonics

    ! Allocate the harmonic coefficients arrays
    allocate( &
         p_nm(nHarmonics+1,nHarmonics+1), &
         g_nm(nHarmonics+1,nHarmonics+1), &
         h_nm(nHarmonics+1,nHarmonics+1), &
         FactRatio1(nHarmonics+1))

    p_nm = 0.0
    g_nm = 0.0
    h_nm = 0.0

  end subroutine read_modified_magnetogram
  !============================================================================
  real function sin_latitude(iTheta)
    integer,intent(in)::iTheta

    sin_latitude = (real(iTheta)+0.5)*dSinTheta-1.0

  end function sin_latitude
  !=================================================================
  real function r_latitude(iTheta)
    integer,intent(in)::iTheta
    !-------------------------------------------------------------- 
    if(UseSinLatitudeGrid)then
       r_latitude = asin(sin_latitude(iTheta))
    else
       r_latitude = (iTheta + 0.50)*dTheta - cPi*0.50
    end if
  end function r_latitude
  !================================================================= 
  real function colatitude(iTheta)
    integer,intent(in)::iTheta
    colatitude = cPi*0.5 - r_latitude(iTheta)
  end function colatitude
  !=================================================================
  subroutine calc_harmonics

    ! This subroutine calculates the spherical harmonics from the
    !  modified - remesed magnetogram data

    use ModReadMagnetogram, ONLY: Br0_II, nTheta0, nPhi0, UseChebyshevNode, &
         iCarringtonRotation, iLong0, UseCosTheta, nThetaorig, nPhiorig, &
         ChebyshevWeightE_I, ChebyshevWeightW_I, dPhi, dTheta, dSinTheta

    integer :: iUnit, iError, iTheta, iPhi, m,inm, nn, mm, n 
    real    :: dThetaChebyshev
    real, allocatable:: Br_II(:,:)

    !-------------------------------------------------------------------------
    write(*,*)'Calculating harmonic coefficients'

    if(.not.UseCosTheta) UseSinLatitudeGrid = .false.

    ! Calculate sqrt(integer) from 1 to 10000::
    do m=1,MaxInt
       Sqrt_I(m) = sqrt(real(m))
    end do

    ! Calculate the ratio sqrt(2m!)/(2^m*m!)::
    factRatio1(:) = 0.0; factRatio1(1) = 1.0
    do m=1,nHarmonics
       factRatio1(m+1) = factRatio1(m)*Sqrt_I(2*m-1)/Sqrt_I(2*m)
    enddo

    ! If Chebyshev is done then remeshed theta & phi are used, else 
    ! original values are used.
    if(UseChebyshevNode) then
       nTheta = nTheta0
       nPhi   = nPhi0
    else 
       nTheta = nThetaorig
       nPhi   = nPhiorig
    endif

    allocate(Br_II(0:nPhi-1,0:nTheta-1))
    Br_II=Br0_II

    iunit=9
    open ( unit = iUnit, &
         file = 'br.dat', &
         form = 'formatted', &
         access = 'sequential')

    do iTheta=0, nTheta-1
       do iPhi=0, nPhi-1
          write(iUnit,'(I3,I3,2f20.10)') iTheta, iPhi, Br_II(iPhi,iTheta)
       enddo
    enddo
    close(iUnit)

    if(UseChebyshevNode) then
       allocate(PNMTheta_III(nHarmonics+1,nHarmonics+1,0:nTheta-1))
       PNMTheta_III = 0.0
       dThetaChebyshev = cPi/(nTheta-1)
       do iTheta=0,nTheta-1
          Theta=cPi-iTheta*dThetaChebyshev
          !write(*,*) Theta
          CosTheta=cos(Theta)
          SinTheta=max(sin(Theta), 1E-9)
          call calc_legendre_polynoms
          PNMTheta_III(:,:,iTheta) = p_nm
       end do
    else
       allocate(PNMTheta_III(nHarmonics+1,nHarmonics+1,0:nTheta-1))
       PNMTheta_III = 0.0
       do iTheta=0,nTheta-1
          Theta=colatitude(iTheta)
          !write(*,*) Theta
          CosTheta=cos(Theta)
          SinTheta=max(sin(Theta), 1E-9)
          call calc_legendre_polynoms
          PNMTheta_III(:,:,iTheta) = p_nm
       end do
    end if
    allocate(CosMPhi_II(0:nPhi-1,0:nHarmonics),&
         SinMPhi_II(0:nPhi-1,0:nHarmonics) )
    !Save sins-cosins

    do iPhi=0,nPhi-1
       Phi=(iPhi)*dPhi
       do m=0,nHarmonics
          CosMPhi_II(iPhi,m) = cos(m*Phi)
          SinMPhi_II(iPhi,m) = sin(m*Phi)
       end do
    end do
    ! Creating an array with the size of SizeOfnm (the total g_nm and h_nm) for
    ! parallel calculation of the coefficients
    SizeOfnm=1
    do iNM=1,nHarmonics
       SizeOfnm=SizeOfnm+(iNM+1)
    enddo

    ! Allocating and initializing arrays
    if(allocated(gArray))deallocate(gArray)
    allocate(gArray(0:SizeOfnm-1))
    if(allocated(hArray))deallocate(hArray)
    allocate(hArray(0:SizeOfnm-1))
    if(allocated(nArray))deallocate(nArray)
    allocate(nArray(0:SizeOfnm-1))
    if(allocated(mArray))deallocate(mArray)
    allocate(mArray(0:SizeOfnm-1))

    gArray=0.0; hArray=0.0
    nArray=0.0; mArray=0.0

    ! Create arrays for n and m for parallelization
    iNM=0
    do nn=0,nHarmonics
       iNM=iNM+nn
       do mm=0,nn
          nArray(iNM+mm)=nn
          mArray(iNM+mm)=mm
       enddo
    end do

    ! Each processor gets part of the array 
    do iNM = 0,SizeOfnm-1

       ! The proper normalization factor is (2n+1)/R_n, where R_n=n+1+n(1/Rs)**(2n+1).
       ! However, in this code the coefficients are normalized only with 2n+1 to reproduce 
       ! the coefficients provided by Stanford. The division by R_n is done after
       ! the coefficients are been read in ModMagnetogram.f90.
       ! R_n=(nArray(iNM)+1.0)+nArray(iNM)*(1.0/Rs_PFSSM)**(2*nArray(iNM)+1)
       NormalizationFactor=(2*nArray(iNM)+1)

       SumArea=0.0

       do iTheta=0,nTheta-1

          if (UseChebyshevNode) then
             ! Use Chebyshev Weight in theta direction
             da = ChebyshevWeightE_I(iTheta)*ChebyshevWeightW_I(iTheta)*dPhi
          else
             Theta = colatitude(iTheta) 
             SinTheta = max(sin(Theta), 0.0)
             if(UseCosTheta)then
                da = dSinTheta*dPhi
             else
                da = SinTheta*dTheta*dPhi
             end if
          end if

          ! Calculate the set of Legendre polynoms (with Schmidt normalization)
          ! for a given CosTheta, SinTheta.
          ! For non-radial magnetogram (LOS), a division in SinTheta is needed.

          gArray(iNM) = gArray(iNM)+&
               sum(Br_II(:,iTheta)*&
               CosMPhi_II(:,mArray(iNM)) )*da*PNMTheta_III(&
               nArray(iNM)+1,mArray(iNM)+1,iTheta)!/SinTheta
          hArray(iNM) = hArray(iNM)+&
               sum(Br_II(:,iTheta)*&
               SinMPhi_II(:,mArray(iNM)) )*da*PNMTheta_III(&
               nArray(iNM)+1,mArray(iNM)+1,iTheta)!/SinTheta
          SumArea = SumArea+da*nPhi
       end do
       gArray(iNM) = NormalizationFactor*gArray(iNM)/SumArea
       hArray(iNM) = NormalizationFactor*hArray(iNM)/SumArea
    end do

    iNM=0
    do nn=0,nHarmonics
       iNM=iNM+nn
       do mm=0,nn
          g_nm(nn+1,mm+1) = gArray(iNM+mm)
          h_nm(nn+1,mm+1) = hArray(iNM+mm) 
       enddo
    end do

    !\
    ! Leave out monopole (n=0) term::
    !/
    g_nm(1,1) = 0.0

    write(*,*)'Done Calculating harmonic coefficients' 

    !\
    ! Writing spherical harmonics file
    !/

    write(*,*)'Writing harmonic coefficients file, named ',NameFileOut
    iUnit = 9
    open ( unit = iUnit, &
         file = NameFileOut, &
         form = 'formatted', &
         access = 'sequential', &
         status = 'replace', iostat = iError )

    if ( iError /= 0 ) then
       write (*,*)' Could not output file ', NameFileOut
       stop
    end if

    write ( iUnit, '(a19,I3,a10,I4.4,a1,i3.3)' ) &
         'Coefficients order=',nHarmonics, &
         ' Central=CT',iCarringtonRotation,':',mod(iLong0+180,360)
    write ( iUnit, '(a)' ) 'Observation time'
    write ( iUnit, '(a45,I3)' ) &
         'B0 angle & Nmax:        0          ',nHarmonics
    write ( iUnit, * )
    write ( iUnit, * )
    write ( iUnit, * )
    write ( iUnit, * )
    write ( iUnit, * )
    write ( iUnit, '(a19,I3,a26)' ) &
         'Max Harmonic Order:',nHarmonics,' Units: Gauss'
    write ( iUnit, '(a)' ) ' '
    write ( iUnit, '(a)' ) '  l   m      g(uT)      h(uT)'
    write ( iUnit, '(a)' ) ' '


    do nn=0, nHarmonics
       do mm = 0, nn
          write(iUnit, '(2I5,2f20.10)') nn,mm,g_nm(nn+1,mm+1),h_nm(nn+1,mm+1)
       enddo
    end do

    close(iUnit)

  end subroutine calc_harmonics

  !===========================================================================

  subroutine calc_legendre_polynoms

    ! Calculate the Legendre polynoms for a particular latitude

    ! Calculate polynomials with appropriate normalization
    ! for Theta_PFSSMa::

    SinThetaM  = 1.0
    SinThetaM1 = 1.0
    p_nm(:,:)  = 0.0

    do m=0,nHarmonics
       if (m == 0) then
          delta_m0 = 1
       else
          delta_m0 = 0
       endif
       !\
       ! Eq.(27) from Altschuler et al. 1976::
       !/
       p_nm(m+1,m+1) = factRatio1(m+1)*Sqrt_I((2-delta_m0)*(2*m+1))* &
            SinThetaM
       !\
       ! Eq.(28) from Altschuler et al. 1976::
       !/
       if (m < nHarmonics) p_nm(m+2,m+1) = p_nm(m+1,m+1)*Sqrt_I(2*m+3)* &
            CosTheta

       SinThetaM1 = SinThetaM
       SinThetaM  = SinThetaM*SinTheta

    enddo
    do m=0,nHarmonics-2; do n=m+2,nHarmonics
       !\
       ! Eq.(29) from Altschuler et al. 1976::
       !/
       stuff1         = Sqrt_I(2*n+1)/Sqrt_I(n**2-m**2)
       stuff2         = Sqrt_I(2*n-1)
       stuff3         = Sqrt_I((n-1)**2-m**2)/Sqrt_I(2*n-3)
       p_nm(n+1,m+1)  = stuff1*(stuff2*CosTheta*p_nm(n,m+1)-  &
            stuff3*p_nm(n-1,m+1))

    enddo; enddo
    !\
    ! Apply Schmidt normalization::
    !/
    do m=0,nHarmonics; do n=m,nHarmonics
       !\
       ! Eq.(33) from Altschuler et al. 1976::
       !/
       stuff1 = 1.0/Sqrt_I(2*n+1)
       !\
       ! Eq.(34) from Altschuler et al. 1976::
       !/
       p_nm(n+1,m+1)  = p_nm(n+1,m+1)*stuff1
    enddo; enddo

  end subroutine calc_legendre_polynoms

end module ModMagHarmonics
