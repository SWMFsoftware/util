!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModMagnetogram

  use ModNumConst
  use ModMpi
  use ModIoUnit,   ONLY: io_unit_new
  use CON_axes,    ONLY: dLongitudeHgrDeg
  use ModUtilities, ONLY: CON_stop
  use ModTiming, ONLY: timing_cpu

  implicit none
  save

  private

  public:: get_pfss_field
  
  ! Rs - radius of outer source surface where field is taken to be radial
  ! Ro - radius of inner boundary for the potential field
  ! H  - height above measurements
  
  real, public :: Rs_PFSSM=2.5, Ro_PFSSM=1.0, H_PFSSM=0.0

  ! Units of the magnetic field in the file including corrections
  ! relative to the magnetic units in the SC (Gauss)
  real, public :: UnitB=1.0, UnitBNew=1.0

  ! Global arrays at the magnetogram grid
  integer, parameter, public :: nDim=3, R_=1, Phi_=2, Theta_=3

  ! Time to which the field in B_DN corresponds

  real, public :: dR=1.0, dPhi=1.0, dSinTheta=1.0, dInv_D(nDim)=1.0
  integer, public :: nThetaPerProc, nRExt=2
  integer, public :: nR=29, nPhi=72, nTheta=29

  ! public available procedures
  public :: read_magnetogram_file, read_new_magnetogram_file

  public :: set_parameters_magnetogram

  public :: get_hlcmm
  ! Read H(eliographic) L(ongitude) of the C(entral) M(eridian) of
  ! the M(ap) from the file header. Assign PhiOffset=HLCMM-180

  public :: get_magnetogram_field

  ! Gives the interpolated values of the Cartesian components of
  ! the macnetic vector in HGR system, input parameters
  ! being the cartesian coordinates in the HGR system

  public :: sin_latitude, r_latitude, colatitude,&
       correct_angles, interpolate_field, PhiOffset

  ! read the potential field source surface solution
  public :: update_magnetogram

  ! Local variables -------------

  ! private variables for MPI
  integer :: iProc, nProc, iComm
  
  ! PFSSM related control parameters

  ! Rotation angle around z-axis, in degrees,
  ! from the coordinate system of the component
  ! towards the coordinate system of the magnetic
  ! map in the positive direction - hence, is positive.
  ! Is set automatically to be equal to the
  ! H(eliographic) L(ongitude) of C(entral) M(eridian)
  ! of the M(ap) minus 180 Deg, if in the input file
  ! PhiOffset is negative.
  real :: PhiOffset=-1.0, NewPhiShift = -1.0

  ! Maximum order of spherical harmonics
  integer, parameter:: nHarmonicsMax=180 ! 90

  ! Weights of the spherical harmonics
  real, dimension(nHarmonicsMax+1,nHarmonicsMax+1):: g_II, h_II

  integer:: N_PFSSM=nHarmonicsMax

  ! Number of header lines in the file
  integer:: iHead_PFSSM=12

  ! Name of the input file
  character (LEN=32):: File_PFSSM='mf.dat', NameNewFile = 'newmf.dat'

  ! Name of output directory
  character(len=32):: NameOutDir

  integer, parameter:: MaxInt=100000
  real, allocatable:: Sqrt_I(:), SqrtRatio_I(:)

  real, allocatable:: B_DN(:,:,:,:), Bnew_DN(:,:,:,:)

  ! write prefix for magnetogram processing
  character(len=*), parameter :: prefix=''

contains
  !============================================================================
  real function sin_latitude(iTheta)
    ! Uniform in sin(latitude) grid enumerated by index iTheta
    ! iTheta varies from 0 to nTheta
    ! dSinTheta = 2.0/(nTheta+1)
    integer,intent(in)::iTheta
    !--------------------------------------------------------------------------
    sin_latitude=(real(iTheta)+0.5)*dSinTheta-1.0
  end function sin_latitude
  !============================================================================
  real function r_latitude(iTheta)
    integer,intent(in)::iTheta
    !--------------------------------------------------------------------------
    r_latitude=asin(sin_latitude(iTheta))
  end function r_latitude
  !============================================================================
  real function colatitude(iTheta)
    integer,intent(in)::iTheta
    !--------------------------------------------------------------------------
    colatitude=0.5*cPi-r_latitude(iTheta)
  end function colatitude
  !============================================================================
  ! SUBROUTINE get_hlcmm
  ! Read H(eliographic) L(ongitude) of the C(entral) M(eridian) of
  ! the M(ap) from the header line. Assign PhiOffset=HLCMM-180
  subroutine get_hlcmm(Head_PFSSM,Shift)
    character (LEN=80),intent(inout):: Head_PFSSM
    real,intent(inout)::Shift
    real::HLCMM     ! Heliographic Longitude of the central meridian of map
    integer::iHLCMM ! The same, but HLCMM is integer at WSO magnetograms
    integer::iErrorRead,iPosition

    !--------------------------------------------------------------------------
    iPosition=index(Head_PFSSM,'Centered')
    if (iPosition>0.and.Shift<0.0)then
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       read(Head_PFSSM,'(i3)',iostat=iErrorRead)iHLCMM
       if(iErrorRead>0)call Con_stop(&
            'Cannot find HLCMM, '//File_PFSSM//&
            ' is not a true WSO magnetogram')
       ! Rotates based on magnetogram central meridian + HGR system
       Shift=modulo(iHLCMM-180-dLongitudeHgrDeg, 360.0)
       if(iProc==0)write(*,*) prefix, ' PhiOffset = ',Shift
       RETURN
    end if
    iPosition=index(Head_PFSSM,'Central')
    if(iPosition>0.and.Shift<0.0)then
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       read(Head_PFSSM,*,iostat=iErrorRead)HLCMM
       if(iErrorRead>0)call CON_stop(&
            'Cannot find HLCMM, '//File_PFSSM//&
            ' is not a true MDI magnetogram')
       ! Rotates based on magnetogram central meridian + HGR system
       Shift=modulo(HLCMM-180-dLongitudeHgrDeg, 360.0)
       if(iProc==0)write(*,*) prefix, ' PhiOffset = ',Shift
       RETURN
    end if

    ! Rotation of HGR system is applied
    if (Shift == 0.)then
       Shift = modulo(Shift-dLongitudeHgrDeg, 360.0)
       if(iProc==0)write(*,*) prefix, ' PhiOffset = ',Shift
       RETURN
    end if
  end subroutine get_hlcmm
  !============================================================================

  subroutine set_parameters_magnetogram(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    character(len=*), parameter:: NameSub = 'set_parameters_magnetogram'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case("#MAGNETOGRAM")
       call read_var('rMagnetogram',        Ro_PFSSM)
       call read_var('rSourceSurface',      Rs_PFSSM)
       call read_var('HeightInnerBc',       H_PFSSM)
       call read_var('NameMagnetogramFile', File_PFSSM)
       call read_var('nHeaderLine',         iHead_PFSSM)
       call read_var('PhiShift',            PhiOffset)
       call read_var('UnitB',               UnitB)

    case("#NEWMAGNETOGRAM")
       call read_var('rMagnetogram',        Ro_PFSSM)
       call read_var('rSourceSurface',      Rs_PFSSM)
       call read_var('HeightInnerBc',       H_PFSSM)
       call read_var('NameNewFile',         NameNewFile)
       call read_var('nHeaderLine',         iHead_PFSSM)
       call read_var('NewPhiShift',         NewPhiShift)
       call read_var('UnitB',               UnitB)

    case("#B0GRID")
       call read_var('nR', nR)

    case default
       call CON_stop(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_magnetogram
  !============================================================================

  subroutine read_magnetogram_file(NamePlotDir, iProcIn, nProcIn, iCommIn)
    character(len=*),intent(in) :: NamePlotDir
    integer,intent(in)          :: iProcIn, nProcIn, iCommIn

    real:: timing0
    !--------------------------------------------------------------------------
    iComm = iCommIn
    nProc = nProcIn
    iProc = iProcIn

    if (iProc==0) then
       write(*,*) prefix, 'nOrder = ',N_PFSSM
       write(*,*) prefix, 'Coefficient file name : ',File_PFSSM
       write(*,*) prefix, 'Number of header lines: ',iHead_PFSSM
    endif
    ! Initialize once g(n+1,m+1) & h(n+1,m+1) by reading a file
    ! created from Web data::
    NameOutDir = NamePlotDir

    call read_harmonics

    timing0 = timing_cpu()
    call set_magnetogram
    if(iProc==0) &
         write(*,*) prefix, 'CPU timing:', timing_cpu() - timing0

  end subroutine read_magnetogram_file
  !============================================================================
  subroutine read_new_magnetogram_file(NamePlotDir, iProcIn, nProcIn, iCommIn)
    character(len=*),intent(in) :: NamePlotDir
    integer, intent(in):: iProcIn, nProcIn, iCommIn
    character(LEN=32):: NameOldFile
    real:: OldPhiShift
    !--------------------------------------------------------------------------
    NameOldFile = File_PFSSM
    File_PFSSM = NameNewFile
    OldPhiShift = PhiOffset
    PhiOffset = NewPhiShift
    call read_magnetogram_file(&
         NamePlotDir=NamePlotDir//'New',&
         iProcIn=iProcIn               ,&
         nProcIn=nProcIn               ,&
         iCommIn=iCommIn)

    ! Allocate the magnetic field array, at the spherical grid.
    if(allocated(BNew_DN))deallocate(BNew_DN)
    allocate(BNew_DN(R_:Theta_,-nRExt:nR,0:nPhi,0:nTheta))

    BNew_DN = B_DN

    PhiOffset = OldPhiShift
    File_PFSSM = NameOldFile

  end subroutine read_new_magnetogram_file
  !============================================================================

  subroutine read_harmonics
    integer :: iUnit,nOrderIn,iPosition,iPosition1,n,m,i,iError
    character (LEN=80):: Head_PFSSM=''
    real::gtemp,htemp,Coef1,Coef
    ! Formats adjusted for wso CR rad coeffs::
    !--------------------------------------------------------------------------
    iUnit = io_unit_new()
    open(iUnit,file=File_PFSSM,status='old',iostat=iError)
    if(iError>0)call CON_stop('Cannot open '//File_PFSSM)
    if (iHead_PFSSM /= 0) then
       do i=1,iHead_PFSSM
          read(iUnit,'(a)') Head_PFSSM
          iPosition=index(Head_PFSSM,'rder')
          iPosition1=0
          if(iPosition>0)&
               iPosition1=max(&
               index(Head_PFSSM(iPosition+4:iPosition+9),'='),&
               index(Head_PFSSM(iPosition+4:iPosition+9),':'))
          if(iPosition1>0)then
             read(Head_PFSSM(iPosition+4+iPosition1:len(Head_PFSSM)),&
                  '(i3)',iostat=iError)nOrderIn
             if(iError>0)call CON_stop('Cannot figure out nOrder')
             if(nOrderIn<N_PFSSM)then
                if(iProc==0)then
                   write(*,*) prefix, 'Reset N_PFSSM=',nOrderIn
                end if
                N_PFSSM=nOrderIn
             end if
          end if
          call get_hlcmm(Head_PFSSM,PhiOffset)
       enddo
       if(PhiOffset<0.0)call CON_stop(&
            'Did not find central meridian longitude')
    endif
    nR=max(nR,N_PFSSM)
    nPhi=max(nPhi,N_PFSSM)
    nTheta=max(nTheta,N_PFSSM)
    nRExt=min(10,1+floor(real(nR)*(0.1*Ro_PFSSM)/(Rs_PFSSM-Ro_PFSSM)))
    if(iProc==0)then
       write(*,*) prefix, &
            'Magnetogram is extended by ',nRExt,' nodes towards the Sun'
    end if
    ! Initialize all coefficient arrays::
    g_II = 0.0; h_II = 0.0
    ! Read file with coefficients, g_II and h_II::
    do
       read(iUnit,*,iostat=iError) n,m,gtemp,htemp
       if (iError /= 0) EXIT
       if (n > N_PFSSM .or. m > N_PFSSM) CYCLE
       g_II(n+1,m+1) = gtemp
       h_II(n+1,m+1) = htemp
    enddo
    close(iUnit)
    ! Add correction factor for radial, not LOS, coefficients::
    ! Note old "coefficients" file are LOS, all new coeffs and
    ! files are radial)
    Coef=Rs_PFSSM
    do n=0,N_PFSSM
       Coef=Coef/Rs_PFSSM**2
       Coef1 = 1.0/(real(n+1)+real(n)*Coef)
       !       Coef1 = 1.0/real(n+1+(n/(Rs_PFSSM**(2*n+1))))
       g_II(n+1,1:n+1) = g_II(n+1,1:n+1)*Coef1
       h_II(n+1,1:n+1) = h_II(n+1,1:n+1)*Coef1
    enddo
    ! Leave out monopole (n=0) term::
    g_II(1,1) = 0.0

    ! Set the Sqrt_I and SqrtRatio_I arrays (n_PFSSM is needed)
    call set_sqrt_arrays(n_PFSSM)

  end subroutine read_harmonics
  !============================================================================
  subroutine set_magnetogram

    ! Reads the file of magnetic field harmonics
    ! and recovers the spatial distribution of the
    ! potential mganetic field at the spherical grid
    ! nR*nPhi*nTheta

    ! This subroutine computes PFSS (Potential Field Source Surface)
    ! field model components in spherical coordinates at user-specified
    ! r (solar radii units, r>1) and theta, phi (both radians).
    ! The subroutine requires a file containing the spherical
    ! harmonic coefficients g(n,m) & h(n.m) obtained from a separate analysis
    ! of a photospheric synoptic map, in a standard ("radial")
    ! format and normalization used by Stanford.
    !
    ! The PFSS field model assumes no currents in the corona and
    ! a pure radial field at the source surface, here R=2.5 Rsun
    !
    ! Get solar coefficients from Todd Hoeksema's files:
    !    1. Go to http://solar.stanford.edu/~wso/forms/prgs.html
    !    2. Fill in name and email as required
    !    3. Chose Carrington rotation (leave default 180 center longitude)
    ! For most requests of integer CRs with order < 20, result will come back
    ! immediately on the web.
    !    4. Count header lines before 1st (0,0) coefficient -this will be asked!

    ! Notes:

    !
    ! The source surface surface radius has been set at Rs=2.5*Ro in
    ! the subroutine. PFSS fields at R>Rs are radial.(br,bthet,bphi)
    ! are the resulting components. Note the units of the B fields
    ! will differ with observatory used for the coefficients. The
    ! computation of the B fields is taken mainly from Altschuler,
    ! Levine, Stix, and Harvey, "High Resolutin Mapping of the
    ! Magnetic Field of the Solar Corona," Solar Physics 51 (1977)
    ! pp. 345-375. It uses Schmidt normalized Legendre polynomials and
    ! the normalization is explained in the paper. The field expansion
    ! in terms of the Schmidt normalized Pnm and dPnm's is best taken
    ! from Todd Hoeksema's notes which can be downloaded from the Web
    ! http://quake.stanford.edu/~wso/Description.ps The expansions
    ! used to get include radial factors to make the field become
    ! purely radial at the source surface. The expans. in Altschuler
    ! et al assumed that the the coefficient g(n,m) and h(n,m) were
    ! the LOS coefficients -- but the g(n,m) and h(n,m) now available
    ! are radial (according to Janet Luhman).  Therefore, she performs
    ! an initial correction to the g(n,m) and h(n,m) to make them
    ! consistent with the the expansion. There is no reference for
    ! this correction.

    integer:: iTheta, iPhi, iR
    real:: Theta,Phi,R_PFSSM

    integer::iBcast, iStart, nSize, iError

    ! True spherical components of the magnetic field and the potential
    real:: Bsph_D(3)
    !--------------------------------------------------------------------------
    ! Introduce a spherical grid with the resolution, depending on the
    ! magnetogram resolution (N_PFSSM)

    dR        = (Rs_PFSSM-Ro_PFSSM)/real(nR)
    dPhi      = cTwoPi/real(nPhi)
    dSinTheta = 2.0/real(nTheta+1)

    dInv_D=1.0/[dR,dPhi,dSinTheta]

    ! Allocate the magnetic field array, at the spherical grid.
    if(allocated(B_DN))deallocate(B_DN)
    allocate(B_DN(R_:Theta_,-nRExt:nR,0:nPhi,0:nTheta))

    B_DN=0.0

    ! Parallel computation of the magnetic field on the grid
    nThetaPerProc = nTheta/nProc + 1

    ! Loop by theta, each processor treats a separate part of the grid
    do iTheta=iProc*nThetaPerProc, (iProc+1)*nThetaPerProc-1

       if(iTheta>nTheta)EXIT ! Some processors may have less amount of work
       Theta = colatitude(iTheta)

       ! Start loop by Phi
       do iPhi=0,nPhi
          Phi=real(iPhi)*dPhi

          ! Loop by radius
          do iR = -nRExt, nR
             R_PFSSM = Ro_PFSSM + dR*iR

             call get_pfss_field(R_PFSSM, Theta, Phi, Bsph_D)

             ! Left handed coordinate system...
             B_DN(R_,iR,iPhi,iTheta)     = Bsph_D(1)
             B_DN(Phi_,iR,iPhi,iTheta)   = Bsph_D(3)
             B_DN(Theta_,iR,iPhi,iTheta) = Bsph_D(2)
          end do
       end do
    end do

    if(nProc>1)then
       ! this should be an mpi_allreduce with in-place
       do iBcast=0,nProc-1
          iStart=iBcast*nThetaPerProc
          if(iStart>nTheta)EXIT
          nSize=min(nThetaPerProc,nTheta+1-iStart)*(nPhi+1)*&
               (nR+1+nRExt)*3
          call MPI_bcast(B_DN(1,-nRExt,0,iStart),nSize,MPI_REAL,iBcast,iComm,iError)
       end do
    end if

    if(iProc==0) call write_Br_plot

  end subroutine set_magnetogram

  !==========================================================================
  subroutine set_sqrt_arrays(nOrder)

    integer, intent(in):: nOrder
    
    ! Calculate square roots and their ratios used in the spherical harmonics
    ! functions up to order nOrder for speeding up the calculations

    integer:: m
    !------------------------------------------------------------------------
    if(allocated(Sqrt_I)) RETURN
    
    allocate(Sqrt_I(MaxInt), SqrtRatio_I(nOrder+1))

    do m = 1, MaxInt
       Sqrt_I(m) = sqrt(real(m))
    end do
    
    ! Calculate the ratio sqrt(2m!)/(2^m*m!) recursively
    SqrtRatio_I(1) = 1.0
    do m = 1, nOrder
       SqrtRatio_I(m+1) = SqrtRatio_I(m)*Sqrt_I(2*m-1)/Sqrt_I(2*m)
    enddo

  end subroutine set_sqrt_arrays

  !============================================================================
  subroutine get_pfss_field(r, Theta, Phi, Bsph_D)
    
    ! Calculate the harmonics based potential field Bsph_D at r, Theta, Phi

    real, intent(in):: r, Theta, Phi ! r in Rs, Theta and Phi in radians
    real, intent(out):: Bsph_D(3)    ! Bsph_D with r, Theta, Phi components
    
    integer:: n, m
    real:: CosTheta, SinTheta

    ! Azimuthal functions sin(m*phi) and cos(m*phi)
    real, allocatable, save:: SinPhi_I(:), CosPhi_I(:)

    ! The spherical components of the magnetic field
    real:: SumR, SumT, SumP

    ! Legendre polynomials and its derivative
    real, save, allocatable:: p_II(:,:), Dp_II(:,:)

    ! Powers of radius ratios
    real, save, allocatable:: rRsPower_I(:), RoRsPower_I(:), RoRPower_I(:)

    real:: Coef1, Coef2, Coef3, Coef4
    !--------------------------------------------------------------------

    ! Calculate the radial part of spherical functions
    call calc_radial_functions(r, N_PFSSM)

    ! Calculate the set of Legendre polynoms for given Theta
    CosTheta = cos(Theta)
    SinTheta = max(sin(Theta), 1E-10)
    call calc_legendre_polynomial(SinTheta, CosTheta, N_PFSSM)

    call calc_azimuthal_functions(Phi, N_PFSSM)

    ! Initialize the values for the sums
    SumR = 0.0; SumT = 0.0; SumP = 0.0

    ! Calculate B from spherical harmonics
    do m = 0, N_PFSSM; do n = m, N_PFSSM

       Coef1 = (n+1)*RoRPower_I(n+2) + RoRsPower_I(n+2)*n*rRsPower_I(n-1)
       Coef3 =       RoRPower_I(n+2) - RoRsPower_I(n+2)*rRsPower_I(n-1)
       Coef2 = g_II(n+1,m+1)*CosPhi_I(m) + h_II(n+1,m+1)*SinPhi_I(m)
       Coef4 = g_II(n+1,m+1)*SinPhi_I(m) - h_II(n+1,m+1)*CosPhi_I(m)

       ! Br = -d(Psi)/dR
       SumR  = SumR + p_II(n+1,m+1)*Coef1*Coef2

       ! Bt = -(1/r)*d(Psi)/dTheta
       SumT  = SumT - Dp_II(n+1,m+1)*Coef2*Coef3

       ! Bp = -(1/r)*d(Psi)/dPhi
       SumP  = SumP + p_II(n+1,m+1)*m/SinTheta * Coef3*Coef4

       ! Potential
       ! SumPsi = SumPsi + r*p_II(n+1,m+1)*Coef2*Coef3

    enddo; enddo

    Bsph_D = [ SumR, SumT, SumP ]

  contains
    !==========================================================================

    subroutine calc_legendre_polynomial(SinTheta, CosTheta, nOrder)

      ! Calculate Legendre polynomials p_II and its derivative Dp_II
      ! with appropriate normalization for Theta

      real, intent(in):: SinTheta, CosTheta
      integer, intent(in):: nOrder

      real:: SinThetaM, SinThetaM1  ! sin(Theta)^m, sin(Theta)^(m-1)
      integer:: m

      integer:: nOrderOld = -10
      real:: SinThetaOld = -10.0, CosThetaOld = -10.0
      !------------------------------------------------------------------------
      ! Cache previous values
      if(nOrder == nOrderOld .and. SinTheta == SinThetaOld &
           .and. CosTheta == CosThetaOld) RETURN

      if(.not.allocated(p_II)) &
           allocate(p_II(nOrder+1,nOrder+1), Dp_II(nOrder+1,nOrder+1))

      nOrderOld = nOrder; SinThetaOld = SinTheta; CosThetaOld = CosTheta

      SinThetaM  = 1.0
      SinThetaM1 = 1.0
      p_II  = 0.0
      Dp_II = 0.0
      
      do m = 0, nOrder
         if (m == 0) then
            Coef1 = Sqrt_I(2*m+1)
         else
            Coef1 = Sqrt_I(2*(2*m+1))
         endif
         ! Eq.(27) from Altschuler et al. 1976::
         p_II(m+1,m+1) = SqrtRatio_I(m+1)*Coef1* SinThetaM
         ! Eq.(28) from Altschuler et al. 1976::
         if (m < nOrder) p_II(m+2,m+1) = p_II(m+1,m+1)*Sqrt_I(2*m+3)* &
              CosTheta
         ! Eq.(30) from Altschuler et al. 1976::
         Dp_II(m+1,m+1) = SqrtRatio_I(m+1)*Coef1*m*CosTheta*SinThetaM1
         ! Eq.(31) from Altschuler et al. 1976::
         if (m < nOrder) &
              Dp_II(m+2,m+1) = Sqrt_I(2*m+3)* &
              (CosTheta*Dp_II(m+1,m+1) - SinTheta*p_II(m+1,m+1))

         SinThetaM1 = SinThetaM
         SinThetaM  = SinThetaM*SinTheta

      enddo

      ! Recursive rules
      do m = 0, nOrder-2; do n = m+2, nOrder
         ! Eq.(29) from Altschuler et al. 1976::
         Coef1 = Sqrt_I(2*n+1)/Sqrt_I(n**2-m**2)
         Coef2 = Sqrt_I(2*n-1)
         Coef3 = Sqrt_I((n-1)**2-m**2)/Sqrt_I(2*n-3)

         p_II(n+1,m+1) = Coef1*(Coef2*CosTheta*p_II(n,m+1)  &
              - Coef3*p_II(n-1,m+1))
         
         ! Eq.(32) from Altschuler et al. 1976::
         Dp_II(n+1,m+1) = Coef1*(Coef2*(CosTheta*Dp_II(n,m+1) &
              - SinTheta*p_II(n,m+1)) - Coef3*Dp_II(n-1,m+1))
      enddo; enddo
      
      ! Apply Schmidt normalization
      do m = 0, nOrder; do n = m, nOrder
         ! Eq.(33) from Altschuler et al. 1976::
         Coef1 = 1.0/Sqrt_I(2*n+1)
         ! Eq.(34) from Altschuler et al. 1976::
         p_II(n+1,m+1)  = p_II(n+1,m+1)*Coef1
         Dp_II(n+1,m+1) = Dp_II(n+1,m+1)*Coef1
      enddo; enddo

    end subroutine calc_legendre_polynomial

    !==========================================================================
    subroutine calc_radial_functions(r, nOrder)

      real,    intent(in):: r
      integer, intent(in):: nOrder

      ! Calculate powers of the ratios of radii up to nOrder

      integer:: m
      real:: RoRs, RoR, rRs

      integer:: nOrderOld = -10
      real::    rOld = -10.0
      !------------------------------------------------------------------------

      if(nOrderOld == nOrder .and. rOld == r) RETURN

      if(nOrderOld /= nOrder .and. allocated(rRsPower_I)) &
           deallocate(rRsPower_I, RoRsPower_I, RoRPower_I)

      nOrderOld = nOrder; rOld = r
      
      if(.not.allocated(rRsPower_I)) allocate( &
           rRsPower_I(-1:nOrder+2), &
           RoRsPower_I(0:nOrder+2), &
           RoRPower_I(0:nOrder+2))

      RoRs = Ro_PFSSM/Rs_PFSSM
      RoR  = Ro_PFSSM/r
      rRs  = r/Rs_PFSSM

      ! Zero and -1 powers
      rRsPower_I(-1) = 1.0/rRs
      rRsPower_I(0)  = 1.0
      RoRsPower_I(0) = 1.0
      RoRPower_I(0)  = 1.0

      ! Recursive: x^m = x^m-1 * x
      do m = 1, nOrder+2
         RoRsPower_I(m) = RoRsPower_I(m-1) * RoRs
         RoRPower_I(m)  = RoRPower_I(m-1)  * RoR
         rRsPower_I(m)  = rRsPower_I(m-1)  * rRs
      end do

    end subroutine calc_radial_functions
    !==========================================================================
    subroutine calc_azimuthal_functions(Phi, nOrder)

      ! Calculate azimuthal harmonics for given Phi

      real,    intent(in):: Phi
      integer, intent(in):: nOrder
      
      integer:: m
      real   :: SinPhi, CosPhi
      
      integer:: nOrderOld = -10
      real   :: PhiOld    = -10.0
      !------------------------------------------------------------------------
      if(nOrderOld == nOrder .and. Phiold == Phi) RETURN

      if(nOrderOld /= nOrder .and. allocated(SinPhi_I)) &
           deallocate(SinPhi_I, CosPhi_I)

      if(.not.allocated(SinPhi_I)) &
           allocate(SinPhi_I(0:nOrder), CosPhi_I(0:nOrder))

      SinPhi_I(0) = 0.0
      CosPhi_I(0) = 1.0
      SinPhi = sin(Phi)
      CosPhi = Cos(Phi)

      ! Recursive: sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
      !            cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
      do m = 1, nOrder
         SinPhi_I(m) = SinPhi*CosPhi_I(m-1) + CosPhi*SinPhi_I(m-1)
         CosPhi_I(m) = CosPhi*CosPhi_I(m-1) - SinPhi*SinPhi_I(m-1)
      end do

    end subroutine calc_azimuthal_functions
  end subroutine get_pfss_field
  !============================================================================
  subroutine write_Br_plot

    use ModPlotFile, ONLY: save_plot_file

    integer :: iError,iPhi,iTheta,iUnit
    real,dimension(2,0:nPhi,0:nTheta):: Coord_DII,State_VII
    character(len=32)                :: FileNameDat
    character(len=32)                :: FileNameOut

    !--------------------------------------------------------------------------
    FileNameDat = trim(NameOutDir)//'PFSSM_Br.dat'
    FileNameOut = trim(NameOutDir)//'PFSSM_Br.out'

    write(*,*) prefix, 'Writing PFSSM_Br output file, named'
    write(*,*) prefix, FileNameDat
    iUnit = io_unit_new()
    open ( unit = iUnit, &
         file = FileNameDat, &
         form = 'formatted', &
         access = 'sequential', &
         status = 'replace', iostat = iError )

    if ( iError /= 0 ) then
       write (*,*) prefix, ' '
       write (*,*) prefix, 'TECPLOT_WRITE_OPEN - Fatal error!'
       write (*,*) prefix, '  Could not open the PFSSM_Br output file.'
       call CON_stop('')
    end if

    write ( iUnit, '(a)' ) 'Title = "'     // trim ('PFSSM_Br') // '"'
    write ( iUnit, '(a)' ) &
         'Variables = ' // trim (&
         '"Longitude [Deg]", "Latitude [Deg]", "Br_0 [G]","Br_SS [G]"')
    write ( iUnit, '(a)' ) ' '
    write ( iUnit, '(a,i6,a,i6,a)' ) 'Zone I = ', nPhi+1, ', J=', nTheta+1,&
         ', F=point'

    do iTheta=0,nTheta
       do iPhi=0,nPhi
          Coord_DII(1,iPhi,iTheta) = real(iPhi)*dPhi/cDegToRad
          Coord_DII(2,iPhi,iTheta) = r_latitude(iTheta)/cDegToRad

          State_VII(1,iPhi,iTheta) = UnitB*B_DN(R_,0,iPhi,iTheta)
          State_VII(2,iPhi,iTheta) = UnitB*B_DN(R_,nR,iPhi,iTheta)

          write ( iUnit, '(4f10.3)' )Coord_DII(:,iPhi,iTheta),&
               State_VII(:,iPhi,iTheta)
       end do
    end do
    close(iUnit)
    call save_plot_file(FileNameOut, &
         StringHeaderIn='Longitude, Latitude [Deg]; Br, BrSS [G]',&
         NameVarIn='Longitude Latitude Br BrSS', &
         CoordIn_DII=Coord_DII, &
         VarIn_VII=State_VII)

  end subroutine Write_Br_plot
  !============================================================================

  ! This subroutine corrects the angles Phi and Theta after every
  ! step of the field-alligned integration
  subroutine correct_angles(TR_D)
    real,dimension(nDim),intent(inout) :: TR_D

    !--------------------------------------------------------------------------
    if(TR_D(Theta_) < 0.0)then
       TR_D(Theta_) = -TR_D(Theta_)
       TR_D(Phi_)=TR_D(Phi_)+cPi
    end if
    if(TR_D(Theta_) > cPi)then
       TR_D(Theta_)=cTwoPi-TR_D(Theta_)
       TR_D(Phi_)=TR_D(Phi_)+cPi
    end if
    TR_D(Phi_)=modulo(TR_D(Phi_),cTwoPi)

  end subroutine correct_angles
  !============================================================================
  ! Calculate B_R, B_Phi, B_Theta in Gauss for given
  ! R, Phi, Theta, Theta being the colatitude.
  subroutine interpolate_field(R_D,BMap_D)
    real,intent(in)::R_D(nDim)
    real,intent(out)::BMap_D(nDim)
    integer::Node_D(nDim)
    real::Res_D(nDim)
    integer::iDim
    real::Weight_III(0:1,0:1,0:1)
    logical::DoCorrection
    real::ReductionCoeff,BRAvr
    !--------------------------------------------------------------------------
    Res_D=R_D
    ! Limit a value of R:
    Res_D(R_)=max(min(Res_D(R_),Rs_PFSSM-cTiny),Ro_PFSSM-nRExt*dR+cTiny)

    Res_D(R_)=Res_D(R_)-Ro_PFSSM

    call correct_angles(Res_D)
    DoCorrection=Res_D(Theta_)/=R_D(Theta_)
    Res_D(Theta_)=cos(Res_D(Theta_)) & ! This is sin(latitude)
         -sin_latitude(0)             ! This is sin(latitude) for iTheta=0
    Res_D=Res_D*dInv_D
    Node_D=floor(Res_D)
    if(Node_D(R_)==nR)Node_D(R_)=Node_D(R_)-1
    Res_D=Res_D-real(Node_D)
    ReductionCoeff=1.0
    BRAvr=0.0
    ! Near poles reduce the Phi and Theta components of the field and
    ! take a weithed average for the R component from using the  value in
    ! the last grid node and that averaged over longitude
    !(iTheta=0 or iTheta=nTheta)
    if(Node_D(Theta_)==nTheta)then
       Node_D(Theta_)=Node_D(Theta_)-1
       ReductionCoeff=sqrt(max(1.0-2.0*Res_D(Theta_),0.0))
       Res_D(Theta_)=1.0
       BRAvr=sum(&
            (1.0-Res_D(R_))*B_DN(R_,Node_D(R_)  ,:,nTheta)&
            +Res_D(R_) *B_DN(R_,Node_D(R_)+1,:,nTheta))/(nPhi+1)
    elseif(Node_D(Theta_)==-1)then
       Node_D(Theta_)= Node_D(Theta_)+1
       ReductionCoeff=sqrt(max(0.0,2.0*Res_D(Theta_)-1.0))
       Res_D(Theta_) = 0.0
       BRAvr=sum(&
            (1.0-Res_D(R_))*B_DN(R_,Node_D(R_)  ,:,     0)&
            +Res_D(R_) *B_DN(R_,Node_D(R_)+1,:,     0))/(nPhi+1)
    end if

    if(Node_D(Phi_)==nPhi)Node_D(Phi_)=0

    Weight_III(0,:,:)=1.0-Res_D(R_)
    Weight_III(1,:,:)=Res_D(R_)
    Weight_III(:,0,:)=Weight_III(:,0,:)*(1.0-Res_D(Phi_))
    Weight_III(:,1,:)=Weight_III(:,1,:)*Res_D(Phi_)
    Weight_III(:,:,0)=Weight_III(:,:,0)*(1.0-Res_D(Theta_))
    Weight_III(:,:,1)=Weight_III(:,:,1)*Res_D(Theta_)
    do iDim=1,nDim
       BMap_D(iDim)=&
            sum(Weight_III*B_DN(iDim,&
            Node_D(R_):Node_D(R_)+1,&
            Node_D(Phi_):Node_D(Phi_)+1,&
            Node_D(Theta_):Node_D(Theta_)+1))
    end do
    if(DoCorrection)then
       BMap_D(Phi_:Theta_)=-ReductionCoeff*BMap_D(Phi_:Theta_)
    else
       BMap_D(Phi_:Theta_)=ReductionCoeff*BMap_D(Phi_:Theta_)
    end if
    BMap_D(R_)=ReductionCoeff*BMap_D(R_)+(1.0-ReductionCoeff)*BRAvr
  end subroutine interpolate_field
  !============================================================================
  subroutine get_magnetogram_field(xInput,yInput,zInput,B0_D)

    ! For Cartesian xInput,yInput,zInput, all normalized per the
    ! solar radius, calculate the magnetic field vector
    ! in Cartesian components, in SI units.

    real, intent(in):: xInput,yInput,zInput
    real, intent(out), dimension(3):: B0_D

    real:: Rin_PFSSM,Theta_PFSSM,Phi_PFSSM
    real:: BMap_D(nDim)
    real:: SinPhi,CosPhi,XY,R_PFSSM
    real:: CosTheta,SinTheta

    ! Variables to obtain FDIPS field below Ro_PFSSM
    real :: BMapo_D(nDim), BMapi_D(nDim), Ri_PFSSM
    !--------------------------------------------------------------------------
    ! Calculate cell-centered spherical coordinates::
    Rin_PFSSM   = sqrt(xInput**2+yInput**2+zInput**2)
    ! Avoid calculating B0 inside a critical radius = 0.9*Rsun
    if(nRExt/=0 .and. Rin_PFSSM < max(Ro_PFSSM-dR*nRExt,0.9*Ro_PFSSM) &
         .or. nRExt==0 .and. Rin_PFSSM < 0.9*Ro_PFSSM)then
       B0_D= 0.0
       RETURN
    end if
    Theta_PFSSM = acos(zInput/Rin_PFSSM)
    Xy          = sqrt(xInput**2+yInput**2+cTiny**2)
    Phi_PFSSM   = atan2(yInput,xInput)
    SinTheta    = Xy/Rin_PFSSM
    CosTheta    = zInput/Rin_PFSSM
    SinPhi      = yInput/Xy
    CosPhi      = xInput/Xy
    ! Set the source surface radius::
    ! The inner boundary in the simulations starts at a height
    ! H_PFSSM above that of the magnetic field measurements!
    R_PFSSM = min(Rin_PFSSM + H_PFSSM, Rs_PFSSM)

    ! Transform Phi_PFSSM from the component's frame to the magnetogram's frame
    ! PhiOffset is the longitudinal shift from the original map +
    ! shift to the cell center.
    Phi_PFSSM = Phi_PFSSM - PhiOffset*cDegToRad

    if(nRExt == 0 .and. R_PFSSM < Ro_PFSSM)then
       ! linearly extrapolate FDIPS field for locations below Ro_PFSSM
       Ri_PFSSM = 2.0*Ro_PFSSM - R_PFSSM

       call interpolate_field([Ri_PFSSM,Phi_PFSSM,Theta_PFSSM],BMapi_D)
       call interpolate_field([Ro_PFSSM,Phi_PFSSM,Theta_PFSSM],BMapo_D)

       BMap_D = 2.0*BMapo_D - BMapi_D
    else
       call interpolate_field([R_PFSSM,Phi_PFSSM,Theta_PFSSM],BMap_D)
    end if
    ! Magnetic field components in global Cartesian coordinates::
    ! Set B0x::
    B0_D(1) = BMap_D(R_)*SinTheta*CosPhi+    &
         BMap_D(Theta_)*CosTheta*CosPhi-&
         BMap_D(Phi_)*SinPhi
    ! Set B0y::
    B0_D(2) =  BMap_D(R_)*SinTheta*SinPhi+    &
         BMap_D(Theta_)*CosTheta*SinPhi+&
         BMap_D(Phi_)*CosPhi
    ! Set B0z::
    B0_D(3) = BMap_D(R_)*CosTheta-           &
         BMap_D(Theta_)*SinTheta
    ! Apply field strength normalization::
    ! UnitB contains the units of the CR file relative to 1 Gauss
    ! and a possible correction factor (e.g. 1.8 or 1.7).
    B0_D=B0_D*UnitB

    if (Rin_PFSSM > Rs_PFSSM) &
         B0_D  =  B0_D*(Rs_PFSSM/Rin_PFSSM)**2

    ! Transform from Gauss to Tesla
    B0_D = B0_D*1.0E-4
  end subroutine get_magnetogram_field
  !============================================================================
  subroutine update_magnetogram(tNow, tMax, tLastUpdate)
    real, intent(in):: tNow, tMax
    real, intent(inout):: tLastUpdate

    !--------------------------------------------------------------------------
    if(tLastUpdate>=tMax)RETURN
    B_DN = ( (tMax - tNow     )/(tMax - tLastUpdate) )*B_DN + &
         ( (tNow - tLastUpdate)/(tMax - tLastUpdate) )*BNew_DN
    tLastUpdate = tNow
  end subroutine update_magnetogram
  !============================================================================
end module ModMagnetogram
!==============================================================================
