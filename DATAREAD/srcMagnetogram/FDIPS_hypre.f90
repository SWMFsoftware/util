module ModHypre

  use ModPotentialField

  implicit none

  private
  public:: hypre_initialize
  public:: hypre_solver
  public:: hypre_preconditioner

  include "HYPREf.h"

  ! This is defined in sstruct_mv/HYPRE_sstruct_mv.h
  integer, parameter:: HYPRE_SSTRUCT_VARIABLE_CELL = 0

  integer, parameter:: nDim = 3
  integer, parameter:: nStencil = 2*nDim+1

  integer, parameter :: Int8_ = selected_int_kind(16)

  ! Hypre uses C type 8-byte integers for pointers
  integer(Int8_):: i8Grid, i8Graph, i8Stencil, i8Precond, i8Solver
  integer(Int8_):: i8A, i8B, i8X, i8ParA, i8ParB, i8ParX

  integer:: nPart = 1 ! grid consists of a single part
  integer:: iPart = 0 ! index of the part
  integer:: nVar  = 1 ! there is only one variable per grid cell
  integer:: iVar  = 0 ! index of the variable
  integer:: iVarType_I(1) = (/ HYPRE_SSTRUCT_VARIABLE_CELL /)
  integer:: iLower_D(nDim), iUpper_D(nDim) ! index limits for each box

  integer:: iLowerBc_D(nDim), iUpperBc_D(nDim) ! index limit for ghost cells
  integer:: jLowerBc_D(nDim), jUpperBc_D(nDim) ! index limit for phys. cells

  ! Mapping between boundaries: same index order and same direction
  integer, parameter:: iIndexMap_D(nDim) = (/0,1,2/) 
  integer, parameter:: iIndexDir_D(nDim) = (/1,1,1/)

  integer:: iStencil
  integer:: iStencil_I(nStencil) = (/ (iStencil, iStencil=0,nStencil-1) /)
  integer:: DiStencil_DI(nDim,nStencil)   ! Stencil description

  integer:: iObjectType

  integer          :: nValue     ! number of matrix elements for 1 box
  real, allocatable:: Value_I(:) ! matrix elements

contains

  !==========================================================================
  subroutine hypre_initialize

    integer:: iValue, iR, iTheta, iPhi, iError

    character(len=*), parameter:: NameSub = 'hypre_initialize'
    !-------------------------------------------------------------------------

    if(DoTestMe)write(*,*) NameSub,' starting'

    if(NamePreconditioner == 'MG')then
       iObjectType = HYPRE_STRUCT
    else
       iObjectType = HYPRE_PARCSR
    end if

    ! Create an empty 3D grid object
    call HYPRE_SStructGridCreate(iComm, nDim, nPart, i8Grid, iError)

    ! Add the box = local domain
    iLower_D = (/  1, iTheta0+1,      iPhi0+1 /)
    iUpper_D = (/ nR, iTheta0+nTheta, iPhi0+nPhi /)
    call HYPRE_SStructGridSetExtents(i8Grid, iPart, iLower_D, iUpper_D, iError)

    ! Single cell centered variable
    call HYPRE_SStructGridSetVariables(i8Grid, iPart, nVar, iVarType_I, iError)

    ! Setup periodic boundaries in Phi direction
    iLowerBc_D = (/  1,         1,    0 /)
    iUpperBc_D = (/ nR, nThetaAll,    0 /)
    jLowerBc_D = (/  1,         1, nPhiAll /)
    jUpperBc_D = (/ nR, nThetaAll, nPhiAll /)

    call HYPRE_SStructGridSetNeighborPart( i8Grid, &
         iPart, iLowerBc_D, iUpperBc_D, &
         iPart, jLowerBc_D, jUpperBc_D, &
         iIndexMap_D, iIndexDir_D, iError)

    iLowerBc_D = (/  1,         1, nPhiAll+1 /)
    iUpperBc_D = (/ nR, nThetaAll, nPhiAll+1 /)
    jLowerBc_D = (/  1,         1,         1 /)
    jUpperBc_D = (/ nR, nThetaAll,         1 /)

    call HYPRE_SStructGridSetNeighborPart( i8Grid, &
         iPart, iLowerBc_D, iUpperBc_D, &
         iPart, jLowerBc_D, jUpperBc_D, &
         iIndexMap_D, iIndexDir_D, iError)

    ! Assemble grid from all processors
    call HYPRE_SStructGridAssemble(i8Grid, iError)
    if(iError/=0)write(*,*)'ERROR: HYPRE_SStructGridAssemble failed'
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructGridAssemble done'

    ! Define index offsets for the 7-point stencil
    DiStencil_DI(:,1) = (/  0, 0, 0 /)
    DiStencil_DI(:,2) = (/ -1, 0, 0 /)
    DiStencil_DI(:,3) = (/ +1, 0, 0 /)
    DiStencil_DI(:,4) = (/  0,-1, 0 /)
    DiStencil_DI(:,5) = (/  0,+1, 0 /)
    DiStencil_DI(:,6) = (/  0, 0,-1 /)
    DiStencil_DI(:,7) = (/  0, 0,+1 /)

    call HYPRE_SStructStencilCreate(nDim, nStencil, i8Stencil, iError)

    do iStencil = 1, nStencil
       call HYPRE_SStructStencilSetEntry(i8Stencil, &
            iStencil-1, DiStencil_DI(1,iStencil), iVar, iError)
    enddo

    ! Create the graph object
    call HYPRE_SStructGraphCreate(iComm, i8Grid, i8Graph, iError)

    ! Tell the graph which stencil to use for each variable on each part 
    call HYPRE_SStructGraphSetStencil(i8Graph, iPart, iVar, i8Stencil, iError)

    ! Since there is only a single part, there is nothing to do, only assemble
    call HYPRE_SStructGraphAssemble(i8Graph, iError)

    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructGraphAssemble done'


    ! Create an empty matrix object
    call HYPRE_SStructMatrixCreate(iComm, i8Graph, i8A, iError)

    ! Set storage type
    call HYPRE_SStructMatrixSetObjectTyp(i8A, iObjectType, iError)

    ! Get ready to set values
    call HYPRE_SStructMatrixInitialize(i8A, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructMatrixInitialize done'

    ! Set non-zero matrix elements
    nValue = nStencil*nR*nTheta*nPhi
    allocate(Value_I(nValue))

    Value_I(1:nValue:nStencil) = -d_I
    Value_I(2:nValue:nStencil) = -e_I
    Value_I(3:nValue:nStencil) = -f_I
    Value_I(4:nValue:nStencil) = -e1_I
    Value_I(5:nValue:nStencil) = -f1_I
    Value_I(6:nValue:nStencil) = -e2_I
    Value_I(7:nValue:nStencil) = -f2_I

    call HYPRE_SStructMatrixSetBoxValues(i8A, iPart, iLower_D, iUpper_D, &
         iVar, nStencil, iStencil_I, Value_I, iError)

    ! Assemble matrix
    call HYPRE_SStructMatrixAssemble(i8A, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructMatrixAssemble done'
    deallocate(Value_I)

    ! Create empty vector objects for RHS and solution
    call HYPRE_SStructVectorCreate(iComm, i8Grid, i8B, iError)
    call HYPRE_SStructVectorCreate(iComm, i8Grid, i8X, iError)

    ! Set object type for the vectors
    call HYPRE_SStructVectorSetObjectTyp(i8B, iObjectType, iError)
    call HYPRE_SStructVectorSetObjectTyp(i8X, iObjectType, iError)

    ! Initialize vectors
    call HYPRE_SStructVectorInitialize(i8B, iError)
    call HYPRE_SStructVectorInitialize(i8X, iError)

    ! Pass matrix to the solvers
    call HYPRE_SStructMatrixGetObject(i8A, i8ParA, iError)

    select case(NamePreconditioner)
    case('MG')
       ! Create the PFMG preconditioner
       call HYPRE_StructPFMGCreate(iComm, i8Precond, iError)

       ! Set PFMG parameters
       call HYPRE_StructPFMGSetMaxIter(i8Precond, 1, iError)
       call HYPRE_StructPFMGSetTol(i8Precond, 0.0, iError)
       call HYPRE_StructPFMGSetZeroGuess(i8Precond, iError)
       call HYPRE_StructPFMGSetNumPreRelax(i8Precond, 2, iError)
       call HYPRE_StructPFMGSetNumPostRelax(i8Precond, 2, iError)
       ! Non-Galerkin coarse grid (more efficient for this problem)
       call HYPRE_StructPFMGSetRAPType(i8Precond, 1, iError)
       ! R/B Gauss-Seidel
       call HYPRE_StructPFMGSetRelaxType(i8Precond, 2, iError)
       ! Skip relaxation on some levels (more efficient for this problem)
       call HYPRE_StructPFMGSetSkipRelax(i8Precond, 1, iError)

       if(DoTestMe)write(*,*) NameSub,' HYPRE_StructPFMGSetSkipRelax done'

    case('AMG')
       ! Create the BoomerAMG as a preconditioner
       call HYPRE_BoomerAMGCreate(i8Precond, iError)

       ! Set BoomerAMG parameters
       call HYPRE_BoomerAMGSetMaxIter(i8Precond, 1, iError)
       call HYPRE_BoomerAMGSetTol(i8Precond, 0.0, iError)

       ! Print AMG solution info
       call HYPRE_BoomerAMGSetPrintLevel(i8Precond, 2, iError)
       call HYPRE_BoomerAMGSetCoarsenType(i8Precond, 6, iError)

       ! Sym G.S./Jacobi hybrid
       call HYPRE_BoomerAMGSetRelaxType(i8Precond, 6, iError)
       call HYPRE_BoomerAMGSetNumSweeps(i8Precond, 1, iError)

       if(UsePreconditioner)then
          ! Setup AMG preconditioner for Krylov solver
          if(UseTiming)write(*,*)NameSub, &
               ' time before BoomerAMGSetup:', MPI_WTIME() - TimeStart
          call HYPRE_BoomerAMGSetup(i8Precond, i8ParA, i8ParB, i8ParX, iError)
          if(UseTiming)write(*,*)NameSub, &
               ' time after BoomerAMGSetup:', MPI_WTIME() - TimeStart
       end if

    end select

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine hypre_initialize

  !==========================================================================

  subroutine hypre_solver

    integer:: iValue, iR, iTheta, iPhi, iError

    character(len=*), parameter:: NameSub = 'hypre_solver'
    !-------------------------------------------------------------------------

    ! Set RHS values
    nValue = nR*nTheta*nPhi
    allocate(Value_I(nValue))

    if(DoTestMe)write(*,*) NameSub,' starting n, maxval, minval(Rhs_C)=', &
         nValue, maxval(Rhs_C), minval(Rhs_C), Rhs_C(1,1,1)

    iValue = 0
    do iPhi = 1, nPhi; do iTheta = 1, nTheta; do iR = 1, nR
       iValue = iValue + 1
       Value_I(iValue) = -Rhs_C(iR,iTheta,iPhi)
    end do; end do; end do

    call HYPRE_SStructVectorSetBoxValues(i8B, iPart, iLower_D, iUpper_D, &
         iVar, Value_I, iError)

    if(DoTestMe)write(*,*) NameSub,' set RHS'

    ! Set initial guess value to zero
    Value_I = 0.0
    call HYPRE_SStructVectorSetBoxValues(i8X, iPart, iLower_D, iUpper_D, &
         iVar, Value_I, iError)

    if(DoTestMe)write(*,*) NameSub,' set X=0'

    ! Assemble vectors
    call HYPRE_SStructVectorAssemble(i8X, iError)
    call HYPRE_SStructVectorAssemble(i8B, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructVectorAssemble done'

    ! Pass vectors to the solvers
    call HYPRE_SStructVectorGetObject(i8B, i8ParB, iError)
    call HYPRE_SStructVectorGetObject(i8X, i8ParX, iError)
    if(DoTestMe)write(*,*) NameSub,' passed vectors to solver'

    select case(NameSolver)
    case('AMG')
       ! Create an empty BoomerAMG solver
       call HYPRE_BoomerAMGCreate(i8Solver, iError)
       if(DoTestMe)write(*,*) NameSub,' HYPRE_BoomerAMGCreate done'
       ! print solve info + parameters
       call HYPRE_BoomerAMGSetPrintLevel(i8Solver, 3, iError)
       ! Falgout coarsening
       call HYPRE_BoomerAMGSetCoarsenType(i8Solver, 6, iError)
       ! G-S/Jacobi hybrid relaxation
       call HYPRE_BoomerAMGSetRelaxType(i8Solver, 3, iError)
       ! Sweeeps on each level
       call HYPRE_BoomerAMGSetNumSweeps(i8Solver, 1, iError)
       ! maximum number of levels
       call HYPRE_BoomerAMGSetMaxLevels(i8Solver, 20, iError)
       ! conv. tolerance
       call HYPRE_BoomerAMGSetTol(i8solver, Tolerance, iError)
       if(DoTestMe)write(*,*) NameSub,' HYPRE_BoomerAMGSetTol done'

       ! Now setup and solve
       if(UseTiming)write(*,*)NameSub, &
            ' time before BoomerAMGSetup:', MPI_WTIME() - TimeStart

       call HYPRE_BoomerAMGSetup(i8Solver, i8ParA, i8ParB, i8ParX, iError)

       if(UseTiming)write(*,*)NameSub, &
            ' time before BoomerAMGSetup:', MPI_WTIME() - TimeStart

       call HYPRE_BoomerAMGSolve(i8Solver, i8ParA, i8ParB, i8ParX, iError)

       if(UseTiming)write(*,*)NameSub, &
            ' time after BoomerAMGSolve:', MPI_WTIME() - TimeStart

       ! Free memory
       call HYPRE_BoomerAMGDestroy(i8Solver, iError)

    case('GMRES')
       select case(NamePreconditioner)
       case('MG')
          ! Create an empty GMRES solver
          call HYPRE_StructGMRESCreate(iComm, i8Solver, iError)

          ! Set GMRES parameters
          call HYPRE_StructGMRESSetTol(i8Solver, Tolerance, iError)
          call HYPRE_StructGMRESSetPrintLevel(i8Solver, 2, iError) !!! 2
          call HYPRE_StructGMRESSetMaxIter(i8Solver, 50, iError)

          if(DoTestMe)write(*,*) NameSub,' HYPRE_StructGMRESSetMaxIter done'

          ! Set preconditioner (PFMG = 1) and solve
          call HYPRE_StructGMRESSetPrecond(i8Solver, 1, i8Precond, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_StructGMRESSetPrecond done'

          if(UseTiming)write(*,*)NameSub, &
               ' time before HYPRE_StructGMRESSetup:', MPI_WTIME() - TimeStart

          call HYPRE_StructGMRESSetup(i8Solver, i8ParA, i8ParB, i8ParX, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_StructGMRESSetup done'

          if(UseTiming)write(*,*)NameSub, &
               ' time after HYPRE_StructGMRESSetup:', MPI_WTIME() - TimeStart

          call HYPRE_StructGMRESSolve(i8Solver, i8ParA, i8ParB, i8ParX, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_StructGMRESSolve done'

          if(UseTiming)write(*,*)NameSub, &
               ' time after HYPRE_StructGMRESSolve:', MPI_WTIME() - TimeStart

          ! Free memory
          call HYPRE_StructGMRESDestroy(i8Solver, iError)
          call HYPRE_StructPFMGDestroy(i8Precond, iError)
       case('AMG')
          ! Create an empty GMRES solver
          call HYPRE_ParCSRGMRESCreate(iComm, i8Solver, iError)

          ! Set GMRES parameters
          call HYPRE_ParCSRGMRESSetTol(i8Solver, Tolerance, iError)
          call HYPRE_ParCSRGMRESSetPrintLevel(i8Solver, 100, iError) !!! 2
          call HYPRE_ParCSRGMRESSetMaxIter(i8Solver, 50, iError)

          if(DoTestMe)write(*,*) NameSub,' HYPRE_ParCSRGMRESSetMaxIter done'

          ! Set preconditioner (BoomerAMG = 2) and solve
          call HYPRE_ParCSRGMRESSetPrecond(i8Solver, 2, i8Precond, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_ParCSRGMRESSetPrecond done'

          if(UseTiming)write(*,*) NameSub, &
               ' time before HYPRE_ParCSRGMRESSetup:', MPI_WTIME() - TimeStart

          call HYPRE_ParCSRGMRESSetup(i8Solver, i8ParA, i8ParB, i8ParX, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_ParCSRGMRESSetup done'

          if(UseTiming)write(*,*)NameSub, &
               ' time after HYPRE_ParCSRGMRESSetup:', MPI_WTIME() - TimeStart

          call HYPRE_ParCSRGMRESSolve(i8Solver, i8ParA, i8ParB, i8ParX, iError)
          if(DoTestMe)write(*,*) NameSub,' HYPRE_ParCSRGMRESSolve done'

          if(UseTiming)write(*,*)NameSub, &
               ' time after HYPRE_ParCSRGMRESSolve:', MPI_WTIME() - TimeStart

          ! Free memory
          call HYPRE_ParCSRGMRESDestroy(i8Solver, iError)
          call HYPRE_BoomerAMGDestroy(i8Precond, iError)
       end select
    end select

    ! Get result
    Potential_C = 0.0
    call HYPRE_SStructVectorGather(i8x, iError);
    call HYPRE_SStructVectorGetBoxValues(i8X, iPart, &
         iLower_D, iUpper_D, iVar, Potential_C, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructVectorGetBoxValues done'

    ! Free memory
    call HYPRE_SStructGridDestroy(i8Grid, iError)
    call HYPRE_SStructStencilDestroy(i8Stencil, iError)
    call HYPRE_SStructGraphDestroy(i8Graph, iError)
    call HYPRE_SStructMatrixDestroy(i8A, iError)
    call HYPRE_SStructVectorDestroy(i8B, iError)
    call HYPRE_SStructVectorDestroy(i8X, iError)

    deallocate(Value_I)

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine hypre_solver

  !===========================================================================

  subroutine hypre_preconditioner(n, y_I)

    integer, intent(in):: n
    real, intent(inout):: y_I(n)

    integer:: iError
    real, allocatable:: Value_I(:)

    logical, parameter:: DoDebug = .false.

    character(len=*), parameter:: NameSub = 'hypre_preconditioner'
    !-------------------------------------------------------------------------

    if(DoDebug)write(*,*) NameSub,' starting n, maxval, minval, y_I(1)=', &
         n, maxval(y_I), minval(y_I), y_I(1)

    ! Preconditioning: y'= AMG.y

    ! Set y_I as the RHS
    call HYPRE_SStructVectorSetBoxValues(i8B, iPart, iLower_D, iUpper_D, &
         iVar, y_I, iError)

    if(DoDebug)write(*,*) NameSub,' set RHS, iLower_D, iUpper_D=', &
         iLower_D, iUpper_D

    ! Set initial guess value to zero
    allocate(Value_I(n))
    Value_I = 0.0
    call HYPRE_SStructVectorSetBoxValues(i8X, iPart, iLower_D, iUpper_D, &
         iVar, Value_I, iError)
    deallocate(Value_I)

    if(DoDebug)write(*,*) NameSub,' set X=0'

    ! Assemble vectors
    call HYPRE_SStructVectorAssemble(i8X, iError)
    call HYPRE_SStructVectorAssemble(i8B, iError)

    if(DoDebug)write(*,*) NameSub,' HYPRE_SStructVectorAssemble done'

    ! Pass the vectors to the solvers
    call HYPRE_SStructVectorGetObject(i8B, i8ParB, iError)
    call HYPRE_SStructVectorGetObject(i8X, i8ParX, iError)

    if(DoDebug)write(*,*) NameSub,' passed vectors to AMG'

    call HYPRE_BoomerAMGSolve(i8Precond, i8ParA, i8ParB, i8ParX, iError)

    if(DoDebug)write(*,*) NameSub,' applied AMG preconditioner'

    ! Get back solution
    call HYPRE_SStructVectorGather(i8x, iError);
    call HYPRE_SStructVectorGetBoxValues(i8X, iPart, &
         iLower_D, iUpper_D, iVar, y_I, iError)

    if(DoDebug)write(*,*) NameSub,' finished'

  end subroutine hypre_preconditioner

end module ModHypre
