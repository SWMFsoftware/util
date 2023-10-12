!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program stand_alone_wsa

  ! Transform spherical harmonics file into 3D lookup table  for B0  field
  use ModMagnetogram, ONLY: read_magnetogram_param, &
       init_magnetogram_lookup_table
  use ModLookupTable, ONLY: read_lookup_table_param
  use ModMpi,         ONLY: MPI_init, MPI_Finalize, MPI_COMM_WORLD
  use ModReadParam,   ONLY: read_file, read_init, read_echo_set, &
       read_line, read_command, lStringLine
  use ModUtilities,   ONLY: CON_stop
  use ModWsa,         ONLY: set_open_threads, read_wsa_param
  implicit none
  character(len=lStringLine) :: NameCommand
  integer :: iError, iProc, nProc, iComm
  character(len=*), parameter:: NameSub = 'stand_alone_wsa'
  !----------------------------------------------------------------------------
  call MPI_init(iError)
  iComm = MPI_COMM_WORLD
  call read_file('WSA.in',iCommIn = iComm)
  call read_init
  call read_echo_set(.false.)

  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE
     select case(NameCommand)
     case("#FACTORB0", "#HARMONICSGRID", "#HARMONICSFILE")
        call read_magnetogram_param(NameCommand)
     case("#LOOKUPTABLE")
        call read_lookup_table_param
     case("#WSAGRID")
        call read_wsa_param
     case("#END")
        EXIT
     case default
        call CON_stop(NameSub//' invalid NameCommand='//NameCommand)
     end select
  end do
  call init_magnetogram_lookup_table(iComm)
  call MPI_comm_rank(iComm, iProc, iError)
  call MPI_comm_size(iComm, nProc, iError)
  call set_open_threads(iProc, nProc, iComm)
  call MPI_finalize(iError)
end program stand_alone_wsa
!==============================================================================
