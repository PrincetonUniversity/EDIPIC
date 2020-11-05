
!-------------------------------------------
! this subroutine is called by every process
!
SUBROUTINE SAVE_CHECKPOINT_MPIIO

  USE CurrentProblemValues
  USE ParallelOperationValues
  USE Diagnostics
  USE Snapshots

  USE rng_wrapper

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER*4 state_ind, func, state_var(624)

  INTEGER ibufer_length
  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER ALLOC_ERR

  INTEGER rbufer_length
  REAL(8), ALLOCATABLE :: rbufer(:)

  INTEGER file_handle

  CHARACTER(20) filename_check         ! Tcntr_TTTTTTTT.check
                                       ! ----x----I----x----I

  INTEGER, ALLOCATABLE :: jbufer(:)
  REAL(8), ALLOCATABLE :: dbufer(:)
  INTEGER s, k, bufsize, pos

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (T_cntr.NE.Save_check_Tcntr) RETURN

  CALL get_rng_state(state_var, state_ind, func)

  ibufer_length = 637

  ALLOCATE(ibufer(1:ibufer_length), STAT = ALLOC_ERR)

  ibufer(1:624) = state_var(1:624)
  ibufer(625) = state_ind
  ibufer(626) = func
  ibufer(627) = T_cntr
  ibufer(628) = Start_diag_Tcntr
  ibufer(629) = current_snap
  ibufer(630) = N_part(1)
  ibufer(631) = Length(1)
  ibufer(632) = N_part(2)
  ibufer(633) = Length(2)
  ibufer(634) = Q_left
  ibufer(635) = Q_right
  ibufer(636) = N_of_saved_records
  ibufer(637) = text_output_counter

!  ibufer(4) = Save_probes_data_T_cntr
  
  rbufer_length = 16

  ALLOCATE(rbufer(1:rbufer_length), STAT = ALLOC_ERR)

  rbufer(1) = full_Q_left
  rbufer(2) = full_Q_right
  rbufer(3) = Energy_full_eV
  rbufer(4) = Init_energy_full_eV
  rbufer(5) = Energy_pot_eV
  rbufer(6) = Energy_heat_eV
  rbufer(7) = Energy_kin_eV(1)
  rbufer(8) = Energy_kin_eV(2)
  rbufer(9) = Energy_wall_eV(1)
  rbufer(10) = Energy_wall_eV(2)
  rbufer(11) = Energy_emit_eV(1)
  rbufer(12) = Energy_emit_eV(2)
  rbufer(13) = Energy_coll_eV(1)
  rbufer(14) = Energy_coll_eV(2)
  rbufer(15) = prev_Q_left
  rbufer(16) = prev_Q_right

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! save everything into a single binary file 

! create filename
  filename_check = 'Tcntr_TTTTTTTT.check'
  filename_check(7:14) = convert_int_to_txt_string(T_cntr, 8)

  CALL MPI_FILE_OPEN( MPI_COMM_WORLD, &
                    & filename_check,  &
                    & MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                    & MPI_INFO_NULL, &
                    & file_handle, &
                    & ierr )
  
  CALL MPI_FILE_WRITE_ORDERED( file_handle, ibufer, ibufer_length, MPI_INTEGER, stattus, ierr )
  
  CALL MPI_FILE_WRITE_ORDERED( file_handle, rbufer, rbufer_length, MPI_DOUBLE_PRECISION, stattus, ierr )

  IF (Rank_of_process.EQ.0) N_part = 0  ! restore it here since it may be modified in collisions procedure 

  DO s = 1, N_spec

     bufsize = 7 * N_part(s)
     IF (ALLOCATED(dbufer)) DEALLOCATE(dbufer, STAT = ALLOC_ERR)
     ALLOCATE(dbufer(1:bufsize), STAT = ALLOC_ERR)
     pos = 1
     DO k = 1, N_part(s)
        dbufer(pos)   = species(s)%part(k)%X
        dbufer(pos+1) = species(s)%part(k)%Z
        dbufer(pos+2) = species(s)%part(k)%VX
        dbufer(pos+3) = species(s)%part(k)%VY
        dbufer(pos+4) = species(s)%part(k)%VZ
        dbufer(pos+5) = species(s)%part(k)%AX
        dbufer(pos+6) = species(s)%part(k)%prev_VX
        pos = pos + 7
     END DO

     CALL MPI_FILE_WRITE_ORDERED( file_handle, dbufer, bufsize, MPI_DOUBLE_PRECISION, stattus, ierr )

     bufsize = N_part(s)
     IF (ALLOCATED(jbufer)) DEALLOCATE(jbufer, STAT = ALLOC_ERR)
     ALLOCATE(jbufer(1:bufsize), STAT = ALLOC_ERR)
     DO k = 1, N_part(s)
        jbufer(k) = species(s)%part(k)%tag
     END DO

     CALL MPI_FILE_WRITE_ORDERED( file_handle, jbufer, bufsize, MPI_INTEGER, stattus, ierr )

  END DO

  CALL MPI_FILE_CLOSE(file_handle, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (ALLOCATED(ibufer)) DEALLOCATE(ibufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(dbufer)) DEALLOCATE(dbufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(jbufer)) DEALLOCATE(jbufer, STAT = ALLOC_ERR)

  Save_check_Tcntr = Save_check_Tcntr + SaveCheck_step

END SUBROUTINE SAVE_CHECKPOINT_MPIIO

!-------------------------------------------
! this subroutine is called by every process
!
SUBROUTINE READ_CHECKPOINT_MPIIO

  USE CurrentProblemValues
  USE ParallelOperationValues
  USE Diagnostics
  USE Snapshots

  USE rng_wrapper

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  LOGICAL exists

  INTEGER*4 state_ind, func, state_var(624)

  INTEGER ibufer_length
  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER ALLOC_ERR

  INTEGER rbufer_length
  REAL(8), ALLOCATABLE :: rbufer(:)

  INTEGER file_handle

  CHARACTER(20) filename_check         ! Tcntr_TTTTTTTT.check
                                       ! ----x----I----x----I

  INTEGER, ALLOCATABLE :: jbufer(:)
  REAL(8), ALLOCATABLE :: dbufer(:)
  INTEGER s, k, bufsize, pos

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

! foolproof: make sure that the unpacked file has proper name

  filename_check = 'Tcntr_TTTTTTTT.check'
  filename_check(7:14) = convert_int_to_txt_string(T_cntr_to_continue, 8)

  INQUIRE (FILE = filename_check, EXIST = exists)
  IF (.NOT.exists) THEN
     PRINT '("Process :: ",i5,", ERROR in READ_CHECKPOINT_MPIIO :: file ",A20," not found program terminated")', &
          & Rank_of_process, filename_check
     STOP
  END IF

! read the binary file 

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  call MPI_FILE_OPEN( MPI_COMM_WORLD, &
                    & filename_check,  &
                    & MPI_MODE_RDONLY, & 
                    & MPI_INFO_NULL, &
                    & file_handle, &
                    & ierr )

  ibufer_length = 637
  rbufer_length = 16

  ALLOCATE(ibufer(1:ibufer_length), STAT = ALLOC_ERR)
  ALLOCATE(rbufer(1:rbufer_length), STAT = ALLOC_ERR)

  CALL MPI_FILE_READ_ORDERED( file_handle, ibufer, ibufer_length, MPI_INTEGER, stattus, ierr)

  CALL MPI_FILE_READ_ORDERED( file_handle, rbufer, rbufer_length, MPI_DOUBLE_PRECISION, stattus, ierr )

  state_var(1:624) = ibufer(1:624)
  state_ind = ibufer(625)
  func      = ibufer(626)
  Start_T_cntr     = ibufer(627)
  Start_diag_Tcntr = ibufer(628)
  current_snap     = ibufer(629)
  N_part(1)    = ibufer(630)
  Length(1)    = ibufer(631)
  N_part(2)    = ibufer(632)
  Length(2)    = ibufer(633)
  Q_left       = ibufer(634)
  Q_right      = ibufer(635)
  N_of_saved_records  = ibufer(636)
  text_output_counter = ibufer(637)

  CALL set_rng_state(state_var, state_ind, func)

  full_Q_left    = rbufer(1)
  full_Q_right   = rbufer(2)
  Energy_full_eV = rbufer(3)
  Init_energy_full_eV = rbufer(4)
  Energy_pot_eV    = rbufer(5)
  Energy_heat_eV   = rbufer(6)
  Energy_kin_eV(1) = rbufer(7)
  Energy_kin_eV(2) = rbufer(8)
  Energy_wall_eV(1) = rbufer(9)
  Energy_wall_eV(2) = rbufer(10)
  Energy_emit_eV(1) = rbufer(11)
  Energy_emit_eV(2) = rbufer(12)
  Energy_coll_eV(1) = rbufer(13)
  Energy_coll_eV(2) = rbufer(14)
  prev_Q_left  = rbufer(15)
  prev_Q_right = rbufer(16)

  CALL CONFIGURE_PARTICLE_DYNAM_ARRAYS

  DO s = 1, N_spec

     bufsize = 7 * N_part(s)
     IF (ALLOCATED(dbufer)) DEALLOCATE(dbufer, STAT = ALLOC_ERR)
     ALLOCATE(dbufer(1:bufsize), STAT = ALLOC_ERR)

print '("proc ",i4," read-3 :: spec ",i2," N_part(s) ",i7," bufsize ",i7)', Rank_of_process, s, N_part(s), bufsize

     CALL MPI_FILE_READ_ORDERED( file_handle, dbufer, bufsize, MPI_DOUBLE_PRECISION, stattus, ierr )

     pos = 1
     DO k = 1, N_part(s)
        species(s)%part(k)%X = dbufer(pos)
        species(s)%part(k)%Z = dbufer(pos+1) 
        species(s)%part(k)%VX = dbufer(pos+2)
        species(s)%part(k)%VY = dbufer(pos+3)
        species(s)%part(k)%VZ = dbufer(pos+4) 
        species(s)%part(k)%AX = dbufer(pos+5)
        species(s)%part(k)%prev_VX = dbufer(pos+6)
        pos = pos + 7
     END DO

     bufsize = N_part(s)
     IF (ALLOCATED(jbufer)) DEALLOCATE(jbufer, STAT = ALLOC_ERR)
     ALLOCATE(jbufer(1:bufsize), STAT = ALLOC_ERR)

     CALL MPI_FILE_READ_ORDERED( file_handle, jbufer, bufsize, MPI_INTEGER, stattus, ierr )

     DO k = 1, N_part(s)
        species(s)%part(k)%tag = jbufer(k)
     END DO

  END DO

  CALL MPI_FILE_CLOSE(file_handle, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (ALLOCATED(ibufer)) DEALLOCATE(ibufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(dbufer)) DEALLOCATE(dbufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(jbufer)) DEALLOCATE(jbufer, STAT = ALLOC_ERR)

END SUBROUTINE READ_CHECKPOINT_MPIIO

!-----------------------------------
! creates a string of length "length_of_string" out of an integer number "int_number"
!
function convert_int_to_txt_string(int_number, length_of_string)

  implicit none

  character*(length_of_string) convert_int_to_txt_string
  integer int_number
  integer length_of_string
  character(5) format_string
  character(2) length2_txt
  character(1) length1_txt

  character*(length_of_string) number_txt
  
  integer blanks_number
  integer i

! produce format string
  if ((length_of_string.gt.0).and.(length_of_string.lt.10)) then
     write (length1_txt, '(i1)') length_of_string
     format_string = '(iN) '
     format_string(3:3) = length1_txt
  else if ((length_of_string.ge.10).and.(length_of_string.lt.100)) then
     write (length2_txt, '(i2)') length_of_string
     format_string = '(iNN)'
     format_string(3:4) = length2_txt
  else
     print *, "ERROR in CONVERT_INT_TO_TXT_STRING:"
     print *, "incorrect string length requested: ", length_of_string
     stop
  end if

  WRITE (number_txt, format_string) int_number
  number_txt = ADJUSTL(TRIM(number_txt))
  blanks_number = length_of_string - LEN_TRIM(number_txt)
  number_txt = ADJUSTR(number_txt)
  do i = 1, blanks_number
     number_txt(i:i) = '0'
  end do

  convert_int_to_txt_string = number_txt

end function convert_int_to_txt_string
