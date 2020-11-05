!===========================================
PROGRAM MainProg

  USE CurrentProblemValues
  USE ParallelOperationValues

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  REAL(8) start, finish

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, Rank_of_process, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_of_processes, ierr)

!  print *, 'my rank is ', Rank_of_process

  CALL INITIATE_PARAMETERS         ! all, with some differences, distribute particles over nodes here

  CALL INITIATE_MC_COLLISIONS
  CALL INITIATE_SE_EMISSION        ! all, with some differences
  CALL INITIATE_ELECTRON_INJECTION
  CALL INITIATE_BEAM_IN_PLASMA
  CALL INITIATE_COULOMB_COLLISIONS

  CALL INITIATE_DIAGNOSTICS

  IF (Rank_of_process.EQ.0) THEN 
     CALL PREPARE_TIME_DEPENDENCE_DATAFILES         ! server only
!restore this later     CALL INITIATE_TEST_PARTICLES                   ! server only
  END IF

  start = MPI_WTIME()

  DO T_cntr = Start_T_cntr, Max_T_cntr

!if (Rank_of_process.eq.0) print '(">>> ",i8," >>> ",i8," >>> ",i8)', T_cntr, N_part(1), N_part(2)

!if (Rank_of_process.eq.0) start = MPI_WTIME()

!print *, T_cntr, Rank_of_process, N_part(1), N_part(2)

     IF (T_cntr.NE.Start_T_cntr) THEN
        CALL CALCULATE_STR_CHARGE_DENSITY                                   ! server and clients, with difference
        CALL CALCULATE_STR_LONG_ELECTR_FIELD                                ! server and clients, with difference   

        IF (Rank_of_process.GT.0) THEN 
           CALL FINAL_PUSH                                                  ! clients
           CALL PROCESS_Z_SHIFT                                             ! clients
           CALL INJECT_ADDITIONAL_EI_PAIRS_LEFT_WALL                        ! clients
           CALL INJECT_ELECTRONS_AT_WALLS                                   ! clients
        ELSE
!restore this later           CALL TEST_PARTICLES_COMBINED_PUSH_AND_SAVE                       ! server
        END IF

        CALL MANAGE_WALL_CHARGE_DENSITIES                                   ! server and clients, with difference

        CALL PROCESS_COLLISIONS_WITH_ATOMS                                  ! server and clients, with difference 

        CALL PROCESS_COLLISIONS_EL_EL_LNGV                                  ! server and clients, with difference

        CALL START_BEAM_IN_PLASMA                                           ! server and clients, with difference, only one time

        IF (Rank_of_process.GT.0) CALL RESIZE_PARTICLE_ARRAYS               ! clients

!        CALL CLEAR_TAGS                                 ! clients

        IF (Rank_of_process.GT.0) CALL CLEAR_ELECTRON_LIST                  ! clients
        IF (Rank_of_process.GT.0) CALL CLEAR_ION_LIST                       ! clients

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)                              ! all, synchronization   !!!! comment me !!!!

     END IF

     CALL DO_DIAGNOSTICS_STEP_1                                             ! server and clients, with differences
     CALL CREATE_SNAPSHOT                                                   ! server and clients, with differences    
     CALL DO_DIAGNOSTICS_STEP_2                                             ! server and clients, with differences

     IF (Rank_of_process.EQ.0) N_part = 0                                   ! server, for correct operation of checkpoints

     CALL SAVE_CHECKPOINT_MPIIO                                             ! server and clients, with differences 

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)                                 ! all, synchronization

     IF (Rank_of_process.GT.0) THEN 
        CALL PREALLOCATE_ELECTRON_LIST                                      ! clients
        CALL PREALLOCATE_ION_LIST                                           ! clients
        CALL PREDICTING_PUSH                                                ! clients
     ELSE
!restore this later        CALL ACTIVATE_TEST_PARTICLE_RUN                                     ! server
     END IF

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)                                 ! all, synchronization

!if (Rank_of_process.eq.0) then
!   finish = MPI_WTIME()
!   if (mod(T_cntr,17).eq.0) PRINT '(2x,"**** cycle ",i8," time is  : ", f11.6," sec")', T_cntr, finish - start
!end if

  END DO

  finish = MPI_WTIME()

  PRINT '(2x,"**** Process ",i3" : Simulation time is  : ", f11.3," sec")', Rank_of_process, finish - start

  IF (Rank_of_process.GT.0) CALL CLEAR_ELECTRON_LIST                        ! clients
  IF (Rank_of_process.GT.0) CALL CLEAR_ION_LIST                             ! clients 

  CALL REMOVE_MESH_ARRAYS                               ! all
  IF (Rank_of_process.GT.0) CALL REMOVE_PARTICLE_DYNAM_ARRAYS               ! clients (server calls it in INITIATE_PARAMETERS : INITIATE_PARTICLE_DYNAM_ARRAYS)

  CALL REMOVE_CRSECT_ARRAYS                              
  CALL REMOVE_COLLSPEC_ARRAYS                              
  CALL REMOVE_COLLISION_ARRAYS
  CALL REMOVE_LANGEVIN_ARRAYS

  CALL FINISH_SNAPSHOTS                                 
  CALL FINISH_DIAGNOSTICS
!restore this later  IF (Rank_of_process.EQ.0) CALL REMOVE_TEST_PARTICLE_ARRAYS                ! server only

  call MPI_FINALIZE(ierr)

END PROGRAM MainProg


