
SUBROUTINE INITIATE_BEAM_IN_PLASMA

  USE ParallelOperationValues
  USE BeamInPlasma
  USE CurrentProblemValues, ONLY : T_e_eV, BC_flag, U_ext_V, N_of_particles_cell, delta_t_s, N_cells
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL exists
  INTEGER left_unused
  INTEGER sum_of_N_ebeam_proc
  INTEGER k 

  CHARACTER (77) buf

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER ALLOC_ERR
 
  INQUIRE (FILE = 'ssc_ebeam.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  OPEN (9, FILE = 'ssc_ebeam.dat')
  IF(exists) THEN

     IF (Rank_of_process.EQ.0) & 
     & PRINT '(2x,"Process ",i3," : The electron beam data file is found. Reading the data file...")', Rank_of_process

     READ (9, '(A77)') buf ! ********************** PERIODICITY OF PARTICLES MOTION **********************")')
     READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
     READ (9, '(A77)') buf ! -------d----- 1 = turned on -------------------------------------------------")')
     READ (9, '(7x,i1)') PeriodicBoundaryFlag

     READ (9, '(A77)') buf ! ****************************** BEAM IN PLASMA *******************************")')
     READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
     READ (9, '(A77)') buf ! ---           1 = monoenergetic beam                                      ---")')
     READ (9, '(A77)') buf ! -------d----- 2 = warm beam -------------------------------------------------")')
     READ (9, '(7x,i1)') eBeamInPlasmaFlag
     READ (9, '(A77)') buf ! --dddddd.ddd- Beam appears at this moment, [ns] (timesteps if < 0) ----------")')
     READ (9, '(2x,f10.3)') Start_ebeam_ns

     READ (9, '(A77)') buf ! ============================= BEAM PARAMETERS ===============================")')
     READ (9, '(A77)') buf ! --dddddd.ddd- Energy, [eV] --------------------------------------------------")')
     READ (9, '(2x,f10.3)') Energy_ebeam_eV
     READ (9, '(A77)') buf ! --dddddd.ddd- Temperature (for warm beams only), [eV] ----------------------")')
     READ (9, '(2x,f10.3)') Te_ebeam_eV
     READ (9, '(A77)') buf ! ------d.dddd- Relative density (0 < alfa < 1) -------------------------------")')
     READ (9, '(6x,f6.4)') Alfa_ebeam
                 
  ELSE

     PeriodicBoundaryFlag = 0           ! periodic motion turned off

     eBeamInPlasmaFlag    = 0           ! the beam is turned off 
     Start_ebeam_ns       = 0.0_8
     Energy_ebeam_eV      = 0.0_8    
     Te_ebeam_eV          = 0.0_8
     Alfa_ebeam           = 0.0_8

     IF (Rank_of_process.EQ.0) THEN

        PRINT '(/2x,"Process ",i3," : File with the name ssc_ebeam.dat not found. Use the default settings ...")', Rank_of_process

        PRINT '(2x,"Process ",i3," : Create ssc_ebeam.dat file . . .")', Rank_of_process

        WRITE (9, '("********************** PERIODICITY OF PARTICLES MOTION **********************")')
        WRITE (9, '("---           0 = turned off                                              ---")')
        WRITE (9, '("-------d----- 1 = turned on -------------------------------------------------")')
        WRITE (9, '(7x,i1)') PeriodicBoundaryFlag

        WRITE (9, '("****************************** BEAM IN PLASMA *******************************")')
        WRITE (9, '("---           0 = turned off                                              ---")')
        WRITE (9, '("---           1 = monoenergetic beam                                      ---")')
        WRITE (9, '("-------d----- 2 = warm beam -------------------------------------------------")')
        WRITE (9, '(7x,i1)') eBeamInPlasmaFlag
        WRITE (9, '("--dddddd.ddd- Beam appears at this moment, [ns] (timesteps if < 0) ----------")')
        WRITE (9, '(2x,f10.3)') Start_ebeam_ns

        WRITE (9, '("============================= BEAM PARAMETERS ===============================")')
        WRITE (9, '("--dddddd.ddd- Energy, [eV] --------------------------------------------------")')
        WRITE (9, '(2x,f10.3)') Energy_ebeam_eV
        WRITE (9, '("--dddddd.ddd- Temperature (for warm beams only), [eV] ----------------------")')
        WRITE (9, '(2x,f10.3)') Te_ebeam_eV
        WRITE (9, '("------d.dddd- Relative density (0 < alfa < 1) -------------------------------")')
        WRITE (9, '(6x,f6.4)') Alfa_ebeam
        
     END IF
 
  END IF

  CLOSE (9, STATUS = 'KEEP')

  IF (PeriodicBoundaryFlag.EQ.1) THEN

     IF (Rank_of_process.EQ.0) PRINT '(/2x,"NOTE, THE MOTION OF PARTICLES IS PERIODIC ...")'     

! check the compatibility of the beam with the potential boundary conditions
     IF ((BC_flag.NE.0).OR.(U_ext_V.NE.0.0_8)) THEN
        IF (Rank_of_process.EQ.0) THEN
           PRINT '(/2x,"Error: the non-periodical field boundary conditions are incompatible with the periodic dynamics of plasma")'
           PRINT '( 2x,"Choose the fixed zero wall potentials...")'
           PRINT '( 2x,"Program will be terminated now, sorry :(")'
        END IF
        STOP
     END IF

  END IF

  IF (eBeamInPlasmaFlag.EQ.0) THEN

     IF (Rank_of_process.EQ.0) PRINT '(/2x,"The electron beam in plasma is turned off ...")'

  ELSE

! check the compatibility of the beam with the potential boundary conditions
!     IF ((BC_flag.NE.0).OR.(U_ext_V.NE.0.0_8)) THEN
!        IF (Rank_of_process.EQ.0) THEN
!           PRINT '(/2x,"Error: the non-periodical boundary conditions are incompatible with the use of the beam in plasma")'
!           PRINT '( 2x,"Choose the fixed zero wall potentials...")'
!           PRINT '( 2x,"Program will be terminated now, sorry :(")'
!        END IF
!        STOP
!     END IF

! check the compatibility of the beam with the conditions of particles dynamics
     IF (PeriodicBoundaryFlag.NE.1) THEN
        IF (Rank_of_process.EQ.0) THEN
           PRINT '(/2x,"Error: the non-periodical dynamics is incompatible with the use of the beam in plasma")'
           PRINT '( 2x,"Choose the periodical particles dynamics ...")'
           PRINT '( 2x,"Program will be terminated now, sorry :(")'
        END IF
        STOP
     END IF

! determine the total number of beam particles     
     N_ebeam_total = Alfa_ebeam * (N_of_particles_cell * N_cells) 
     N_ebeam_proc  = 0

     IF (Rank_of_process.GT.0) THEN
! for each CLIENT process determine the number of beam macroparticles to be managed by this process
        left_unused = N_ebeam_total  
        DO
           IF ((left_unused - Rank_of_process).GE.0) THEN
              N_ebeam_proc = N_ebeam_proc + 1
              left_unused  = left_unused - N_of_processes + 1
           ELSE
              EXIT
           END IF
        END DO      
! send the number of beam macroparticles to the server process
        CALL MPI_SEND(N_ebeam_proc, 1, MPI_INTEGER, 0, 30, MPI_COMM_WORLD, ierr)

        PRINT '(2x,"Process ",i3," will manage ",i5," macroparticles of the electron beam")', Rank_of_process, N_ebeam_proc

     ELSE

! for the SERVER process print out the total number of macroparticles to be injected
        PRINT '(2x,"Total number of electron beam macroparticles is : ",i7)', N_ebeam_total

        IF (.NOT.ALLOCATED(N_ebeam_proc_array)) THEN
           ALLOCATE(N_ebeam_proc_array(1:(N_of_processes-1)), STAT=ALLOC_ERR)
           IF(ALLOC_ERR.NE.0)THEN
              PRINT '(/2x,"Process ",i3," : Error in ALLOCATE N_ebeam_proc_array !!!")', Rank_of_process
              PRINT  '(2x,"Program will be terminated now :(")'
              STOP
           END IF
        END IF
        N_ebeam_proc_array = 0
        sum_of_N_ebeam_proc = 0
        DO k = 1, N_of_processes - 1            
! get the number of beam macroparticles from each process
           CALL MPI_PROBE(k, 30, MPI_COMM_WORLD, stattus, ierr)
           CALL MPI_RECV(N_ebeam_proc_array(k), 1, MPI_INTEGER, k, 30, MPI_COMM_WORLD, stattus, ierr)
           PRINT '(2x,"Process ",i3," reported about ",i5," macroparticles of the electron beam")', k, N_ebeam_proc_array(k)
           sum_of_N_ebeam_proc = sum_of_N_ebeam_proc + N_ebeam_proc_array(k)
        END DO

        IF (sum_of_N_ebeam_proc.NE.N_ebeam_total) THEN
           PRINT '(/2x,"Process ",i3," : Error in INITIATE_BEAM_IN_PLASMA !!!")', Rank_of_process
           PRINT  '(2x,"Client processes use ",i7," beam particles, which is NOT the desired value ",i7)', sum_of_N_ebeam_proc, N_ebeam_total
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        
     END IF

! if the delay of the start of injection is negative then it is given in timesteps, assumed to be integer !!!
     IF (Start_ebeam_ns.LT.0.0_8) THEN  
        IF (Rank_of_process.EQ.0) PRINT &
             &'(/2x,"The moment of the electron beam appearance is given in timesteps (negative Start_ebeam_ns was given)")'
        Start_ebeam_T_cntr = - Start_ebeam_ns
        Start_ebeam_ns = Start_ebeam_T_cntr * delta_t_s * 1.0e9
        
        IF (Rank_of_process.EQ.0) PRINT '(/2x,"The recalculated moment of the electron beam appearance is : ", f10.3)', Start_ebeam_ns
     
     ELSE
        Start_ebeam_T_cntr = 1.0e-9 * Start_ebeam_ns / delta_t_s                               ! Number of time steps
     END IF

  END IF


END SUBROUTINE INITIATE_BEAM_IN_PLASMA

!------------------------------
!
SUBROUTINE ALLOCATE_BEAM_ARRAYS

  USE ParallelOperationValues
  USE BeamInPlasma
  IMPLICIT NONE
  INTEGER ALLOC_ERR

! for the SERVER process 
  IF (Rank_of_process.EQ.0) THEN  

! allocate arrays using the total number of beam particles
     IF (.NOT.ALLOCATED(X_ebeam_part)) THEN
        ALLOCATE(X_ebeam_part(1:N_ebeam_total), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE VX_ebeam_part !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
     
     IF (.NOT.ALLOCATED(VX_ebeam_part)) THEN
        ALLOCATE(VX_ebeam_part(1:N_ebeam_total), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE VX_ebeam_part !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF    

! for the CLIENT processes
  ELSE

! allocate arrays using the local number of beam particles
     IF (.NOT.ALLOCATED(X_ebeam_part)) THEN
        ALLOCATE(X_ebeam_part(1:N_ebeam_proc), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE VX_ebeam_part !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
     
     IF (.NOT.ALLOCATED(VX_ebeam_part)) THEN
        ALLOCATE(VX_ebeam_part(1:N_ebeam_proc), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE VX_ebeam_part !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF    

  END IF

END SUBROUTINE ALLOCATE_BEAM_ARRAYS

!----------------------------------
!
SUBROUTINE DEALLOCATE_BEAM_ARRAYS

  USE ParallelOperationValues
  USE BeamInPlasma
  IMPLICIT NONE
  INTEGER DEALLOC_ERR

  IF (ALLOCATED(X_ebeam_part)) THEN
     DEALLOCATE(X_ebeam_part, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE X_ebeam_part !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(VX_ebeam_part)) THEN
     DEALLOCATE(VX_ebeam_part, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE VX_ebeam_part !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(N_ebeam_proc_array)) THEN
     DEALLOCATE(N_ebeam_proc_array, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE N_ebeam_proc_array !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

END SUBROUTINE DEALLOCATE_BEAM_ARRAYS


!-------------------------------------------------------------------------------------------------
! the subroutine works only once during the program execution - 
! it creates the beam particles and place them into the list of particles to be injected
!
SUBROUTINE START_BEAM_IN_PLASMA

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BeamInPlasma
  USE Diagnostics, ONLY : Rate_energy_leftemit, Rate_number_leftemit
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER k_start, k_fin
  INTEGER i, j
  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER ALLOC_ERR
  REAL(8) v2

!  real(8) bufer(1:1000000)

  IF ((eBeamInPlasmaFlag.EQ.0).OR.(T_cntr.NE.Start_ebeam_T_cntr)) RETURN

  CALL ALLOCATE_BEAM_ARRAYS

! for the SERVER process 
  IF (Rank_of_process.EQ.0) THEN  

! initiate beam particles
     CALL BEAM_DISTRIB_IN_SPACE
     CALL BEAM_DISTRIB_IN_VELOCITY

     DO i = 1, 3
        CALL SHUFFLE_BEAM            ! use this later, when the warm beam will be implemented
     END DO
! transmit data to the client processes
     k_start = 1
     k_fin   = 0
     DO i = 1, N_of_processes - 1

        k_fin   = k_fin + N_ebeam_proc_array(i)

!        bufer(1:N_ebeam_proc_array(i)) = X_ebeam_part(k_start:k_fin)
        CALL MPI_SEND(X_ebeam_part(k_start:k_fin), N_ebeam_proc_array(i), MPI_DOUBLE_PRECISION, i, 31, MPI_COMM_WORLD, ierr) 

!        bufer(1:N_ebeam_proc_array(i)) = VX_ebeam_part(k_start:k_fin)
        CALL MPI_SEND(VX_ebeam_part(k_start:k_fin), N_ebeam_proc_array(i), MPI_DOUBLE_PRECISION, i, 32, MPI_COMM_WORLD, ierr) 

        k_start = k_fin + 1

     END DO

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ELSE                                                        ! client process >>>>

     CALL MPI_PROBE(0, 31, MPI_COMM_WORLD, stattus, ierr)
        
     CALL MPI_RECV(X_ebeam_part, N_ebeam_proc, MPI_DOUBLE_PRECISION, 0, 31, MPI_COMM_WORLD, stattus, ierr)
!        X_of_spec(s)%part(1:N_ebeam_proc)   = bufer(1:N_ebeam_proc)
     
     CALL MPI_RECV(VX_ebeam_part, N_ebeam_proc, MPI_DOUBLE_PRECISION, 0, 32, MPI_COMM_WORLD, stattus, ierr)
!        VX_ebeam_part(1:N_ebeam_proc)  = bufer(1:N_ebeam_proc)

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! place beam particles into the list of particles to be injected
     DO j = 1, N_ebeam_proc

! calculate the squared absolute velocity
        v2 = VX_ebeam_part(j) * VX_ebeam_part(j)
! Save the parameters of the refluxed electron in the linked list
        IF (ASSOCIATED(Current_electron)) THEN
           N_inject(1)            = N_inject(1) + 1
!           electron_reemit_count  = electron_reemit_count + 1
!           electron_reemit_energy = electron_reemit_energy + v2
           Rate_number_leftemit(1) = Rate_number_leftemit(1) + 1             !
           Rate_energy_leftemit(1) = Rate_energy_leftemit(1) + v2            !
           Current_electron%X       =  X_ebeam_part(j)
           Current_electron%VX      = VX_ebeam_part(j)
           Current_electron%VY      = 0.0_8
           Current_electron%VZ      = 0.0_8
!           Current_electron%AX      = 0.0_8
           Current_electron%Tag     = eTag_Emit_Left                         ! 
!           CALL ADD_EMITTED_E_TO_LEFT_DF(vx, vy, vz)
           ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
           IF (ALLOC_ERR.NE.0) THEN
              PRINT '(/2x,"Process ",i3," : Error in START_BEAM_IN_PLASMA:")', Rank_of_process
              PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
              PRINT  '(2x,"Program will be terminated now :(")'
              STOP
           END IF
           Current_electron => Current_electron%next
           NULLIFY(Current_electron%next)
        ELSE
           PRINT '(/2x,"Process ",i3," : Error in INJECT_BEAM_ELECTRONS:")', Rank_of_process
           PRINT  '(2x,"Current_electron is NOT associated!")'
           PRINT  '(2x,"The program will be terminated now :(")'
           STOP
        END IF

! add the same number of zero velocity ions to the linked list for injected ions
        IF (ASSOCIATED(Current_ion)) THEN
           N_inject(2)            = N_inject(2) + 1
!           electron_reemit_count  = electron_reemit_count + 1
!           electron_reemit_energy = electron_reemit_energy + v2
           Rate_number_leftemit(2) = Rate_number_leftemit(2) + 1             !
           Rate_energy_leftemit(2) = Rate_energy_leftemit(2) !+ v2            !
           Current_ion%X       =  X_ebeam_part(j)
           Current_ion%VX      = 0.0_8
           Current_ion%VY      = 0.0_8
           Current_ion%VZ      = 0.0_8
!           Current_ion%AX      = 0.0_8
           Current_ion%Tag     = 0 !eTag_Emit_Left                         ! 
!           CALL ADD_EMITTED_E_TO_LEFT_DF(vx, vy, vz)
           ALLOCATE(Current_ion%next, STAT = ALLOC_ERR)
           IF (ALLOC_ERR.NE.0) THEN
              PRINT '(/2x,"Process ",i3," : Error in START_BEAM_IN_PLASMA:")', Rank_of_process
              PRINT  '(2x,"Error in ALLOCATE Current_ion%next !!!")'
              PRINT  '(2x,"Program will be terminated now :(")'
              STOP
           END IF
           Current_ion => Current_ion%next
           NULLIFY(Current_ion%next)
        ELSE
           PRINT '(/2x,"Process ",i3," : Error in INJECT_BEAM_IONS:")', Rank_of_process
           PRINT  '(2x,"Current_ion is NOT associated!")'
           PRINT  '(2x,"The program will be terminated now :(")'
           STOP
        END IF

     END DO

  END IF

  CALL DEALLOCATE_BEAM_ARRAYS

END SUBROUTINE START_BEAM_IN_PLASMA

!=======================================================================
! This subroutine uniformly distributed the beam electrons in space 
! 
SUBROUTINE BEAM_DISTRIB_IN_SPACE

  USE CurrentProblemValues, ONLY : N_cells
  USE BeamInPlasma
  IMPLICIT NONE

  REAL(8) x_left
  INTEGER j
  REAL(8) Gap_betw_part              ! distance between neighbour beam particles

  IF (N_ebeam_total.GT.0) THEN
     Gap_betw_part = DBLE(N_cells) / N_ebeam_total
  ELSE
     Gap_betw_part = 0.0_8
  END IF

  x_left = 0.0_8              ! we start distribution at the left boundary

  DO j = 1, N_ebeam_total
     X_ebeam_part(j) = x_left + Gap_betw_part * (j-1)
  END DO

END SUBROUTINE BEAM_DISTRIB_IN_SPACE

!=======================================================================
! This subroutine distributes the beam electrons in velocity
!
SUBROUTINE BEAM_DISTRIB_IN_VELOCITY

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : T_e_eV, N_box_vel
  USE BeamInPlasma
  IMPLICIT NONE

  REAL(8) vx
  INTEGER n                          ! the counter of the particles
  REAL(8) VX_ebeam                   ! dimensionless initial flow velocity of electron beam

  VX_ebeam = SQRT(Energy_ebeam_eV / T_e_eV) * N_box_vel

  DO n = 1, N_ebeam_total

     IF (eBeamInPlasmaFlag.EQ.1) THEN
        vx = VX_ebeam                        ! cold beam
     ELSE
     ! in fact this is not a good idea, better do not use warm beams yet
        CALL GetMaxwellVelocity(vx)        ! warm beam, get the new x-velocity according to the maxwell distribution 
        vx = SQRT(VX_ebeam**2 + vx * ABS(vx) * (Te_ebeam_eV / T_e_eV))    ! renormalize velocity according to the beam's temperature
        
     END IF

     VX_ebeam_part(n) = vx

  END DO

END SUBROUTINE BEAM_DISTRIB_IN_VELOCITY


!==========================================================
! This subroutine MUST be called AFTER the subroutines
!  BEAM_DISTRIB_IN_SPACE
!  BEAM_DISTRIB_IN_VELOCITY
! in order to avoid the correlation
SUBROUTINE SHUFFLE_BEAM
  
  USE ParallelOperationValues
  USE BeamInPlasma

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) tmp_value
  INTEGER i, random_i

!!!return     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  PRINT '(6x,"Process ",i3," : Shuffling the beam particles arrays ...")', Rank_of_process

! First of all we shuffle the coordinates
  DO i = 1, N_ebeam_total  
     tmp_value  = X_ebeam_part(i) 
! We use the unconditional DO loop here in order to avoid the situation when
! random_i < 1 or random_i > N_part
     DO
! The new random value will be retrieved until it will satisfy the desired range
! 1 <= random_i <= N_of_part_initial
        random_i = well_random_number() * N_ebeam_total
        IF ((random_i.GE.1).AND.(random_i.LE.N_ebeam_total)) THEN
           X_ebeam_part(i) = X_ebeam_part(random_i)
           X_ebeam_part(random_i) = tmp_value
           EXIT
        END IF
     END DO
  END DO
! Shuffle the X-velocity
  DO i = 1, N_ebeam_total  
     tmp_value = VX_ebeam_part(i) 
     DO
        random_i = well_random_number() * N_ebeam_total
        IF ((random_i.GE.1).AND.(random_i.LE.N_ebeam_total)) THEN
           VX_ebeam_part(i)        = VX_ebeam_part(random_i)
           VX_ebeam_part(random_i) = tmp_value
           EXIT
        END IF
     END DO
  END DO

END SUBROUTINE SHUFFLE_BEAM
