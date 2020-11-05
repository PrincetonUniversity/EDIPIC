
SUBROUTINE INITIATE_ELECTRON_INJECTION

  USE ParallelOperationValues
  USE ELectronInjection
  USE CurrentProblemValues, ONLY : T_e_eV, BC_flag, N_box_vel, N_max_vel, N_of_particles_cell, delta_t_s, e_Cl, N_plasma_m3, delta_x_m
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL exists
  REAL(8) length_of_injection
  INTEGER left_for_injection

  CHARACTER (77) buf
  INTEGER ierr

  REAL(8) temp
 
  INQUIRE (FILE = 'ssc_einject.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  OPEN (9, FILE = 'ssc_einject.dat')
  IF(exists) THEN

     IF (Rank_of_process.EQ.0) & 
     & PRINT '(2x,"Process ",i3," : The electron injection data file is found. Reading the data file...")', Rank_of_process

     READ (9, '(A77)') buf ! ****************** LEFT WALL, CONSTANT ELECTRON INJECTION *******************")')
     READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
     READ (9, '(A77)') buf ! ---           1 = monoenergetic beam                                      ---")')
     READ (9, '(A77)') buf ! -------d----- 2 = warm maxwellian source ------------------------------------")')
     READ (9, '(7x,i1)') BeamInjectFlag_left
     READ (9, '(A77)') buf ! --dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
     READ (9, '(2x,f10.3)') Delay_of_injection_ns_left

     READ (9, '(A77)') buf ! ====================== ELECTRON INJECTION, PARAMETERS =======================")')
     READ (9, '(A77)') buf ! --dddddd.ddd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
     READ (9, '(2x,f10.3)') Beam_energy_eV_left
     READ (9, '(A77)') buf ! --dddddd.ddd- Injection current [A/m^2] -------------------------------------")')
     READ (9, '(2x,f10.3)') Beam_J_Am2_left

     READ (9, '(A77)') buf ! ****************** RIGHT WALL, CONSTANT ELECTRON INJECTION ******************")')
     READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
     READ (9, '(A77)') buf ! ---           1 = monoenergetic beam                                      ---")')
     READ (9, '(A77)') buf ! -------d----- 2 = warm maxwellian source ------------------------------------")')
     READ (9, '(7x,i1)') BeamInjectFlag_right
     READ (9, '(A77)') buf ! --dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
     READ (9, '(2x,f10.3)') Delay_of_injection_ns_right

     READ (9, '(A77)') buf ! ====================== ELECTRON INJECTION, PARAMETERS =======================")')
     READ (9, '(A77)') buf ! --dddddd.ddd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
     READ (9, '(2x,f10.3)') Beam_energy_eV_right
     READ (9, '(A77)') buf ! --dddddd.ddd- Injection current [A/m^2] -------------------------------------")')
     READ (9, '(2x,f10.3)') Beam_J_Am2_right

     READ (9, '(A77)') buf ! ******************************* SMART TAGS **********************************")')
     READ (9, '(A77)') buf ! ---           0 = turned off (changing velocity sign does not change tag) ---")')
     READ (9, '(A77)') buf ! -------d----- 1 = turned on -------------------------------------------------")')
     READ (9, '(7x,i1)') UseSmartTagsFlag
                 
  ELSE

     BeamInjectFlag_left         = 0           ! injection at the left wall is turned off 
     Delay_of_injection_ns_left  = 0.0_8
     Beam_energy_eV_left         = 0.0_8    
     Beam_J_Am2_left             = 0.0_8

     BeamInjectFlag_right        = 0           ! injection at the right wall is turned off 
     Delay_of_injection_ns_right = 0.0_8
     Beam_energy_eV_right        = 0.0_8    
     Beam_J_Am2_right            = 0.0_8

     UseSmartTagsFlag = 1                      ! smart tags are on

     IF (Rank_of_process.EQ.0) THEN

        PRINT '(/2x,"Process ",i3," : File with the name ssc_einject.dat not found. Use the default settings ...")'

        PRINT '(2x,"Process ",i3," : Create ssc_beaminject.dat file . . .")', Rank_of_process

        WRITE (9, '("****************** LEFT WALL, CONSTANT ELECTRON INJECTION *******************")')
        WRITE (9, '("---           0 = turned off                                              ---")')
        WRITE (9, '("---           1 = monoenergetic beam                                      ---")')
        WRITE (9, '("-------d----- 2 = warm maxwellian source ------------------------------------")')
        WRITE (9, '(7x,i1)') BeamInjectFlag_left
        WRITE (9, '("--dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
        WRITE (9, '(2x,f10.3)') Delay_of_injection_ns_left

        WRITE (9, '("====================== ELECTRON INJECTION, PARAMETERS =======================")')
        WRITE (9, '("--dddddd.ddd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
        WRITE (9, '(2x,f10.3)') Beam_energy_eV_left
        WRITE (9, '("--dddddd.ddd- Injection current [A/m^2] -------------------------------------")')
        WRITE (9, '(2x,f10.3)') Beam_J_Am2_left

        WRITE (9, '("****************** RIGHT WALL, CONSTANT ELECTRON INJECTION ******************")')
        WRITE (9, '("---           0 = turned off                                              ---")')
        WRITE (9, '("---           1 = monoenergetic beam                                      ---")')
        WRITE (9, '("-------d----- 2 = warm maxwellian source ------------------------------------")')
        WRITE (9, '(7x,i1)') BeamInjectFlag_right
        WRITE (9, '("--dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
        WRITE (9, '(2x,f10.3)') Delay_of_injection_ns_right

        WRITE (9, '("====================== ELECTRON INJECTION, PARAMETERS =======================")')
        WRITE (9, '("--dddddd.ddd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
        WRITE (9, '(2x,f10.3)') Beam_energy_eV_right
        WRITE (9, '("--dddddd.ddd- Injection current [A/m^2] -------------------------------------")')
        WRITE (9, '(2x,f10.3)') Beam_J_Am2_right

        WRITE (9, '("******************************* SMART TAGS **********************************")')
        WRITE (9, '("---           0 = turned off (changing velocity sign does not change tag) ---")')
        WRITE (9, '("-------d----- 1 = turned on -------------------------------------------------")')
        WRITE (9, '(7x,i1)') UseSmartTagsFlag
        
     END IF
 
  END IF

  CLOSE (9, STATUS = 'KEEP')

! for the left wall ----------------------------

  IF (BeamInjectFlag_left.EQ.0) THEN

     IF (Rank_of_process.EQ.0) PRINT '(/2x,"The electron injection at the left wall is turned off ...")'
     inject_every_this_many_timesteps_left = 1

  ELSE

     temp = (Beam_J_Am2_left / e_Cl) / (N_plasma_m3 * delta_x_m / (delta_t_s * DBLE(N_of_particles_cell)))

     N_to_inject_total_left = INT(temp)

     IF (N_to_inject_total_left.GT.0) THEN

        const_N_to_inject_by_proc_left = N_to_inject_total_left / (N_of_processes-1)     ! each node will emit at least this number of electrons
        variable_N_to_inject_left = MOD(N_to_inject_total_left, N_of_processes-1)        ! leftover, will be distributed between all client nodes
        inject_every_this_many_timesteps_left = 1

     ELSE IF (temp.GE.0.001_8) THEN  ! let the minimal permitted current is when a macroparticle is injected every 1000 timesteps

        const_N_to_inject_by_proc_left = 0                               ! only one macroparticle will be injected by one (each time different) process 
        variable_N_to_inject_left = 1                                    ! every inject_every_this_many_timesteps_left timesteps
        inject_every_this_many_timesteps_left = INT(1.0_8 / temp) + 1    ! note that skip_timesteps_beam_inject_left >= 2

     ELSE

        const_N_to_inject_by_proc_left = 0          ! if the beam was requested but its current is too small
        variable_N_to_inject_left = 0               ! nothing will be injected
        inject_every_this_many_timesteps_left = 1   !

     END IF

     IF (Rank_of_process.GT.0) THEN
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
        PRINT '(2x,"Process ",i3," will inject at least ",i5," electron macroparticles at the LEFT boundary EACH timestep")', &
             & Rank_of_process, const_N_to_inject_by_proc_left
     ELSE
! for the server process print out the total number of macroparticles to be injected
        PRINT '(/2x,"Constant REQUESTED flux of electron injection from the LEFT WALL is ",e14.7," m^-2 s^-1")', (Beam_J_Am2_left / e_Cl)
        PRINT '(2x,"TOTAL INTEGER number of electron macroparticles to be injected from the LEFT boundary at each timestep is : ",i7)', N_to_inject_total_left
        PRINT '(2x,"TOTAL REAL    number of electron macroparticles to be injected from the LEFT boundary at each timestep is : ",e12.5)', temp
        PRINT '(2x,"COMMON injection leftover distributed between all nodes is : ",i5)', variable_N_to_inject_left
        PRINT '(2x,"The leftover macroparticle(s) will be injected every ",i5," timestep(s)")', inject_every_this_many_timesteps_left
        PRINT '("-------")'
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)   ! place the synchronization point here to print server's output before the clients'
     END IF

     VX_e_beam_left = SQRT(Beam_energy_eV_left / T_e_eV) * N_box_vel
     
     length_of_injection = SQRT(Beam_energy_eV_left / T_e_eV) / N_max_vel
     
     IF (N_to_inject_total_left.GT.0) THEN
        Gap_betw_part_left = length_of_injection / N_to_inject_total_left
     ELSE
        Gap_betw_part_left = 0.0_8
     END IF

! if the delay of the start of injection is negative then it is given in timesteps, assumed to be integer !!!
     IF (Delay_of_injection_ns_left.LT.0.0_8) THEN  
        IF (Rank_of_process.EQ.0) PRINT &
             &'(/2x,"The delay of electron injection is given in timesteps (negative value of Delay_of_injection_ns_left was read from the file)")'
        Delay_of_injection_T_cntr_left = - Delay_of_injection_ns_left
        Delay_of_injection_ns_left = Delay_of_injection_T_cntr_left * delta_t_s * 1.0e9
        
        IF (Rank_of_process.EQ.0) PRINT '(/2x,"The recalculated delay of electron injection is : ", f10.3)', Delay_of_injection_ns_left
     
     ELSE
        Delay_of_injection_T_cntr_left  = 1.0d-9 * Delay_of_injection_ns_left / delta_t_s                               ! Number of time steps
     END IF

  END IF

! for the right wall -------------------------

  IF (BeamInjectFlag_right.EQ.0) THEN

     IF (Rank_of_process.EQ.0) PRINT '(/2x,"The electron injection at the right wall is turned off ...")'
     inject_every_this_many_timesteps_right = 1

  ELSE

     temp = (Beam_J_Am2_right / e_Cl) / (N_plasma_m3 * delta_x_m / (delta_t_s * DBLE(N_of_particles_cell)))

     N_to_inject_total_right = INT(temp)

     IF (N_to_inject_total_right.GT.0) THEN

        const_N_to_inject_by_proc_right = N_to_inject_total_right / (N_of_processes-1)     ! each node will emit at least this number of electrons
        variable_N_to_inject_right = MOD(N_to_inject_total_right, N_of_processes-1)        ! leftover, will be distributed between all client nodes
        inject_every_this_many_timesteps_right = 1

     ELSE IF (temp.GE.0.001_8) THEN  ! let the minimal permitted current is when a macroparticle is injected every 1000 timesteps

        const_N_to_inject_by_proc_right = 0                               ! only one macroparticle will be injected by one (each time different) process 
        variable_N_to_inject_right = 1                                    ! every inject_every_this_many_timesteps_right timesteps
        inject_every_this_many_timesteps_right = INT(1.0_8 / temp) + 1    ! note that skip_timesteps_beam_inject_right >= 2

     ELSE

        const_N_to_inject_by_proc_right = 0          ! if the beam was requested but its current is too small
        variable_N_to_inject_right = 0               ! nothing will be injected
        inject_every_this_many_timesteps_right = 1   !

     END IF

     IF (Rank_of_process.GT.0) THEN
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
        PRINT '(2x,"Process ",i3," will inject at least ",i5," electron macroparticles at the RIGHT boundary EACH timestep")', &
             & Rank_of_process, const_N_to_inject_by_proc_right
     ELSE
! for the server process print out the total number of macroparticles to be injected
        PRINT '(/2x,"Constant REQUESTED flux of electron injection from the RIGHT WALL is ",e14.7," m^-2 s^-1")', (Beam_J_Am2_right / e_Cl)
        PRINT '(2x,"TOTAL INTEGER number of electron macroparticles to be injected from the RIGHT boundary at each timestep is : ",i7)', N_to_inject_total_right
        PRINT '(2x,"TOTAL REAL    number of electron macroparticles to be injected from the RIGHT boundary at each timestep is : ",e12.5)', temp
        PRINT '(2x,"COMMON injection leftover distributed between all nodes is : ",i5)', variable_N_to_inject_right
        PRINT '(2x,"The leftover macroparticle(s) will be injected every ",i5," timestep(s)")', inject_every_this_many_timesteps_right
        PRINT '("-------")'
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)   ! place the synchronization point here to print server's output before the clients'
     END IF

     VX_e_beam_right = SQRT(Beam_energy_eV_right / T_e_eV) * N_box_vel
     
     length_of_injection = SQRT(Beam_energy_eV_right / T_e_eV) / N_max_vel
   
     IF (N_to_inject_total_right.GT.0) THEN
        Gap_betw_part_right = length_of_injection / N_to_inject_total_right
     ELSE
        Gap_betw_part_right = 0.0_8
     END IF

! if the delay of the start of injection is negative then it is given in timesteps, assumed to be integer !!!
     IF (Delay_of_injection_ns_right.LT.0.0_8) THEN  
        IF (Rank_of_process.EQ.0) PRINT &
             &'(/2x,"The delay of electron injection is given in timesteps (negative value of Delay_of_injection_ns_right was read from the file)")'
        Delay_of_injection_T_cntr_right = - Delay_of_injection_ns_right
        Delay_of_injection_ns_right = Delay_of_injection_T_cntr_right * delta_t_s * 1.0e9
        
        IF (Rank_of_process.EQ.0) PRINT '(/2x,"The recalculated delay of electron injection is : ", f10.3)', Delay_of_injection_ns_right
     
     ELSE
        Delay_of_injection_T_cntr_right  = 1.0e-9 * Delay_of_injection_ns_right / delta_t_s                               ! Number of time steps
     END IF

  END IF

END SUBROUTINE INITIATE_ELECTRON_INJECTION 

!-------------------------------------------------------------------------------------------------
! the subroutine injects electrons with the constant rate at both walls
! this subroutine must be called after the final push before rearranging of the particle arays.
! therefore it is not necessary to modify the mesh charge density here
!
SUBROUTINE INJECT_ELECTRONS_AT_WALLS

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : T_cntr
  USE ElectronInjection

  IMPLICIT NONE

  INTEGER var_add, N_clients, T_cntr_flashed, k, T_inject_k

! left wall -------------------------------

  N_clients = (N_of_processes-1) * inject_every_this_many_timesteps_left

  T_cntr_flashed = MOD(T_cntr, N_clients) + 1  ! periodically grows from 1 to N_clients

  IF ((BeamInjectFlag_left.NE.0).AND.(T_cntr.GE.Delay_of_injection_T_cntr_left)) THEN

     var_add = 0
     
     DO k = 0, variable_N_to_inject_left-1                ! this cycle works only if variable_N_to_inject_left>0

! find "timestep" suitable for injection for given client        
        T_inject_k = Rank_of_process * inject_every_this_many_timesteps_left + k
        IF (T_inject_k.GT.N_clients) T_inject_k = T_inject_k - N_clients
        
        IF (T_inject_k.EQ.T_cntr_flashed) THEN
           var_add=1
           EXIT
        END IF
     END DO

     CALL INJECT_N_ELECTRONS_LEFT_WALL(const_N_to_inject_by_proc_left + var_add)

  END IF

! right wall -------------------------------

  N_clients = (N_of_processes-1) * inject_every_this_many_timesteps_right

  T_cntr_flashed = MOD(T_cntr, N_clients) + 1  ! periodically grows from 1 to N_clients

  IF ((BeamInjectFlag_right.NE.0).AND.(T_cntr.GE.Delay_of_injection_T_cntr_right)) THEN

     var_add = 0
     
     DO k = 0, variable_N_to_inject_right-1                ! this cycle works only if variable_N_to_inject_right>0

! find "timestep" suitable for injection for given client        
        T_inject_k = Rank_of_process * inject_every_this_many_timesteps_right + k
        IF (T_inject_k.GT.N_clients) T_inject_k = T_inject_k - N_clients
        
        IF (T_inject_k.EQ.T_cntr_flashed) THEN
           var_add=1
           EXIT
        END IF
     END DO

     CALL INJECT_N_ELECTRONS_RIGHT_WALL(const_N_to_inject_by_proc_right + var_add)

  END IF

END SUBROUTINE INJECT_ELECTRONS_AT_WALLS

!-------------------------------------------------------------------------------------------------
! the subroutine injects electrons with the constant rate at the left wall
! this subroutine must be called after the final push before rearranging of the particle arays.
! therefore it is not necessary to modify the mesh charge density here
!
SUBROUTINE INJECT_N_ELECTRONS_LEFT_WALL(N_to_inject)

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronInjection
  USE Diagnostics, ONLY : Rate_energy_leftemit, Rate_number_leftemit

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER N_to_inject

  REAL(8) x             ! coordinate of refluxed particle (after refluxing)
  REAL(8) vx, vy, vz    ! velocity components of refluxed particle
  REAL(8) v2            ! absolute value of velocity (squared)

  INTEGER ALLOC_ERR, j
  INTEGER left_node, right_node

! left wall -------------------------------

  DO j = 1, N_to_inject

     Q_left = Q_left - Qs(1)

! get the new velocities for electrons -----------------
     vy = 0.0_8      
     vz = 0.0_8          
     IF (BeamInjectFlag_left.EQ.1) THEN
        vx = VX_e_beam_left                        ! cold beam
! calculate the new coordinate 
!
! note, if only one particle must be injected by all processes
! the coordinate of this particle will be jumping depending on the rank of the processes
! which is bad
!  
!  basically, at this point the cold beam is not working properly #############
!
!        x = (Rank_of_process + (N_of_processes - 1) * (j - 1) - 0.5_8) * Gap_betw_part_left
!
        x = well_random_number() * vx * KVx
     ELSE
        CALL GetInjMaxwellVelocity(vx)        ! warm beam, get the new x-velocity according to the maxwell distribution 
        vx = vx * SQRT(Beam_energy_eV_left / T_e_eV)    ! renormalize velocity according to the beam's temperature
! calculate the new coordinate 
        x = 0.0_8
!####        x = well_random_number() * vx * KVx
     END IF
! calculate the squared absolute velocity
     v2 = vx * vx + vy * vy + vz * vz
! Save the parameters of the refluxed electron in the linked list
     IF (ASSOCIATED(Current_electron)) THEN
        N_inject(1)            = N_inject(1) + 1
        Rate_number_leftemit(1) = Rate_number_leftemit(1) + 1             !
        Rate_energy_leftemit(1) = Rate_energy_leftemit(1) + v2            !
        Current_electron%X       = x
        Current_electron%VX      = vx
        Current_electron%VY      = vy
        Current_electron%VZ      = vz
!        Current_electron%AX      = 0.0_8
        Current_electron%Tag     = eTag_Emit_Left                         ! 
        CALL ADD_EMITTED_E_TO_LEFT_DF(vx, vy, vz)
        ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
        IF (ALLOC_ERR.NE.0) THEN
           PRINT '(/2x,"Process ",i3," : Error in INJECT_N_ELECTRONS_LEFT_WALL:")', Rank_of_process
           PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        Current_electron => Current_electron%next
        NULLIFY(Current_electron%next)
     ELSE
        PRINT '(/2x,"Process ",i3," : Error in INJECT_N_ELECTRONS_LEFT_WALL:")', Rank_of_process
        PRINT  '(2x,"Current_electron is NOT associated!")'
        PRINT  '(2x,"The program will be terminated now :(")'
        STOP
     END IF

!       left_node  = INT(x)               
!       right_node = left_node + 1
!       Q_strm_spec(left_node, 1)  = Q_strm_spec(left_node, 1)  + right_node - x    ! s = 1
!       Q_strm_spec(right_node, 1) = Q_strm_spec(right_node, 1) + x - left_node     ! s = 1
     Q_strm_spec(0, 1) = Q_strm_spec(0, 1)  + 1.0_8

  END DO

END SUBROUTINE INJECT_N_ELECTRONS_LEFT_WALL

!-------------------------------------------------------------------------------------------------
! the subroutine injects electrons with the constant rate at the right wall
! this subroutine must be called after the final push before rearranging of the particle arays.
! therefore it is not necessary to modify the mesh charge density here
!
SUBROUTINE INJECT_N_ELECTRONS_RIGHT_WALL(N_to_inject)

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronInjection
  USE Diagnostics, ONLY : Rate_energy_rightemit, Rate_number_rightemit

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER N_to_inject

  REAL(8) x             ! coordinate of refluxed particle (after refluxing)
  REAL(8) vx, vy, vz    ! velocity components of refluxed particle
  REAL(8) v2            ! absolute value of velocity (squared)

  INTEGER ALLOC_ERR, j
  INTEGER left_node, right_node

! right wall -------------------------------

  DO j = 1, N_to_inject

     Q_right =  Q_right - Qs(1)   

! get the new velocities for electrons -----------------
     vy = 0.0_8      
     vz = 0.0_8          
     IF (BeamInjectFlag_right.EQ.1) THEN
        vx = - VX_e_beam_right                        ! cold beam
! calculate the new coordinate 
!
! note, if only one particle must be injected by all processes
! the coordinate of this particle will be jumping depending on the rank of the processes
! which is bad
!  
!  basically, at this point the cold beam is not working properly #############
!
!        x = N_cells - (Rank_of_process + (N_of_processes - 1) * (j - 1) - 0.5_8) * Gap_betw_part_right
        x = N_cells + well_random_number() * vx * KVx
     ELSE
        CALL GetInjMaxwellVelocity(vx)        ! warm beam, get the new x-velocity according to the maxwell distribution 
        vx = - vx * SQRT(Beam_energy_eV_right / T_e_eV)    ! renormalize velocity according to the beam's temperature
! calculate the new coordinate 
        x = DBLE(N_cells)
!####        x = N_cells + well_random_number() * vx * KVx
     END IF
! calculate the squared absolute velocity
     v2 = vx * vx + vy * vy + vz * vz
! Save the parameters of the refluxed electron in the linked list
     IF (ASSOCIATED(Current_electron)) THEN
        N_inject(1)            = N_inject(1) + 1
        Rate_number_rightemit(1) = Rate_number_rightemit(1) + 1             !
        Rate_energy_rightemit(1) = Rate_energy_rightemit(1) + v2            !
        Current_electron%X       = x
        Current_electron%VX      = vx
        Current_electron%VY      = vy
        Current_electron%VZ      = vz
!        Current_electron%AX      = 0.0_8
        Current_electron%Tag     = eTag_Emit_Right                         ! 
        CALL ADD_EMITTED_E_TO_RIGHT_DF(vx, vy, vz)
        ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
        IF (ALLOC_ERR.NE.0) THEN
           PRINT '(/2x,"Process ",i3," : Error in INJECT_N_ELECTRONS_RIGHT_WALL:")', Rank_of_process
           PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        Current_electron => Current_electron%next
        NULLIFY(Current_electron%next)
     ELSE
        PRINT '(/2x,"Process ",i3," : Error in INJECT_N_ELECTRONS_RIGHT_WALL:")', Rank_of_process
        PRINT  '(2x,"Current_electron is NOT associated!")'
        PRINT  '(2x,"The program will be terminated now :(")'
        STOP
     END IF

!! account for the density of emitted particle
!     left_node  = MIN(INT(x), N_cells-1)               
!     right_node = left_node + 1
!     Q_strm_spec(left_node, 1)  = Q_strm_spec(left_node, 1)  + right_node - x    ! s = 1
!     Q_strm_spec(right_node, 1) = Q_strm_spec(right_node, 1) + x - left_node     ! s = 1

     Q_strm_spec(N_cells, 1) = Q_strm_spec(N_cells, 1) + 1.0_8

  END DO

END SUBROUTINE INJECT_N_ELECTRONS_RIGHT_WALL
