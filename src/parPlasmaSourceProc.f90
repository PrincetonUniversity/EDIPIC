!--------------------------------------------------------------------
! This subrputine is called in subroutine INITIATE_SE_EMISSION
!
SUBROUTINE INITIATE_ADDITIONAL_CONSTANT_INJECTION

  USE ParallelOperationValues
  USE SEEmission, ONLY : const_add_N_to_inject_by_proc, variable_add_N_to_inject
  USE CurrentProblemValues, ONLY : N_distrib_m3, e_Cl, Tx_e_eV, Tz_e_eV, T_i_eV, M_i_amu, delta_x_m, delta_t_s, N_plasma_m3, N_of_particles_cell
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  REAL(8) G_Bohm_m2s

  INTEGER total_to_inject
 
!  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! calculate the Bohm ion flux for the initial plasma bulk parameters

  G_Bohm_m2s = 0.61_8 * N_distrib_m3 * SQRT(e_Cl * (Tx_e_eV + T_i_eV) / (M_i_amu * 1.660565d-27))   ! amu_kg=1.660565d-27

  total_to_inject = G_Bohm_m2s / (N_plasma_m3 * delta_x_m / (delta_t_s * DBLE(N_of_particles_cell)))

  const_add_N_to_inject_by_proc = total_to_inject / (N_of_processes-1)                                       ! each node will emit at least this number of particles

  variable_add_N_to_inject = MOD(total_to_inject, N_of_processes-1)           ! leftover, will be distributed between all client nodes
  
  IF (Rank_of_process.GT.0) THEN
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
     PRINT '(2x,"Process ",i3," will inject at least ",i5," e-i pairs at the LEFT boundary each timestep")', Rank_of_process, const_add_N_to_inject_by_proc
  ELSE
! for the server process print out the total number of macroparticles to be injected
     PRINT '(/2x,"Constant Bohm ion flux is ",e12.5," m^-2 s^-1")', G_Bohm_m2s
     PRINT '(2x,"TOTAL number of electron-ion macroparticle pairs to be injected from the LEFT boundary at each timestep is : ",i7)', total_to_inject
     PRINT '(2x,"COMMON injection leftover distributed between all nodes is : ",i5)', variable_add_N_to_inject
     PRINT '("-------")'
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)   ! place the synchronization point here to print server's output before the clients'
  END IF

END SUBROUTINE INITIATE_ADDITIONAL_CONSTANT_INJECTION

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE LEFT_WALL_REFLUX_THERMAL(x_after, v_after, s)

  USE ParallelOperationValues
  USE SEEmission
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_emit, Rate_energy_leftemit, Rate_number_leftemit
  IMPLICIT NONE

  REAL(8) x_after       ! coordinate after collision, we must restore the initial particle position
  REAL(8) v_after       ! velocity after collision, we must restore the initial particle position
  INTEGER s             ! species

!  REAL RAN

  REAL(8) x             ! coordinate of refluxed particle (after refluxing)
  REAL(8) vx, vy, vz    ! velocity components of refluxed particle
  REAL(8) v2            ! absolute value of velocity (squared)

  INTEGER ALLOC_ERR
  INTEGER left_node, right_node

! restore the surface charge density, because it is modified in the dynamic procedures
  Q_left = 0    

  IF (s.EQ.1) THEN                        ! for ELECTRONs
! get the new velocities according to the maxwell distribution
     CALL GetInjMaxwellVelocity(vx) 
     CALL GetMaxwellVelocity(vy)  
     CALL GetMaxwellVelocity(vz)     
! introduce the anisotropy
     vx = vx * SQRT(Tx_e_eV / T_e_eV) + Vx_e_drift
     vy = vy * SQRT(Tz_e_eV / T_e_eV) + Vy_e_drift
     vz = vz * SQRT(Tz_e_eV / T_e_eV)
! calculate the squared absolute velocity
     v2 = vx * vx + vy * vy + vz * vz
! calculate the new coordinate                  x_after = x + v_after * KVx
!     x = x_after - v_after * KVx !RAN(I_random_seed) * vx * KVx         ! #################### #################### ###################### ################
!@@@     x = RAN(I_random_seed) * vx * KVx
     x = 0.0_8

     IF (x.LT.0.0_8) THEN
        PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLUX_THERMAL:")', Rank_of_process
        PRINT  '(2x,"Wrong coordinate of the refluxed particle :",e12.5," , species ",i2)', x, s
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF  

! Save the parameters of the refluxed electron in the linked list
     IF (ASSOCIATED(Current_electron)) THEN
        N_inject(1)            = N_inject(1) + 1
        electron_reflux_count  = electron_reflux_count + 1
        electron_reflux_energy = electron_reflux_energy + v2
        Rate_number_leftemit(1) = Rate_number_leftemit(1) + 1             !
        Rate_energy_leftemit(1) = Rate_energy_leftemit(1) + v2            !
        Current_electron%X       = x
        Current_electron%VX      = vx
        Current_electron%VY      = vy
        Current_electron%VZ      = vz
!        Current_electron%AX      = 0.0_8
        Current_electron%Tag     = 0 !eTag_Emit_Left                ! mark the refluxed electron
        CALL ADD_EMITTED_E_TO_LEFT_DF(vx, vy, vz)
        ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
        IF (ALLOC_ERR.NE.0) THEN
           PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLUX_THERMAL:")', Rank_of_process
           PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        Current_electron => Current_electron%next
        NULLIFY(Current_electron%next)
     ELSE
        PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLUX_THERMAL:")', Rank_of_process
        PRINT  '(2x,"Current_electron is NOT associated!")'
        PRINT  '(2x,"The program will be terminated now :(")'
        STOP
     END IF

  ELSE IF (s.EQ.2) THEN                     ! for IONs
! get the new velocities according to the maxwell distribution
     CALL GetInjMaxwellVelocity(vx) 
     CALL GetMaxwellVelocity(vy)  
     CALL GetMaxwellVelocity(vz)     
! renormalize velocity for ions
     vx = vx * VT(2) / VT(1) + Vx_i_drift
     vy = vy * VT(2) / VT(1) + Vy_i_drift
     vz = vz * VT(2) / VT(1)
! calculate the squared absolute velocity
     v2 = vx * vx + vy * vy + vz * vz
! calculate the new coordinate
!     x = 0.0_8 !RAN(I_random_seed) * vx * KVx                           ! #################### #################### ###################### ################
!     x = x_after - v_after * KVx !RAN(I_random_seed) * vx * KVx         ! #################### #################### ###################### ################
!@@@     x = RAN(I_random_seed) * vx * KVx
     x = 0.0_8

     IF (x.LT.0.0_8) THEN
        PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLUX_THERMAL:")', Rank_of_process
        PRINT  '(2x,"Wrong coordinate of the refluxed particle :",e12.5," , species ",i2)', x, s
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF

! Save the parameters of the refluxed ion in the linked list
     IF (ASSOCIATED(Current_ion)) THEN
        N_inject(2)       = N_inject(2) + 1
        ion_reflux_count  = ion_reflux_count + 1
        ion_reflux_energy = ion_reflux_energy + v2
        Rate_number_leftemit(2) = Rate_number_leftemit(2) + 1                !
        Rate_energy_leftemit(2) = Rate_energy_leftemit(2) + v2               !
        Current_ion%X       = x
        Current_ion%VX      = vx
        Current_ion%VY      = vy
        Current_ion%VZ      = vz
!        Current_ion%AX      = 0.0_8
        Current_ion%Tag     = 0            ! later this value can be modified to mark the refluxed ion
        ALLOCATE(Current_ion%next, STAT = ALLOC_ERR)
        IF (ALLOC_ERR.NE.0) THEN
           PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLUX_THERMAL:")', Rank_of_process
           PRINT  '(2x,"Error in ALLOCATE Current_ion%next !!!")'
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        Current_ion => Current_ion%next
        NULLIFY(Current_ion%next)
     ELSE
        PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLUX_THERMAL:")', Rank_of_process
        PRINT  '(2x,"Current_ion is NOT associated!")'
        PRINT  '(2x,"The program will be terminated now :(")'
        STOP
     END IF

  END IF

! Below we modify the new array of streaming density "Q_strm_new" for the first group (either electron or ion) near the left wall.
! The charge of the particle, which was refluxed, is not accounted in the calling dynamic procedure PRE_PUSH_ACTIVE_BLOCKS_LNGV. 
! This (primary) particle was treated like the particle colliding with the left wall and producing only the surface charge on the left wall as usual. 
! However, no surface charge can be accumulated at the plasma-plasma boundary, therefore  the surface charge at the left wall was set to zero (see above). 
! As far as the particle after reefluxing returned into the simulated area, its charge must be accounted as the volume charge, 
! i.e. distributed among the two enclosing nodes. 
! If the refluxing event is requested in the dynamic procedure FINAL_PUSH_FOR_ACTIVE_BLOCKS, the modification of "Q_strm_new" and the related
! modification of "Q_strm_old" does not affect the resulting charge distribution, which is used for the field calculation. 
! We do not distinguish the two possible calls of the LEFT_WALL_REFLUX_THERMAL for simplicity.

  left_node  = INT(x)
  right_node = left_node + 1
  Q_strm_spec(left_node, s)  = Q_strm_spec(left_node, s)  + right_node - x
  Q_strm_spec(right_node, s) = Q_strm_spec(right_node, s) + x - left_node  

END SUBROUTINE LEFT_WALL_REFLUX_THERMAL

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE LEFT_WALL_REFLUX_SPECULAR(x, vx, vy, vz, v2, s)

  USE ParallelOperationValues
  USE SEEmission
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_emit, Rate_energy_leftemit, Rate_number_leftemit
  IMPLICIT NONE

  REAL(8) x             ! coordinate of incident particle, must be beyond the plasma boundaries 
  REAL(8) vx, vy, vz    ! velocity components of incident particle
  REAL(8) v2            ! squared absolute value of velocity
  INTEGER s             ! species

  REAL(8) x_new

  INTEGER ALLOC_ERR

  INTEGER left_node     ! left node of the electron after reflection (to be injected)
  INTEGER right_node    ! right node - " -

! restore the surface charge density, because it is modified in the dynamic procedures
  Q_left = 0    

!@@@  x_new = ABS(x)
  x_new = 0.0_8

! electrons
  IF (s.EQ.1) THEN

! Save the parameters of the injected (elastically reflected) electron in the linked list
     IF (ASSOCIATED(Current_electron)) THEN
        N_inject(1)            = N_inject(1) + 1
        electron_reflux_count  = electron_reflux_count + 1
        electron_reflux_energy = electron_reflux_energy + v2
        Rate_number_leftemit(1) = Rate_number_leftemit(1) + 1           !
        Rate_energy_leftemit(1) = Rate_energy_leftemit(1) + v2          ! for diagnostics
        Current_electron%X  = x_new
        Current_electron%VX = ABS(vx)
        Current_electron%VY = vy
        Current_electron%VZ = vz
!        Current_electron%AX = 0.0_8
        Current_electron%Tag = 0 ! eTag_Emit_Left                               ! mark the electron
        CALL ADD_EMITTED_E_TO_LEFT_DF(ABS(vx), vy, vz)
        ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
        IF (ALLOC_ERR.NE.0) THEN
           PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLUX_SPECULAR:")', Rank_of_process
           PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        Current_electron => Current_electron%next
        NULLIFY(Current_electron%next)
     ELSE
        PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLUX_SPECULAR:")', Rank_of_process
        PRINT  '(2x,"Current_electron is NOT associated!")'
        PRINT  '(2x,"The program will be terminated now :(")'
        STOP
     END IF

! ions
  ELSE IF (s.EQ.2) THEN

! Save the parameters of the injected (elastically reflected) ion in the linked list
     IF (ASSOCIATED(Current_ion)) THEN
        N_inject(2)       = N_inject(2) + 1
        ion_reflux_count  = ion_reflux_count + 1
        ion_reflux_energy = ion_reflux_energy + v2
        Rate_number_leftemit(2) = Rate_number_leftemit(2) + 1                !
        Rate_energy_leftemit(2) = Rate_energy_leftemit(2) + v2               !
        Current_ion%X  = x_new
        Current_ion%VX = ABS(vx)
        Current_ion%VY = vy
        Current_ion%VZ = vz
!        Current_ion%AX = 0.0_8
        Current_ion%Tag = 0            ! later the tag can be modified here in order to mark the collided ion
        ALLOCATE(Current_ion%next, STAT = ALLOC_ERR)
        IF (ALLOC_ERR.NE.0) THEN
           PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLUX_SPECULAR:")', Rank_of_process
           PRINT  '(2x,"Error in ALLOCATE Current_ion%next !!!")'
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        Current_ion => Current_ion%next
        NULLIFY(Current_ion%next)
     ELSE
        PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLUX_SPECULAR:")', Rank_of_process
        PRINT  '(2x,"Current_ion is NOT associated!")'
        PRINT  '(2x,"The program will be terminated now :(")'
        STOP
     END IF

  END IF

! account for the density of emitted particle
  left_node  = INT(ABS(x))
  right_node = left_node + 1
  IF (left_node.LT.N_cells) THEN   
     Q_strm_spec(left_node, s)  = Q_strm_spec(left_node, s)  + right_node - x_new    ! s = 2
     Q_strm_spec(right_node, s) = Q_strm_spec(right_node, s) + x_new - left_node     ! s = 2
  ELSE   ! if particle was placed exactly on the wall
     PRINT '(/2x,"Process ",i3," : LEFT_WALL_REFLUX_SPECULAR: ERROR! emitted particle of species ",i2," was placed NOT inside the plasma!")', Rank_of_process, s
     PRINT '(2x, "The program will be terminated now!")'
     STOP
  END IF

END SUBROUTINE LEFT_WALL_REFLUX_SPECULAR

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE INJECT_ADDITIONAL_EI_PAIRS_LEFT_WALL

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : T_cntr
  USE SEEmission,           ONLY : AddInjectionFlag, N_of_lost_ions, const_add_N_to_inject_by_proc, variable_add_N_to_inject

  IMPLICIT NONE

  INTEGER var_add, N_clients, T_cntr_flashed, k, T_inject_k

! additional injection turned off
  IF (AddInjectionFlag.EQ.0) THEN
     N_of_lost_ions = 0     
     RETURN
  END IF

! skip initial burst of wall flux
  IF (T_cntr.LT.100) THEN
     N_of_lost_ions = 0
     RETURN
  END IF

! feedback between the ion flux to the right wall and the additional injected flux
  IF (AddInjectionFlag.EQ.1) THEN

     CALL INJECT_N_EI_PAIRS_LEFT_WALL(N_of_lost_ions)

! constant additional injected flux
  ELSE IF (AddInjectionFlag.EQ.2) THEN

     var_add = 0

     N_clients = N_of_processes-1

     T_cntr_flashed = MOD(T_cntr, N_clients) + 1  ! periodically grows from 1 to N_clients
     
     DO k = 0, variable_add_N_to_inject-1                ! this cycle works only if variable_add_N_to_inject>0

! find "timestep" suitable for injection for given client        
        T_inject_k = Rank_of_process + k
        IF (T_inject_k.GT.N_clients) T_inject_k = T_inject_k - N_clients
        
        IF (T_inject_k.EQ.T_cntr_flashed) THEN
           var_add=1
           EXIT
        END IF
     END DO

     CALL INJECT_N_EI_PAIRS_LEFT_WALL(const_add_N_to_inject_by_proc + var_add)

  END IF

  N_of_lost_ions = 0       ! clear the counter before the next pre-push

END SUBROUTINE INJECT_ADDITIONAL_EI_PAIRS_LEFT_WALL

!-------------------------------------------------------------------------------------------------
! this subroutine must be called after the final push before rearranging of the particle arays.
! therefore it is not necessary to modify the mesh charge density here
!
SUBROUTINE INJECT_N_EI_PAIRS_LEFT_WALL(N_to_inject)

  USE ParallelOperationValues
  USE SEEmission 
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_emit, Rate_energy_leftemit, Rate_number_leftemit

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N_to_inject

  REAL(8) x             ! coordinate of refluxed particle (after refluxing)
  REAL(8) vx, vy, vz    ! velocity components of refluxed particle
  REAL(8) v2            ! absolute value of velocity (squared)

  INTEGER ALLOC_ERR, j

  DO j = 1, N_to_inject !N_of_lost_ions

! get the new velocities for electrons -----------------
     CALL GetInjMaxwellVelocity(vx) 
     CALL GetMaxwellVelocity(vy)  
     CALL GetMaxwellVelocity(vz)     
! introduce the anisotropy
     vx = vx * SQRT(Tx_e_eV / T_e_eV) + Vx_e_drift
     vy = vy * SQRT(Tz_e_eV / T_e_eV) + Vy_e_drift
     vz = vz * SQRT(Tz_e_eV / T_e_eV)

! calculate the new coordinate 
     x = well_random_number() * vx * KVx
     
! calculate the squared absolute velocity
     v2 = vx * vx + vy * vy + vz * vz
! Save the parameters of the refluxed electron in the linked list
     IF (ASSOCIATED(Current_electron)) THEN
        N_inject(1)            = N_inject(1) + 1
        electron_reemit_count  = electron_reemit_count + 1
        electron_reemit_energy = electron_reemit_energy + v2
        Rate_number_leftemit(1) = Rate_number_leftemit(1) + 1             !
        Rate_energy_leftemit(1) = Rate_energy_leftemit(1) + v2            !
        Current_electron%X       = x
        Current_electron%VX      = vx
        Current_electron%VY      = vy
        Current_electron%VZ      = vz
!        Current_electron%AX      = 0.0_8
        Current_electron%Tag     = eTag_Emit_Left                         ! mark the reemitted electron !!!
        CALL ADD_EMITTED_E_TO_LEFT_DF(vx, vy, vz)
        ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
        IF (ALLOC_ERR.NE.0) THEN
           PRINT '(/2x,"Process ",i3," : Error in INJECT_N_EI_PAIRS_LEFT_WALL:")', Rank_of_process
           PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        Current_electron => Current_electron%next
        NULLIFY(Current_electron%next)
     ELSE
        PRINT '(/2x,"Process ",i3," : Error in INJECT_N_EI_PAIRS_LEFT_WALL:")', Rank_of_process
        PRINT  '(2x,"Current_electron is NOT associated!")'
        PRINT  '(2x,"The program will be terminated now :(")'
        STOP
     END IF

! get the new velocities for ions ----------------------
     CALL GetInjMaxwellVelocity(vx) 
     CALL GetMaxwellVelocity(vy)  
     CALL GetMaxwellVelocity(vz)     
! renormalize velocity for ions
     vx = vx * VT(2) / VT(1) + Vx_i_drift
     vy = vy * VT(2) / VT(1) + Vy_i_drift
     vz = vz * VT(2) / VT(1)

! calculate the squared absolute velocity
     v2 = vx * vx + vy * vy + vz * vz

! calculate the new coordinate  
     x = x !RAN(I_random_seed) * vx * KVx        !############### ##################### ######################

! Save the parameters of the refluxed ion in the linked list
     IF (ASSOCIATED(Current_ion)) THEN
        N_inject(2)       = N_inject(2) + 1
        ion_reemit_count  = ion_reemit_count + 1
        ion_reemit_energy = ion_reemit_energy + v2
        Rate_number_leftemit(2) = Rate_number_leftemit(2) + 1                !
        Rate_energy_leftemit(2) = Rate_energy_leftemit(2) + v2               !
        Current_ion%X       = x
        Current_ion%VX      = vx
        Current_ion%VY      = vy
        Current_ion%VZ      = vz
!        Current_ion%AX      = 0.0_8
        Current_ion%Tag     = 0                     ! later this value can be modified to mark the reemitted ion !!!!!!!!!!
        ALLOCATE(Current_ion%next, STAT = ALLOC_ERR)
        IF (ALLOC_ERR.NE.0) THEN
           PRINT '(/2x,"Process ",i3," : Error in INJECT_N_EI_PAIRS_LEFT_WALL:")', Rank_of_process
           PRINT  '(2x,"Error in ALLOCATE Current_ion%next !!!")'
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        Current_ion => Current_ion%next
        NULLIFY(Current_ion%next)
     ELSE
        PRINT '(/2x,"Process ",i3," : Error in INJECT_N_EI_PAIRS_LEFT_WALL:")', Rank_of_process
        PRINT  '(2x,"Current_ion is NOT associated!")'
        PRINT  '(2x,"The program will be terminated now :(")'
        STOP
     END IF

  END DO

END SUBROUTINE INJECT_N_EI_PAIRS_LEFT_WALL
