!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%>>>>>>>>>>> LEFT WALL >>>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE PROCESS_ION_COLLISION_WITH_LEFT_WALL(x, vx, vy, vz, v2, s)

  USE SEEmission, ONLY : N_of_lost_ions, PlasmaSourceFlag, Ion_interac_model
  USE Diagnostics, ONLY : Factor_energy_eV
  USE IonInducedSEEmission
  USE CurrentProblemValues

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) x            ! coordinate of particle
  REAL(8) vx           ! VX of colliding particle
  REAL(8) vy           ! VY of colliding particle
  REAL(8) vz           ! VZ of colliding particle
  REAL(8) v2           ! v2=vx^2+vy^2+vz^2
  INTEGER s            ! species

  REAL(8) R                  ! random number 

  REAL(8) energy_inc_eV   ! energy of incident particle

! if necessary, perform refluxing at the left wall
  IF (x.LT.0.0_8) THEN
     IF (PlasmaSourceFlag.EQ.1) THEN
        CALL LEFT_WALL_REFLUX_THERMAL(x, vx, s)
        RETURN
     ELSE IF (PlasmaSourceFlag.EQ.2) THEN
        CALL LEFT_WALL_REFLUX_SPECULAR(x, vx, vy, vz, v2, s)
        RETURN
     END IF
  END IF

  N_of_lost_ions = N_of_lost_ions + 1

  IF (Ion_interac_model.EQ.0) RETURN

  IF (Ion_interac_model.EQ.1) THEN
     CALL LEFT_WALL_REFLECT_ION(x, vx, vy, vz, v2, s) 
     RETURN
  END IF

! calculate the energy of incident ion
  energy_inc_eV = v2 * Ms(2) * Factor_energy_eV

  IF (energy_inc_eV.LT.minE_ionsee_eV) RETURN

! try to inject a secondary electron
  R = well_random_number()
  IF (R.LT.setD_ionsee_true) CALL LEFT_WALL_INJECT_ION_INDUCED_SECONDARY

END SUBROUTINE PROCESS_ION_COLLISION_WITH_LEFT_WALL

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE LEFT_WALL_REFLECT_ION(x, vx, vy, vz, v2, s)

  USE ParallelOperationValues
  USE SEEmission
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_emit, Rate_energy_leftemit, Rate_number_leftemit
  IMPLICIT NONE

  REAL(8) x             ! coordinate of incident particle, must be beyond the plasma boundaries 
  REAL(8) vx, vy, vz    ! velocity components of incident particle
  REAL(8) v2            ! squared absolute value of velocity
  INTEGER s             ! type of species

  INTEGER ALLOC_ERR

  REAL(8) x_new         ! x-coordinate of the electron after reflection (to be injected)
  INTEGER left_node     ! left node  - " -
  INTEGER right_node    ! right node - " -

! restore the surface charge density
  Q_left = Q_left - Qs(s)
! calculate the new coordinate
!  x_new = -x
  x_new = 0.0_8                            ! calculate the new coordinate
  Q_strm_spec(0, s) = Q_strm_spec(0, s) + 1.0_8 

  IF ((x_new.LT.0.0_8).OR.(x_new.GT.N_cells)) THEN
     PRINT '(/2x,"Process ",i3," : LEFT_WALL_REFLECT_ION: ERROR! emitted particle was placed NOT inside the plasma!")', Rank_of_process
     PRINT '(2x, "The program will be terminated now!")'
     STOP
  END IF

! Save the parameters of the injected (elastically reflected) ion in the linked list
  IF (ASSOCIATED(Current_ion)) THEN
     N_inject(s)         = N_inject(s) + 1

     Rate_energy_emit(s)     = Rate_energy_emit(s) + v2         !
     Rate_number_leftemit(s) = Rate_number_leftemit(s) + 1      !
     Rate_energy_leftemit(s) = Rate_energy_leftemit(s) + v2     ! for diagnostics
     ion_left_reflect_count  = ion_left_reflect_count + 1       !
     ion_left_reflect_energy = ion_left_reflect_energy + v2     !

     Current_ion%X   = x_new
     Current_ion%VX  = -vx
     Current_ion%VY  = vy
     Current_ion%VZ  = vz
!     Current_ion%AX  = 0.0_8
     Current_ion%Tag = 0            ! later the tag can be modified here in order to mark the collided ion

     ALLOCATE(Current_ion%next, STAT = ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLECT_ION:")', Rank_of_process
        PRINT  '(2x,"Error in ALLOCATE Current_ion%next !!!")'
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
     Current_ion => Current_ion%next
     NULLIFY(Current_ion%next)
  ELSE
     PRINT '(/2x,"Process ",i3," : Error in LEFT_WALL_REFLECT_ION:")', Rank_of_process
     PRINT  '(2x,"Current_ion is NOT associated!")'
     PRINT  '(2x,"The program will be terminated now :(")'
     STOP
  END IF

!! account for the density of emitted particle
!! account for an [unlikely] situation when the ion is so fast that the reflected ion is not in the very first cell
!  left_node  = INT(x_new)
!  right_node = left_node + 1
!  Q_strm_spec(left_node, s)  = Q_strm_spec(left_node, s)  + right_node - x_new    ! s = 2
!  Q_strm_spec(right_node, s) = Q_strm_spec(right_node, s) + x_new - left_node     ! s = 2

! decrease the counter of lost ions
  N_of_lost_ions = N_of_lost_ions - 1    

END SUBROUTINE LEFT_WALL_REFLECT_ION

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE LEFT_WALL_INJECT_ION_INDUCED_SECONDARY

  USE ParallelOperationValues
  USE SEEmission, ONLY : T_see_true_eV
  USE IonInducedSEEmission
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_leftemit, Rate_number_leftemit

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER ALLOC_ERR

  REAL(8) energy        ! energy of injected (true secondary emitted) electron
  REAL(8) teta, fi      ! angles of injection (emission)
  REAL(8) v             ! absolute velocity   - " -
  REAL(8) vx, vy, vz    ! velocity components - " -
  REAL(8) x_new         ! x-coordinate        - " -
  INTEGER left_node     ! left node           - " -
  INTEGER right_node    ! right node          - " -

! get the energy of a secondary electron
  CALL GetTrueSecondaryEnergy(energy)                ! note that this function is tuned for the electron-induced SEE energy
  energy = energy * (T_ionsee_eV / T_see_true_eV)

! get the absolute velocity of a secondary electron
  v = SQRT(2.0_8 * energy)

! get the angles of backscattering
  teta = ASIN(well_random_number())                ! assume that the distribution over the teta angle is f(teta) = COS(teta)   
  fi   = 6.283185307_8 * well_random_number()      

! get the velocity components
  vx = v * COS(teta)
  vy = v * SIN(teta) * SIN(fi)
  vz = v * SIN(teta) * COS(fi)

! decrease the surface charge density
  Q_left = Q_left - Qs(1)

! set the new coordinate
  x_new = 0.0_8

!  left_node  = 0
!  right_node = 1

! Save the parameters of the injected (inelastically backscattered) electron in the linked list
  IF (ASSOCIATED(Current_electron)) THEN
     N_inject(1) = N_inject(1) + 1

     Rate_number_leftemit(1) = Rate_number_leftemit(1) + 1               !
     Rate_energy_leftemit(1) = Rate_energy_leftemit(1) + 0.5_8 * energy  !
     ionsee_left_count  = ionsee_left_count  + 1                         ! for diagnostics
     ionsee_left_energy = ionsee_left_energy + 0.5_8 * energy            !
     CALL ADD_EMITTED_E_TO_LEFT_DF(vx, vy, vz)

     Current_electron%X   = x_new
     Current_electron%VX  = vx
     Current_electron%VY  = vy
     Current_electron%VZ  = vz
!     Current_electron%AX  = 0.0_8
     Current_electron%Tag = eTag_Emit_Left                               ! mark the electron

     ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in INJECT_ION_INDUCED_SECONDARY_LEFT :")', Rank_of_process
        PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
     Current_electron => Current_electron%next
     NULLIFY(Current_electron%next)
  ELSE
     PRINT '(/2x,"Process ",i3," : Error in INJECT_ION_INDUCED_SECONDARY_LEFT :")', Rank_of_process
     PRINT  '(2x,"Current_electron is NOT associated!")'
     PRINT  '(2x,"The program will be terminated now :(")'
     STOP
  END IF

! account for the density of emitted particle
  Q_strm_spec(0, 1) = Q_strm_spec(0, 1) + 1.0_8
!  Q_strm_spec(left_node, 1)  = Q_strm_spec(left_node, 1)  + right_node - x_new    ! s = 1
!  Q_strm_spec(right_node, 1) = Q_strm_spec(right_node, 1) + x_new - left_node     ! s = 1

END SUBROUTINE LEFT_WALL_INJECT_ION_INDUCED_SECONDARY

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%>>>>>>>>>>> RIGHT WALL >>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE PROCESS_ION_COLLISION_WITH_RIGHT_WALL(x, vx, vy, vz, v2, s)

  USE SEEmission, ONLY : N_of_lost_ions, PlasmaSourceFlag, Ion_interac_model
  USE Diagnostics, ONLY : Factor_energy_eV
  USE IonInducedSEEmission
  USE CurrentProblemValues

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) x            ! coordinate of particle
  REAL(8) vx           ! VX of colliding particle
  REAL(8) vy           ! VY of colliding particle
  REAL(8) vz           ! VZ of colliding particle
  REAL(8) v2           ! v2=vx^2+vy^2+vz^2
  INTEGER s            ! species

  REAL(8) R               ! random number 

  REAL(8) energy_inc_eV   ! energy of incident particle

  N_of_lost_ions = N_of_lost_ions + 1

  IF (Ion_interac_model.EQ.0) RETURN

  IF (Ion_interac_model.EQ.1) THEN
     CALL RIGHT_WALL_REFLECT_ION(x, vx, vy, vz, v2, s) 
     RETURN
  END IF

! calculate the energy of incident ion
  energy_inc_eV = v2 * Ms(2) * Factor_energy_eV

  IF (energy_inc_eV.LT.minE_ionsee_eV) RETURN

! try to inject a secondary electron
  R = well_random_number()
  IF (R.LT.setD_ionsee_true) CALL RIGHT_WALL_INJECT_ION_INDUCED_SECONDARY

END SUBROUTINE PROCESS_ION_COLLISION_WITH_RIGHT_WALL

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE RIGHT_WALL_REFLECT_ION(x, vx, vy, vz, v2, s)

  USE ParallelOperationValues
  USE SEEmission
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_emit, Rate_energy_rightemit, Rate_number_rightemit
  IMPLICIT NONE

  REAL(8) x             ! coordinate of incident particle, must be beyond the plasma boundaries 
  REAL(8) vx, vy, vz    ! velocity components of incident particle
  REAL(8) v2            ! squared absolute value of velocity
  INTEGER s             ! type of species

  INTEGER ALLOC_ERR

  REAL(8) x_new         ! x-coordinate of the electron after reflection (to be injected)
  INTEGER left_node     ! left node  - " -
  INTEGER right_node    ! right node - " -

! restore the surface charge density
  Q_right =  Q_right - Qs(s)
! calculate the new coordinate
!  x_new = 2.0_8 * N_cells - x
  x_new = dble(N_cells)                  ! calculate the new coordinate
  Q_strm_spec(N_cells, s) = Q_strm_spec(N_cells, s) + 1.0_8 

  IF ((x_new.LT.0.0_8).OR.(x_new.GT.N_cells)) THEN
     PRINT '(/2x,"Process ",i3," : RIGHT_WALL_REFLECT_ION: ERROR! emitted particle was placed NOT inside the plasma!")', Rank_of_process
     PRINT '(2x, "The program will be terminated now!")'
     STOP
  END IF

! Save the parameters of the injected (elastically reflected) ion in the linked list
  IF (ASSOCIATED(Current_ion)) THEN
     N_inject(s)         = N_inject(s) + 1

     Rate_energy_emit(s)      = Rate_energy_emit(s) + v2           !
     Rate_number_rightemit(s) = Rate_number_rightemit(s) + 1       !
     Rate_energy_rightemit(s) = Rate_energy_rightemit(s) + v2      ! for diagnostics
     ion_right_reflect_count  = ion_right_reflect_count + 1        !
     ion_right_reflect_energy = ion_right_reflect_energy + v2      !

     Current_ion%X   = x_new
     Current_ion%VX  = -vx
     Current_ion%VY  = vy
     Current_ion%VZ  = vz
!     Current_ion%AX  = 0.0_8
     Current_ion%Tag = 0            ! later the tag can be modified here in order to mark the collided ion

     ALLOCATE(Current_ion%next, STAT = ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in RIGHT_WALL_REFLECT_ION:")', Rank_of_process
        PRINT  '(2x,"Error in ALLOCATE Current_ion%next !!!")'
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
     Current_ion => Current_ion%next
     NULLIFY(Current_ion%next)
  ELSE
     PRINT '(/2x,"Process ",i3," : Error in RIGHT_WALL_REFLECT_ION:")', Rank_of_process
     PRINT  '(2x,"Current_ion is NOT associated!")'
     PRINT  '(2x,"The program will be terminated now :(")'
     STOP
  END IF

!! account for the density of emitted particle
!! account for an [unlikely] situation when the ion is so fast that the reflected ion is not in the very first cell
!  left_node  = INT(x_new)
!  right_node = left_node + 1
!! account for a [very unlikely] situation when the ion is put exactly at the right boundary x=N_cells
!  IF (left_node.EQ.N_cells) THEN
!     left_node = N_cells-1
!     right_node = N_cells
!  END IF
!  Q_strm_spec(left_node, s)  = Q_strm_spec(left_node, s)  + right_node - x_new    ! s = 2
!  Q_strm_spec(right_node, s) = Q_strm_spec(right_node, s) + x_new - left_node     ! s = 2

! decrease the counter of lost ions
  N_of_lost_ions = N_of_lost_ions - 1    

END SUBROUTINE RIGHT_WALL_REFLECT_ION

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE RIGHT_WALL_INJECT_ION_INDUCED_SECONDARY

  USE ParallelOperationValues
  USE SEEmission, ONLY : T_see_true_eV
  USE IonInducedSEEmission
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_rightemit, Rate_number_rightemit

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER ALLOC_ERR

  REAL(8) energy        ! energy of injected (true secondary emitted) electron
  REAL(8) teta, fi      ! angles of injection (emission)
  REAL(8) v             ! absolute velocity   - " -
  REAL(8) vx, vy, vz    ! velocity components - " -
  REAL(8) x_new         ! x-coordinate        - " -
  INTEGER left_node     ! left node           - " -
  INTEGER right_node    ! right node          - " -

! get the energy of a secondary electron
  CALL GetTrueSecondaryEnergy(energy)                ! note that this function is tuned for the electron-induced SEE energy
  energy = energy * (T_ionsee_eV / T_see_true_eV) 

! get the absolute velocity of a secondary electron
  v = SQRT(2.0_8 * energy)

! get the angles of backscattering
  teta = ASIN(well_random_number())                ! assume that the distribution over the teta angle is f(teta) = COS(teta)   
  fi   = 6.283185307_8 * well_random_number()      

! get the velocity components
  vx =-v * COS(teta)
  vy = v * SIN(teta) * SIN(fi)
  vz = v * SIN(teta) * COS(fi)

! decrease the surface charge density
  Q_right =  Q_right - Qs(1)

! set the new coordinate
  x_new = DBLE(N_cells)

!  left_node = N_cells - 1
!  right_node = N_cells

! Save the parameters of the injected (inelastically backscattered) electron in the linked list
  IF (ASSOCIATED(Current_electron)) THEN
     N_inject(1) = N_inject(1) + 1

     Rate_number_rightemit(1) = Rate_number_rightemit(1) + 1                ! 
     Rate_energy_rightemit(1) = Rate_energy_rightemit(1) + 0.5_8 * energy   !
     ionsee_right_count  = ionsee_right_count  + 1                          ! for diagnostics
     ionsee_right_energy = ionsee_right_energy + 0.5_8 * energy             !
     CALL ADD_EMITTED_E_TO_RIGHT_DF(vx, vy, vz)

     Current_electron%X   = x_new
     Current_electron%VX  = vx
     Current_electron%VY  = vy
     Current_electron%VZ  = vz
!     Current_electron%AX  = 0.0_8
     Current_electron%Tag = eTag_Emit_Right                              ! mark the electron

     ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in INJECT_ION_INDUCED_SECONDARY_RIGHT :")', Rank_of_process
        PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
     Current_electron => Current_electron%next
     NULLIFY(Current_electron%next)
  ELSE
     PRINT '(/2x,"Process ",i3," : Error in INJECT_ION_INDUCED_SECONDARY_RIGHT :")', Rank_of_process
     PRINT  '(2x,"Current_electron is NOT associated!")'
     PRINT  '(2x,"The program will be terminated now :(")'
     STOP
  END IF

! account for the density of emitted particle
  Q_strm_spec(N_cells, 1) = Q_strm_spec(N_cells, 1) + 1.0_8
!  Q_strm_spec(left_node, 1)  = Q_strm_spec(left_node, 1)  + right_node - x_new    ! s = 1
!  Q_strm_spec(right_node, 1) = Q_strm_spec(right_node, 1) + x_new - left_node     ! s = 1

END SUBROUTINE RIGHT_WALL_INJECT_ION_INDUCED_SECONDARY
