!==================================================================================
! This subroutine performs the final push
!
SUBROUTINE FINAL_PUSH 

  USE CurrentProblemValues
  USE Diagnostics
  USE SEEmission
  USE ElectronInjection
  USE BeamInPlasma

  IMPLICIT NONE
  
  INTEGER s            ! species type
  INTEGER k            ! particle
! temporary values  
  INTEGER n_str_left     ! node to the left of the particle at the streaming position
  INTEGER n_left         ! node to the left of the particle at the final position
  REAL(8) Xtmp           ! X-coordinate
  REAL(8) Xstr           ! streaming coordinate before correction is applied
  REAL(8) EXstr          ! streaming X-electric field

  REAL(8) v_2            ! squared absolute value of velocity

  Q_strm_spec = 0.0_8  ! Clean the densities because we will need non-streaming densities
                       ! in the near-wall nodes

! Cycle over electrons
  s = 1
  k = 0

  DO WHILE (k.LT.N_part(s))

     k = k + 1

     IF (species(s)%part(k)%X.GT.N_cells) THEN     !
        EXstr = EX(N_cells)                        ! for particles that crossed the plasma borders during the predicting push
     ELSE IF (species(s)%part(k)%X.LT.0.0_8) THEN  !
        EXstr = EX(0)                              !
     ELSE
        n_str_left   = INT(species(s)%part(k)%X)           ! obtained after pre-push, corresponds to streaming position 
        IF (n_str_left.EQ.N_cells) n_str_left = N_cells - 1
        EXstr = EX(n_str_left) + GradEX(n_str_left) * (species(s)%part(k)%X - n_str_left)    
     END IF

! velocity correction
     species(s)%part(k)%VX = species(s)%part(k)%VX + KvEx(s) * EXstr
     species(s)%part(k)%VY = species(s)%part(k)%VY + KvEy(s) * EXstr
     species(s)%part(k)%VZ = species(s)%part(k)%VZ + KvEz(s) * EXstr

     IF (UseSmartTagsFlag.EQ.1) THEN
! if tracing of whether a particle changed its direction of propagation along the x-axis is requested
! clear the tag of the particle if the x-velocity changes sign after x-acceleration (i.e. particle turns back or stops)
        IF (species(s)%part(k)%prev_VX.GE.0.0_8) THEN                                          !
           IF (species(s)%part(k)%VX.LE.0.0_8) species(s)%part(k)%Tag = 0     !
        ELSE                                                                  ! diagnostics
           IF (species(s)%part(k)%VX.GE.0.0_8) species(s)%part(k)%Tag = 0     !
        END IF                                                                !
     END IF

! x-move
     Xstr = species(s)%part(k)%X
     Xtmp = species(s)%part(k)%X + KxEx(s) * EXstr
     species(s)%part(k)%X = Xtmp

! z-move (if requested)
     IF (trace_z_shift) THEN
        species(s)%part(k)%Z = species(s)%part(k)%Z + KVx * species(s)%part(k)%VZ
     END IF

     IF (Xtmp.GT.N_cells) THEN
! if crossed the right boundary
        v_2 = species(s)%part(k)%VX**2 + species(s)%part(k)%VY**2 + species(s)%part(k)%VZ**2       !
        Rate_energy_rightwall(s) = Rate_energy_rightwall(s) + v_2                                  ! diagnostics
        Rate_number_rightwall(s) = Rate_number_rightwall(s) + 1                                    !
        IF (PeriodicBoundaryFlag.EQ.1) THEN
           species(s)%part(k)%X = Xtmp - N_cells   ! move particle to the left end of the system
           Rate_energy_leftemit(s) = Rate_energy_leftemit(s) + v_2                                  ! diagnostics
           Rate_number_leftemit(s) = Rate_number_leftemit(s) + 1                                    !
        ELSE
           IF (Xstr.LE.N_cells) THEN
! since we are here, collision occurred NOT at the pre-push, but after the final correction was applied
! so we must update the wall charge 
              Q_right = Q_right + Qs(s)
           END IF
           SELECT CASE (species(s)%part(k)%Tag)                                              !
              CASE (eTag_Emit_Left)                                                          !
                 prie_right_from_left_count   = prie_right_from_left_count   + 1             !
                 prie_right_from_left_energy  = prie_right_from_left_energy  + v_2           ! diagnostics
              CASE (eTag_Coll_Neutral)                                                       !
                 prie_right_after_coll_count  = prie_right_after_coll_count  + 1             !
                 prie_right_after_coll_energy = prie_right_after_coll_energy + v_2           !
           END SELECT                                                                        !
           CALL ADD_PRIMARY_E_TO_RIGHT_DF(species(s)%part(k)%VX, species(s)%part(k)%VY, species(s)%part(k)%VZ)      ! diagnostics
           CALL PROCESS_ELECTRON_COLLISION_WITH_WALL(s, Xtmp, k)    
           CALL SUBSTITUTE_LEAVING_PARTICLE(s, k)
           k = k - 1
           N_part(s) = N_part(s) - 1
           CYCLE
        END IF
     ELSE IF (Xtmp.LT.0.0_8) THEN
! if crossed the left boundary
        v_2 = species(s)%part(k)%VX**2 + species(s)%part(k)%VY**2 + species(s)%part(k)%VZ**2       !
        Rate_energy_leftwall(s)  = Rate_energy_leftwall(s)  + v_2                                  ! diagnostics
        Rate_number_leftwall(s)  = Rate_number_leftwall(s)  + 1                                    !
        IF (PeriodicBoundaryFlag.EQ.1) THEN
           species(s)%part(k)%X = N_cells + Xtmp    ! move particle to the right end of the system
           Rate_energy_rightemit(s) = Rate_energy_rightemit(s) + v_2                                  ! diagnostics
           Rate_number_rightemit(s) = Rate_number_rightemit(s) + 1                                    !
        ELSE
           IF (Xstr.GE.0.0_8) then
! since we are here, collision occurred NOT at the pre-push, but after the final correction was applied
! so we must update the wall charge
              Q_left = Q_left + Qs(s)
           END IF
           SELECT CASE (species(s)%part(k)%Tag)                                        !
           CASE (eTag_Emit_Right)                                                   !
              prie_left_from_right_count  = prie_left_from_right_count  + 1         !
              prie_left_from_right_energy = prie_left_from_right_energy + v_2       ! diagnostics
           CASE (eTag_Coll_Neutral)                                                 !
              prie_left_after_coll_count  = prie_left_after_coll_count  + 1         !
              prie_left_after_coll_energy = prie_left_after_coll_energy + v_2       !
           END SELECT                                                                  !
           CALL ADD_PRIMARY_E_TO_LEFT_DF(species(s)%part(k)%VX, species(s)%part(k)%VY, species(s)%part(k)%VZ)    ! diagnostics
           CALL PROCESS_ELECTRON_COLLISION_WITH_WALL(s, Xtmp, k)   
           CALL SUBSTITUTE_LEAVING_PARTICLE(s, k)
           k = k - 1
           N_part(s) = N_part(s) - 1
           CYCLE
        END IF
     ELSE IF (Xstr.GT.N_cells) THEN
! crossed the right boundary at the pre-push but after correction returned to the simulation area
! so we must remove contribution to the wall charge
        Q_right = Q_right - Qs(s)
     ELSE IF (Xstr.LT.0.0_8) THEN
! crossed the left boundary at the pre-push but after correction returned to the simulation area
! so we must remove contribution to the wall charge
        Q_left = Q_left - Qs(s)
     END IF

! calculate the non-streaming density
     n_left = INT(Xtmp) 
     IF (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
     Q_strm_spec(n_left,   s)  = Q_strm_spec(n_left, s) + (n_left + 1 - Xtmp)
     Q_strm_spec(n_left+1, s)  = Q_strm_spec(n_left+1, s) + (Xtmp - n_left)

! note that for periodic boundary conditions, accelerations of electrons crossing the boundaries are updated here as well        
     species(s)%part(k)%AX = 0.5_8 * (species(s)%part(k)%AX + QMs(s) * EXstr)
        
  END DO ! End of cycle over electrons

  IF (N_spec.EQ.2) THEN

! Cycle over ions
     s = 2
     k = 0
     DO WHILE (k.LT.N_part(s))

        k = k + 1

        IF (species(s)%part(k)%X.GT.N_cells) THEN     !
           EXstr = EX(N_cells)                        ! for particles that crossed the plasma borders during the predicting push
        ELSE IF (species(s)%part(k)%X.LT.0.0_8) THEN  !
           EXstr = EX(0)                              !
        ELSE
           n_str_left   = INT(species(s)%part(k)%X)           ! obtained after pre-push, corresponds to streaming position 
           IF (n_str_left.EQ.N_cells) n_str_left = N_cells - 1
           EXstr        = EX(n_str_left) + GradEX(n_str_left) * (species(s)%part(k)%X - n_str_left)     
        END IF

! velocity correction
        species(s)%part(k)%VX = species(s)%part(k)%VX + KvEx(s) * EXstr
        species(s)%part(k)%VY = species(s)%part(k)%VY + KvEy(s) * EXstr
        species(s)%part(k)%VZ = species(s)%part(k)%VZ + KvEz(s) * EXstr

! x-move
        Xstr = species(s)%part(k)%X
        Xtmp = species(s)%part(k)%X + KxEx(s) * EXstr
        species(s)%part(k)%X = Xtmp

! z-move (if requested) 
        IF (trace_z_shift_ion) THEN
           species(s)%part(k)%Z = species(s)%part(k)%Z + KVx * species(s)%part(k)%VZ
        END IF

        IF (Xtmp.GT.N_cells) THEN               
! if crossed the right boundary
           v_2 = species(s)%part(k)%VX**2 + species(s)%part(k)%VY**2 + species(s)%part(k)%VZ**2       !
           Rate_energy_rightwall(s) = Rate_energy_rightwall(s) + v_2                                  ! diagnostics
           Rate_number_rightwall(s) = Rate_number_rightwall(s) + 1                                    !
           IF (PeriodicBoundaryFlag.EQ.1) THEN
              species(s)%part(k)%X = Xtmp - N_cells   ! move particle to the left end of the system
              Rate_energy_leftemit(s) = Rate_energy_leftemit(s) + v_2                                  ! diagnostics 
              Rate_number_leftemit(s) = Rate_number_leftemit(s) + 1                                    !
           ELSE
              IF (Xstr.LE.N_cells) THEN 
! since we are here, collision occurred NOT at the pre-push, but after the final correction was applied
! so we must update the wall charge
                 Q_right = Q_right + Qs(s)
              END IF
              CALL PROCESS_ION_COLLISION_WITH_RIGHT_WALL(Xtmp, species(s)%part(k)%VX, species(s)%part(k)%VY, species(s)%part(k)%VZ, v_2, s)    
              CALL SUBSTITUTE_LEAVING_PARTICLE(s, k)
              k = k - 1
              N_part(s) = N_part(s) - 1
              CYCLE
           END IF
        ELSE IF (Xtmp.LT.0.0_8) THEN                 
! if crossed the left boundary
           v_2 = species(s)%part(k)%VX**2 + species(s)%part(k)%VY**2 + species(s)%part(k)%VZ**2     !
           Rate_energy_leftwall(s) = Rate_energy_leftwall(s) + v_2                                  ! diagnostics
           Rate_number_leftwall(s) = Rate_number_leftwall(s) + 1                                    !
           IF (PeriodicBoundaryFlag.EQ.1) THEN
              species(s)%part(k)%X = N_cells + Xtmp    ! move particle to the right end of the system
              Rate_energy_rightemit(s) = Rate_energy_rightemit(s) + v_2                                ! diagnostics
              Rate_number_rightemit(s) = Rate_number_rightemit(s) + 1                                  !
           ELSE
              IF (Xstr.GE.0.0_8) THEN
! since we are here, collision occurred NOT at the pre-push, but after the final correction was applied
! so we must update the wall charge
                 Q_left = Q_left + Qs(s)  ! only if it was not updated at the pre-push
              END IF
              CALL PROCESS_ION_COLLISION_WITH_LEFT_WALL(Xtmp, species(s)%part(k)%VX, species(s)%part(k)%VY, species(s)%part(k)%VZ, v_2, s)   
              CALL SUBSTITUTE_LEAVING_PARTICLE(s, k)
              k = k - 1
              N_part(s) = N_part(s) - 1
              CYCLE
           END IF
        ELSE IF (Xstr.GT.N_cells) THEN
! crossed the right boundary at the pre-push but after correction returned to the simulation area
! so we must remove contribution to the wall charge
           Q_right = Q_right - Qs(s)
        ELSE IF (Xstr.LT.0.0_8) THEN
! crossed the left boundary at the pre-push but after correction returned to the simulation area
! so we must remove contribution to the wall charge
           Q_left = Q_left - Qs(s)
        END IF

! calculate the non-streaming density
        n_left = INT(Xtmp) 
        IF (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
        Q_strm_spec(n_left,   s)  = Q_strm_spec(n_left, s) + (n_left + 1 - Xtmp)
        Q_strm_spec(n_left+1, s)  = Q_strm_spec(n_left+1, s) + (Xtmp - n_left)

! note that for periodic boundary conditions, accelerations of ions crossing the boundaries are updated here as well        
        species(s)%part(k)%AX = 0.5_8 * (species(s)%part(k)%AX + QMs(s) * EXstr)      

     END DO ! End of cycle over ions

  END IF

! this should not be necessary, just a paranoidal precaution...
  IF (PeriodicBoundaryFlag.EQ.1) THEN    
! for periodic boundaries prevent charge accumulation
     Q_right = 0
     Q_left = 0
  END IF
  
END SUBROUTINE FINAL_PUSH 


!===============================================================================
! This subroutine performs the predicting push 
!
SUBROUTINE PREDICTING_PUSH

  USE CurrentProblemValues
  USE Diagnostics
  USE SEEmission
  USE ElectronInjection
  USE BeamInPlasma

  IMPLICIT NONE
  
  INTEGER s            ! species type
  INTEGER k            ! particle

  REAL(8) new_Vx, new_Vy

  INTEGER n_left, n_right            ! bounding nodes

  REAL(8) Xstr

! if necessary set flags for accumulation of distribution functions at the walls
  CALL REQUEST_WALL_DF_FOR_SNAPSHOT

  N_inject = 0         ! Clean the counters
  Q_strm_spec = 0.0_8  ! Clean the densities

! clean wall charge density accumulators if we simulate a situation with given wall potential or external circuit
! because in this case the final wall charge density at the end of previous time step is calculated from the Gauss law
  IF (BC_flag.NE.1) THEN
     Q_left = 0
     Q_right = 0
  END IF

  IF (PeriodicBoundaryFlag.EQ.1) THEN
! periodic boundary ---------------------------------------

! electrons
     s = 1
! cycle over particles
     DO k = 1, N_part(s)
        new_Vx = K11(s) * species(s)%part(k)%VX + &
               & K12(s) * species(s)%part(k)%VY - &
               & K13(s) * species(s)%part(k)%VZ + &
               & A11(s) * species(s)%part(k)%AX - &
               & A13(s)
! y-acceleration
        new_Vy = K12(s) * species(s)%part(k)%VX + &
               & K22(s) * species(s)%part(k)%VY + &
               & K23(s) * species(s)%part(k)%VZ + &
               & A21(s) * species(s)%part(k)%AX + &
               & A23(s)
! z-acceleration
        species(s)%part(k)%VZ = K13(s) * species(s)%part(k)%VX - &
                              & K23(s) * species(s)%part(k)%VY + &
                              & K33(s) * species(s)%part(k)%VZ + &
                              & A31(s) * species(s)%part(k)%AX + &
                              & A33(s)
        species(s)%part(k)%prev_VX = species(s)%part(k)%VX
        species(s)%part(k)%VX = new_Vx
        species(s)%part(k)%VY = new_Vy
! x-move
        Xstr = species(s)%part(k)%X + KVx * species(s)%part(k)%VX

        IF (Xstr.GT.N_cells) THEN
! crossed the right boundary
           Xstr = Xstr - N_cells   ! move particle to the left end of the system
        ELSE IF (Xstr.LT.0.0_8) THEN
! crossed the left boundary
           Xstr = Xstr + N_cells   ! move particle to the right end of the system
        END IF
        species(s)%part(k)%X = Xstr  

        n_left = INT(Xstr) 
        IF (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
        n_right = n_left + 1
     
        Q_strm_spec(n_left, s)  = Q_strm_spec(n_left, s)  + n_right - Xstr
        Q_strm_spec(n_right, s) = Q_strm_spec(n_right, s) + Xstr - n_left  
     END DO
  
     IF (N_spec.EQ.2) THEN
! ions
        s = 2

        IF (ions_magnetized) THEN

! cycle over particles
           DO k = 1, N_part(s)
              new_Vx = K11(s) * species(s)%part(k)%VX + &
                     & K12(s) * species(s)%part(k)%VY - &
                     & K13(s) * species(s)%part(k)%VZ + &
                     & A11(s) * species(s)%part(k)%AX - &
                     & A13(s)
! y-acceleration
              new_Vy = K12(s) * species(s)%part(k)%VX + &
                     & K22(s) * species(s)%part(k)%VY + &
                     & K23(s) * species(s)%part(k)%VZ + &
                     & A21(s) * species(s)%part(k)%AX + &
                     & A23(s)
! z-acceleration
              species(s)%part(k)%VZ = K13(s) * species(s)%part(k)%VX - &
                                    & K23(s) * species(s)%part(k)%VY + &
                                    & K33(s) * species(s)%part(k)%VZ + &
                                    & A31(s) * species(s)%part(k)%AX + &
                                    & A33(s)
              species(s)%part(k)%prev_VX = species(s)%part(k)%VX
              species(s)%part(k)%VX = new_Vx
              species(s)%part(k)%VY = new_Vy
! x-move
              Xstr = species(s)%part(k)%X + KVx * species(s)%part(k)%VX

              IF (Xstr.GT.N_cells) THEN
! crossed the right boundary
                 Xstr = Xstr - N_cells   ! move particle to the left end of the system
              ELSE IF (Xstr.LT.0.0_8) THEN
! crossed the left boundary
                 Xstr = Xstr + N_cells   ! move particle to the right end of the system
              END IF
              species(s)%part(k)%X = Xstr  

              n_left = INT(Xstr) 
              IF (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
              n_right = n_left + 1
     
              Q_strm_spec(n_left, s)  = Q_strm_spec(n_left, s)  + n_right - Xstr
              Q_strm_spec(n_right, s) = Q_strm_spec(n_right, s) + Xstr - n_left  
           END DO

        ELSE
! ions are not magnetized

! cycle over particles
           DO k = 1, N_part(s)
! x-acceleration
              species(s)%part(k)%prev_VX = species(s)%part(k)%VX
              species(s)%part(k)%VX      = species(s)%part(k)%VX + A11(s) * species(s)%part(k)%AX     
! no y-acceleration
! z-acceleration if requested (omit_ion_z_acceleration=.false.)
              species(s)%part(k)%VZ = species(s)%part(k)%VZ + A33(s)
! x-move
              Xstr = species(s)%part(k)%X + KVx * species(s)%part(k)%VX

              IF (Xstr.GT.N_cells) THEN
! crossed the right boundary
                 Xstr = Xstr - N_cells   ! move particle to the left end of the system
              ELSE IF (Xstr.LT.0.0_8) THEN
! crossed the left boundary
                 Xstr = Xstr + N_cells   ! move particle to the right end of the system
              END IF
              species(s)%part(k)%X = Xstr  

              n_left = INT(Xstr) 
              IF (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
              n_right = n_left + 1
     
              Q_strm_spec(n_left, s)  = Q_strm_spec(n_left, s)  + n_right - Xstr
              Q_strm_spec(n_right, s) = Q_strm_spec(n_right, s) + Xstr - n_left  
           END DO

        END IF       ! IF (ions_magnetized) THEN

     END IF          ! IF (N_spec.EQ.2) THEN

  ELSE
! non-periodic boundary ------------------------------------------

! electrons
     s = 1
! cycle over particles
     DO k = 1, N_part(s)
! x-acceleration
        new_Vx = K11(s) * species(s)%part(k)%VX + &
               & K12(s) * species(s)%part(k)%VY - &
               & K13(s) * species(s)%part(k)%VZ + &
               & A11(s) * species(s)%part(k)%AX - &
               & A13(s)
! y-acceleration
        new_Vy = K12(s) * species(s)%part(k)%VX + &
               & K22(s) * species(s)%part(k)%VY + &
               & K23(s) * species(s)%part(k)%VZ + &
               & A21(s) * species(s)%part(k)%AX + &
               & A23(s)
! z-acceleration
        species(s)%part(k)%VZ = K13(s) * species(s)%part(k)%VX - &
                              & K23(s) * species(s)%part(k)%VY + &
                              & K33(s) * species(s)%part(k)%VZ + &
                              & A31(s) * species(s)%part(k)%AX + &
                              & A33(s)
        species(s)%part(k)%prev_VX = species(s)%part(k)%VX
        species(s)%part(k)%VX = new_Vx
        species(s)%part(k)%VY = new_Vy
! x-move
        Xstr = species(s)%part(k)%X + KVx * species(s)%part(k)%VX
        species(s)%part(k)%X = Xstr  

        IF (Xstr.GT.N_cells) THEN
! crossed the right boundary
           Q_right = Q_right + Qs(s)
        ELSE IF (Xstr.LT.0.0_8) THEN
! crossed the left boundary
           Q_left = Q_left + Qs(s)
        ELSE
! stays inside the system
           n_left = INT(Xstr) 
           IF (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
           n_right = n_left + 1
     
           Q_strm_spec(n_left, s)  = Q_strm_spec(n_left, s)  + n_right - Xstr
           Q_strm_spec(n_right, s) = Q_strm_spec(n_right, s) + Xstr - n_left  
        END IF
     END DO
     
     IF (N_spec.EQ.2) THEN
! ions
        s = 2

        IF (ions_magnetized) THEN

! cycle over particles
           DO k = 1, N_part(s)
! x-acceleration
              new_Vx = K11(s) * species(s)%part(k)%VX + &
                     & K12(s) * species(s)%part(k)%VY - &
                     & K13(s) * species(s)%part(k)%VZ + &
                     & A11(s) * species(s)%part(k)%AX - &
                     & A13(s)
! y-acceleration
              new_Vy = K12(s) * species(s)%part(k)%VX + &
                     & K22(s) * species(s)%part(k)%VY + &
                     & K23(s) * species(s)%part(k)%VZ + &
                     & A21(s) * species(s)%part(k)%AX + &
                     & A23(s)
! z-acceleration
              species(s)%part(k)%VZ = K13(s) * species(s)%part(k)%VX - &
                                    & K23(s) * species(s)%part(k)%VY + &
                                    & K33(s) * species(s)%part(k)%VZ + &
                                    & A31(s) * species(s)%part(k)%AX + &
                                    & A33(s)
              species(s)%part(k)%prev_VX = species(s)%part(k)%VX
              species(s)%part(k)%VX = new_Vx
              species(s)%part(k)%VY = new_Vy
! x-move
              Xstr = species(s)%part(k)%X + KVx * species(s)%part(k)%VX
              species(s)%part(k)%X = Xstr  

              IF (Xstr.GT.N_cells) THEN
! crossed the right boundary
                 Q_right = Q_right + Qs(s)
              ELSE IF (Xstr.LT.0.0_8) THEN
! crossed the left boundary
                 Q_left = Q_left + Qs(s)
              ELSE
! stays inside the system
                 n_left = INT(Xstr) 
                 IF (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
                 n_right = n_left + 1
     
                 Q_strm_spec(n_left, s)  = Q_strm_spec(n_left, s)  + n_right - Xstr
                 Q_strm_spec(n_right, s) = Q_strm_spec(n_right, s) + Xstr - n_left  
              END IF
           END DO

        ELSE
! ions are not magnetized

! cycle over particles
           DO k = 1, N_part(s)
! x-acceleration
              species(s)%part(k)%prev_VX = species(s)%part(k)%VX
              species(s)%part(k)%VX      = species(s)%part(k)%VX + A11(s) * species(s)%part(k)%AX     
! no y-acceleration
! z-acceleration if requested (omit_ion_z_acceleration=.false.)
              species(s)%part(k)%VZ = species(s)%part(k)%VZ + A33(s)
! x-move
              Xstr = species(s)%part(k)%X + KVx * species(s)%part(k)%VX
              species(s)%part(k)%X = Xstr  

              IF (Xstr.GT.N_cells) THEN
! crossed the right boundary
                 Q_right = Q_right + Qs(s)
              ELSE IF (Xstr.LT.0.0_8) THEN
! crossed the left boundary
                 Q_left = Q_left + Qs(s)
              ELSE
! stays inside the system
                 n_left = INT(Xstr) 
                 IF (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
                 n_right = n_left + 1
     
                 Q_strm_spec(n_left, s)  = Q_strm_spec(n_left, s)  + n_right - Xstr
                 Q_strm_spec(n_right, s) = Q_strm_spec(n_right, s) + Xstr - n_left  
              END IF
           END DO
           
        END IF      ! IF (ions_magnetized) THEN

     END IF         ! IF (N_spec.EQ.2) THEN

  END IF            ! IF (PeriodicBoundaryFlag.EQ.1) THEN

END SUBROUTINE PREDICTING_PUSH

!---------------------------------------------------------------------------------
! Subroutine removes particle, which must leave the plasma, 
!
SUBROUTINE SUBSTITUTE_LEAVING_PARTICLE(s, k)

  USE CurrentProblemValues, ONLY : N_part, species 
  IMPLICIT NONE
  
  INTEGER s            ! species type
  INTEGER k            ! index of particle to be removed
     
  IF (k.LT.N_part(s)) THEN                     ! if the data are not the last active one
     species(s)%part(k)%X       = species(s)%part(N_part(s))%X        ! 
     species(s)%part(k)%Z       = species(s)%part(N_part(s))%Z        ! 
     species(s)%part(k)%VX      = species(s)%part(N_part(s))%VX       ! SET the data 
     species(s)%part(k)%VY      = species(s)%part(N_part(s))%VY       ! equal to the
     species(s)%part(k)%VZ      = species(s)%part(N_part(s))%VZ       ! data from LAST ACTIVE position
     species(s)%part(k)%AX      = species(s)%part(N_part(s))%AX       !
     species(s)%part(k)%Tag     = species(s)%part(N_part(s))%Tag      !
     species(s)%part(k)%prev_VX = species(s)%part(N_part(s))%prev_VX  !
  END IF                                                ! 

  species(s)%part(N_part(s))%X       = 0.0_8 !  
  species(s)%part(N_part(s))%Z       = 0.0_8 !  
  species(s)%part(N_part(s))%VX      = 0.0_8 ! CLEAR 
  species(s)%part(N_part(s))%VY      = 0.0_8 ! the LAST ACTIVE data
  species(s)%part(N_part(s))%VZ      = 0.0_8 ! 
  species(s)%part(N_part(s))%AX      = 0.0_8 ! 
  species(s)%part(N_part(s))%Tag     = 0     ! 
  species(s)%part(N_part(s))%prev_VX = 0.0_8 !

END SUBROUTINE SUBSTITUTE_LEAVING_PARTICLE

!------------------------------------------
!
SUBROUTINE MANAGE_WALL_CHARGE_DENSITIES
    
  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(8) rbufer(1:2), rbufer2(1:2)
  INTEGER ierr

  IF (BC_flag.EQ.1) RETURN

  rbufer  = 0.0_8
  rbufer2 = 0.0_8

  IF (Rank_of_process.GT.0) THEN
! send volume charge densities at the walls to server 

     rbufer(1) = Qs(1) * Q_strm_spec(0, 1)       + Qs(2) * Q_strm_spec(0, 2)
     rbufer(2) = Qs(1) * Q_strm_spec(N_cells, 1) + Qs(2) * Q_strm_spec(N_cells, 2)

     CALL MPI_REDUCE(rbufer, rbufer2, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     
  ELSE
! receive from clients volume charge densities assigned to wall nodes

     CALL MPI_REDUCE(rbufer2, rbufer, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     full_Q_left  = (F(0) - F(1)) / K_Q - rbufer(1)
     full_Q_right = (F(N_cells) - F(N_cells-1)) / K_Q - rbufer(2)

  END IF
  
END SUBROUTINE MANAGE_WALL_CHARGE_DENSITIES

!--------------------------------------------
!
SUBROUTINE PROCESS_Z_SHIFT

  USE CurrentProblemValues
  IMPLICIT NONE

  REAL(8) factor_x, factor_z
  INTEGER s, k
  REAL(8) vx, vy, vz

  REAL(8) factor_ion

!integer count

  IF (trace_z_shift) THEN
     factor_x = SQRT(Tx_e_eV / T_e_eV)
     factor_z = SQRT(Tz_e_eV / T_e_eV)
     s = 1   ! ### for electrons only
     DO k = 1, N_part(s)
        IF (species(s)%part(k)%Z.LT.(-max_Z_shift)) THEN
           CALL GetMaxwellVelocity(vx)  
           CALL GetMaxwellVelocity(vy)  
           CALL GetMaxwellVelocity(vz)  
           species(s)%part(k)%Z = 0.0_8
           species(s)%part(k)%VX = vx * factor_x + Vx_e_drift
           species(s)%part(k)%VY = vy * factor_z + Vy_e_drift
           species(s)%part(k)%VZ = vz * factor_z
           species(s)%part(k)%AX = 0.0_8     
        END IF
     END DO
  END IF

!count=0

  IF (trace_z_shift_ion) THEN
     factor_ion = VT(2) / VT(1)
     s = 2   ! ### for ions only
     DO k = 1, N_part(s)
        IF (species(s)%part(k)%Z.GT.max_Z_shift_ion) THEN
!count=count+1
           CALL GetMaxwellVelocity(vx)  
           CALL GetMaxwellVelocity(vy)  
           CALL GetMaxwellVelocity(vz)  
           species(s)%part(k)%Z = 0.0_8
           species(s)%part(k)%VX = vx * factor_ion + Vx_i_drift
           species(s)%part(k)%VY = vy * factor_ion + Vy_i_drift
           species(s)%part(k)%VZ = vz * factor_ion + VZi_0
           species(s)%part(k)%AX = 0.0_8      
        END IF
     END DO
  END IF

!if (count.gt.0) print '("thermalized ions via z-shift :: ",i8)', count

END SUBROUTINE PROCESS_Z_SHIFT
