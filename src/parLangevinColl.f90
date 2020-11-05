!--------------------------------------------------------------------------
!
SUBROUTINE INITIATE_COULOMB_COLLISIONS

  USE ParallelOperationValues 
  USE LangevinCollisions
  USE CurrentProblemValues

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL exists
  CHARACTER (77) buf
  INTEGER ierr

! read / write the data file 
  INQUIRE (FILE = 'ssc_langevin.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  OPEN (9, FILE = 'ssc_langevin.dat')
  IF (exists) THEN
     
     IF (Rank_of_process.EQ.0) THEN
        PRINT '(2x,"Process ",i3," : Coulomb collisions control data file is found. Reading the data file...")', Rank_of_process
     END IF

     READ (9, '(A77)') buf ! '("************************ COULOMB COLLISIONS CONTROL *************************")')
     READ (9, '(A77)') buf ! '("-------d--- Langevin model {off / e-e / e-e and e-i} = {0/1/2} --------------")')
     READ (9, '(7x,i1)') Accounted_Langevin
     READ (9, '(A77)') buf ! '("----dddd--- Number of timesteps to skip between e-e collisions (>=0) --------")')
     READ (9, '(4x,i4)') T_skip_Lang
     READ (9, '(A77)') buf ! '("----dddd--- Number of additional cells to calculate Langevin coeffs (>=0) ---")')
     READ (9, '(4x,i4)') added_cells
     READ (9, '(A77)') buf ! '("------dd--- Maximal velocity for EVDF (in units of scale thermal velocity) --")')
     READ (9, '(4x,i4)') W_max_Lang

  ELSE

     Accounted_Langevin  = 0
     T_skip_Lang         = 0
     added_cells         = 0
     W_max_Lang          = 0

     IF (Rank_of_process.EQ.0) THEN

        PRINT '(/2x,"Process ",i3," : File with the name ssc_langevin.dat not found. Use the default settings ...")', &
                                                                                                   & Rank_of_process
        PRINT '(2x,"Process ",i3," : Create ssc_langevin.dat file ...")', Rank_of_process

        WRITE (9, '("************************ COULOMB COLLISIONS CONTROL *************************")')
        WRITE (9, '("-------d--- Langevin model {off / e-e / e-e and e-i} = {0/1/2} --------------")')
        WRITE (9, '(7x,i1)') Accounted_Langevin
        WRITE (9, '("----dddd--- Number of timesteps to skip between e-e collisions (>=0) --------")')
        WRITE (9, '(4x,i4)') T_skip_Lang
        WRITE (9, '("----dddd--- Number of additional cells to calculate Langevin coeffs (>=0) ---")')
        WRITE (9, '(4x,i4)') added_cells
        WRITE (9, '("------dd--- Maximal velocity for EVDF (in units of scale thermal velocity) --")')
        WRITE (9, '(6x,i2)') W_max_Lang

     END IF

  END IF

  CLOSE (9, STATUS = 'KEEP')

! report about the general status of e-e collisions
  IF (Accounted_Langevin.EQ.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '(/2x,"Electron-electron collisions are turned OFF ...")'
     RETURN
  END IF

! if collisions are accounted for ...

  IF (Rank_of_process.EQ.0) THEN 
     IF (Accounted_Langevin.EQ.1) PRINT '(/2x,"Only electron-electron collisions are turned ON ...")'
     IF (Accounted_Langevin.EQ.2) PRINT '(/2x,"Both electron-electron and electron-ion collisions are turned ON ...")'
  END IF

  T_skip_Lang = T_skip_Lang + 1       ! if initially T_skip_Lang = 0 then no timesteps will be skipped

  IF (Rank_of_process.EQ.0) THEN 
     PRINT '( 2x,"For Langevin model the EVDF in each cell will be calculated including up to ",i6," neighbor cells on each side")', added_cells
     PRINT '( 2x,"The Langevin model will be applied each ",i5," timestep")', T_skip_Lang
  END IF

  N_box_w_Lang = W_max_Lang * N_box_vel
     
  CALL ALLOCATE_LANGEVIN_ARRAYS

  Factor_F_drag = 2.0_8 * (DBLE(T_skip_Lang) * (DBLE(N_box_vel))**3) / &
                & (DBLE(N_max_vel) * DBLE(N_of_cells_debye) * DBLE(N_of_particles_cell) * N_plasma_m3 * L_debye_m**3)

  Factor_D_1_sqrt = SQRT(0.333333333333_8 * Factor_F_drag) / (DBLE(N_box_vel))  ! divide by N_box_vel here
                                                                                ! to account (after taking the square root) 
                                                                                ! the normalization of velocity
                                                                                ! produced by GetMaxwellVelocity
  Factor_D_3_sqrt = 1.414213562_8 * Factor_D_1_sqrt

END SUBROUTINE INITIATE_COULOMB_COLLISIONS

!--------------------------------------------------------------------------
! Note, this subroutine is called in INITIATE_COULOMB_COLLISIONS
!
SUBROUTINE ALLOCATE_LANGEVIN_ARRAYS

  USE LangevinCollisions
  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : N_cells

  IMPLICIT NONE

  INTEGER ALLOC_ERR

  ALLOCATE(    Fd_w(1:N_box_w_Lang, 0:N_cells-1), STAT=ALLOC_ERR)    
  ALLOCATE(  F_drag(0:N_box_w_Lang, 0:N_cells-1), STAT=ALLOC_ERR)    
  ALLOCATE(D_1_sqrt(0:N_box_w_Lang, 0:N_cells-1), STAT=ALLOC_ERR)    
  ALLOCATE(D_3_sqrt(0:N_box_w_Lang, 0:N_cells-1), STAT=ALLOC_ERR)    

  ALLOCATE(V_ion_threshold(0:N_cells-1), STAT=ALLOC_ERR)
  ALLOCATE(         Vx_avg(0:N_cells-1), STAT=ALLOC_ERR)
  ALLOCATE(         Vy_avg(0:N_cells-1), STAT=ALLOC_ERR)
  ALLOCATE(         Vz_avg(0:N_cells-1), STAT=ALLOC_ERR)
  ALLOCATE(  energy_before(0:N_cells-1), STAT=ALLOC_ERR)
  ALLOCATE(   energy_after(0:N_cells-1), STAT=ALLOC_ERR)
  ALLOCATE(       Log_coul(0:N_cells-1), STAT=ALLOC_ERR)

  IF (Rank_of_process.GT.0) THEN        ! for client processes only

     ALLOCATE(   V_ion_threshold_2(0:N_cells-1), STAT=ALLOC_ERR)
     ALLOCATE(V_ion_threshold_sqrt(0:N_cells-1), STAT=ALLOC_ERR)
     ALLOCATE(      Density_e_sqrt(0:N_cells-1), STAT=ALLOC_ERR)    

  END IF

END SUBROUTINE ALLOCATE_LANGEVIN_ARRAYS

!--------------------------------------------------------------------------
!
SUBROUTINE PROCESS_COLLISIONS_EL_EL_LNGV

  USE ParallelOperationValues
  USE LangevinCollisions
  USE CurrentProblemValues 
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER k                                       ! index of particle
  INTEGER n_of_loc                                ! index of location (left node of cell containing the particle)
  REAL(8) wx, wy, wz                              ! velocity components of an electron in the electron flow frame
  REAL(8) abs_w                                   ! absolute value of the velocity of an electron in the electron flow frame
  REAL(8) perp_w                                  ! electron velocity in the electron flow frame, perpendicular to the relative velocity
  INTEGER indx_w                                  ! index of absolute velocity
  REAL(8), ALLOCATABLE :: alpha(:)                ! ratio of energy before to energy after
  REAL(8) SinTeta, CosTeta
  REAL(8) SinFi, CosFi

  REAL(8) drag_force, diffus_1, diffus_3          ! values of drag force and velocity diff-n coeff-s for interpolation
  REAL(8) dvx_drag, dvy_drag, dvz_drag            ! increment of electron velocity due to drag force
  REAL(8) Q1, Q2, Q3                              ! vector, distributed randomly according to Maxwellian distribution
  REAL(8) dvx_dif, dvy_dif, dvz_dif               ! increment of electron velocity due to diffusion in velocity space

  REAL(8), ALLOCATABLE :: rbufer(:), rbufer2(:)
  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER ierr

  IF (Accounted_Langevin.EQ.0) RETURN             ! if the Langevin model is turned off

  IF (MOD(T_cntr,T_skip_Lang).NE.0) RETURN   ! we apply the Langevin collisions model NOT each timestep,
                                             ! but with greater time intervals

  CALL CALCULATE_COEFFICIENTS_LNGV           ! the flow velocity and the electron energy before collisions are
                                             ! calculated here

  ALLOCATE(rbufer(1:N_cells),  STAT=ALLOC_ERR)
  ALLOCATE(rbufer2(1:N_cells), STAT=ALLOC_ERR)
  ALLOCATE(alpha(0:N_cells-1), STAT=ALLOC_ERR)

! clear the final energy storages
  energy_after  = 0.0_8

  IF (Rank_of_process.GT.0) THEN  !------------------------------------------------------------------ CLIENT

! cycle over all electrons of this process
     DO k = 1, N_part(1)
! define index of location
        n_of_loc  = INT(species(1)%part(k)%X)
! calculate the electron velocity in the electron flow frame 
        wx = species(1)%part(k)%VX - Vx_avg(n_of_loc)
        wy = species(1)%part(k)%VY - Vy_avg(n_of_loc)
        wz = species(1)%part(k)%VZ - Vz_avg(n_of_loc)
        abs_w = SQRT(wx*wx + wy*wy + wz*wz)      
! interpolate drag force coefficient and velocity diffusion coefficients
        indx_w = INT(abs_w)
        IF (indx_w.LT.N_box_w_Lang) THEN
           drag_force =   F_drag(indx_w,n_of_loc) + (  F_drag(indx_w+1,n_of_loc) -   F_drag(indx_w,n_of_loc)) * (abs_w - indx_w)
           diffus_1   = D_1_sqrt(indx_w,n_of_loc) + (D_1_sqrt(indx_w+1,n_of_loc) - D_1_sqrt(indx_w,n_of_loc)) * (abs_w - indx_w)
           diffus_3   = D_3_sqrt(indx_w,n_of_loc) + (D_3_sqrt(indx_w+1,n_of_loc) - D_3_sqrt(indx_w,n_of_loc)) * (abs_w - indx_w)
        ELSE
           drag_force =   F_drag(N_box_w_Lang,n_of_loc)
           diffus_1   = D_1_sqrt(N_box_w_Lang,n_of_loc)
           diffus_3   = D_3_sqrt(N_box_w_Lang,n_of_loc)
        END IF
! Note, the time factor is already included in Langevin coefficients, see CALCULATE_COEFFICIENTS_LNGV

! ### future work: since the coefficients below are local, it is reasonable to include the densities into the coefficients

        drag_force = drag_force * Density_e_sqrt(INT(species(1)%part(k)%X))**2   ! we store the root of the density -
                                                                                 ! one extra "*" is faster than three "sqrt"
        diffus_1   = diffus_1   * Density_e_sqrt(INT(species(1)%part(k)%X))
        diffus_3   = diffus_3   * Density_e_sqrt(INT(species(1)%part(k)%X))

! process the possible zero-velocity particle correctly
        IF (abs_w.EQ.0.0_8) THEN
! calculate the velocity change due to drag force
           dvx_drag = 0.0_8 
           dvy_drag = 0.0_8
           dvz_drag = 0.0_8
! Calculate the trigonometric functions 
           SinTeta = 0.0_8
           CosTeta = 1.0_8
        ELSE
! calculate the velocity change due to drag force
           dvx_drag = drag_force * wx / abs_w
           dvy_drag = drag_force * wy / abs_w
           dvz_drag = drag_force * wz / abs_w
! Calculate the trigonometric functions 
           SinTeta = SQRT(wx*wx + wy*wy) / abs_w
           CosTeta = wz / abs_w
        END IF

        perp_w = SQRT(wx*wx + wy*wy)
        IF (perp_w.EQ.0.0_8) THEN
           SinFi = 0.0_8
           CosFi = 1.0_8
        ELSE 
           SinFi =  wx / perp_w
           CosFi = -wy / perp_w
        END IF
! obtain random velocities according to Gaussian (Maxwell) distribution
! note, SUBROUTINE PrepareMaxwellDistribIntegral must be called before once, now it is called in  SUBROUTINE INITIATE_MC_COLLISIONS 
! Also, here Q1-3 are obtained as the values normalized by [v_th_e / N_box_vel], while they are supposed to be in units of v_th_e
! this extra factor should be accounted for in CALCULATE_COEFFICIENTS_LNGV
        CALL GetMaxwellVelocity(Q1) 
        CALL GetMaxwellVelocity(Q2) 
        CALL GetMaxwellVelocity(Q3) 

! transform the random gaussian velocity to the velocity increment (scale and turn), 
! the factor of 2 (inside the roots) is included in CALCULATE_COEFFICIENTS_LNGV
        Q1 = Q1 * diffus_1 !SQRT(diffus_1)
        Q2 = Q2 * diffus_1 !SQRT(diffus_1) !we will calculate and interpolate already the square roots
        Q3 = Q3 * diffus_3 !SQRT(diffus_3)

        dvx_dif = CosFi * Q1 - CosTeta * SinFi * Q2 + SinTeta * SinFi * Q3
        dvy_dif = SinFi * Q1 + CosTeta * CosFi * Q2 - SinTeta * CosFi * Q3
        dvz_dif =              SinTeta         * Q2 + CosTeta         * Q3

! process the electron-ion collisions if they are turned on --------------------------------------------------------------->>>
        IF (Accounted_Langevin.EQ.2) THEN
           
! calculate the absolute electron velocity (assume that ions are at rest for now) 
           wx = species(1)%part(k)%VX !- Vx_avg
           wy = species(1)%part(k)%VY !- Vy_avg
           wz = species(1)%part(k)%VZ !- Vz_avg
           abs_w = SQRT(wx * wx + wy * wy + wz * wz)      
! we use electron density for simplicity (although in future it should be changed for the ion density)
! below 0.039788736 = 1/(8*pi)
!       0.488602512 = SQRT{3/(4*pi)}
           drag_force = - 0.039788736_8 * Factor_F_drag   * Log_coul(n_of_loc) * Density_e_sqrt(INT(species(1)%part(k)%X))**2  
           diffus_1   =   0.488602512_8 * Factor_D_1_sqrt * SQRT(Log_coul(n_of_loc)) * Density_e_sqrt(INT(species(1)%part(k)%X))
! below the threshold velocity the langevin coefficient are set constant
! for now we neglect the modification of Langevin coefficients by the ion EVDF at low
! electron energy
           IF (abs_w.LT.V_ion_threshold(n_of_loc)) THEN
              drag_force = drag_force / V_ion_threshold_2(n_of_loc)
              diffus_1   =   diffus_1 / V_ion_threshold_sqrt(n_of_loc)
           ELSE
              drag_force = drag_force / abs_w**2
              diffus_1   =   diffus_1 / SQRT(abs_w)
           END IF

! process the possible zero-velocity particle correctly
           IF (abs_w.EQ.0.0_8) THEN
! the velocity change due to drag force is zero
!              dvx_drag = dvx_drag + 0.0_8 
!              dvy_drag = dvy_drag + 0.0_8
!              dvz_drag = dvz_drag + 0.0_8
! Calculate the trigonometric functions 
              SinTeta = 0.0_8
              CosTeta = 1.0_8
           ELSE
! calculate the velocity change due to drag force
              dvx_drag = dvx_drag + drag_force * wx / abs_w
              dvy_drag = dvy_drag + drag_force * wy / abs_w
              dvz_drag = dvz_drag + drag_force * wz / abs_w
! Calculate the trigonometric functions 
              SinTeta = SQRT(wx*wx + wy*wy) / abs_w
              CosTeta = wz / abs_w
           END IF

           perp_w = SQRT(wx*wx + wy*wy)
           IF (perp_w.EQ.0.0_8) THEN
              SinFi = 0.0_8
              CosFi = 1.0_8
           ELSE 
              SinFi =  wx / perp_w
              CosFi = -wy / perp_w
           END IF

           CALL GetMaxwellVelocity(Q1) 
           CALL GetMaxwellVelocity(Q2) 
           Q1 = Q1 * diffus_1 
           Q2 = Q2 * diffus_1 

           dvx_dif = dvx_dif + CosFi * Q1 - CosTeta * SinFi * Q2 
           dvy_dif = dvy_dif + SinFi * Q1 + CosTeta * CosFi * Q2 
           dvz_dif = dvz_dif +              SinTeta         * Q2

        END IF !-----------------------------------------------------------------------------------------------------------<<<

! calculate the velocity of particle after collision
        species(1)%part(k)%VX = species(1)%part(k)%VX + dvx_drag + dvx_dif
        species(1)%part(k)%VY = species(1)%part(k)%VY + dvy_drag + dvy_dif
        species(1)%part(k)%VZ = species(1)%part(k)%VZ + dvz_drag + dvz_dif

! collect the electron energy after collision IN THIS CELL ONLY, NO NEIGHBORS
        energy_after(n_of_loc) = energy_after(n_of_loc) + (species(1)%part(k)%VX**2 + species(1)%part(k)%VY**2 + species(1)%part(k)%VZ**2)

     END DO

! send "energy after collision" collected in this process to the server node
     rbufer(1:N_cells)  = energy_after(0:N_cells-1)
     rbufer2 = 0.0_8
     CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! receive the normalization coefficient
     rbufer = 0.0_8
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     alpha(0:N_cells-1) = rbufer(1:N_cells)

! renormalize all electron velocities to ensure the energy conservation
     DO k = 1, N_part(1)
        n_of_loc  = INT(species(1)%part(k)%X)
        species(1)%part(k)%VX = species(1)%part(k)%VX * alpha(n_of_loc)
        species(1)%part(k)%VY = species(1)%part(k)%VY * alpha(n_of_loc)
        species(1)%part(k)%VZ = species(1)%part(k)%VZ * alpha(n_of_loc)
     END DO

  ELSE  !----------------------------------------------------------------------------- SERVER

! collect "energy after collision" sent by client processes
     rbufer  = 0.0_8
     rbufer2 = 0.0_8
     CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     energy_after(0:N_cells-1) = rbufer(1:N_cells)
! calculate the normalization coefficient
     alpha = 1.0_8
     DO n_of_loc = 0, N_cells-1
        IF (energy_after(n_of_loc).GT.0.0_8) alpha(n_of_loc) = SQRT(energy_before(n_of_loc) / energy_after(n_of_loc))
     END DO
! send the normalization coefficient
     rbufer(1:N_cells) = alpha(0:N_cells-1)
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  END IF

  DEALLOCATE(rbufer, STAT=DEALLOC_ERR)
  DEALLOCATE(rbufer2, STAT=DEALLOC_ERR)
  DEALLOCATE(alpha, STAT=DEALLOC_ERR)

END SUBROUTINE PROCESS_COLLISIONS_EL_EL_LNGV 

!--------------------------------------------------------------------------
!
SUBROUTINE CALCULATE_COEFFICIENTS_LNGV

  USE LangevinCollisions
  USE CurrentProblemValues
  USE ParallelOperationValues
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER k      ! particle index
  INTEGER j
  INTEGER n      ! location index

  REAL(8) Fd_w_integral
  INTEGER n_of_loc                      ! node to the left of particle
  INTEGER N_cells_in_loc                ! number of cells associated with given locations
                                        ! usually it is 1+2*added_cells but near the edges of the simulation area
                                        ! it gradually decreases to 1+added_cells

!  INTEGER Fd_wx(0:(N_box_w_Lang-1))                              ! non-normalized distribution functions over velocity components,
!  INTEGER Fd_wy(0:(N_box_w_Lang-1))                              ! averaged over positive and negative velocity,
!  INTEGER Fd_wz(0:(N_box_w_Lang-1))                              ! to be averaged to produce the final EVDF

  INTEGER, ALLOCATABLE :: N_part_in_loc(:)                        ! number of particles associated with the given location
  INTEGER, ALLOCATABLE :: Fd_wx(:,:)                              ! non-normalized distribution functions over velocity components,
  INTEGER, ALLOCATABLE :: Fd_wy(:,:)                              ! averaged over positive and negative velocity,
  INTEGER, ALLOCATABLE :: Fd_wz(:,:)                              ! to be averaged to produce the final EVDF

  INTEGER indx_wx, indx_wy, indx_wz

  REAL(8) Int_1(0:N_box_w_Lang)                                  !
  REAL(8) Int_2(0:N_box_w_Lang)                                  ! integrals of distribution function
  REAL(8) Int_4(0:N_box_w_Lang)                                  ! 

  REAL(8) average_energy_eV
  REAL(8) average_density_m3
  REAL(8) debye_length_m
  REAL(8) w
  REAL(8) max_rho

  REAL(8), ALLOCATABLE :: rbufer(:), rbufer2(:)
  INTEGER, ALLOCATABLE :: ibufer(:), ibufer2(:)
  INTEGER buf_dim
  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER ierr

  IF (Accounted_Langevin.EQ.0) RETURN

! determine the maximal necessary buffer size
  buf_dim = MAX(N_cells, N_box_w_Lang+1)

! allocate buffers
  ALLOCATE(rbufer( 1:buf_dim), STAT=ALLOC_ERR)
  ALLOCATE(rbufer2(1:buf_dim), STAT=ALLOC_ERR)

  ALLOCATE(ibufer( 1:N_cells), STAT=ALLOC_ERR)
  ALLOCATE(ibufer2(1:N_cells), STAT=ALLOC_ERR)

  ALLOCATE(N_part_in_loc(0:N_cells-1), STAT=ALLOC_ERR)

! clear the initial energy and flow velocity storages
  energy_before = 0.0_8
  Vx_avg = 0.0_8  
  Vy_avg = 0.0_8  
  Vz_avg = 0.0_8   
  N_part_in_loc = 0

  IF (Rank_of_process.GT.0) THEN    !---------------------------------------------------------------------------- CLIENT PROCESS

     ALLOCATE(Fd_wx(0:(N_box_w_Lang-1),0:N_cells-1), STAT=ALLOC_ERR)
     ALLOCATE(Fd_wy(0:(N_box_w_Lang-1),0:N_cells-1), STAT=ALLOC_ERR)
     ALLOCATE(Fd_wz(0:(N_box_w_Lang-1),0:N_cells-1), STAT=ALLOC_ERR)

     DO k = 1, N_part(1)
! define location of particle
        n_of_loc  = INT(species(1)%part(k)%X)
! collect the electron energy before collision IN THIS CELL ONLY, NO NEIGHBORS !!!
        energy_before(n_of_loc) = energy_before(n_of_loc) + (species(1)%part(k)%VX**2 + species(1)%part(k)%VY**2 + species(1)%part(k)%VZ**2)
! collect the electron velocity and number of particles over several neighbor cells (to reduce noise)
        DO n = MAX(0,n_of_loc-added_cells), MIN(n_of_loc+added_cells,N_cells-1)
           Vx_avg(n) = Vx_avg(n) + species(1)%part(k)%VX
           Vy_avg(n) = Vy_avg(n) + species(1)%part(k)%VY      
           Vz_avg(n) = Vz_avg(n) + species(1)%part(k)%VZ
           N_part_in_loc(n) = N_part_in_loc(n)+1
        END DO
     END DO

! send collected velocity and energy to the server process
! Vx_avg
     rbufer(1:N_cells) = Vx_avg(0:N_cells-1)
     rbufer2 = 0.0_8
     CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
! Vy_avg
     rbufer(1:N_cells) = Vy_avg(0:N_cells-1)
     rbufer2 = 0.0_8
     CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
! Vz_avg
     rbufer(1:N_cells) = Vz_avg(0:N_cells-1)
     rbufer2 = 0.0_8
     CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
! energy_before
     rbufer(1:N_cells) = energy_before(0:N_cells-1)
     rbufer2 = 0.0_8
     CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! send the number of particles to the server process
     ibufer(1:N_cells) = N_part_in_loc(0:N_cells-1)
     ibufer2 = 0
     CALL MPI_REDUCE(ibufer, ibufer2, N_cells, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! receive the electron flow velocity, the coulomb logarithm and the ion threshold velocity from the server process
! Vx_avg
     rbufer = 0.0_8
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     Vx_avg(0:N_cells-1) = rbufer(1:N_cells)
! Vy_avg
     rbufer = 0.0_8
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     Vy_avg(0:N_cells-1) = rbufer(1:N_cells)
! Vz_avg
     rbufer = 0.0_8
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     Vz_avg(0:N_cells-1) = rbufer(1:N_cells)
! Log_coul
     rbufer = 0.0_8
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     Log_coul(0:N_cells-1) = rbufer(1:N_cells)
! V_ion_threshold
     rbufer = 0.0_8
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     V_ion_threshold(0:N_cells-1) = rbufer(1:N_cells)

     V_ion_threshold_2    = V_ion_threshold**2  
     V_ion_threshold_sqrt = SQRT(V_ion_threshold)

! build the local (for this process) EVDF      
     Fd_wx = 0
     Fd_wy = 0
     Fd_wz = 0
     DO k = 1, N_part(1)            
        n_of_loc = INT(species(1)%part(k)%X)

        IF ((n_of_loc.LT.0).OR.(n_of_loc.GE.N_cells)) CYCLE
!        IF ((left_node.LT.N_of_left_node_Lang).OR.(left_node.GE.N_of_right_node_Lang)) CYCLE

        indx_wx = MIN( INT( ABS(species(1)%part(k)%VX - Vx_avg(n_of_loc)) ), N_box_w_Lang-1 )
        indx_wy = MIN( INT( ABS(species(1)%part(k)%VY - Vy_avg(n_of_loc)) ), N_box_w_Lang-1 )
        indx_wz = MIN( INT( ABS(species(1)%part(k)%VZ - Vz_avg(n_of_loc)) ), N_box_w_Lang-1 )

!        indx_wx = INT( ABS(VX_of_spec(1)%part(k) - Vx_avg(n_of_loc)) )
!        indx_wy = INT( ABS(VY_of_spec(1)%part(k) - Vy_avg(n_of_loc)) )
!        indx_wz = INT( ABS(VZ_of_spec(1)%part(k) - Vz_avg(n_of_loc)) )
!        IF ( (indx_wx.GE.N_box_w_Lang).OR.(indx_wy.GE.N_box_w_Lang).OR.(indx_wz.GE.N_box_w_Lang) ) CYCLE

        DO n = MAX(0,n_of_loc-added_cells), MIN(n_of_loc+added_cells,N_cells-1)
           Fd_wx(indx_wx, n) = Fd_wx(indx_wx, n) + 1
           Fd_wy(indx_wy, n) = Fd_wy(indx_wy, n) + 1
           Fd_wz(indx_wz, n) = Fd_wz(indx_wz, n) + 1
        END DO
     END DO

     DO n_of_loc = 0, N_cells-1
        DO j = 0, N_box_w_Lang-1 
           Fd_w(j+1, n_of_loc) = DBLE(Fd_wx(j, n_of_loc) + Fd_wy(j, n_of_loc) + Fd_wz(j, n_of_loc))
        END DO
     END DO

! send the local EVDF to the server node
     DO n_of_loc = 0, N_cells-1
        rbufer(1:N_box_w_Lang) = Fd_w(1:N_box_w_Lang, n_of_loc)
        rbufer2 = 0.0_8
        CALL MPI_REDUCE(rbufer, rbufer2, N_box_w_Lang, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     END DO

! RECEIVE the tabulated Langevin coefficients

! F_drag
     DO n_of_loc = 0, N_cells-1
        rbufer = 0.0_8
        CALL MPI_BCAST(rbufer, (N_box_w_Lang+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        F_drag(0:N_box_w_Lang, n_of_loc) = rbufer(1:(N_box_w_Lang+1))
     END DO

! D_1_sqrt
     DO n_of_loc = 0, N_cells-1
        rbufer = 0.0_8
        CALL MPI_BCAST(rbufer, (N_box_w_Lang+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        D_1_sqrt(0:N_box_w_Lang, n_of_loc) = rbufer(1:(N_box_w_Lang+1))
     END DO

! D_3_sqrt
     DO n_of_loc = 0, N_cells-1
        rbufer = 0.0_8
        CALL MPI_BCAST(rbufer, (N_box_w_Lang+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        D_3_sqrt(0:N_box_w_Lang, n_of_loc) = rbufer(1:(N_box_w_Lang+1))
     END DO

! receive the array of the square root of streaming density in the middles of cells
     rbufer = 0.0_8
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     Density_e_sqrt(0:(N_cells-1)) = rbufer(1:N_cells)

     DEALLOCATE(Fd_wx, STAT=DEALLOC_ERR)
     DEALLOCATE(Fd_wy, STAT=DEALLOC_ERR)
     DEALLOCATE(Fd_wz, STAT=DEALLOC_ERR)

  ELSE    ! ------------------------------------------------------------------------------------- SERVER
! receive collected velocity and energy from client processes
! Vx_avg
     rbufer  = 0.0_8
     rbufer2 = 0.0_8
     CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     Vx_avg(0:N_cells-1) = rbufer(1:N_cells)
! Vy_avg
     rbufer  = 0.0_8
     rbufer2 = 0.0_8
     CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     Vy_avg(0:N_cells-1) = rbufer(1:N_cells)
! Vz_avg
     rbufer  = 0.0_8
     rbufer2 = 0.0_8
     CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     Vz_avg(0:N_cells-1) = rbufer(1:N_cells)
! energy_before
     rbufer  = 0.0_8
     rbufer2 = 0.0_8
     CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     energy_before(0:N_cells-1) = rbufer(1:N_cells)

! receive the number of particles from client processe
     ibufer  = 0
     ibufer2 = 0
     CALL MPI_REDUCE(ibufer2, ibufer, N_cells, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     N_part_in_loc(0:N_cells-1) = ibufer(1:N_cells)

! for every cell
     DO n_of_loc = 0, N_cells-1

! prepare for calculation of average energy and density for the coulomb logarithm
        N_cells_in_loc = 0
        average_energy_eV = 0.0_8
        DO n = MAX(0,n_of_loc-added_cells), MIN(n_of_loc+added_cells,N_cells-1)
           average_energy_eV = average_energy_eV +  energy_before(n)
! calculate number of cells associated with the given location (to find average density)
           N_cells_in_loc = N_cells_in_loc + 1
        END DO

! calculate the electron flow velocity and coulomb logarithm
        IF (N_part_in_loc(n_of_loc).GT.0) THEN
           Vx_avg(n_of_loc) = Vx_avg(n_of_loc) / N_part_in_loc(n_of_loc)
           Vy_avg(n_of_loc) = Vy_avg(n_of_loc) / N_part_in_loc(n_of_loc)
           Vz_avg(n_of_loc) = Vz_avg(n_of_loc) / N_part_in_loc(n_of_loc)
! calculate the Coulomb Logarithm
           average_energy_eV  = average_energy_eV * T_e_eV / ( N_box_vel * N_box_vel * DBLE(N_part_in_loc(n_of_loc)) )
           average_density_m3 = N_plasma_m3 * DBLE(N_part_in_loc(n_of_loc)) / DBLE(N_cells_in_loc*N_of_particles_cell)
           debye_length_m     = SQRT(2.0_8 * average_energy_eV * eps_0_Fm / (e_Cl * average_density_m3))
           Log_coul(n_of_loc) = LOG(3.141592654_8 * eps_0_Fm * average_energy_eV * debye_length_m / e_Cl) 
        ELSE
           Vx_avg(n_of_loc) = 0.0_8
           Vy_avg(n_of_loc) = 0.0_8
           Vz_avg(n_of_loc) = 0.0_8
           Log_coul(n_of_loc) = 10.0_8
        END IF
! calculate the threshold velocity for ion-electron collisions
! we calculate the threshold as the velocity value when the correction due to the e-i drag force
! becomes comparable with the velocity itself
!     max_rho = 0.0_8
!     DO j = 0, N_cells
!        IF (Q_strm_spec(j,1).GT.max_rho) max_rho = Q_strm_spec(j,1)
!     END DO
!     V_ion_threshold      = (0.039788736_8 * Factor_F_drag * Log_coul * max_rho)**0.3333333_8
        max_rho = 0.5_8 * (Q_strm_spec(n_of_loc,1) + Q_strm_spec(n_of_loc+1,1))
        V_ion_threshold(n_of_loc) = (0.039788736_8 * Factor_F_drag * Log_coul(n_of_loc) * max_rho)**0.3333333_8
!print *, '***********', V_ion_threshold
!        V_ion_threshold_2    = V_ion_threshold**2  
!        V_ion_threshold_sqrt = SQRT(V_ion_threshold)
     END DO

! send the electron flow velocity, the coulomb logarithm and the threshold velocity to client processes
! Vx_avg
     rbufer(1:N_cells) = Vx_avg(0:N_cells-1) 
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
! Vy_avg
     rbufer(1:N_cells) = Vy_avg(0:N_cells-1) 
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
! Vz_avg
     rbufer(1:N_cells) = Vz_avg(0:N_cells-1) 
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
! Log_coul
     rbufer(1:N_cells) = Log_coul(0:N_cells-1)
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
! V_ion_threshold
     rbufer(1:N_cells) = V_ion_threshold(0:N_cells-1)
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! receive the local EVDFs from client processes
     DO n_of_loc = 0, N_cells-1
        rbufer  = 0.0_8
        rbufer2 = 0.0_8
        CALL MPI_REDUCE(rbufer2, rbufer, N_box_w_Lang, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        Fd_w(1:N_box_w_Lang, n_of_loc) = rbufer(1:N_box_w_Lang)
     END DO

! for every cell
     DO n_of_loc = 0, N_cells-1
! normalize the collected EVDF
        Fd_w_integral = 0.0_8
        DO j = 1, N_box_w_Lang   
           Fd_w_integral = Fd_w_integral + (DBLE(j)-0.5_8)**2 * Fd_w(j, n_of_loc) 
        END DO                                                                    ! 
        Fd_w(1:N_box_w_Lang, n_of_loc) = Fd_w(1:N_box_w_Lang, n_of_loc) / (Fd_w_integral * 12.56637061_8)  ! (4*pi*integral = 1), here dw = 1
! tabulate integrals #1,2,4
        Int_1 = 0.0_8 
        Int_2 = 0.0_8                                                             
        Int_4 = 0.0_8 
        DO j = 1, N_box_w_Lang 
           w = DBLE(j)-0.5_8
           Int_1(j) = Int_1(j-1) + w    * Fd_w(j, n_of_loc)
           Int_2(j) = Int_2(j-1) + w**2 * Fd_w(j, n_of_loc)
           Int_4(j) = Int_4(j-1) + w**4 * Fd_w(j, n_of_loc)
        END DO
        DO j = 0, N_box_w_Lang 
           Int_1(j) = Int_1(N_box_w_Lang) - Int_1(j)  
        END DO
! calculate the tabulated Langevin coefficients
        F_drag(0, n_of_loc) = 0.0_8
        D_1_sqrt(0, n_of_loc)    = Factor_D_1_sqrt * SQRT(Log_coul(n_of_loc) * 2.0_8 * Int_1(0))
        D_3_sqrt(0, n_of_loc)    = Factor_D_3_sqrt * SQRT(Log_coul(n_of_loc) * Int_1(0))
        DO j = 1, N_box_w_Lang
           w = DBLE(j)-0.5_8 
           F_drag(j, n_of_loc) = - Factor_F_drag * Log_coul(n_of_loc) * Int_2(j) / (w**2)      
           D_1_sqrt(j, n_of_loc) = Factor_D_1_sqrt * SQRT(Log_coul(n_of_loc) * ( 3.0_8 * Int_2(j) / w - Int_4(j) / (w**3) + 2.0_8 * Int_1(j)))
           D_3_sqrt(j, n_of_loc) = Factor_D_3_sqrt * SQRT(Log_coul(n_of_loc) * (Int_4(j) / (w**3) + Int_1(j)))
        END DO
     END DO

! SEND the tabulated Langevin coefficients

! F_drag
     DO n_of_loc = 0, N_cells-1
        rbufer(1:(N_box_w_Lang+1)) = F_drag(0:N_box_w_Lang, n_of_loc)
        CALL MPI_BCAST(rbufer, (N_box_w_Lang+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     END DO

! D_1_sqrt
     DO n_of_loc = 0, N_cells-1
        rbufer(1:(N_box_w_Lang+1)) = D_1_sqrt(0:N_box_w_Lang, n_of_loc)
        CALL MPI_BCAST(rbufer, (N_box_w_Lang+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     END DO

! D_3_sqrt
     DO n_of_loc = 0, N_cells-1
        rbufer(1:(N_box_w_Lang+1)) = D_3_sqrt(0:N_box_w_Lang, n_of_loc)
        CALL MPI_BCAST(rbufer, (N_box_w_Lang+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     END DO

! prepare the array of the square root of streaming density in the middles of cells
     rbufer(1) = SQRT(0.5_8 * (2.0_8 * Q_strm_spec(0,1) + Q_strm_spec(1,1)))
     DO j = 1, N_cells - 2
        rbufer(j+1) = SQRT(0.5_8 * (Q_strm_spec(j,1) + Q_strm_spec(j+1,1)))
     END DO
     rbufer(N_cells) = SQRT(0.5_8 * (2.0_8 * Q_strm_spec(N_cells,1) + Q_strm_spec(N_cells-1,1)))
! receive the streaming density
     CALL MPI_BCAST(rbufer, N_cells, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  END IF

! deallocate buffers
  DEALLOCATE(rbufer, STAT=DEALLOC_ERR)
  DEALLOCATE(rbufer2, STAT=DEALLOC_ERR)
  DEALLOCATE(ibufer, STAT=DEALLOC_ERR)
  DEALLOCATE(ibufer2, STAT=DEALLOC_ERR)

  DEALLOCATE(N_part_in_loc, STAT=DEALLOC_ERR) 

END SUBROUTINE CALCULATE_COEFFICIENTS_LNGV

!--------------------------------------------------------------------------
! Note, this subroutine is called in INITIATE_COULOMB_COLLISIONS
!
SUBROUTINE REMOVE_LANGEVIN_ARRAYS

  USE LangevinCollisions
  USE ParallelOperationValues
  IMPLICIT NONE
  INTEGER DEALLOC_ERR

  IF (Accounted_Langevin.EQ.0) RETURN 

  DEALLOCATE(    Fd_w, STAT=DEALLOC_ERR)    
  DEALLOCATE(  F_drag, STAT=DEALLOC_ERR)    
  DEALLOCATE(D_1_sqrt, STAT=DEALLOC_ERR)    
  DEALLOCATE(D_3_sqrt, STAT=DEALLOC_ERR)    

  DEALLOCATE(V_ion_threshold, STAT=DEALLOC_ERR)
  DEALLOCATE(         Vx_avg, STAT=DEALLOC_ERR)
  DEALLOCATE(         Vy_avg, STAT=DEALLOC_ERR)
  DEALLOCATE(         Vz_avg, STAT=DEALLOC_ERR)
  DEALLOCATE(  energy_before, STAT=DEALLOC_ERR)
  DEALLOCATE(   energy_after, STAT=DEALLOC_ERR)
  DEALLOCATE(       Log_coul, STAT=DEALLOC_ERR)

  IF (ALLOCATED(V_ion_threshold_2)) THEN
     DEALLOCATE(V_ion_threshold_2, STAT=DEALLOC_ERR)    
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE V_ion_threshold_2 !!!")', Rank_of_process
        PRINT '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(V_ion_threshold_sqrt)) THEN
     DEALLOCATE(V_ion_threshold_sqrt, STAT=DEALLOC_ERR)    
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE V_ion_threshold_sqrt !!!")', Rank_of_process
        PRINT '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(Density_e_sqrt)) THEN
     DEALLOCATE(Density_e_sqrt, STAT=DEALLOC_ERR)    
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Density_e_sqrt !!!")', Rank_of_process
        PRINT '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

END SUBROUTINE REMOVE_LANGEVIN_ARRAYS
