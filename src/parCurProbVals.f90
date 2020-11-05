!===================================================================================================
!
SUBROUTINE INITIATE_PARAMETERS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE Diagnostics
  USE Snapshots

  USE rng_wrapper

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL exists
  INTEGER i, s, j
  INTEGER itmp
  INTEGER flag_ions_magnetized
  INTEGER flag_ions_sense_Ez
  INTEGER flag_trace_z_shift
  INTEGER flag_trace_z_shift_ion
  INTEGER flag_use_ExB_for_electrons
  INTEGER flag_use_ExB_for_ions
  REAL(8) Vx_drift        ! dim-less ExB-drift x-velocity, -Ez/By
  REAL(8) Vy_drift        ! dim-less ExB-drift y-velocity,  Ez/Bx

  REAL(8) L_plasma_deb                    ! plasma length in debye lengths                  ! ####### ######### ######## was INTEGER ######## ###### ######

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)

  INTEGER init_random_seed
  INTEGER my_random_seed
  REAL(8) myran

  CHARACTER (77) buf
  CHARACTER(3)   procnumber_txt
  INTEGER proc_blnks

  INTEGER dummy_i
  REAL(8) dummy_r
  REAL(8) Energy_ebeam_eV

  INQUIRE (FILE = 'ssc_initial.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  OPEN (9, FILE = 'ssc_initial.dat')

  IF(exists) THEN

     IF (Rank_of_process.EQ.0) THEN
        PRINT '(2x,"Process ",i3," : ssc_initial.dat is found. Reading the data file...")', Rank_of_process
     END IF

     READ (9, '(A77)') buf !============================ SYSTEM CONFIGURATION ===========================")')
     READ (9, '(A77)') buf !-------d----- BC: Left wall potential is fixed/floating/ext.circuit {0/1/2}--")')
     READ (9, '(7x,i1)') BC_flag
     READ (9, '(A77)') buf !--dddddd.ddd- Resistance of the resistor in the external circuit (ohm) ------")')
     READ (9, '(2x,f10.3)') R_ext_ohm
     READ (9, '(A77)') buf !--dddddd.ddd- Electrode area (cm^2) -----------------------------------------")')
     READ (9, '(2x,f10.3)') S_electrode_cm2
     READ (9, '(A77)') buf !--dddddd.ddd- External z-electric field parallel to walls (V/m) -------------")')
     READ (9, '(2x,f10.3)') E_z_ext_Vm
     READ (9, '(A77)') buf !--dddddd.ddd- External x-magnetic field normal to walls (Gauss) -------------")')
     READ (9, '(2x,f10.3)') B_x_ext_Gs
     READ (9, '(A77)') buf !--dddddd.ddd- External y-magnetic field parallel to walls (Gauss) -----------")')
     READ (9, '(2x,f10.3)') B_y_ext_Gs
     READ (9, '(A77)') buf !--dddddd.ddd- Left wall potential/battery voltage (V) -----------------------")')
     READ (9, '(2x,f10.3)') U_ext_V
     READ (9, '(A77)') buf !--dddd.ddddd- Plasma layer width (m, beam wavelengths if negative) ----------")')
     READ (9, '(2x,f10.5)') L_plasma_m
     READ (9, '(A77)') buf !--dddddd.ddd- Ion mass (a.m.u.) ---------------------------------------------")')
     READ (9, '(2x,f10.3)') M_i_amu
     READ (9, '(A77)') buf !--dddddd.ddd- Ion temperature (eV) ------------------------------------------")')
     READ (9, '(2x,f10.3)') T_i_eV
     READ (9, '(A77)') buf !-------d----- Ions sense the external MAGNETIC field (1/0=Yes/No) -----------")')
     READ (9, '(7x,i1)') flag_ions_magnetized
     READ (9, '(A77)') buf !-------d----- Ions sense the external ELECTRIC field (1/0=Yes/No) -----------")')
     READ (9, '(7x,i1)') flag_ions_sense_Ez
     READ (9, '(A77)') buf !-------d----- Thermalize electron if Z-shift exceeds limit (1/0=Yes/No) -----")')
     READ (9, '(7x,i1)') flag_trace_z_shift
     READ (9, '(A77)') buf !--dddd.ddddd- Limit for the electron Z-shift (m) ----------------------------")')
     READ (9, '(2x,f10.5)') max_Z_shift
     READ (9, '(A77)') buf !-------d----- Thermalize ion if Z-shift exceeds limit (1/0=Yes/No) ----------")')
     READ (9, '(7x,i1)') flag_trace_z_shift_ion
     READ (9, '(A77)') buf !--dddd.ddddd- Limit for the ion Z-shift (m) ---------------------------------")')
     READ (9, '(2x,f10.5)') max_Z_shift_ion
     READ (9, '(A77)') buf !======================= PARAMETERS DEFINING SCALE VALUES ====================")')
     READ (9, '(A77)') buf !--#d.dddE#dd- Plasma density (m^-3) -----------------------------------------")')
     READ (9, '(2x,e10.3)') N_plasma_m3
     READ (9, '(A77)') buf !--dddddd.ddd- Electron temperature (eV) -------------------------------------")')
     READ (9, '(2x,f10.3)') T_e_eV
     READ (9, '(A77)') buf !--dddddd----- Number of macroparticles per cell -----------------------------")')
     READ (9, '(2x,i6)') N_of_particles_cell
     READ (9, '(A77)') buf !--dddddd----- Number of cells per Debye length ------------------------------")')
     READ (9, '(2x,i6)') N_of_cells_debye
     READ (9, '(A77)') buf !------dd----- Maximal expected velocity (in V_therm_e) ----------------------")')
     READ (9, '(6x,i2)') N_max_vel
     READ (9, '(A77)') buf !--dddddd----- Number of velocity boxes per unit of V_therm ------------------")')
     READ (9, '(2x,i6)') N_box_vel
     READ (9, '(A77)') buf !============================ SIMULATION CONTROL =============================")')
     READ (9, '(A77)') buf !--dddddd.ddd- Duration of simulation (ns) -----------------------------------")')
     READ (9, '(2x,f10.3)') T_sim_ns   
     READ (9, '(A77)') buf !--dddddd----- Step for saving checkpoints (timesteps, skip if <=0) ----------")')
     READ (9, '(2x,i6)') SaveCheck_step
     READ (9, '(A77)') buf !--dddddd----- Seed for random numbers generator -----------------------------")')
     READ (9, '(2x,i6)') init_random_seed
     READ (9, '(A77)') buf !------------- Method of initialization --------------------------------------")')
     READ (9, '(A77)') buf !--            0 = start a new run, ordinary                                --")')
     READ (9, '(A77)') buf !--            1 = continue the old run, start at a checkpoint              --")')
     READ (9, '(A77)') buf !-------d----- 2 = start a new run, take particles from a checkpoint ---------")')
     READ (9, '(7x,i1)') Restore_from_checkpoint
     READ (9, '(A77)') buf !---dddddddd-- Time step when the checkpoint to be used was created ----------")')
     READ (9, '(3x,i8)') T_cntr_to_continue
           
  ELSE                                   ! if file does not exist

!======================= System configuration =====================
     BC_flag         = 0 
     R_ext_ohm       = 0.0_8
     S_electrode_cm2 = 400.0_8
     E_z_ext_Vm      = 10000.0_8
     B_x_ext_Gs      = 100.0_8
     B_y_ext_Gs      = 0.0_8
     U_ext_V         = 0.0_8      ! 
     L_plasma_m      = 0.02_8
     M_i_amu         = 131.0_8    ! Xenon
     T_i_eV          = 1.0_8
     flag_ions_magnetized = 1
     flag_ions_sense_Ez   = 1
     flag_trace_z_shift = 0
     max_Z_shift = 1.0_8
     flag_trace_z_shift_ion = 0
     max_Z_shift_ion = 1.0_8
!==================== Parameters defining scale values =============
     N_plasma_m3          = 1.0d17
     T_e_eV               = 40.0_8
     N_of_particles_cell  = 500
     N_of_cells_debye     = 16
     N_max_vel            = 4
     N_box_vel            = 20
!======================= Simulation control ========================
     T_sim_ns                = 1.0_8
     SaveCheck_step          = 200000
     init_random_seed        = 1234
     Restore_from_checkpoint = 0
     T_cntr_to_continue      = 1

     IF (Rank_of_process.EQ.0) THEN

        PRINT '(/2x,"Process ",i3," : File with the name ssc_initial.dat not found. Use the default settings ...")', Rank_of_process

        PRINT '(2x,"Process ",i3," : Create ssc_initial.dat file . . .")', Rank_of_process

        WRITE (9, '("============================ SYSTEM CONFIGURATION ===========================")')
        WRITE (9, '("-------d----- BC: Left wall potential is fixed/floating/ext.circuit {0/1/2}--")')
        WRITE (9, '(7x,i1)') BC_flag
        WRITE (9, '("--dddddd.ddd- Resistor in the external circuit (ohm) ------------------------")')
        WRITE (9, '(2x,f10.3)') R_ext_ohm
        WRITE (9, '("--dddddd.ddd- Electrode area (cm^2) -----------------------------------------")')
        WRITE (9, '(2x,f10.3)') S_electrode_cm2
        WRITE (9, '("--dddddd.ddd- External z-electric field parallel to walls (V/m) -------------")')
        WRITE (9, '(2x,f10.3)') E_z_ext_Vm
        WRITE (9, '("--dddddd.ddd- External x-magnetic field normal to walls (Gauss) -------------")')
        WRITE (9, '(2x,f10.3)') B_x_ext_Gs
        WRITE (9, '("--dddddd.ddd- External y-magnetic field parallel to walls (Gauss) -----------")')
        WRITE (9, '(2x,f10.3)') B_y_ext_Gs
        WRITE (9, '("--dddddd.ddd- Left wall potential/battery voltage (V) -----------------------")')
        WRITE (9, '(2x,f10.3)') U_ext_V
        WRITE (9, '("--dddd.ddddd- Plasma layer width (m, beam wavelengths if negative) ----------")')
        WRITE (9, '(2x,f10.5)') L_plasma_m
        WRITE (9, '("--dddddd.ddd- Ion mass (a.m.u., m_e if negative) ----------------------------")')
        WRITE (9, '(2x,f10.3)') M_i_amu
        WRITE (9, '("--dddddd.ddd- Ion temperature (eV) ------------------------------------------")')
        WRITE (9, '(2x,f10.3)') T_i_eV
        WRITE (9, '("-------d----- Ions sense the external MAGNETIC field (1/0=Yes/No) -----------")')
        WRITE (9, '(7x,i1)') flag_ions_magnetized
        WRITE (9, '("-------d----- Ions sense the external ELECTRIC field (1/0=Yes/No) -----------")')
        WRITE (9, '(7x,i1)') flag_ions_sense_Ez
        WRITE (9, '("-------d----- Thermalize electron if Z-shift exceeds limit (1/0=Yes/No) -----")')
        WRITE (9, '(7x,i1)') flag_trace_z_shift
        WRITE (9, '("--dddd.ddddd- Limit for the electron Z-shift (m) ----------------------------")')
        WRITE (9, '(2x,f10.5)') max_Z_shift
        WRITE (9, '("-------d----- Thermalize ion if Z-shift exceeds limit (1/0=Yes/No) ----------")')
        WRITE (9, '(7x,i1)') flag_trace_z_shift_ion
        WRITE (9, '("--dddd.ddddd- Limit for the ion Z-shift (m) ---------------------------------")')
        WRITE (9, '(2x,f10.5)') max_Z_shift_ion
        WRITE (9, '("======================= PARAMETERS DEFINING SCALE VALUES ====================")')
        WRITE (9, '("--#d.dddE#dd- Plasma density (m^-3) -----------------------------------------")')
        WRITE (9, '(2x,e10.3)') N_plasma_m3
        WRITE (9, '("--dddddd.ddd- Electron temperature (eV) -------------------------------------")')
        WRITE (9, '(2x,f10.3)') T_e_eV
        WRITE (9, '("--dddddd----- Number of macroparticles per cell -----------------------------")')
        WRITE (9, '(2x,i6)') N_of_particles_cell
        WRITE (9, '("--dddddd----- Number of cells per Debye length ------------------------------")')
        WRITE (9, '(2x,i6)') N_of_cells_debye
        WRITE (9, '("------dd----- Maximal expected velocity (in V_therm_e) ----------------------")')
        WRITE (9, '(6x,i2)') N_max_vel
        WRITE (9, '("--dddddd----- Number of velocity boxes per unit of V_therm ------------------")')
        WRITE (9, '(2x,i6)') N_box_vel
        WRITE (9, '("============================ SIMULATION CONTROL =============================")')
        WRITE (9, '("--dddddd.ddd- Duration of simulation (ns) -----------------------------------")')
        WRITE (9, '(2x,f10.3)') T_sim_ns   
        WRITE (9, '("--dddddd----- Step for saving checkpoints (timesteps, skip if <=0) ----------")')
        WRITE (9, '(2x,i6)') SaveCheck_step
        WRITE (9, '("--dddddd----- Seed for random numbers generator -----------------------------")')
        WRITE (9, '(2x,i6)') init_random_seed
        WRITE (9, '("------------- Method of initialization --------------------------------------")')
        WRITE (9, '("--            0 = start a new run, ordinary                                --")')
        WRITE (9, '("--            1 = continue the old run, start at the last checkpoint       --")')
        WRITE (9, '("-------d----- 2 = start a new run, take particles from a checkpoint ---------")')
        WRITE (9, '(7x,i1)') Restore_from_checkpoint
        WRITE (9, '("---dddddddd-- Time step when the checkpoint to be used was created ----------")')
        WRITE (9, '(3x,i8)') T_cntr_to_continue

     END IF

  END IF

  CLOSE (9, STATUS = 'KEEP')

! fool-proof
  IF ((BC_flag.EQ.2).AND.(R_ext_ohm.LE.0.0_8)) THEN
     IF (Rank_of_process.EQ.0) THEN
        PRINT '("INITIATE_PARAMETERS: ERROR: external contour boundary condition requires positive resistance")'
        PRINT '("the present value of R_ext_ohm is ",e12.5," ohm")', R_ext_ohm
        PRINT '("terminating the program")'
     END IF
     STOP
  END IF

  IF ((BC_flag.EQ.2).AND.(S_electrode_cm2.LE.0.0_8)) THEN
     IF (Rank_of_process.EQ.0) THEN
        PRINT '("INITIATE_PARAMETERS: ERROR: external contour boundary condition requires positive electrode area")'
        PRINT '("the present value of S_electrode_cm2 is ",e12.5," cm^2")', R_ext_ohm
        PRINT '("terminating the program")'
     END IF
     STOP
  END IF

  IF (flag_ions_magnetized.EQ.0) THEN
     ions_magnetized = .FALSE.
  ELSE
     ions_magnetized = .TRUE.
  END IF

  IF (flag_ions_sense_Ez.EQ.0) THEN
     ions_sense_Ez = .FALSE.
  ELSE
     ions_sense_Ez = .TRUE.
  END IF

  IF (flag_trace_z_shift.EQ.0) THEN
     trace_z_shift = .FALSE.
  ELSE
     trace_z_shift = .TRUE.
  END IF

  IF (flag_trace_z_shift_ion.EQ.0) THEN
     trace_z_shift_ion = .FALSE.
  ELSE
     trace_z_shift_ion = .TRUE.
  END IF

! fool proof
  IF (Rank_of_process.EQ.0) THEN
     IF ((trace_z_shift_ion).AND.(.NOT.ions_sense_Ez)) THEN
        PRINT '("ERROR :: Incompatible options in ssc_initial.dat :: tracing ion z shift requires ions sensing the Ez field")'
        STOP
     END IF
  END IF
  IF ((trace_z_shift_ion).AND.(.NOT.ions_sense_Ez)) STOP

  INQUIRE (FILE = 'ssc_anisotropy.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  OPEN (9, FILE = 'ssc_anisotropy.dat')

  IF (exists) THEN

     IF (Rank_of_process.EQ.0) THEN
        PRINT '(2x,"Process ",i3," : ssc_anisotropy.dat is found. Reading the data file...")', Rank_of_process
     END IF

     READ (9, '(A77)') buf !============== Parameters of anisotropic electron distribution ==============")')
     READ (9, '(A77)') buf !--dddddd.ddd- X-temperature (eV) --------------------------------------------")')
     READ (9, '(2x,f10.3)') Tx_e_eV
     READ (9, '(A77)') buf !--dddddd.ddd- Z-temperature (eV) --------------------------------------------")')
     READ (9, '(2x,f10.3)') Tz_e_eV
     READ (9, '(A77)') buf !------dd----- Maximal velocity for initial distribution (in V_therm_e) ------")')
     READ (9, '(6x,i2)') N_max_vel_distrib
     READ (9, '(A77)') buf !--#d.dddE#dd- Plasma density (m^-3), does not affect scaling ----------------")')
     READ (9, '(2x,e10.3)') N_distrib_m3
     READ (9, '(A77)') buf !-------d----- Account for ExB drift due to external fields (1/0=Yes/No) -----")')
     READ (9, '(7x,i1)') flag_use_ExB_for_electrons
     READ (9, '(A77)') buf !=============== Parameters of ion distribution ==============================")')
     READ (9, '(A77)') buf !-------d----- Account for ExB drift due to external fields (1/0=Yes/No) -----")')
     READ (9, '(7x,i1)') flag_use_ExB_for_ions
           
  ELSE                                   ! if file does not exist

     Tx_e_eV           = T_e_eV 
     Tz_e_eV           = T_e_eV
     N_max_vel_distrib = 3
     N_distrib_m3      = N_plasma_m3
     flag_use_ExB_for_electrons = 1
     flag_use_ExB_for_ions = 1

     IF (Rank_of_process.EQ.0) THEN

        PRINT '(/2x,"Process ",i3," : File with the name ssc_anisotropy.dat not found. Use the default settings ...")', Rank_of_process

        PRINT '(2x,"Process ",i3," : Create ssc_anisotropy.dat file . . .")', Rank_of_process

        WRITE (9, '("============== Parameters of anisotropic electron distribution ==============")')
        WRITE (9, '("--dddddd.ddd- X-temperature (eV) --------------------------------------------")')
        WRITE (9, '(2x,f10.3)') Tx_e_eV
        WRITE (9, '("--dddddd.ddd- Z-temperature (eV) --------------------------------------------")')
        WRITE (9, '(2x,f10.3)') Tz_e_eV
        WRITE (9, '("------dd----- Maximal velocity for initial distribution (in V_therm_e) ------")')
        WRITE (9, '(6x,i2)') N_max_vel_distrib
        WRITE (9, '("--#d.dddE#dd- Plasma density (m^-3), does not affect scaling ----------------")')
        WRITE (9, '(2x,e10.3)') N_distrib_m3
        WRITE (9, '("-------d----- Account for ExB drift due to external fields (1/0=Yes/No) -----")')
        WRITE (9, '(7x,i1)') flag_use_ExB_for_electrons
        WRITE (9, '("=============== Parameters of ion distribution ==============================")')
        WRITE (9, '("-------d----- Account for ExB drift due to external fields (1/0=Yes/No) -----")')
        WRITE (9, '(7x,i1)') flag_use_ExB_for_ions

     END IF

  END IF

  CLOSE (9, STATUS = 'KEEP')

  IF (flag_use_ExB_for_electrons.EQ.0) THEN
     use_ExB_for_electrons = .FALSE.
  ELSE
     use_ExB_for_electrons = .TRUE.
  END IF

  IF (flag_use_ExB_for_ions.EQ.0) THEN
     use_ExB_for_ions = .FALSE.
  ELSE
     use_ExB_for_ions = .TRUE.
  END IF

! no ExB drift for ions if the magnetic field is omitted
  IF (.NOT.ions_magnetized) use_ExB_for_ions = .FALSE.

  v_Te_ms     = SQRT(2.0_8 * T_e_eV * e_Cl / m_e_kg)
  W_plasma_s1 = SQRT(N_plasma_m3 * e_Cl**2 / (eps_0_Fm * m_e_kg))
  W_cycl_x_s1 = e_Cl * B_x_ext_Gs * 0.0001_8 / m_e_kg
  W_cycl_y_s1 = e_Cl * B_y_ext_Gs * 0.0001_8 / m_e_kg
  L_debye_m   = v_Te_ms / W_plasma_s1
  delta_x_m   = L_debye_m / N_of_cells_debye    
  delta_t_s   = delta_x_m / (N_max_vel * v_Te_ms)

! make the maximal Z-shift dimensionless
  max_Z_shift = max_Z_shift / delta_x_m
  max_Z_shift_ion = max_Z_shift_ion / delta_x_m
  
  IF ((R_ext_ohm.GT.0.0_8).AND.(S_electrode_cm2.GT.0.0_8)) THEN
     factor_SR = (delta_x_m * delta_t_s) / (1.0d-4 * S_electrode_cm2 * R_ext_ohm * eps_0_Fm)
  ELSE
     factor_SR = 0.0_8
  END IF

  E_scl_Vm    = m_e_kg * v_Te_ms * W_plasma_s1 / e_Cl 
  U_scl_V     = E_scl_Vm * delta_x_m
  N_scl_m3    = N_plasma_m3 / DBLE(N_of_particles_cell)  
  V_scl_ms    = v_Te_ms / DBLE(N_box_vel)

  U_ext = U_ext_V / U_scl_V

  Max_T_cntr  = 1.0d-9 * T_sim_ns / delta_t_s                               ! Number of time steps

! if the plasma length is negative then it is given in beam wavelengths, assumed to be integer !!!
  IF (L_plasma_m.LT.0.0_8) THEN  
     IF (Rank_of_process.EQ.0) PRINT '(/2x,"The plasma length was given in units of electron beam wavelength")'
! to determine the length of the system as integer number of the wavelengths we will need the beam energy. Read it from the beam file.
     INQUIRE (FILE = 'ssc_ebeam.dat', EXIST = exists)
     IF(exists) THEN
        OPEN (9, FILE = 'ssc_ebeam.dat')
        READ (9, '(A77)') buf ! ********************** PERIODICITY OF PARTICLES MOTION **********************")')
        READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
        READ (9, '(A77)') buf ! -------d----- 1 = turned on -------------------------------------------------")')
        READ (9, '(7x,i1)') dummy_i
        READ (9, '(A77)') buf ! ****************************** BEAM IN PLASMA *******************************")')
        READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
        READ (9, '(A77)') buf ! ---           1 = monoenergetic beam                                      ---")')
        READ (9, '(A77)') buf ! -------d----- 2 = warm beam -------------------------------------------------")')
        READ (9, '(7x,i1)') dummy_i
        READ (9, '(A77)') buf ! --dddddd.ddd- Beam appears at this moment, [ns] (timesteps if < 0) ----------")')
        READ (9, '(2x,f10.3)') dummy_r
        READ (9, '(A77)') buf ! ============================= BEAM PARAMETERS ===============================")')
        READ (9, '(A77)') buf ! --dddddd.ddd- Energy, [eV] --------------------------------------------------")')
        READ (9, '(2x,f10.3)') Energy_ebeam_eV
        CLOSE (9, STATUS = 'keep')   !## this line was added on 2008-10-11
     ELSE
        IF  (Rank_of_process.EQ.0) THEN
           PRINT '(/2x,"ERROR: File with the name ssc_ebeam.dat not found ...")'
           PRINT  '(2x,"The system lengths was given in the wavelengthes, but the beam energy is unknown !")'
           PRINT  '(2x,"Program will be terminated now :(")'
        END IF
        STOP
     END IF
!     L_plasma_deb = - L_plasma_m
     L_plasma_deb = - L_plasma_m * 6.2831853_8 * SQRT(Energy_ebeam_eV / T_e_eV)  ! ####### ####### was like in the previous line ##### #####        
                                                ! ####### ####### now we assume that "-L_plasma_m" is the number of wavelengths ##### #####
  ELSE
     L_plasma_deb = L_plasma_m / L_debye_m
  END IF

!  N_cells  = L_plasma_deb * N_of_cells_debye      ! Number of cells is calculated    
  N_cells  = L_plasma_deb * DBLE(N_of_cells_debye)      ! Number of cells is calculated      ! ####### ####### was like in the previous line ##### #####
  N_nodes  = N_cells + 1                          ! Number of nodes is calculated

  L_plasma_m = N_cells * delta_x_m

  IF (M_i_amu.LT.0.0_8) THEN
     IF (Rank_of_process.EQ.0) THEN 
        PRINT '(/2x,"Ion mass is ",f10.3," electron masses or ",e12.5," kg")', Rank_of_process, ABS(M_i_amu), ABS(M_i_amu) * m_e_kg / amu_kg
     END IF
     M_i_amu = ABS(M_i_amu) * m_e_kg / amu_kg
  END IF

  CALL SETVALUES_SPECIES_ARRAYS

! set the drift velocity for the initial distribution
  IF ((E_z_ext_Vm.EQ.0.0_8).OR.(B_y_ext_Gs.EQ.0.0_8)) THEN
     Vx_drift = 0.0_8
  ELSE
     Vx_drift = -N_box_vel * (10000.0_8 * E_z_ext_Vm * B_y_ext_Gs / (B_x_ext_Gs**2 + B_y_ext_Gs**2) ) / v_Te_ms
  END IF
  IF ((E_z_ext_Vm.EQ.0.0_8).OR.(B_x_ext_Gs.EQ.0.0_8)) THEN
     Vy_drift = 0.0_8
  ELSE
     Vy_drift = N_box_vel * (10000.0_8 * E_z_ext_Vm * B_x_ext_Gs / (B_x_ext_Gs**2 + B_y_ext_Gs**2) ) / v_Te_ms
  END IF

!###Vx_drift = -N_box_vel * (1.0d6) / v_Te_ms   !############## introduced to simulated Buneman-type instability
!###Vy_drift = 0.0_8                            !############## due to reative flow of electrons and ions
!###                                            !############## is expected to be used with zero external fields
!###                                            !############## must be commented out / removed for ordinary simulations

  IF (use_ExB_for_electrons) THEN
     Vx_e_drift = Vx_drift
     Vy_e_drift = Vy_drift
  ELSE
     Vx_e_drift = 0.0_8
     Vy_e_drift = 0.0_8
  END IF
     
  IF (use_ExB_for_ions) THEN
     Vx_i_drift = Vx_drift
     Vy_i_drift = Vy_drift
  ELSE
     Vx_i_drift = 0.0_8
     Vy_i_drift = 0.0_8
  END IF

  IF (Rank_of_process.EQ.0) THEN 
     PRINT '(/2x,"Initial electron drift X-velocity is ",e12.5," m/s")', Vx_e_drift * v_Te_ms / N_box_vel
     PRINT  '(2x,"Initial electron drift Y-velocity is ",e12.5," m/s")', Vy_e_drift * v_Te_ms / N_box_vel
     PRINT  '(2x,"Initial      ion drift X-velocity is ",e12.5," m/s")', Vx_i_drift * v_Te_ms / N_box_vel
     PRINT  '(2x,"Initial      ion drift Y-velocity is ",e12.5," m/s")', Vy_i_drift * v_Te_ms / N_box_vel
  END IF

! set velocity of ion drift in the Z-direction at z=0, to be used if thermalization of ions that travelled too far along z is requested
  IF (N_spec.GT.1) VZi_0 = N_box_vel / SQRT(2.0_8 * Ms(2))

! REAL(8) KVx          ! Coefficients for equations of motion (X-PRE-FIN-moving)                         
  KVx = 1.0_8 / (DBLE(N_max_vel) * DBLE(N_box_vel))

! REAL(8) K_Q          ! Coefficient used in Poisson's equation
  K_Q = 1.0_8 / (DBLE(N_of_cells_debye) * DBLE(N_of_particles_cell))

  CALL CONFIGURE_MESH_ARRAYS
  CALL SETZEROS_MESH_ARRAYS

  CALL PrepareMaxwellDistribIntegral          ! now we can use GetMaxwellVelocity once the random number generator is ready

! initialize particles
  IF (Restore_from_checkpoint.EQ.0) THEN    ! for ordinary initialization

! initialize random numbers generators
     IF (Rank_of_process.EQ.0) THEN 
        PRINT '(/2x,"Process ",i3," : Seed for well_random_seed: ",i12)', Rank_of_process, init_random_seed

        CALL well_random_seed(init_random_seed)

        DO j = 1, 1000000
           myran = well_random_number()
        END DO

        DO i = 1, N_of_processes - 1
           itmp = INT(2000000000.0_8 * well_random_number())                       
           CALL MPI_SEND(itmp, 1, MPI_INTEGER, i, 101, MPI_COMM_WORLD, ierr)
        END DO
     ELSE
        CALL MPI_RECV(my_random_seed, 1, MPI_INTEGER, 0, 101, MPI_COMM_WORLD, stattus, ierr)

        PRINT '(/2x,"Process ",i3," : Seed for well_random_seed: ",i12)', Rank_of_process, my_random_seed

        CALL well_random_seed(my_random_seed)

        DO j = 1, 1000000
           myran = well_random_number()
        END DO

     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL DETERMINE_NUMBER_OF_PARTICLES
!  N_of_part_initial = N_of_particles_cell * N_cells   

! determine approximate number of particles for a client
     do s = 1, N_spec
        N_part(s) = N_of_part_initial / (N_of_processes-1)
     end do

! Length(1:N_spec)     Lengthes of allocated arrays in blocks (they are not allocated yet!!!)
     DO s = 1, N_spec
        Length(s) = N_part(s) + INT(0.01 * REAL(N_part(s)))
     END DO

     if (Rank_of_process.eq.0) then
        N_part = 0
        Length = 0
     end if

     CALL CONFIGURE_PARTICLE_DYNAM_ARRAYS
     CALL INITIATE_PARTICLE_DYNAM_ARRAYS_LOCAL

     Q_left  = 0
     Q_right = 0
     Start_T_cntr  = 0

     full_Q_left  = 0.0_8
     full_Q_right = 0.0_8
    
     prev_Q_left  = 0.0_8   !
     prev_Q_right = 0.0_8   !
     
  ELSE ! if we must read data from checkpoint files

     CALL READ_CHECKPOINT_MPIIO

     IF (Restore_from_checkpoint.EQ.2) Start_T_cntr = 0 

  END IF

  IF (Rank_of_process.EQ.0) THEN 
     PRINT '(/2x,"Electron Langmuir frequency                      : ",e10.3, " s^-1")', W_plasma_s1
     PRINT  '(2x,"     Ion Langmuir frequency                      : ",e10.3, " s^-1")', W_plasma_s1 * SQRT(m_e_kg / (M_i_amu * amu_kg))
     PRINT  '(2x,"Electron cyclotron frequency                     : ",e10.3, " s^-1")', SQRT(W_cycl_x_s1**2 + W_cycl_y_s1**2)
     PRINT  '(2x,"     Ion cyclotron frequency                     : ",e10.3, " s^-1")', SQRT(W_cycl_x_s1**2 + W_cycl_y_s1**2) * m_e_kg / (M_i_amu * amu_kg)
     PRINT  '(2x,"The electron thermal velocity                    : ",e10.3, " m/s")',  v_Te_ms
     PRINT  '(2x,"The electron Debye length                        : ",e10.3, " m")',    L_debye_m 

     PRINT '(/2x,"After mesh distribution the plasma length is     : ",f9.3," cm")', L_plasma_m * 100.0_8
     PRINT  '(2x,"The new plasma length is                         : ",i6," (debye len.)")', (N_cells / N_of_cells_debye)
     PRINT  '(2x,"Total number of plasma nodes is                  : ",i9)', N_nodes
     PRINT  '(2x,"The cell size is                                 : ",f9.3," mkm")', delta_x_m * 1.0e6
     IF (Restore_from_checkpoint.EQ.0) THEN
        PRINT  '(2x,"Initial total number of particles of one species : ",i9)', N_of_part_initial
     ELSE
        PRINT  '(2x,"Initial total number of particles (e)            : ",i9)', N_part(1)
        PRINT  '(2x,"Initial total number of particles (i)            : ",i9)', N_part(2)        
        PRINT  '(2x,"Charge assigned to the left wall                 : ",i9," charges of one particle")', Q_left
        PRINT  '(2x,"Charge assigned to the right wall                : ",i9," charges of one particle")', Q_right
     END IF
     PRINT  '(2x,"Number of time steps                             : ",i9)', Max_T_cntr
     PRINT  '(2x,"The time step is                                 : ",f9.3," ps")', delta_t_s * 1.0e12
     PRINT '(/2x,"Ratio of electron plasma period to the time step : ",f9.2)', 6.28318530718_8 / (W_plasma_s1 * delta_t_s)

     IF ((W_cycl_x_s1.NE.0.0_8).OR.(W_cycl_y_s1.NE.0.0_8)) THEN
        PRINT  '(2x,"Ratio of electron gyro period to the time step   : ",f9.2)', 6.28318530718_8 / (SQRT(W_cycl_x_s1**2 + W_cycl_y_s1**2) * delta_t_s)
     END IF
  END IF

END SUBROUTINE INITIATE_PARAMETERS 
 
!===================================================================================================
! Note, this subroutine is called in the INITIATE_PARAMETERS  
SUBROUTINE SETVALUES_SPECIES_ARRAYS 

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER s       ! species
  REAL(8) alfa_x, alfa_y, alfa_x2, alfa_y2, theta2

! particle parameters -------------------------
  Ms  = 0.0_8
  Qs  = 0
  QMs = 0.0_8
  VT  = 0.0_8

  Ms(1)  =  1.0_8
  Qs(1)  = -1
  QMs(1) = -1.0_8
  VT(1)  =  N_box_vel
     
  DO s = 2, N_spec
     Ms(s)  = M_i_amu * amu_kg / m_e_kg
     Qs(s)  = 1          
     QMs(s) = Qs(s) / Ms(s)
     VT(s)  = SQRT( T_i_eV / (Ms(s) * T_e_eV) ) * N_box_vel
  END DO

! precalculated coefficients -----------------

  DO s = 1, N_spec

     alfa_x = 0.5_8 * QMs(s) * W_cycl_x_s1 / (W_plasma_s1 * DBLE(N_max_vel) * DBLE(N_of_cells_debye))
     alfa_y = 0.5_8 * QMs(s) * W_cycl_y_s1 / (W_plasma_s1 * DBLE(N_max_vel) * DBLE(N_of_cells_debye))
     alfa_x2 = alfa_x**2
     alfa_y2 = alfa_y**2
     theta2 = alfa_x2 + alfa_y2

     K_Xi(s) = (0.25_8 * QMs(s) * Qs(s) / (DBLE(N_of_particles_cell) * (DBLE(N_of_cells_debye) * DBLE(N_max_vel))**2)) * (1.0_8 + alfa_x2) / (1.0_8 + theta2) 

     K11(s) = (1.0_8 - theta2 + 2.0_8 * alfa_x2) / (1.0_8 + theta2)

     K12(s) =            2.0_8 * alfa_x * alfa_y / (1.0_8 + theta2)

     K13(s) =                     2.0_8 * alfa_y / (1.0_8 + theta2)

     K22(s) = (1.0_8 - theta2 + 2.0_8 * alfa_y2) / (1.0_8 + theta2)

     K23(s) =                     2.0_8 * alfa_x / (1.0_8 + theta2)

     K33(s) =                   (1.0_8 - theta2) / (1.0_8 + theta2)

     A11(s) =  0.5_8 * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))) * ((1.0_8 + alfa_x2) / (1.0_8 + theta2))

     A21(s) =  0.5_8 * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))) *   (alfa_x * alfa_y / (1.0_8 + theta2))

     A31(s) =  0.5_8 * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))) *            (alfa_y / (1.0_8 + theta2))
     
     A13(s) = QMs(s) * (E_z_ext_Vm / E_scl_Vm) * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))) * (alfa_y / (1.0_8 + theta2))

     A23(s) = QMs(s) * (E_z_ext_Vm / E_scl_Vm) * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))) * (alfa_x / (1.0_8 + theta2))
     
     A33(s) = QMs(s) * (E_z_ext_Vm / E_scl_Vm) * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye)))           / (1.0_8 + theta2)

     KvEx(s) = ((1.0_8 + alfa_x2) / (1.0_8 + theta2)) * 0.5_8 * QMs(s) * DBLE(N_box_vel) / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))

     KvEy(s) =   (alfa_x * alfa_y / (1.0_8 + theta2)) * 0.5_8 * QMs(s) * DBLE(N_box_vel) / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))

     KvEz(s) =            (alfa_y / (1.0_8 + theta2)) * 0.5_8 * QMs(s) * DBLE(N_box_vel) / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))

     KxEx(s) = KvEx(s) / (DBLE(N_max_vel) * DBLE(N_box_vel))    !0.5_8 * QMs(s) / (DBLE(N_of_cells_debye) * (DBLE(N_max_vel))**2) 

  END DO

  IF (N_spec.EQ.2) THEN

     s = 2

     IF (.NOT.ions_magnetized) THEN

        alfa_x = 0.0_8 ! 0.5_8 * QMs(s) * W_cycl_x_s1 / (W_plasma_s1 * DBLE(N_max_vel) * DBLE(N_of_cells_debye))
        alfa_y = 0.0_8 ! 0.5_8 * QMs(s) * W_cycl_y_s1 / (W_plasma_s1 * DBLE(N_max_vel) * DBLE(N_of_cells_debye))
        alfa_x2 = 0.0_8 !alfa_x**2
        alfa_y2 = 0.0_8 !alfa_y**2
        theta2 = 0.0_8 !alfa_x2 + alfa_y2

        K_Xi(s) = (0.25_8 * QMs(s) * Qs(s) / (DBLE(N_of_particles_cell) * (DBLE(N_of_cells_debye) * DBLE(N_max_vel))**2)) !* (1.0_8 + alfa_x2) / (1.0_8 + theta2) 

        K11(s) = 1.0_8 !(1.0_8 - theta2 + 2.0_8 * alfa_x2) / (1.0_8 + theta2)

        K12(s) = 0.0_8 !           2.0_8 * alfa_x * alfa_y / (1.0_8 + theta2)

        K13(s) = 0.0_8 !                    2.0_8 * alfa_y / (1.0_8 + theta2)

        K22(s) = 1.0_8 !(1.0_8 - theta2 + 2.0_8 * alfa_y2) / (1.0_8 + theta2)

        K23(s) = 0.0_8 !                    2.0_8 * alfa_x / (1.0_8 + theta2)

        K33(s) = 1.0_8 !                  (1.0_8 - theta2) / (1.0_8 + theta2)

        A11(s) =  0.5_8 * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))) !* ((1.0_8 + alfa_x2) / (1.0_8 + theta2))

        A21(s) =  0.0_8 !0.5_8 * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))) *   (alfa_x * alfa_y / (1.0_8 + theta2))

        A31(s) =  0.0_8 !0.5_8 * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))) *            (alfa_y / (1.0_8 + theta2))
     
        A13(s) = 0.0_8 !QMs(s) * (E_z_ext_Vm / E_scl_Vm) * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))) * (alfa_y / (1.0_8 + theta2))

        A23(s) = 0.0_8 !QMs(s) * (E_z_ext_Vm / E_scl_Vm) * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))) * (alfa_x / (1.0_8 + theta2))
     
        A33(s) = QMs(s) * (E_z_ext_Vm / E_scl_Vm) * (N_box_vel / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))) !           / (1.0_8 + theta2)

        IF (.NOT.ions_sense_Ez) A33(s) = 0.0_8  ! this is the only place where ions_sense_Ez does something

        KvEx(s) = 0.5_8 * QMs(s) * DBLE(N_box_vel) / (DBLE(N_max_vel) * DBLE(N_of_cells_debye)) ! ((1.0_8 + alfa_x2) / (1.0_8 + theta2)) * 

        KvEy(s) = 0.0_8 !  (alfa_x * alfa_y / (1.0_8 + theta2)) * 0.5_8 * QMs(s) * DBLE(N_box_vel) / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))

        KvEz(s) = 0.0_8 !           (alfa_y / (1.0_8 + theta2)) * 0.5_8 * QMs(s) * DBLE(N_box_vel) / (DBLE(N_max_vel) * DBLE(N_of_cells_debye))

        KxEx(s) = KvEx(s) / (DBLE(N_max_vel) * DBLE(N_box_vel))    !0.5_8 * QMs(s) / (DBLE(N_of_cells_debye) * (DBLE(N_max_vel))**2) 

     END IF

  END IF
  
!  N_inject = 0  ! Number (counter) of injected particles of one species (1 = e, 2 = i) due to ionization or SEE

END SUBROUTINE SETVALUES_SPECIES_ARRAYS

!===================================================================================================
! Note, this subroutine is called in the INITIATE_PARAMETERS  
! Note, the new numbering assumes that node # 0       is the left plasma boundary  x = 0
!                                  and node # N_cells is the right plasma boundart x = L_plasma
SUBROUTINE CONFIGURE_MESH_ARRAYS 

  USE CurrentProblemValues
  IMPLICIT NONE

  INTEGER ALLOC_ERR

  ALLOCATE(EX(0:N_cells), STAT=ALLOC_ERR)
  ALLOCATE(GradEX(0:(N_cells-1)), STAT=ALLOC_ERR)
  ALLOCATE(F(0:N_cells), STAT=ALLOC_ERR)  
  ALLOCATE(Xi(0:N_cells), STAT=ALLOC_ERR)
  ALLOCATE(Q_stream(0:N_cells), STAT=ALLOC_ERR)
  ALLOCATE( Q_strm_spec(0:N_cells,1:N_spec), STAT=ALLOC_ERR)

END SUBROUTINE CONFIGURE_MESH_ARRAYS 

!===================================================================================================
!  Note, this subroutine is called in the INITIATE_PARAMETERS  
SUBROUTINE SETZEROS_MESH_ARRAYS

  USE CurrentProblemValues
  IMPLICIT NONE

  EX          = 0.0_8
  GradEX      = 0.0_8
  F           = 0.0_8
  Xi          = 0.0_8
  Q_stream    = 0.0_8
  Q_strm_spec = 0.0_8

END SUBROUTINE SETZEROS_MESH_ARRAYS

!=============================================================================================
! This subroutine calculates the number of particles required to 
! reproduce the desired density profile.
! 
SUBROUTINE DETERMINE_NUMBER_OF_PARTICLES

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE

  INTEGER(8) n_int
  REAL(8) x_left
  REAL(8) x_right
  REAL(8) dx

  REAL(8) integral_of_n
  INTEGER(8) i

! function
  REAL(8) density

  n_int   = N_of_particles_cell * N_cells * 20
  x_left  = 0.5_8 / N_of_particles_cell !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! step aside from the walls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  x_right = N_cells - x_left
  dx = (x_right - x_left) / n_int     

! Calculate the integral of the function "density(z)", which describes the density profile.
  integral_of_n = 0.0_8
  DO i = 1_8, n_int 
     integral_of_n = integral_of_n + density(x_left + i * dx - 0.5_8 * dx)
  END DO

  N_of_part_initial = N_of_particles_cell * (integral_of_n * dx)

END SUBROUTINE DETERMINE_NUMBER_OF_PARTICLES

!===================================================================================================
SUBROUTINE CONFIGURE_PARTICLE_DYNAM_ARRAYS 

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER ALLOC_ERR
  INTEGER s 

  if (Rank_of_process.eq.0) return

  ALLOCATE(species(1:N_spec), STAT=ALLOC_ERR)
  DO s = 1, N_spec
     ALLOCATE(species(s)%part(1:Length(s)), STAT=ALLOC_ERR)
  END DO

END SUBROUTINE CONFIGURE_PARTICLE_DYNAM_ARRAYS 

!===================================================================================================
SUBROUTINE INITIATE_PARTICLE_DYNAM_ARRAYS_LOCAL 

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER i, s

  REAL(8) factor_x, factor_z

! function
  REAL(8) ion_velocity_z

  PRINT '(4x,"Process ",i3," : Initiating particle dynamic arrays with estimated N_part ", i8,2x,i8)', Rank_of_process, N_part(1), N_part(2)

  IF (Rank_of_process.EQ.0) RETURN

  CALL INITIAL_DISTRIB_IN_SPACE_LOCAL   ! N_part may be modified here

  CALL INITIAL_DISTRIB_IN_VELOCITY_LOCAL

  species(1)%part(1:N_part(1))%Z = 0.0_8

  IF (N_spec.EQ.2) THEN
     IF (trace_z_shift_ion) THEN
        CALL INITIAL_ION_DISTRIBUTION_ALONG_Z_LOCAL
     ELSE
        species(2)%part(1:N_part(2))%Z = 0.0_8
     END IF
  END IF

! introduce anisotropy and drift for the electrons
  factor_x = SQRT(Tx_e_eV / T_e_eV)
  factor_z = SQRT(Tz_e_eV / T_e_eV)
  DO i = 1, N_part(1)
     species(1)%part(i)%VX = species(1)%part(i)%VX * factor_x + Vx_e_drift
     species(1)%part(i)%VY = species(1)%part(i)%VY * factor_z + Vy_e_drift
     species(1)%part(i)%VZ = species(1)%part(i)%VZ * factor_z
  END DO

! introduce drift for the ions
  IF (N_spec.EQ.2) THEN
     DO i = 1, N_part(2)
        species(2)%part(i)%VX = species(2)%part(i)%VX + Vx_i_drift
        species(2)%part(i)%VY = species(2)%part(i)%VY + Vy_i_drift
        species(2)%part(i)%VZ = species(2)%part(i)%VZ + ion_velocity_z(species(2)%part(i)%Z)  ! ion_velocity_z=0 if .NOT.trace_z_shift_ion
     END DO
  END IF

  DO s = 1, N_spec
     DO i = 1, N_part(s)
        species(s)%part(i)%AX = 0.0_8
        species(s)%part(i)%tag = 0
     END DO
  END DO

  PRINT '(2x,"Finished :)")'

END SUBROUTINE INITIATE_PARTICLE_DYNAM_ARRAYS_LOCAL

!=============================================================================================
! This subroutine distributed the electrons and ions in space according to the density profile.
! 
SUBROUTINE INITIAL_DISTRIB_IN_SPACE_LOCAL

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE

  INTEGER(8) n_int
  REAL(8) x_left
  REAL(8) x_right
  REAL(8) dx

  REAL(8) integral_of_n
  INTEGER(8) i 

  REAL(8) normalization_factor

  INTEGER(8) n_global
  INTEGER n_local

  INTEGER(8) i_prev     ! These integer variables are used to trace the following situation: 
  INTEGER(8) i_new      ! a number ABC{D}.**** increases to ABC{D+1}.**** (ABCD - decimal digits) 

  INTEGER(8) n_to_use

  REAL(8) integral_of_n_normalized

  INTEGER ii

! function
  REAL(8) density

  PRINT '(6x,"Process ",i3," : Performing initial distribution of particles in space ...")', Rank_of_process 
  
  n_int   = N_of_part_initial * 20
  x_left  = 0.5_8 / N_of_particles_cell !+ 20.0_8 !0.0_8  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! step aside from the walls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  x_right = N_cells - x_left
  dx = (x_right - x_left) / n_int   

  integral_of_n = 0.0_8
  DO i = 1_8, n_int 
     integral_of_n = integral_of_n + density(x_left + i * dx - 0.5_8 * dx)
  END DO

! Normalize the values of the integral
  normalization_factor = DBLE(N_of_part_initial - 1_8) / integral_of_n

! Note, particle #1 is distributed to x_left, therefore integral distributes (N_of_part_initial - 1) particles

! Distribute the particles
  n_global = 1_8                       ! particles counter
  n_local=0
  i_prev = 0_8
  integral_of_n = 0.0_8

  n_to_use = Rank_of_process

  if (n_global.eq.n_to_use) then
     n_local = n_local+1
     species(1)%part(n_local)%X = x_left
     n_to_use = n_to_use + N_of_processes - 1
  end if

  DO i = 1_8, n_int-1_8   

     integral_of_n = integral_of_n + density(x_left + i * dx - 0.5_8 * dx)

     integral_of_n_normalized = integral_of_n * normalization_factor

     i_new = integral_of_n_normalized
! If the integral value increases by 1,2,3,etc. we place 1,2,3,etc. particles (electron and ion) at this point.

     do ii = 1, i_new - i_prev

        n_global = n_global + 1_8

        if (n_global.eq.n_to_use) then
           n_local = n_local+1
! check that size of array is not exceeded
           if (n_local.gt.Length(1)) then
              print '("Process ",i4," :: error in INITIAL_DISTRIB_IN_SPACE_LOCAL :: too many particles assigned ",i12," vs ",i12)', Rank_of_process, n_local, Length(1)
              stop
           end if
           species(1)%part(n_local)%X = x_left + i * dx
           n_to_use = n_to_use + N_of_processes - 1
        end if

     end do

     i_prev = i_new

  END DO

! there must be a particle in the very last point of the integration interval.
! so we just place it there and assign it to process rank 1
! thus we assume that there is no sharp density jumps near the end which may cause distribution of more than one particle at the same point

  IF (Rank_of_process.EQ.1) THEN
     n_local = n_local+1
     species(1)%part(n_local)%X = x_right
  END IF

  N_part(1) = n_local
  if (N_spec.eq.2) then
     N_part(2) = n_local
     do ii = 1, n_local
        species(2)%part(ii)%X = species(1)%part(ii)%X
     end do
     print '(/4x,"Process ",i4," : ",i12," electron particles were distributed in space")', Rank_of_process, n_local
  else
     print '(/4x,"Process ",i4," : ",i12," particles of each species were distributed in space")', Rank_of_process, n_local
  end if

  PRINT '(2x,"Finished :)")'

END SUBROUTINE INITIAL_DISTRIB_IN_SPACE_LOCAL

!=======================================================================
! This subroutine sets particle velocities according to Maxwellian distribution with proper temperatures
!
SUBROUTINE INITIAL_DISTRIB_IN_VELOCITY_LOCAL

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE

  INTEGER s, i

  s=1
  DO i = 1, N_part(1)
     CALL GetMaxwellVelocity(species(1)%part(i)%VX)
     CALL GetMaxwellVelocity(species(1)%part(i)%VY)
     CALL GetMaxwellVelocity(species(1)%part(i)%VZ)
  END DO
  
  IF (N_spec.EQ.2) THEN
     DO i = 1, N_part(2)
        CALL GetMaxwellVelocity(species(2)%part(i)%VX)
        CALL GetMaxwellVelocity(species(2)%part(i)%VY)
        CALL GetMaxwellVelocity(species(2)%part(i)%VZ)

        species(2)%part(i)%VX = species(2)%part(i)%VX * VT(2) / VT(1)
        species(2)%part(i)%VY = species(2)%part(i)%VY * VT(2) / VT(1)
        species(2)%part(i)%VZ = species(2)%part(i)%VZ * VT(2) / VT(1)

     END DO
  END IF
 
END SUBROUTINE INITIAL_DISTRIB_IN_VELOCITY_LOCAL

!==========================================================
REAL(8) FUNCTION density(x)

  USE CurrentProblemValues, ONLY : N_cells, N_distrib_m3, N_plasma_m3

  IMPLICIT NONE
  
  REAL(8) x
! We assume the constant profile of density distribution

!  density = 1.0_8

  IF (N_distrib_m3.EQ.N_plasma_m3) THEN
     density = 1.0_8
  ELSE   
     density = N_distrib_m3 / N_plasma_m3
!     IF (x.GT.(0.95_8*N_cells)) density = 0.0_8
!     IF (x.LT.(0.05_8*N_cells)) density = 0.0_8
   END IF
 
!  density = 0.01_8

!  density = 0.25_8

!  density = -6.0_8 * (x / N_cells - 0.5_8)**2 + 2.0_8


!density  = 0.0_8  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if ((x.ge.0.0_8).and.(x.le.100.0_8)) density = 1.0_8    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!if((x.lt.100.0_8).or.(x.gt.505.0_8)) density = 0.1_8
!if((x.gt.600.0_8)) density = 0.1_8

END FUNCTION density

!=============================================================================================
! This subroutine distributes IONS along the Z-direction to ensure ION current continuity
! in constant EZ electric field
! 
SUBROUTINE INITIAL_ION_DISTRIBUTION_ALONG_Z_LOCAL

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE

  INTEGER i

  IF (N_spec.NE.2) RETURN

  PRINT '(6x,"Process ",i3," : Performing initial distribution of particles along Z ...")', Rank_of_process 

  CALL PrepareIonZDistribIntegral

  DO i = 1, N_part(2)
     CALL GetIonZ(species(2)%part(i)%Z)
  END DO

  PRINT '(2x,"INITIAL_ION_DISTRIBUTION_ALONG_Z-LOCAL :: Finished")'

END SUBROUTINE INITIAL_ION_DISTRIBUTION_ALONG_Z_LOCAL

!*******************************************************************************************
! This module is used by two subroutines, calculating ion particle z coordinate
! according to the density profile corresponding to constant ion flux along z (that is n(z)*vz(z)=const)
MODULE IonDistributionOverZ

  integer, parameter :: N_of_reference_points = 10000
  REAL(8) z_ion(0:N_of_reference_points)

END MODULE IonDistributionOverZ

!-------------------------------------------------------------------------------------------
! Prepares the tabulated values of integral of the ion density profile along z when the ions are 
! accelerated by a constant electric field EZ
!  
SUBROUTINE PrepareIonZDistribIntegral

  USE ParallelOperationValues
  USE IonDistributionOverZ
  USE CurrentProblemValues, ONLY : max_Z_shift_ion
  IMPLICIT NONE

  integer, parameter :: zoom_factor = 1000
  integer, parameter :: roundoff_threshold = 5
  integer N_pnts

  REAL(8) dz
  INTEGER i
  REAL(8) integral_of_n

  REAL(8) converting_factor

  LOGICAL not_finished

  INTEGER actual_N_pnts

  REAL(8) prev_integral_of_n
  INTEGER count
  INTEGER n

! function
  REAL(8) density_z

  N_pnts = zoom_factor * N_of_reference_points

  dz = max_Z_shift_ion / N_pnts

  integral_of_n = 0.0_8
  DO i = 1, N_pnts
     integral_of_n = integral_of_n + dz * density_z((DBLE(i)-0.5_8) * dz)
  END DO

  converting_factor = N_of_reference_points / integral_of_n

! due to roundoff errors, the integral calculated below most likely will be slightly different 
! from the desired value of N_of_reference_points after N_pnts steps (that is at z = max_Z_shift_ion)
! we find the difference of integration steps between theoretical N_pnts and the actual number
! which gives the closest approximation of INT(integral) = N_of_reference_points
! if this difference is small compared to typical number of steps required to increase the integral by 1
! then we integrate only till N_of_reference_points-1 are obtained,
! the last point z(N_of_reference_points)=max_Z_shift_ion is simply imposed. 

! find how many integration steps is required now to reach the desired integral value
  integral_of_n = 0.0_8
  not_finished = .TRUE.
  DO i = 1, N_pnts+zoom_factor
     integral_of_n = integral_of_n + converting_factor * dz * density_z((DBLE(i)-0.5_8) * dz)
     IF (INT(integral_of_n).EQ.N_of_reference_points) THEN
        actual_N_pnts = i
        not_finished = .FALSE.
        EXIT
     END IF
  END DO

  IF (not_finished) THEN
     PRINT '("Error-1 in PrepareIonZDistribIntegral ",2x,e16.9,2x,i9,2x,i9)', integral_of_n, INT(integral_of_n), N_of_reference_points
     STOP
  END IF

  IF (ABS(actual_N_pnts-N_pnts).GT.roundoff_threshold) THEN
     PRINT '("Error-2 in PrepareIonZDistribIntegral ",2x,i9,2x,i9)', N_pnts, actual_N_pnts 
     STOP
  END IF

! since we are here everything looks shiny so far

! repeat integration

  integral_of_n = 0.0_8
  prev_integral_of_n = integral_of_n
  count = 0
  z_ion(0) = 0.0_8
  i = 0
  not_finished = .TRUE.

  DO WHILE (not_finished)
     i = i+1
     integral_of_n = integral_of_n + converting_factor * dz * density_z((DBLE(i)-0.5_8) * dz)

     DO n = 1, INT(integral_of_n)-INT(prev_integral_of_n)
        count = count + 1
        z_ion(count) = DBLE(i) * dz
        IF (count.EQ.(N_of_reference_points-1)) THEN
           not_finished = .FALSE.
           EXIT
        END IF
     END DO

     prev_integral_of_n = integral_of_n

     IF (i.GT.(N_pnts+zoom_factor)) THEN
        PRINT '("Error-3 in PrepareIonZDistribIntegral ",2x,i9,2x,e16.9,2x,i9,2x,i9)', i, integral_of_n, N_pnts, N_of_reference_points
        STOP
     END IF
  END DO

  z_ion(N_of_reference_points) = max_Z_shift_ion

END SUBROUTINE PrepareIonZDistribIntegral

!-------------------------------------------------------------------------------------------
!  
SUBROUTINE GetIonZ(z) 

  USE IonDistributionOverZ

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8), intent(out) :: z

  REAL(8) R
  INTEGER indx
  
  R = N_of_reference_points * well_random_number()

  indx = INT(R)

  IF (indx.LT.N_of_reference_points) THEN
     z = z_ion(indx) + (R - indx) * (z_ion(indx+1) - z_ion(indx))
  ELSE
     z = z_ion(N_of_reference_points)
  END IF
  RETURN
  
END SUBROUTINE GetIonZ

!==========================================================
REAL(8) FUNCTION density_z(z)

  IMPLICIT NONE
  
  REAL(8) z
! function
  REAL(8) ion_velocity_z

  density_z = 1.0_8 / ion_velocity_z(z)

END FUNCTION density_z

!==========================================================
!
REAL(8) FUNCTION ion_velocity_z(z)

  USE CurrentProblemValues, ONLY : trace_z_shift_ion, VZi_0, N_box_vel, E_z_ext_Vm, delta_x_m, T_e_eV, Ms, max_Z_shift_ion

  IMPLICIT NONE

  REAL(8) z
  REAL(8) VZi_max2_div_by_zmax

  IF (.NOT.trace_z_shift_ion) THEN
     ion_velocity_z = 0.0_8
     RETURN
  END IF

  IF (z.LT.0.0_8) THEN
! do it for completeness, should not be here normally
     ion_velocity_z = VZi_0
     RETURN
  END IF

!  VZi_max2_div_by_zmax = (N_box_vel / v_Te_ms)**2 * (2.0_8 * e_Cl * E_z_ext_Vm * delta_x_m) / (M_i_amu * amu_kg)        ! use max_Z_shift_ion / max_Z_shift_ion = 1

  VZi_max2_div_by_zmax = N_box_vel * N_box_vel * E_z_ext_Vm * delta_x_m / (T_e_eV * Ms(2))

!  IF (z.GT.max_Z_shift_ion) THEN
!! do it for completeness, should not be here normally
!     ion_velocity_z = SQRT(VZi_0**2 + VZi_max2_div_by_zmax * max_Z_shift_ion)
!     RETURN
!  END IF

  ion_velocity_z = SQRT(VZi_0**2 + VZi_max2_div_by_zmax * z)

END FUNCTION ion_velocity_z

!==========================================================
REAL(8) FUNCTION maxwell(v, vt)

  IMPLICIT NONE
  REAL(8) v, vt
! We assume the isothermic maxwell plasma

!maxwell = 1.0_8
!if (abs(v).gt.vt) maxwell = 0.0_8

  maxwell = DEXP(-v*v/(vt*vt))   ! 1.0_8                          !!!!!!!!&^*^%&$%&^%R(&^)*(&*(&)(^

!  maxwell = DEXP(-4.0_8 * v*v/(vt*vt))   ! 1.0_8         ! ####### ####### we use the 4 times smaller temperature for distribution

!  maxwell = DEXP(-2.0_8 * v*v/(vt*vt))   ! 1.0_8         ! ####### ####### we use the 2 times smaller temperature for distribution

! if (abs(v).gt.vt) maxwell = 0.0_8


END FUNCTION maxwell

!=============================================================================================
! preallocating the linked list for electrons
!
SUBROUTINE PREALLOCATE_ELECTRON_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER ALLOC_ERR

  NULLIFY(Inj_electron)

  IF (.NOT.ASSOCIATED(Inj_electron)) THEN
     ALLOCATE(Inj_electron, STAT=ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Inj_electron !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  Current_electron => Inj_electron
  NULLIFY(Current_electron%next)

END SUBROUTINE PREALLOCATE_ELECTRON_LIST

!---------------------------------------------------------------------------------------------
! clear the linked list for electrons
!
SUBROUTINE CLEAR_ELECTRON_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER DEALLOC_ERR

  TYPE (injected_particle), POINTER :: current, temp

  current => Inj_electron

  DO WHILE(ASSOCIATED(current))
     temp => current
     current => current%next
     DEALLOCATE(temp, STAT=DEALLOC_ERR)
     IF (DEALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Inj_electron !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END DO

END SUBROUTINE CLEAR_ELECTRON_LIST

!=============================================================================================
! preallocating the linked list for ions
!
SUBROUTINE PREALLOCATE_ION_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER ALLOC_ERR

  NULLIFY(Inj_ion)

  IF (.NOT.ASSOCIATED(Inj_ion)) THEN
     ALLOCATE(Inj_ion, STAT=ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Inj_ion !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  Current_ion => Inj_ion
  NULLIFY(Current_ion%next)

END SUBROUTINE PREALLOCATE_ION_LIST

!---------------------------------------------------------------------------------------------
! clear the linked list for ions
!
SUBROUTINE CLEAR_ION_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER DEALLOC_ERR

  TYPE (injected_particle), POINTER :: current, temp

  current => Inj_ion

  DO WHILE(ASSOCIATED(current))
     temp => current
     current => current%next
     DEALLOCATE(temp, STAT=DEALLOC_ERR)
     IF (DEALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Inj_ion !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END DO

END SUBROUTINE CLEAR_ION_LIST

!===================================================================================================
SUBROUTINE REMOVE_MESH_ARRAYS

  USE CurrentProblemValues
  IMPLICIT NONE

  INTEGER DEALLOC_ERR

  DEALLOCATE(EX, STAT=DEALLOC_ERR)
  DEALLOCATE(GradEX, STAT=DEALLOC_ERR)
  DEALLOCATE(F, STAT=DEALLOC_ERR)
  DEALLOCATE(Xi, STAT=DEALLOC_ERR)
  DEALLOCATE(Q_stream, STAT=DEALLOC_ERR)
  DEALLOCATE(Q_strm_spec, STAT=DEALLOC_ERR)

END SUBROUTINE REMOVE_MESH_ARRAYS

!===================================================================================================
SUBROUTINE REMOVE_PARTICLE_DYNAM_ARRAYS

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE

  INTEGER s
  INTEGER DEALLOC_ERR

!  PRINT '("Deleting particle dynamic arrays ...")'

  DO s = 1, N_spec
     IF (ALLOCATED(species(s)%part)) DEALLOCATE(species(s)%part, STAT = DEALLOC_ERR)
  END DO
  IF (ALLOCATED(species)) DEALLOCATE(species, STAT = DEALLOC_ERR)

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE REMOVE_PARTICLE_DYNAM_ARRAYS
