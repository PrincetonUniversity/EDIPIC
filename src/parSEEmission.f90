
!===================================================================================================

SUBROUTINE INITIATE_SE_EMISSION

  USE ParallelOperationValues
  USE SEEmission
  USE IonInducedSEEmission
  USE CurrentProblemValues, ONLY : T_e_eV, BC_flag, N_box_vel
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(8) angle_rad, angle_deg ! for testing purpose, initial_angle_of_incidence_deg expressed in radians

  LOGICAL exists
  INTEGER j !, s, k
  REAL(8) energy, tot_see_coef
  REAL(8) Coeff_SEE_Elastic
  REAL(8) Coeff_SEE_Inelastic
  REAL(8) Coeff_SEE_True

  CHARACTER (77) buf
  INTEGER ierr
 
  INQUIRE (FILE = 'ssc_secondary.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  OPEN (9, FILE = 'ssc_secondary.dat')
  IF(exists) THEN

     IF (Rank_of_process.EQ.0) THEN 
        PRINT '(/2x,"Process ",i3," : Secondary electron emission data file is found. Reading the data file...")', Rank_of_process
     END IF

     READ (9, '(A77)') buf ! **************************** LEFT PLASMA BOUNDARY ***************************")')
     READ (9, '(A77)') buf ! --- Refluxing of both electrons and ions:                                 ---")')
     READ (9, '(A77)') buf ! ---           0 = off, solid wall                                         ---")')
     READ (9, '(A77)') buf ! ---           1 = thermalization with bulk Te and Ti                      ---")')
     READ (9, '(A77)') buf ! -------d----- 2 = specular reflection ---------------------------------------")')
     READ (9, '(7x,i1)') PlasmaSourceFlag
     READ (9, '(A77)') buf ! --- Additional injection of electron-ion pairs (requires 1 or 2 above):   ---")')
     READ (9, '(A77)') buf ! ---           0 = off                                                     ---")')
     READ (9, '(A77)') buf ! ---           1 = additional flux equals to ion flux to the right wall    ---")')
     READ (9, '(A77)') buf ! -------d----- 2 = constant (Maxwellian source with Te and Ti of the bulk) ---")')
     READ (9, '(7x,i1)') AddInjectionFlag     

     READ (9, '(A77)') buf ! ***************************** ELASTIC REFLECTION ****************************")')
     READ (9, '(A77)') buf ! --- Select the model:                                                     ---")')
     READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
     READ (9, '(A77)') buf ! ---           1 = constant, nonzero for E_min < energy < E_max            ---")')
     READ (9, '(A77)') buf ! -------d----- 2 = function 1-1 of incident electron energy ------------------")')
     READ (9, '(7x,i1)') Emitted_model(1)
     READ (9, '(A77)') buf ! ================== ELASTIC REFLECTION, TYPE OF REFLECTION ===================")')
     READ (9, '(A77)') buf ! --- Select the type:                                                      ---")')
     READ (9, '(A77)') buf ! ---           0 = specular                                                ---")')
     READ (9, '(A77)') buf ! -------d----- 1 = random                                                  ---")')
     READ (9, '(7x,i1)') Elast_refl_type
     READ (9, '(A77)') buf ! ====================== ELASTIC, MODEL 1, PARAMETERS =========================")')
     READ (9, '(A77)') buf ! -------d.ddd- Value of the coefficient (from 0 to 1) ------------------------")')
     READ (9, '(7x,f5.3)') setD_see_elastic
     READ (9, '(A77)') buf ! --dddddd.ddd- LOWER energy boundary E_min, [eV] -----------------------------")')
     READ (9, '(2x,f10.3)') minE_see_elastic_eV
     READ (9, '(A77)') buf ! --dddddd.ddd- UPPER energy boundary E_max, [eV] -----------------------------")')
     READ (9, '(2x,f10.3)') maxE_see_elastic_eV
     READ (9, '(A77)') buf ! ====================== ELASTIC, MODEL 2, PARAMETERS =========================")')
     READ (9, '(A77)') buf ! --dddddd.ddd- Threshold energy, [eV] ----------------------------------------")')
     READ (9, '(2x,f10.3)') E_elast_0_eV
     READ (9, '(A77)') buf ! --dddddd.ddd- Maximal yield energy, [eV] ------------------------------------")')
     READ (9, '(2x,f10.3)') E_elast_max_eV
     READ (9, '(A77)') buf ! -------d.ddd- Maximal yield at low energies, [eV] ---------------------------")')
     READ (9, '(7x,f5.3)') maxD_elast
     READ (9, '(A77)') buf ! --dddddd.ddd- Half-width of decaying part after the maximum, [eV] -----------")')
     READ (9, '(2x,f10.3)') dE_elast_eV
     READ (9, '(A77)') buf ! -------d.ddd- Yield, high energy (fraction of total yield, 0 if not used) ---")')
     READ (9, '(7x,f5.3)') Frac_elast_highenergy

     READ (9, '(A77)') buf ! *************************** INELASTIC BACKSCATTERING ************************")')
     READ (9, '(A77)') buf ! --- Select the model:                                                     ---")')
     READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
     READ (9, '(A77)') buf ! ---           1 = constant, nonzero for E_min < energy < E_max            ---")')
     READ (9, '(A77)') buf ! -------d----- 2 = function 2-1 of incident electron energy ------------------")')
     READ (9, '(7x,i1)') Emitted_model(2)
     READ (9, '(A77)') buf ! ====================== INELASTIC, MODEL 1, PARAMETERS =======================")')
     READ (9, '(A77)') buf ! -------d.ddd- Value of the coefficient (from 0 to 1) ------------------------")')
     READ (9, '(7x,f5.3)') setD_see_inelastic
     READ (9, '(A77)') buf ! --dddddd.ddd- LOWER energy boundary E_min, [eV] -----------------------------")')
     READ (9, '(2x,f10.3)') minE_see_inelastic_eV
     READ (9, '(A77)') buf ! --dddddd.ddd- UPPER energy boundary E_max, [eV] -----------------------------")')
     READ (9, '(2x,f10.3)') maxE_see_inelastic_eV
     READ (9, '(A77)') buf ! ====================== INELASTIC, MODEL 2, PARAMETERS =======================")')
     READ (9, '(A77)') buf ! -------d.ddd- Yield (fraction of total yield, 0 if not used) ----------------")')
     READ (9, '(7x,f5.3)') Frac_inelastic
     
     READ (9, '(A77)') buf ! ************************** TRUE SECONDARY EMISSION **************************")')
     READ (9, '(A77)') buf ! --- Select the model:                                                     ---")')
     READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
     READ (9, '(A77)') buf ! ---           1 = constant, nonzero for E_min < energy < E_max            ---")')
     READ (9, '(A77)') buf ! -------d----- 2 = function 3-1 of incident electron energy ------------------")')
     READ (9, '(7x,i1)') Emitted_model(3)
     READ (9, '(A77)') buf ! ================= TRUE SECONDARY, MODEL 1, PARAMETERS =======================")')
     READ (9, '(A77)') buf ! -------d.ddd- Value of the coefficient --------------------------------------")')
     READ (9, '(7x,f5.3)') setD_see_true
     READ (9, '(A77)') buf ! --dddddd.ddd- LOWER energy boundary E_min, [eV] -----------------------------")')
     READ (9, '(2x,f10.3)') minE_see_true_eV
     READ (9, '(A77)') buf ! --dddddd.ddd- UPPER energy boundary E_max, [eV] -----------------------------")')
     READ (9, '(2x,f10.3)') maxE_see_true_eV
     READ (9, '(A77)') buf ! ============= TRUE SECONDARY, MODEL 2 (AND CLASSIC), PARAMETERS =============")')
     READ (9, '(A77)') buf ! --dddddd.ddd- Threshold energy, [eV] ----------------------------------------")')
     READ (9, '(2x,f10.3)') E_see_0_eV
     READ (9, '(A77)') buf ! --dddddd.ddd- Maximal emission energy, [eV] ---------------------------------")')
     READ (9, '(2x,f10.3)') E_see_max_eV
     READ (9, '(A77)') buf ! ------dd.ddd- Maximal emission coefficient (normal to the surface) ----------")')
     READ (9, '(6x,f6.3)') maxD_see_classic
     READ (9, '(A77)') buf ! -------d.ddd- Smoothness factor (0 = very rough, 2 = polished) --------------")')
     READ (9, '(7x,f5.3)') k_smooth
     READ (9, '(A77)') buf ! ============= PARAMETERS OF INJECTED TRUE SECONDARY ELECTRONS ===============")')
     READ (9, '(A77)') buf ! --dddddd.ddd- Temperature, [eV] ---------------------------------------------")')
     READ (9, '(2x,f10.3)') T_see_true_eV 

     READ (9, '(A77)') buf ! **************************** INTERACTION OF IONS ****************************")')
     READ (9, '(A77)') buf ! --- Select the model:                                                     ---")')
     READ (9, '(A77)') buf ! ---           0 = 100% adsorption                                         ---")')
     READ (9, '(A77)') buf ! ---           1 = 100% specular reflection                                ---")')
     READ (9, '(A77)') buf ! -------d----- 2 = 100% ion adsorption, secondary electrons may be emitted ---")')
     READ (9, '(7x,i1)') Ion_interac_model
     READ (9, '(A77)') buf ! ============ ION INDUCED SECONDARY ELECTRON EMISSION, PARAMETERS ============")')
     READ (9, '(A77)') buf ! -------d.ddd- Value of the coefficient (<1) ---------------------------------")')
     READ (9, '(7x,f5.3)') setD_ionsee_true
     READ (9, '(A77)') buf ! --dddddd.ddd- Minimal ion energy to start the emission, [eV] ----------------")')
     READ (9, '(2x,f10.3)') minE_ionsee_eV
     READ (9, '(A77)') buf ! --dddddd.ddd- Temperature of injected secondary electron, [eV] --------------")')
     READ (9, '(2x,f10.3)') T_ionsee_eV
   
                
  ELSE

     PlasmaSourceFlag      = 0                ! at the left the plasma is bounded by the dielectric wall
     AddInjectionFlag      = 0                ! and there is no additional injection

     Emitted_model(1)      = 2                ! elastic backscattering is turned on
     Elast_refl_type       = 1                ! type of reflection set to random
! elastic, model 1
     setD_see_elastic      = 1.0_8
     minE_see_elastic_eV   = 0.0_8
     maxE_see_elastic_eV   = 10000.0_8
! elastic, model 2
     E_elast_0_eV          =  2.0_8            
     E_elast_max_eV        = 10.0_8
     maxD_elast            =  0.55_8
     dE_elast_eV           = 14.0_8
     Frac_elast_highenergy =  0.03_8

     Emitted_model(2)      = 2                ! inelastic backscattering is turned on
! inelastic, model 1
     setD_see_inelastic    = 1.0_8
     minE_see_inelastic_eV = 0.0_8
     maxE_see_inelastic_eV = 10000.0_8
! inelastic, model 2
     Frac_inelastic        = 0.07_8

     Emitted_model(3)      = 2                ! true secondary emission is turned on
! true, model 1
     setD_see_true         = 0.3_8
     minE_see_true_eV      = 4.0_8
     maxE_see_true_eV      = 10000.0_8
! true, model 2
     E_see_0_eV            = 13.0_8
     E_see_max_eV          = 500.0_8
     maxD_see_classic      = 3.0_8
     k_smooth              = 1.0_8

     T_see_true_eV         = 2.0_8

     Ion_interac_model     = 0                 ! all ions are adsorbed at the walls

     setD_ionsee_true      = 0.1_8
     minE_ionsee_eV        = 100.0_8
     T_ionsee_eV           = 1.0_8

     IF (Rank_of_process.EQ.0) THEN

        PRINT '(/2x,"Process ",i3," : File with the name ssc_secondary.dat not found. Use the default settings ...")', Rank_of_process

        PRINT '(/2x,"Process ",i3," : Create ssc_secondary.dat file . . .")', Rank_of_process

        WRITE (9, '("**************************** LEFT PLASMA BOUNDARY ***************************")')
        WRITE (9, '("--- Refluxing of both electrons and ions:                                 ---")')
        WRITE (9, '("---           0 = off, solid wall                                         ---")')
        WRITE (9, '("---           1 = thermalization with bulk Te and Ti                      ---")')
        WRITE (9, '("-------d----- 2 = specular reflection ---------------------------------------")')
        WRITE (9, '(7x,i1)') PlasmaSourceFlag
        WRITE (9, '("--- Additional injection of electron-ion pairs (requires 1 or 2 above):   ---")')
        WRITE (9, '("---           0 = off                                                     ---")')
        WRITE (9, '("---           1 = additional flux equals to ion flux to the right wall    ---")')
        WRITE (9, '("-------d----- 2 = constant (Maxwellian source with Te and Ti of the bulk) ---")')
        WRITE (9, '(7x,i1)') AddInjectionFlag     

        WRITE (9, '("***************************** ELASTIC REFLECTION ****************************")')
        WRITE (9, '("--- Select the model:                                                     ---")')
        WRITE (9, '("---           0 = turned off                                              ---")')
        WRITE (9, '("---           1 = constant, nonzero for E_min < energy < E_max            ---")')
        WRITE (9, '("-------d----- 2 = function 1-1 of incident electron energy ------------------")')
        WRITE (9, '(7x,i1)') Emitted_model(1)
        WRITE (9, '("================== ELASTIC REFLECTION, TYPE OF REFLECTION ===================")')
        WRITE (9, '("--- Select the type:                                                      ---")')
        WRITE (9, '("---           0 = specular                                                ---")')
        WRITE (9, '("-------d----- 1 = random                                                  ---")')
        WRITE (9, '(7x,i1)') Elast_refl_type
        WRITE (9, '("====================== ELASTIC, MODEL 1, PARAMETERS =========================")')
        WRITE (9, '("-------d.ddd- Value of the coefficient (from 0 to 1) ------------------------")')
        WRITE (9, '(7x,f5.3)') setD_see_elastic
        WRITE (9, '("--dddddd.ddd- LOWER energy boundary E_min, [eV] -----------------------------")')
        WRITE (9, '(2x,f10.3)') minE_see_elastic_eV
        WRITE (9, '("--dddddd.ddd- UPPER energy boundary E_max, [eV] -----------------------------")')
        WRITE (9, '(2x,f10.3)') maxE_see_elastic_eV
        WRITE (9, '("====================== ELASTIC, MODEL 2, PARAMETERS =========================")')
        WRITE (9, '("--dddddd.ddd- Threshold energy, [eV] ----------------------------------------")')
        WRITE (9, '(2x,f10.3)') E_elast_0_eV
        WRITE (9, '("--dddddd.ddd- Maximal yield energy, [eV] ------------------------------------")')
        WRITE (9, '(2x,f10.3)') E_elast_max_eV
        WRITE (9, '("-------d.ddd- Maximal yield at low energies, [eV] ---------------------------")')
        WRITE (9, '(7x,f5.3)') maxD_elast
        WRITE (9, '("--dddddd.ddd- Half-width of decaying part after the maximum, [eV] -----------")')
        WRITE (9, '(2x,f10.3)') dE_elast_eV
        WRITE (9, '("-------d.ddd- Yield, high energy (fraction of total yield, 0 if not used) ---")')
        WRITE (9, '(7x,f5.3)') Frac_elast_highenergy

        WRITE (9, '("*************************** INELASTIC BACKSCATTERING ************************")')
        WRITE (9, '("--- Select the model:                                                     ---")')
        WRITE (9, '("---           0 = turned off                                              ---")')
        WRITE (9, '("---           1 = constant, nonzero for E_min < energy < E_max            ---")')
        WRITE (9, '("-------d----- 2 = function 2-1 of incident electron energy ------------------")')
        WRITE (9, '(7x,i1)') Emitted_model(2)
        WRITE (9, '("====================== INELASTIC, MODEL 1, PARAMETERS =======================")')
        WRITE (9, '("-------d.ddd- Value of the coefficient (from 0 to 1) ------------------------")')
        WRITE (9, '(7x,f5.3)') setD_see_inelastic
        WRITE (9, '("--dddddd.ddd- LOWER energy boundary E_min, [eV] -----------------------------")')
        WRITE (9, '(2x,f10.3)') minE_see_inelastic_eV
        WRITE (9, '("--dddddd.ddd- UPPER energy boundary E_max, [eV] -----------------------------")')
        WRITE (9, '(2x,f10.3)') maxE_see_inelastic_eV
        WRITE (9, '("====================== INELASTIC, MODEL 2, PARAMETERS =======================")')
        WRITE (9, '("-------d.ddd- Yield (fraction of total yield, 0 if not used) ----------------")')
        WRITE (9, '(7x,f5.3)') Frac_inelastic

        WRITE (9, '("************************** TRUE SECONDARY EMISSION **************************")')
        WRITE (9, '("--- Select the model:                                                     ---")')
        WRITE (9, '("---           0 = turned off                                              ---")')
        WRITE (9, '("---           1 = constant, nonzero for E_min < energy < E_max            ---")')
        WRITE (9, '("-------d----- 2 = function 3-1 of incident electron energy ------------------")')
        WRITE (9, '(7x,i1)') Emitted_model(3)
        WRITE (9, '("================= TRUE SECONDARY, MODEL 1, PARAMETERS =======================")')
        WRITE (9, '("-------d.ddd- Value of the coefficient --------------------------------------")')
        WRITE (9, '(7x,f5.3)') setD_see_true
        WRITE (9, '("--dddddd.ddd- LOWER energy boundary E_min, [eV] -----------------------------")')
        WRITE (9, '(2x,f10.3)') minE_see_true_eV
        WRITE (9, '("--dddddd.ddd- UPPER energy boundary E_max, [eV] -----------------------------")')
        WRITE (9, '(2x,f10.3)') maxE_see_true_eV
        WRITE (9, '("============= TRUE SECONDARY, MODEL 2 (AND CLASSIC), PARAMETERS =============")')
        WRITE (9, '("--dddddd.ddd- Threshold energy, [eV] ----------------------------------------")')
        WRITE (9, '(2x,f10.3)') E_see_0_eV
        WRITE (9, '("--dddddd.ddd- Maximal emission energy, [eV] ---------------------------------")')
        WRITE (9, '(2x,f10.3)') E_see_max_eV
        WRITE (9, '("------dd.ddd- Maximal emission coefficient (normal to the surface) ----------")')
        WRITE (9, '(6x,f6.3)') maxD_see_classic
        WRITE (9, '("-------d.ddd- Smoothness factor (0 = very rough, 2 = polished) --------------")')
        WRITE (9, '(7x,f5.3)') k_smooth
        WRITE (9, '("============= PARAMETERS OF INJECTED TRUE SECONDARY ELECTRONS ===============")')
        WRITE (9, '("--dddddd.ddd- Temperature, [eV] ---------------------------------------------")')
        WRITE (9, '(2x,f10.3)') T_see_true_eV 

        WRITE (9, '("**************************** INTERACTION OF IONS ****************************")')
        WRITE (9, '("--- Select the model:                                                     ---")')
        WRITE (9, '("---           0 = 100% adsorption                                         ---")')
        WRITE (9, '("---           1 = 100% specular reflection                                ---")')
        WRITE (9, '("-------d----- 2 = 100% ion adsorption, secondary electrons may be emitted ---")')
        WRITE (9, '(7x,i1)') Ion_interac_model
        WRITE (9, '("============ ION INDUCED SECONDARY ELECTRON EMISSION, PARAMETERS ============")')
        WRITE (9, '("-------d.ddd- Value of the coefficient (<1) ---------------------------------")')
        WRITE (9, '(7x,f5.3)') setD_ionsee_true
        WRITE (9, '("--dddddd.ddd- Minimal ion energy to start the emission, [eV] ----------------")')
        WRITE (9, '(2x,f10.3)') minE_ionsee_eV
        WRITE (9, '("--dddddd.ddd- Temperature of injected secondary electron, [eV] --------------")')
        WRITE (9, '(2x,f10.3)') T_ionsee_eV
        
     END IF
 
  END IF

  CLOSE (9, STATUS = 'KEEP')

  minE_see_elastic = 0.5_8 * minE_see_elastic_eV * N_box_vel**2 / T_e_eV
  maxE_see_elastic = 0.5_8 * maxE_see_elastic_eV * N_box_vel**2 / T_e_eV

  E_elast_0   = 0.5_8 * E_elast_0_eV   * N_box_vel**2 / T_e_eV
  E_elast_max = 0.5_8 * E_elast_max_eV * N_box_vel**2 / T_e_eV
  dE_elast    = 0.5_8 * dE_elast_eV    * N_box_vel**2 / T_e_eV
  
  minE_see_inelastic = 0.5_8 * minE_see_inelastic_eV * N_box_vel**2 / T_e_eV
  maxE_see_inelastic = 0.5_8 * maxE_see_inelastic_eV * N_box_vel**2 / T_e_eV
  
  minE_see_true = 0.5_8 * minE_see_true_eV * N_box_vel**2 / T_e_eV
  maxE_see_true = 0.5_8 * maxE_see_true_eV * N_box_vel**2 / T_e_eV
  
  E_see_0   = 0.5_8 * E_see_0_eV   * N_box_vel**2 / T_e_eV        
  E_see_max = 0.5_8 * E_see_max_eV * N_box_vel**2 / T_e_eV

  IF ((PlasmaSourceFlag.NE.0).AND.(AddInjectionFlag.EQ.2)) THEN
     CALL INITIATE_ADDITIONAL_CONSTANT_INJECTION                    ! both server and clients
  END IF

! for the server process only:
  IF (Rank_of_process.EQ.0) THEN

! write the file of the see coefficient vs incident energy at normal incidence (teta = 0.0): total / elast./ inelast./ true
     OPEN (40, FILE = 'dim_seecoef_vs_energy.dat')
     angle_deg = 0.0_8                                  ! set the angle of incidence in degrees (for testing only)
     angle_rad = angle_deg * 3.141592653_8 / 180.0_8    ! conver angle of incidence from degrees to radians
     DO j = 0, 1000
        energy = 0.5_8 * DBLE(j) * N_box_vel**2 / T_e_eV     
        tot_see_coef = Coeff_SEE_Elastic(energy, angle_rad) + Coeff_SEE_Inelastic(energy, angle_rad) + &
                     & Coeff_SEE_True(energy, angle_rad) 
        WRITE (40, '(5(2x,f10.4))') DBLE(j), tot_see_coef, & 
                                           & Coeff_SEE_Elastic(energy, angle_rad), & 
                                           & Coeff_SEE_Inelastic(energy, angle_rad), & 
                                           & Coeff_SEE_True(energy, angle_rad)
     END DO
     CLOSE (40, STATUS = 'KEEP')

     IF (PlasmaSourceFlag.EQ.0) THEN
        PRINT '(/2x,"The plasma is bounded by two solid walls")'
        IF (AddInjectionFlag.NE.0) THEN
           PRINT '(/2x,"Error: A solid wall should not be combined with additional injection of e-i pairs!")'
           PRINT '( 2x,"Program will be terminated now, sorry :(")'
           STOP
        END IF
     ELSE IF (PlasmaSourceFlag.EQ.1) THEN
        PRINT '(/2x,"At the LEFT the plasma borders with the PLASMA SOURCE")'
        PRINT '(2x, "Bulk particles leaving the plasma here are returned after THERMALIZATION")'
     ELSE IF (PlasmaSourceFlag.EQ.2) THEN
        PRINT '(/2x,"At the LEFT the plasma borders with the PLASMA SOURCE")'
        PRINT '(2x, "Bulk particles leaving the plasma here reflect SPECULARLY from the left boundary")'
     ELSE
        PRINT '(/2x,"Error: unknown plasma source flag specified at the left wall: ",i4)', PlasmaSourceFlag
        PRINT '( 2x,"Program will be terminated now, sorry :(")'
        STOP
     END IF

     IF (PlasmaSourceFlag.NE.0) THEN
        IF (BC_flag.NE.1) THEN
           PRINT '(/2x,"Error: the potential boundary conditions are incompatible with the plasma source")'
           PRINT '( 2x,"Choose the floating potentials (zero field inside the plasma source)...")'
           PRINT '( 2x,"Program will be terminated now, sorry :(")'
           STOP
        END IF
        IF (AddInjectionFlag.EQ.0) THEN
           PRINT '(2x, "Additional injection of electron-ion pairs turned off")'
        ELSE IF (AddInjectionFlag.EQ.1) THEN
           PRINT '(2x, "The flux of additionally injected electron-ion pairs is set equal to the flux of ions at the right wall")'
        ELSE IF (AddInjectionFlag.EQ.2) THEN
           PRINT '(2x, "The flux of additionally injected electron-ion pairs is constant and equal to initial ion Bohm flux of bulk plasma")'
        ELSE
           PRINT '(/2x,"Error: unknown additional injection flag specified at the left wall: ",i4)', AddInjectionFlag
           PRINT '( 2x,"Program will be terminated now, sorry :(")'
           STOP
        END IF
     END IF
     
     PRINT '(/2x,"Elastic electron backscattering:")'
     SELECT CASE (Emitted_model(1))
     CASE (0)
        PRINT '(4x,"turned off :)")'
     CASE (1)
        PRINT '(4x,"constant (delta = ",f5.2,") for ",f7.1," (eV) < energy < ",f7.1," (eV)")', setD_see_elastic,    &
             & minE_see_elastic_eV, &
             & maxE_see_elastic_eV
     CASE (2) 
        PRINT '(4x,"function 1-1 of incident electron energy")'
     END SELECT
     
     PRINT '(/2x,"Inelastic electron backscattering:")'
     SELECT CASE (Emitted_model(2))
     CASE (0)
        PRINT '(4x,"turned off :)")'
     CASE (1)
        PRINT '(4x,"constant (delta = ",f5.2,") for ",f7.1," (eV) < energy < ",f7.1," (eV)")', setD_see_inelastic,    &
             & minE_see_inelastic_eV, &
             & maxE_see_inelastic_eV
     CASE (2) 
        PRINT '(4x,"function 2-1 of incident electron energy")'
     END SELECT
     
     PRINT '(/2x,"True secondary electron emission:")'
     SELECT CASE (Emitted_model(3))
     CASE (0)
        PRINT '(4x,"turned off :)")'
     CASE (1)
        PRINT '(4x,"constant (delta = ",f5.2,") for ",f7.1," (eV) < energy < ",f7.1," (eV)")', setD_see_true,    &
             & minE_see_true_eV, &
             & maxE_see_true_eV
     CASE (2) 
        PRINT '(4x,"function 3-1 of incident energy")'
     END SELECT
     
     PRINT '(/2x,"Ions at the walls are:")'
     SELECT CASE (Ion_interac_model)
     CASE (0)
        PRINT '(4x,"100 % adsorbed :)")'
     CASE (1)
        PRINT '(4x,"100 % elastically reflected")'
     CASE (2)
        PRINT '(4x,"100 % adsorbed, ion induced secondary electron emission is ON")'
     END SELECT
     
  END IF
  
  CALL PrepareEnergyDistribIntegral

  electron_reflux_count  = 0;  electron_reflux_energy = 0.0_8   
  ion_reflux_count       = 0;  ion_reflux_energy      = 0.0_8
  electron_reemit_count  = 0;  electron_reemit_energy = 0.0_8   
  ion_reemit_count       = 0;  ion_reemit_energy      = 0.0_8

  N_of_lost_ions = 0

  ion_left_reflect_count  = 0; ion_left_reflect_energy  = 0.0_8
  ion_right_reflect_count = 0; ion_right_reflect_energy = 0.0_8

  see_left_elastic_count   = 0; see_left_elastic_energy   = 0.0_8 
  see_left_inelastic_count = 0; see_left_inelastic_energy = 0.0_8
  see_left_true_count      = 0; see_left_true_energy      = 0.0_8

  see_right_elastic_count   = 0; see_right_elastic_energy   = 0.0_8 
  see_right_inelastic_count = 0; see_right_inelastic_energy = 0.0_8
  see_right_true_count      = 0; see_right_true_energy      = 0.0_8

  prie_left_from_right_count = 0; prie_left_from_right_energy = 0.0_8
  prie_left_after_coll_count = 0; prie_left_after_coll_energy = 0.0_8

  prie_right_from_left_count  = 0; prie_right_from_left_energy  = 0.0_8
  prie_right_after_coll_count = 0; prie_right_after_coll_energy = 0.0_8

  sece_left_from_right_count = 0; sece_left_after_coll_count = 0
  sece_right_from_left_count = 0; sece_right_after_coll_count = 0

  ionsee_left_count  = 0; ionsee_left_energy = 0.0_8
  ionsee_right_count = 0; ionsee_right_energy = 0.0_8

END SUBROUTINE INITIATE_SE_EMISSION 

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLLISION_WITH_WALL(s, x, k)

  USE ParallelOperationValues
  USE SEEmission
  USE CurrentProblemValues

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER s            ! species
  REAL(8) x            ! coordinate of particle, must be out of the plasma limits
                       ! depending on the moment of detection of the collision with the walls,
                       ! x can be either streaming or corrected coordinate
  INTEGER k            ! index of colliding particle

  REAL(8) Coeff_SEE_Elastic
  REAL(8) Coeff_SEE_Inelastic
  REAL(8) Coeff_SEE_True
  
  REAL(8) R                  ! random number 

  REAL(8) vx, vy, vz, v   ! velocity components of incident particle
  INTEGER tag             ! tag of incident particle
  REAL(8) energy_inc      ! energy of incident particle
  REAL(8) teta_inc        ! angle of incidence with respect to wall normal 

  REAL(8) coef_elastic    !
  REAL(8) coef_inelastic  ! corresponding emission coefficients at current energy and angle of incidence
  REAL(8) coef_true       !
  REAL(8) coef_total      !


! take the particle
  vx = species(s)%part(k)%VX
  vy = species(s)%part(k)%VY
  vz = species(s)%part(k)%VZ
  tag = species(s)%part(k)%Tag

  v = SQRT(vx**2 + vy**2 + vz**2)

! if necessary, perform refluxing at the left wall
  IF (x.LT.0.0_8) THEN
     IF (PlasmaSourceFlag.EQ.1) THEN
        CALL LEFT_WALL_REFLUX_THERMAL(x, vx, s)
        RETURN
     ELSE IF (PlasmaSourceFlag.EQ.2) THEN
        CALL LEFT_WALL_REFLUX_SPECULAR(x, vx, vy, vz, v*v, s)
        RETURN
     END IF
  END IF

! calculate the dim-less energy of incident electron
  energy_inc = 0.5_8 * v * v

! calculate the angle of incidence
  teta_inc = ACOS(ABS(vx)/v)  !!!!  0.0_8  !!!!   

! calculate the coefficients of elastic/inelastic scattering and true secondary emission
  coef_elastic   = Coeff_SEE_Elastic(energy_inc, teta_inc)
  coef_inelastic = Coeff_SEE_Inelastic(energy_inc, teta_inc)
  IF ((coef_elastic + coef_inelastic).GT.1.0_8) THEN
     PRINT '("Process ",i3," : Error in secondary emission model!")', Rank_of_process
     PRINT '("The total coefficient of elastic/inelastic backscattering is greater than 1 !!!")'
     PRINT '("elastic: ",f5.2," inelastic: ",f5.2)', coef_elastic, coef_inelastic
     PRINT '("Program will be terminated now :(")'
     STOP
  END IF  
  coef_true      = Coeff_SEE_True(energy_inc, teta_inc)

! then we combine the rest (fractional part) of the true secondary emission coefficient 
! with the coefficients of elastic/inelastic backscattering to complete the total emission
  coef_total = coef_elastic + coef_inelastic + coef_true

! quit subroutine if SEE is already exhausted / turned off
  IF (coef_total.EQ.0.0_8) RETURN

! We can always treat coefficients (relative yields) of elastic/inelastic backscattering as probabilities.
! If the total coefficient is less or equal than 1, we simply treat ALL coefficients as probabilities.
! If at this stage the total coefficient is greater than 1, this means that, on one hand, 
! not all electrons will produce true secondary electrons, 
! on the other hand, those incident electrons, which produce the true secondary electrons, 
! will (in average) produce more than 1 secondary electron per 1 incident electron. 
! We cannot treat the true secondary coefficient as the probability yet,
! (the difference from the other coefficients is in that only true secondary emission 
! can produce more than 1 electron per incident electron).
! Also, we cannot start here with the injection of the true secondaries, because afterwards it is 
! possible that the elastically or inelastically reflected particle will be injected.
! Therefore we must start from the attempts to inject a reflected (elastically or inelastically) particle.

! take a random number
  R = well_random_number()

! try to inject the inelasticaly backscattered electron
  IF ((R.GT.coef_elastic).AND.(R.LE.(coef_elastic + coef_inelastic))) THEN
     CALL INJECT_INELASTIC_BACKSCATTERED(x, vx, v, tag)
     RETURN
  END IF

! try  to inject the elastically reflected electron
  IF (R.LE.coef_elastic) THEN
     CALL INJECT_ELASTIC_REFLECTED(x, vx, vy, vz, v, tag)       
     RETURN
  END IF

  IF (coef_total.LE.1.0_8) THEN   ! in this case we can treat the true secondary yield as the probability
! try to inject the true secondary electron
     IF ((R.GT.(coef_elastic + coef_inelastic)).AND.(R.LE.coef_total)) THEN
        CALL INJECT_TRUE_SECONDARY(x, vx, tag)
     END IF
! if nothing was injected (R.GT.coef_total), the particle is attached to the wall and we leave the subroutine
     RETURN
  END IF

! if we are here the total SEE yield is greater than 1 (coef_total.GT.1) and at least one true secondary electron must be injected
! we calculate the ratio of the emitted true secondary flux 
! to the portion of the incident flux, which does not produce elastically/inelastically reflected electrons:
  coef_true = coef_true / (1.0_8 - (coef_elastic + coef_inelastic))   ! note, we cannot be here if coef_elastic + coef_inelastic = 1

! then we inject true secondary electrons until we exhaust the integer part of the ratio
  DO WHILE (coef_true.GE.1.0_8)
     CALL INJECT_TRUE_SECONDARY(x, vx, tag)
     coef_true = coef_true - 1.0_8
  END DO
! the remaining part gives us the probability for statistical injection of fractional part of true secondary electron

! take a random number
  R = well_random_number()
  IF (R.LT.coef_true) CALL INJECT_TRUE_SECONDARY(x, vx, tag)

END SUBROUTINE PROCESS_ELECTRON_COLLISION_WITH_WALL

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE INJECT_ELASTIC_REFLECTED(x, vx, vy, vz, v, tag)

  USE ParallelOperationValues
  USE SEEmission
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_emit, Rate_energy_leftemit, Rate_energy_rightemit, &
                        & Rate_number_leftemit, Rate_number_rightemit

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) x             ! coordinate of incident particle, must be beyond the plasma boundaries 
  REAL(8) vx, vy, vz    ! velocity components of incident particle
  REAL(8) v             ! absolute value of velocity
  INTEGER tag           ! tag of incident particle

  INTEGER ALLOC_ERR

  REAL(8) teta, fi                 ! angles of backscattering

  REAL(8) x_new                    ! x-coordinate of the electron after reflection (to be injected)
!  REAL(8) v_new                    ! absolute velocity   - " -
  REAL(8) vx_new, vy_new, vz_new   ! velocity components - " -

  INTEGER left_node                ! left node  - " -
  INTEGER right_node               ! right node - " -

!!!  IF (Elast_refl_type.NE.0) Elast_refl_type = 0       

  IF (Elast_refl_type.EQ.0) THEN    ! if elastic reflection occurs specularly

     IF (x.GT.DBLE(N_cells)) THEN              ! if the particle hits the RIGHT WALL
        Q_right =  Q_right - Qs(1)             ! restore the surface charge density
!        x_new = 2.0_8 * N_cells - x            ! calculate the new coordinate
        x_new = dble(N_cells)                  ! calculate the new coordinate
        Q_strm_spec(N_cells, 1) = Q_strm_spec(N_cells, 1) + 1.0_8 
     END IF

     IF (x.LT.0.0_8) THEN                      ! if the particle hits the LEFT WALL
        Q_left = Q_left - Qs(1)                ! restore the surface charge density 
!        x_new = - x                            ! calculate the new coordinate
        x_new = 0.0_8                            ! calculate the new coordinate
        Q_strm_spec(0, 1) = Q_strm_spec(0, 1) + 1.0_8 
     END IF

     vx_new = -vx
     vy_new =  vy
     vz_new =  vz
  
  ELSE                             ! if elastic reflection occurs at random angle

!     teta = 1.570796327_8 * RAN(I_random_seed)                  ! get the angles of reflection
     teta = ASIN(MIN(well_random_number(),1.0_8))                ! assume that the distribution over the teta angle is f(teta) = COS(teta)   
     fi   = 6.283185307_8 * well_random_number()

     vy_new = v * SIN(teta) * SIN(fi)          ! get the Y,Z-velocity components
     vz_new = v * SIN(teta) * COS(fi)

     IF (x.GT.DBLE(N_cells)) THEN                                   ! if the particle hits the RIGHT WALL
        Q_right =  Q_right - Qs(1)                                  ! restore the surface charge density
        vx_new  = - v * COS(teta)                                   ! calculate the new X-velocity
!        x_new   = N_cells - (x - N_cells) * ABS(vx_new / vx)        ! calculate the new coordinate
        x_new = dble(N_cells)                  ! calculate the new coordinate
        Q_strm_spec(N_cells, 1) = Q_strm_spec(N_cells, 1) + 1.0_8 
      END IF

     IF (x.LT.0.0_8) THEN                      ! if the particle hits the LEFT WALL
        Q_left = Q_left - Qs(1)                ! restore the surface charge density
        vx_new = v * COS(teta)                 ! calculate the new X-velocity
!        x_new = ABS(x * vx_new / vx)           ! calculate the new coordinate
        x_new = 0.0_8                            ! calculate the new coordinate
        Q_strm_spec(0, 1) = Q_strm_spec(0, 1) + 1.0_8 
     END IF
      
  END IF

! Save the parameters of the injected (elastically reflected) electron in the linked list
  IF (ASSOCIATED(Current_electron)) THEN
     N_inject(1)         = N_inject(1) + 1
     IF (x_new.LT.(N_cells/2)) THEN                                         ! ------------------------------ if hits the LEFT wall
        Rate_number_leftemit(1) = Rate_number_leftemit(1) + 1               !
        Rate_energy_leftemit(1) = Rate_energy_leftemit(1) + (v**2)          !
        see_left_elastic_count  = see_left_elastic_count  + 1               ! for diagnostics
        see_left_elastic_energy = see_left_elastic_energy + (v**2)          !    
        SELECT CASE (tag)                                                   
           CASE (eTag_Emit_Right)                                           !
              sece_left_from_right_count = sece_left_from_right_count + 1   ! if primary electron was emitted from the right wall
           CASE (eTag_Coll_Neutral)                                         !
              sece_left_after_coll_count = sece_left_after_coll_count + 1   ! if primary electron collided with neutral atom
        END SELECT
        Current_electron%Tag = eTag_Emit_Left                               ! mark the electron
        CALL ADD_EMITTED_E_TO_LEFT_DF(vx_new, vy_new, vz_new)
     ELSE                                                                   ! ------------------------------ if hits the RIGHT wall
        Rate_number_rightemit(1) = Rate_number_rightemit(1) + 1             ! 
        Rate_energy_rightemit(1) = Rate_energy_rightemit(1) + (v**2)        !
        see_right_elastic_count  = see_right_elastic_count  + 1             ! for diagnostics
        see_right_elastic_energy = see_right_elastic_energy + (v**2)        !    
        SELECT CASE (tag)
           CASE (eTag_Emit_Left)                                            ! 
              sece_right_from_left_count = sece_right_from_left_count + 1   ! if primary electron was emitted from the left wall
           CASE (eTag_Coll_Neutral)                                         !
              sece_right_after_coll_count = sece_right_after_coll_count + 1 ! if primary electron collided with neutral atom
        END SELECT
        Current_electron%Tag = eTag_Emit_Right                              ! mark the electron
        CALL ADD_EMITTED_E_TO_RIGHT_DF(vx_new, vy_new, vz_new)
     END IF
     Current_electron%X       = x_new
     Current_electron%VX      = vx_new     ! -vx
     Current_electron%VY      = vy_new
     Current_electron%VZ      = vz_new
!     Current_electron%AX      = 0.0_8
     ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in INJECT_ELASTIC_REFLECTED:")', Rank_of_process
        PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
     Current_electron => Current_electron%next
     NULLIFY(Current_electron%next)
  ELSE
     PRINT '(/2x,"Process ",i3," : Error in INJECT_ELASTIC_REFLECTED:")', Rank_of_process
     PRINT  '(2x,"Current_electron is NOT associated!")'
     PRINT  '(2x,"The program will be terminated now :(")'
     STOP
  END IF
     
!! account for the density of emitted particle
!  left_node  = INT(x_new)               
!  right_node = left_node + 1
!  IF (left_node.LT.N_cells) THEN
!     Q_strm_spec(left_node, 1)  = Q_strm_spec(left_node, 1)  + right_node - x_new    ! s = 1
!     Q_strm_spec(right_node, 1) = Q_strm_spec(right_node, 1) + x_new - left_node     ! s = 1
!  ELSE   ! if particle was placed exactly on the wall
!     PRINT '(/2x,"Process ",i3," : INJECT_ELASTIC_REFLECTED: ERROR! emitted particle was placed NOT inside the plasma!")', &
!                                                                                                           & Rank_of_process
!     PRINT '(2x, "The program will be terminated now!")'
!     STOP
!  END IF

END SUBROUTINE INJECT_ELASTIC_REFLECTED

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE INJECT_INELASTIC_BACKSCATTERED(x, vx, v, tag) 

  USE ParallelOperationValues
  USE SEEmission
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_emit, Rate_energy_leftemit, Rate_energy_rightemit, &
                        & Rate_number_leftemit, Rate_number_rightemit

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) x             ! coordinate of incident particle, must be beyond the plasma boundaries 
  REAL(8) vx            ! x-velocity component of incident particle
  REAL(8) v             ! absolute velocity of incident electron
  INTEGER tag           ! tag of incident electron

  INTEGER ALLOC_ERR

  REAL(8) teta, fi                 ! angles of backscattering

  REAL(8) x_new                    ! x-coordinate of the electron after reflection (to be injected)
  REAL(8) v_new                    ! absolute velocity   - " -
  REAL(8) vx_new, vy_new, vz_new   ! velocity components - " -
  INTEGER left_node                ! left node           - " -
  INTEGER right_node               ! right node          - " -

! get the angles of backscattering
!  teta = 1.570796327_8 * RAN(I_random_seed) 
  teta = ASIN(MIN(1.0_8,well_random_number()))                ! assume that the distribution over the teta angle is f(teta) = COS(teta)   
  fi   = 6.283185307_8 * well_random_number() 
 
! get the absolute velocity of backscattered
  v_new = v * SQRT(well_random_number())           ! assume that the energy (NOT the velocity) is uniformly distributed between 0 and initial_energy

! get the Y,Z-velocity components
  vy_new = v_new * SIN(teta) * SIN(fi)
  vz_new = v_new * SIN(teta) * COS(fi)

! if the particle hits the RIGHT WALL
  IF (x.GT.DBLE(N_cells)) THEN
! restore the surface charge density
     Q_right =  Q_right - Qs(1)
! calculate the new X-velocity
     vx_new = - v_new * COS(teta)
! calculate the new coordinate
!     x_new = N_cells - (x - N_cells) * ABS(vx_new / vx)
     x_new = dble(N_cells)                  ! calculate the new coordinate
     Q_strm_spec(N_cells, 1) = Q_strm_spec(N_cells, 1) + 1.0_8 
  END IF

! if the particle hits the LEFT WALL
  IF (x.LT.0.0_8) THEN
! restore the surface charge density
     Q_left = Q_left - Qs(1)
! calculate the new X-velocity
     vx_new = v_new * COS(teta)
! calculate the new coordinate
!     x_new = ABS(x * vx_new / vx)
     x_new = 0.0_8                            ! calculate the new coordinate
     Q_strm_spec(0, 1) = Q_strm_spec(0, 1) + 1.0_8 
  END IF

!  left_node  = INT(x_new)
!  IF (left_node.GE.N_cells) THEN   
!     left_node = N_cells - 1
!     x_new = N_cells
!  END IF
!  right_node = left_node + 1

! Save the parameters of the injected (inelastically backscattered) electron in the linked list
  IF (ASSOCIATED(Current_electron)) THEN
     N_inject(1)         = N_inject(1) + 1
     IF (x_new.LT.(N_cells/2)) THEN                                             ! ---------------------------- if hits the LEFT wall
        Rate_number_leftemit(1) = Rate_number_leftemit(1) + 1                   !
        Rate_energy_leftemit(1) = Rate_energy_leftemit(1) + (v_new**2)          !
        see_left_inelastic_count  = see_left_inelastic_count + 1                ! for diagnostics
        see_left_inelastic_energy = see_left_inelastic_energy + (v_new**2)      !
        SELECT CASE (tag)                                                   
           CASE (eTag_Emit_Right)                                               !
              sece_left_from_right_count = sece_left_from_right_count + 1       ! if primary electron was emitted from the right wall
           CASE (eTag_Coll_Neutral)                                             !
              sece_left_after_coll_count = sece_left_after_coll_count + 1       ! if primary electron collided with neutral atom
        END SELECT
        Current_electron%Tag = eTag_Emit_Left                                   ! mark the electron
        CALL ADD_EMITTED_E_TO_LEFT_DF(vx_new, vy_new, vz_new)
     ELSE                                                                       ! ---------------------------- if hits the RIGHT wall
        Rate_number_rightemit(1) = Rate_number_rightemit(1) + 1                 ! 
        Rate_energy_rightemit(1) = Rate_energy_rightemit(1) + (v_new**2)        !
        see_right_inelastic_count  = see_right_inelastic_count + 1              ! for diagnostics
        see_right_inelastic_energy = see_right_inelastic_energy + (v_new**2)    !
        SELECT CASE (tag)
           CASE (eTag_Emit_Left)                                                ! 
              sece_right_from_left_count = sece_right_from_left_count + 1       ! if primary electron was emitted from the left wall
           CASE (eTag_Coll_Neutral)                                             !
              sece_right_after_coll_count = sece_right_after_coll_count + 1     ! if primary electron collided with neutral atom
        END SELECT
        Current_electron%Tag = eTag_Emit_Right                                  ! mark the electron
        CALL ADD_EMITTED_E_TO_RIGHT_DF(vx_new, vy_new, vz_new)
     END IF
     Current_electron%X       = x_new
     Current_electron%VX      = vx_new
     Current_electron%VY      = vy_new
     Current_electron%VZ      = vz_new
!     Current_electron%AX      = 0.0_8
     ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in INJECT_INELASTIC_BACKSCATTERED:")', Rank_of_process
        PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
     Current_electron => Current_electron%next
     NULLIFY(Current_electron%next)
  ELSE
     PRINT '(/2x,"Process ",i3," : Error in INJECT_INELASTIC_BACKSCATTERED:")', Rank_of_process
     PRINT  '(2x,"Current_electron is NOT associated!")'
     PRINT  '(2x,"The program will be terminated now :(")'
     STOP
  END IF

!! account for the density of emitted particle
!  Q_strm_spec(left_node, 1)  = Q_strm_spec(left_node, 1)  + right_node - x_new    ! s = 1
!  Q_strm_spec(right_node, 1) = Q_strm_spec(right_node, 1) + x_new - left_node     ! s = 1

END SUBROUTINE INJECT_INELASTIC_BACKSCATTERED

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE INJECT_TRUE_SECONDARY(x, vx_inc, tag)

  USE ParallelOperationValues
  USE SEEmission
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_emit, Rate_energy_leftemit, Rate_energy_rightemit, &
                        & Rate_number_leftemit, Rate_number_rightemit

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) x             ! coordinate of incident electron, must be beyond the plasma boundaries 
  REAL(8) vx_inc        ! x-velocity of the electron before the collision
  INTEGER tag           ! tag of the incident electron

  INTEGER ALLOC_ERR

  REAL(8) energy        ! energy of injected (true secondary emitted) electron
  REAL(8) teta, fi      ! angles of injection (emission)
  REAL(8) v             ! absolute velocity   - " -
  REAL(8) vx, vy, vz    ! velocity components - " -
  REAL(8) x_new         ! x-coordinate        - " -
  INTEGER left_node     ! left node           - " -
  INTEGER right_node    ! right node          - " -

if (vx_inc.eq.0.0_8) then
   print '(/2x,"Process ",i3," : Warning! Or maybe Error!!!")', Rank_of_process
   print *, 'the x-velocity of particle hitting the wall is zero!'
   print *, 'the program will be terminated now :('
   stop
end if

! get the energy of true secondary electron
  CALL GetTrueSecondaryEnergy(energy) 

! get the absolute velocity of true secondary electron
  v = SQRT(2.0_8 * energy)

! get the angles of backscattering
  teta = ASIN(MIN(1.0_8,well_random_number()))                ! assume that the distribution over the teta angle is f(teta) = COS(teta)   
  fi   = 6.283185307_8 * well_random_number()     

! get the Y,Z-velocity components
  vy = v * SIN(teta) * SIN(fi)
  vz = v * SIN(teta) * COS(fi)

! if the particle hits the RIGHT WALL
  IF (x.GT.N_cells) THEN
! decrease the surface charge density
     Q_right =  Q_right - Qs(1)
! calculate the new X-velocity
     vx = - v * COS(teta)
! set the new coordinate
! calculate the new coordinate
!     x_new = DBLE(N_cells) - (x - DBLE(N_cells)) * ABS(vx / vx_inc)
!!!!!!!!     x_new = X_node(N_nodes)
! set the left / right nodes
!!!!!!!!     left_node  = N_nodes - 1
!!!!!!!!     right_node = N_nodes
     x_new = dble(N_cells)                  ! calculate the new coordinate
     Q_strm_spec(N_cells, 1) = Q_strm_spec(N_cells, 1) + 1.0_8 
  END IF

! if the particle hits the LEFT WALL
  IF (x.LT.0.0_8) THEN
! decrease the surface charge density
     Q_left = Q_left - Qs(1)
! calculate the new X-velocity
     vx = v * COS(teta)
! set the new coordinate
!     x_new = ABS(x * vx / vx_inc)
!!!!!!!!!!!!!     x_new = 0.0_8
! set the left / right nodes
!!!!!!!!!!!!!     left_node = 1
!!!!!!!!!!!!!     right_node = 2
     x_new = 0.0_8                            ! calculate the new coordinate
     Q_strm_spec(0, 1) = Q_strm_spec(0, 1) + 1.0_8 
  END IF

!  left_node  = INT(x_new)
!  IF (left_node.EQ.N_cells) THEN
!     x_new = N_cells
!     left_node = N_cells - 1
!  END IF
!  right_node = left_node + 1

! Save the parameters of the injected (inelastically backscattered) electron in the linked list
  IF (ASSOCIATED(Current_electron)) THEN
     N_inject(1) = N_inject(1) + 1
     IF (x_new.LT.(N_cells/2)) THEN                                         ! --------------------------------- if hits the LEFT wall
        Rate_number_leftemit(1) = Rate_number_leftemit(1) + 1               !
        Rate_energy_leftemit(1) = Rate_energy_leftemit(1) + (v**2)          !
        see_left_true_count  = see_left_true_count  + 1                     ! for diagnostics
        see_left_true_energy = see_left_true_energy + (v**2)                !
        SELECT CASE (tag)                                                   
           CASE (eTag_Emit_Right)                                           !
              sece_left_from_right_count = sece_left_from_right_count + 1   ! if primary electron was emitted from the right wall
           CASE (eTag_Coll_Neutral)                                         !
              sece_left_after_coll_count = sece_left_after_coll_count + 1   ! if primary electron collided with neutral atom
        END SELECT
        Current_electron%Tag = eTag_Emit_Left                               ! mark the electron
        CALL ADD_EMITTED_E_TO_LEFT_DF(vx, vy, vz)
     ELSE                                                                   ! --------------------------------- if hits the RIGHT wall
        Rate_number_rightemit(1) = Rate_number_rightemit(1) + 1             ! 
        Rate_energy_rightemit(1) = Rate_energy_rightemit(1) + (v**2)        !
        see_right_true_count  = see_right_true_count  + 1                   ! for diagnostics
        see_right_true_energy = see_right_true_energy + (v**2)              !
        SELECT CASE (tag)
           CASE (eTag_Emit_Left)                                            ! 
              sece_right_from_left_count = sece_right_from_left_count + 1   ! if primary electron was emitted from the left wall
           CASE (eTag_Coll_Neutral)                                         !
              sece_right_after_coll_count = sece_right_after_coll_count + 1 ! if primary electron collided with neutral atom
        END SELECT
        Current_electron%Tag = eTag_Emit_Right                              ! mark the electron
        CALL ADD_EMITTED_E_TO_RIGHT_DF(vx, vy, vz)
     END IF
     Current_electron%X       = x_new
     Current_electron%VX      = vx
     Current_electron%VY      = vy
     Current_electron%VZ      = vz
!     Current_electron%AX      = 0.0_8
     ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in INJECT_TRUE_SECONDARY :")', Rank_of_process
        PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
     Current_electron => Current_electron%next
     NULLIFY(Current_electron%next)
  ELSE
     PRINT '(/2x,"Process ",i3," : Error in INJECT_TRUE_SECONDARY :")', Rank_of_process
     PRINT  '(2x,"Current_electron is NOT associated!")'
     PRINT  '(2x,"The program will be terminated now :(")'
     STOP
  END IF

!! account for the density of emitted particle
!  Q_strm_spec(left_node, 1)  = Q_strm_spec(left_node, 1)  + right_node - x_new    ! s = 1
!  Q_strm_spec(right_node, 1) = Q_strm_spec(right_node, 1) + x_new - left_node     ! s = 1

END SUBROUTINE INJECT_TRUE_SECONDARY

!-------------------------------------------------------------------------------------------------
!
REAL(8) FUNCTION Coeff_SEE_Elastic(energy, teta)

  USE SEEmission
  IMPLICIT NONE

  REAL(8) energy       ! dim-less energy of incident electron
  REAL(8) teta         ! angle of incidence 

  REAL(8) coeff_max_teta
  REAL(8) Coeff_SEE_Classic

  SELECT CASE (Emitted_model(1))

     CASE (0)
        Coeff_SEE_Elastic = 0.0_8

     CASE (1)
        IF ((energy.GT.minE_see_elastic).AND.(energy.LE.maxE_see_elastic)) THEN
           Coeff_SEE_Elastic = setD_see_elastic
        END IF

     CASE (2) 
        IF (energy.GE.E_elast_max) THEN
           coeff_max_teta    = maxD_elast !!!* (1.0_8 + 0.159154943_8 * k_smooth * teta**2) ! 0.159154943 = 1 / 2*pi
           Coeff_SEE_Elastic = coeff_max_teta * &
                             & EXP( - (energy - E_elast_max) / dE_elast) * &
                             & (energy - E_elast_max + dE_elast) / dE_elast + &
                             & Frac_elast_highenergy * Coeff_SEE_Classic(energy, teta)

        ELSE IF (energy.GT.E_elast_0) THEN
           coeff_max_teta    = maxD_elast !!!* (1.0_8 + 0.159154943_8 * k_smooth * teta**2) ! 0.159154943 = 1 / 2*pi         
           Coeff_SEE_Elastic = coeff_max_teta * & 
                             & EXP(-(energy - E_elast_max) / (E_elast_max - E_elast_0) ) * &
                             & (energy - E_elast_0) / (E_elast_max - E_elast_0) + & 
                             & Frac_elast_highenergy * Coeff_SEE_Classic(energy, teta)

        ELSE
           Coeff_SEE_Elastic = 0.0_8

        END IF
        
  END SELECT

END FUNCTION Coeff_SEE_Elastic

!-------------------------------------------------------------------------------------------------
!
REAL(8) FUNCTION Coeff_SEE_Inelastic(energy, teta)

  USE SEEmission
  IMPLICIT NONE

  REAL(8) energy       ! dim-less energy of incident electron
  REAL(8) teta         ! angle of incidence 

  REAL(8) Coeff_SEE_Classic

  SELECT CASE (Emitted_model(2))

     CASE (0)
        Coeff_SEE_Inelastic = 0.0_8

     CASE (1)
        IF ((energy.GT.minE_see_inelastic).AND.(energy.LE.maxE_see_inelastic)) THEN
           Coeff_SEE_Inelastic = setD_see_inelastic
        END IF

     CASE (2) 
        Coeff_SEE_Inelastic = Frac_inelastic * Coeff_SEE_Classic(energy, teta)
        
  END SELECT

END FUNCTION Coeff_SEE_Inelastic


!-------------------------------------------------------------------------------------------------
!
REAL(8) FUNCTION Coeff_SEE_True(energy, teta)

  USE SEEmission
  IMPLICIT NONE

  REAL(8) energy       ! dim-less energy of incident electron
  REAL(8) teta         ! angle of incidence, with respect to the surface normal  (0 = normal incidence)

  REAL(8) Coeff_SEE_Classic

  SELECT CASE (Emitted_model(3))

     CASE (0)
        Coeff_SEE_True = 0.0_8

     CASE (1)
        IF ((energy.GT.minE_see_true).AND.(energy.LE.maxE_see_true)) THEN
           Coeff_SEE_True = setD_see_true
        END IF

     CASE (2) 
        Coeff_SEE_True = Coeff_SEE_Classic(energy, teta) * (1.0_8 - Frac_elast_highenergy - Frac_inelastic)
        IF (Coeff_SEE_True.LT.0.0_8) Coeff_SEE_True = 0.0_8        
        
  END SELECT

END FUNCTION Coeff_SEE_True

!-------------------------------------------------------------------------------------------------
!
REAL(8) FUNCTION Coeff_SEE_Classic(energy, teta)

  USE SEEmission
  IMPLICIT NONE

  REAL(8) energy       ! dim-less energy of incident electron
  REAL(8) teta         ! angle of incidence, with respect to the surface normal  (0 = normal incidence)
  REAL(8) v
  REAL k
  REAL(8) energy_max_teta
  REAL(8) coeff_max_teta

  IF (energy.GT.E_see_0) THEN

     energy_max_teta = E_see_max * (1.0_8 + 0.318309886_8 * k_smooth * teta**2)     ! 0.318309886 = 1 / pi
     v               = (energy - E_see_0) / (energy_max_teta - E_see_0) 
     coeff_max_teta  = maxD_see_classic * (1.0_8 + 0.159154943_8 * k_smooth * teta**2) ! 0.159154943 = 1 / 2*pi
     
     IF (v.LT.1.0_8) THEN 
        k = 0.62
     ELSE
        k = 0.25
     END IF
     
     Coeff_SEE_Classic = coeff_max_teta * (v * EXP(1.0_8 - v)) ** k
     
  ELSE
     
     Coeff_SEE_Classic = 0.0_8

  END IF


END FUNCTION Coeff_SEE_Classic

!*******************************************************************************************
! This module is used by two subroutines, calculating the arbitrary energy
! according to the maxwell energy distribution
MODULE MaxwellEnergy
  REAL(8) E(0:200)
  REAL(8) :: E_max = 8.0_8
  INTEGER :: R_max = 200
END MODULE MaxwellEnergy

!-------------------------------------------------------------------------------------------
! Prepares the tabulated values of integral of the maxwell distribution function
!  
SUBROUTINE PrepareEnergyDistribIntegral

  USE ParallelOperationValues
  USE MaxwellEnergy
  USE SEEmission, ONLY : T_see_true_eV
  USE CurrentProblemValues, ONLY : T_e_eV, N_box_vel
  IMPLICIT NONE

  INTEGER i
  INTEGER :: N_pnts = 20000
  INTEGER count
  REAL(8) F(0:20003)            ! to be sure that we overcome E_max
  REAL(8) temp
  REAL(8) dE
  
  LOGICAL :: check = .FALSE.

  dE = E_max / N_pnts 

  F = 0.0_8 
  DO i = 1, N_pnts + 3
     temp = (REAL(i) - 0.5_8) * dE
!     F(i) = F(i-1) + temp * EXP( -temp )     
     F(i) = F(i-1) + SQRT(temp) * EXP( -temp )     
  END DO

  temp = F(N_pnts)
  F = F * R_max / temp   ! normalize integral such that F(N_pnts) = R_max

  E(0) = 0.0_8
  count = 0
  DO i = 1, N_pnts + 3
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        E(count) = i * dE
        IF (count.EQ.R_max) THEN
           check = .TRUE.
           EXIT
        END IF
     END IF
  END DO

  E = E * 0.5_8 * T_see_true_eV * N_box_vel**2 / T_e_eV    !??????????????????????????????????????????????????????????

  IF (check) THEN
!     PRINT &
!& '(2x,"Process ",i3," : Integral for producing energy distribution of secondary electrons is successfully obtained ...")', &
!                                                                                                            & Rank_of_process
  ELSE
     PRINT '(2x,"Process ",i3," : ERROR in PrepareEnergyDistribIntegral !!!")', Rank_of_process
     PRINT '(2x,"The initialization in PrepareMaxwellDistribIntegral is not performed !!!")'
     PRINT '(2x,"The program will be terminated now :(")'
     STOP
  END IF

END SUBROUTINE PrepareEnergyDistribIntegral

!-------------------------------------------------------------------------------------------
!  
SUBROUTINE GetTrueSecondaryEnergy(W) 

  USE MaxwellEnergy
  USE SEEmission, ONLY : T_see_true_eV
  USE CurrentProblemValues, ONLY : T_e_Ev

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) W

  REAL(8) R
  INTEGER indx 
  
  R = R_max * well_random_number()

  indx = INT(R)

  IF (indx.LT.R_max) THEN
     W = E(indx) + (R - indx) * (E(indx+1) - E(indx))
  ELSE
     W = E(R_max)
  END IF
  RETURN
  
END SUBROUTINE GetTrueSecondaryEnergy

