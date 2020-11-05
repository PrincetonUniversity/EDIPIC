
!===================================================================================================
SUBROUTINE INITIATE_MC_COLLISIONS

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_max_vel, T_e_eV, T_i_eV, N_spec, W_cycl_x_s1, Ms, N_box_vel

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL exists
  INTEGER s, k

  CHARACTER (67) buf
  INTEGER ierr

  maxNcolkind_spec(1) = 4            ! Currently we have 4 kinds of e-n collisions         !@#$
  maxNcolkind_spec(2) = 3            ! Currently we have 3 kinds of i-n collisions         !@#$
 
  INQUIRE (FILE = 'ssc_partcolls.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  OPEN (9, FILE = 'ssc_partcolls.dat')

  IF(exists) THEN

     PRINT '(2x,"Process ",i3," : ssc_partcolls.dat file is found. Reading the data file...")', Rank_of_process

     READ (9, '(A67)') buf !================= NEUTRAL COMPONENT PARAMETERS ====================")')
     READ (9, '(A67)') buf !--dddddd.ddd- Mass (a.m.u.) ---------------------------------------")')
     READ (9, '(2x,f10.3)') M_neutral_amu 
     READ (9, '(A67)') buf !--#d.dddE#dd- Density (m^-3) --------------------------------------")')
     READ (9, '(2x,e10.3)') N_neutral_m3 
     READ (9, '(A67)') buf !--dddddd.ddd- Temperature (eV) ------------------------------------")')
     READ (9, '(2x,f10.3)') T_neutral_eV     
     READ (9, '(A67)') buf !============ ELECTRON - NEUTRAL COLLISIONS, ACTIVATION ============")')
     READ (9, '(A67)') buf !-------d----- Elastic-1 (1 = yes, 0 = no) -------------------------")')
     READ (9, '(7x,i1)') Colflag_kind_spec(1,1)
     READ (9, '(A67)') buf !-------d----- Excitation-1 (1 = yes, 0 = no) ----------------------")')
     READ (9, '(7x,i1)') Colflag_kind_spec(2,1)
     READ (9, '(A67)') buf !-------d----- Ionization-1 (1 = yes, 0 = no) ----------------------")')
     READ (9, '(7x,i1)') Colflag_kind_spec(3,1)
     READ (9, '(A67)') buf !============== ELECTRON - NEUTRAL COLLISIONS, CONTROL =============")')
     READ (9, '(A67)') buf !--dddddd.ddd- Maximal electron energy (eV), default if negative ---")')
     READ (9, '(2x,f10.3)') maxEnergy_eV_spec(1) 
     READ (9, '(A67)') buf !--dddddd----- Number of energy values (>0) ------------------------")')
     READ (9, '(2x,i6)') Nenergyval_spec(1)
     READ (9, '(A67)') buf !========== ELECTRON - TURBULENCE, ACTIVATION and CONTROL ==========")')         !@#$
     READ (9, '(A67)') buf !-------d----- Turbulence-1 (1 = yes, 0 = no) ----------------------")')         !@#$
     READ (9, '(7x,i1)') Colflag_kind_spec(4,1)                                                            !@#$
     READ (9, '(A67)') buf !--#d.dddE#dd- Frequency, model-1 (s^-1) ---------------------------")')         !@#$      
     READ (9, '(2x,e10.3)') Freq_turb_e_1_s1                                                                  !@#$
     READ (9, '(A67)') buf !============== ION - NEUTRAL COLLISIONS, ACTIVATION ===============")')
     READ (9, '(A67)') buf !-------d----- Elastic-1 (1 = yes, 0 = no) -------------------------")')
     READ (9, '(7x,i1)') Colflag_kind_spec(1,2)
     READ (9, '(A67)') buf !-------d----- Charge exchange-1 (1 = yes, 0 = no) -----------------")')
     READ (9, '(7x,i1)') Colflag_kind_spec(2,2)
     READ (9, '(A67)') buf !================ ION - NEUTRAL COLLISIONS, CONTROL ================")')
     READ (9, '(A67)') buf !--dddddd.ddd- Maximal ion energy (eV), default if negative --------")')
     READ (9, '(2x,f10.3)') maxEnergy_eV_spec(2)
     READ (9, '(A67)') buf !--dddddd----- Number of energy values (>0) ------------------------")')
     READ (9, '(2x,i6)') Nenergyval_spec(2)
     READ (9, '(A67)') buf !============ ION - TURBULENCE, ACTIVATION and CONTROL =============")')         !@#$
     READ (9, '(A67)') buf !-------d----- Turbulence-1 (1 = yes, 0 = no) ----------------------")')         !@#$
     READ (9, '(7x,i1)') Colflag_kind_spec(3,2)                                                            !@#$
     READ (9, '(A67)') buf !--#d.dddE#dd- Frequency, model-1 (s^-1) ---------------------------")')         !@#$      
     READ (9, '(2x,e10.3)') Freq_turb_i_1_s1                                                                  !@#$
     
  ELSE

     M_neutral_amu          = 131.0_8
     N_neutral_m3           = 2.0e18
     T_neutral_eV           = 0.1_8 
     Colflag_kind_spec(1,1) = 1           ! e-n, elastic-1,       yes
     Colflag_kind_spec(2,1) = 1           ! e-n, excitation-1,    yes
     Colflag_kind_spec(3,1) = 1           ! e-n, ionization-1,    yes
     maxEnergy_eV_spec(1)   = 1000.0_8    ! 
     Nenergyval_spec(1)     = 20001
     Colflag_kind_spec(4,1) = 1           ! e-turbulence-1,       no                                !@#$
     Freq_turb_e_1_s1       = -80.0_8                                                                 !@#$
     Colflag_kind_spec(1,2) = 0           ! i-n, elastic-1,         no
     Colflag_kind_spec(2,2) = 0           ! i-n, charge exchange-1, no
     maxEnergy_eV_spec(2)   = -1.0_8      ! default value will be calculated
     Nenergyval_spec(2)     = 101
     Colflag_kind_spec(3,2) = 0           ! i-turbulence-1,       no                                !@#$
     Freq_turb_i_1_s1       = 0.0_8                                                                 !@#$
     Colflag_kind_spec(4,2) = 0           ! empty                                                   !@#$

     PRINT '(2x,"Process ",i3," : File with the name ssc_partcolls.dat not found. Use the default settings ...")', &
                                                                                                    & Rank_of_process

     IF (Rank_of_process.EQ.0) THEN

        PRINT '(/2x,"Process ",i3," : Create ssc_partcolls.dat file . . .")', Rank_of_process

        WRITE (9, '("================= NEUTRAL COMPONENT PARAMETERS ====================")')
        WRITE (9, '("--dddddd.ddd- Mass (a.m.u.) ---------------------------------------")')
        WRITE (9, '(2x,f10.3)') M_neutral_amu 
        WRITE (9, '("--#d.dddE#dd- Density (m^-3) --------------------------------------")')
        WRITE (9, '(2x,e10.3)') N_neutral_m3 
        WRITE (9, '("--dddddd.ddd- Temperature (eV) ------------------------------------")')
        WRITE (9, '(2x,f10.3)') T_neutral_eV     
        WRITE (9, '("============ ELECTRON - NEUTRAL COLLISIONS, ACTIVATION ============")')
        WRITE (9, '("-------d----- Elastic-1 (1 = yes, 0 = no) -------------------------")')
        WRITE (9, '(7x,i1)') Colflag_kind_spec(1,1)
        WRITE (9, '("-------d----- Excitation-1 (1 = yes, 0 = no) ----------------------")')
        WRITE (9, '(7x,i1)') Colflag_kind_spec(2,1)
        WRITE (9, '("-------d----- Ionization-1 (1 = yes, 0 = no) ----------------------")')
        WRITE (9, '(7x,i1)') Colflag_kind_spec(3,1)
        WRITE (9, '("============== ELECTRON - NEUTRAL COLLISIONS, CONTROL =============")')
        WRITE (9, '("--dddddd.ddd- Maximal electron energy (eV), default if negative ---")')
        WRITE (9, '(2x,f10.3)') maxEnergy_eV_spec(1) 
        WRITE (9, '("--dddddd----- Number of energy values (>0) ------------------------")')
        WRITE (9, '(2x,i6)') Nenergyval_spec(1)
        WRITE (9, '("========== ELECTRON - TURBULENCE, ACTIVATION and CONTROL ==========")')         !@#$
        WRITE (9, '("-------d----- Turbulence-1 (1 = yes, 0 = no) ----------------------")')         !@#$
        WRITE (9, '(7x,i1)') Colflag_kind_spec(4,1)                                                            !@#$
        WRITE (9, '("--#d.dddE#dd- Frequency, model-1 (s^-1) ---------------------------")')         !@#$      
        WRITE (9, '(2x,e10.3)') Freq_turb_e_1_s1                                                                  !@#$
        WRITE (9, '("============== ION - NEUTRAL COLLISIONS, ACTIVATION ===============")')
        WRITE (9, '("-------d----- Elastic-1 (1 = yes, 0 = no) -------------------------")')
        WRITE (9, '(7x,i1)') Colflag_kind_spec(1,2)
        WRITE (9, '("-------d----- Charge exchange-1 (1 = yes, 0 = no) -----------------")')
        WRITE (9, '(7x,i1)') Colflag_kind_spec(2,2)
        WRITE (9, '("================ ION - NEUTRAL COLLISIONS, CONTROL ================")')
        WRITE (9, '("--dddddd.ddd- Maximal ion energy (eV), default if negative --------")')
        WRITE (9, '(2x,f10.3)') maxEnergy_eV_spec(2)
        WRITE (9, '("--dddddd----- Number of energy values (>0) ------------------------")')
        WRITE (9, '(2x,i6)') Nenergyval_spec(2)
        WRITE (9, '("============ ION - TURBULENCE, ACTIVATION and CONTROL =============")')         !@#$
        WRITE (9, '("-------d----- Turbulence-1 (1 = yes, 0 = no) ----------------------")')         !@#$
        WRITE (9, '(7x,i1)') Colflag_kind_spec(3,2)                                                            !@#$
        WRITE (9, '("--#d.dddE#dd- Frequency, model-1 (s^-1) ---------------------------")')         !@#$      
        WRITE (9, '(2x,e10.3)') Freq_turb_i_1_s1                                                                  !@#$
     END IF

  END IF

  CLOSE (9, STATUS = 'KEEP')

  IF (maxEnergy_eV_spec(1).LT.0.0_8) THEN                                 ! 
     maxEnergy_eV_spec(1) = 2.0_8 * N_max_vel * N_max_vel * T_e_eV        ! default value of the maximal energy for electrons
  END IF                                                                  !

  IF (maxEnergy_eV_spec(2).LT.0.0_8) THEN                                 ! 
     maxEnergy_eV_spec(2) = 2.0_8 * N_max_vel * N_max_vel * T_i_eV        ! default value of the maximal energy for ions
  END IF                                                                  !
  
  IF (Freq_turb_e_1_s1.LT.0.0_8) THEN                                     ! frequency of electron turbulence collisions is defined
     Freq_turb_e_1_s1 = ABS(W_cycl_x_s1 / Freq_turb_e_1_s1)                 ! by the electron cyclotron frequency
  END IF                                                                  ! (if the input data is negative)
!print *, '1 <- ', Rank_of_process
  CALL CONFIG_READ_CRSECT_ARRAYS
!print *, '2 <- ', Rank_of_process
  CALL CONFIGURE_COLLSPEC_ARRAYS
  CALL SETVALUES_COLLSPEC_ARRAYS
!print *, '3 <- ', Rank_of_process
  CALL CONFIGURE_COLLISION_ARRAYS
  CALL SETVALUES_COLLISION_ARRAYS
!print *, '4 <- ', Rank_of_process

!####  CALL PrepareMaxwellDistribIntegral   !### NEW :: moved this call to INITIATE_PARAMETERS to use GetMaxwellVelocity

  IF (Rank_of_process.EQ.0) THEN

     DO s = 1, N_spec
        SELECT CASE (s)
           CASE (1)
              PRINT '(/2x,"Activated electron-neutral/turbulence collisions:")'
              IF (Ncolkind_spec(s).EQ.0) THEN
                 PRINT '(4x,"None :)")'
                 CYCLE
              END IF
              DO k = 1, Ncolkind_spec(s)
                 SELECT CASE (Colkind_of_spec(s)%activated(k))
                    CASE(1)  
                       PRINT '(4x,"Elastic, model 1")'
                    CASE(2)
                       PRINT '(4x,"Excitation, model 1")'
                    CASE(3)
                       PRINT '(4x,"Ionization, model 1")'
                    CASE(4)                                                                                                       !@#$
                       PRINT '(4x,"Turbulence, model 1, frequency is ",e10.3," s^-1")', Freq_turb_e_1_s1                          !@#$
                 END SELECT
              END DO

           CASE (2)
              PRINT '(/2x,"Activated ion-neutral/turbulence collisions:")'
              IF (Ncolkind_spec(s).EQ.0) THEN
                 PRINT '(4x,"None :)")'
                 CYCLE
              END IF
              DO k = 1, Ncolkind_spec(s)
                 SELECT CASE (Colkind_of_spec(s)%activated(k))
                    CASE(1)
                       PRINT '(4x,"Elastic, model 1")'
                    CASE(2)
                       PRINT '(4x,"Charge exchange, model 1")'
                    CASE(3)                                                                                                       !@#$
                       PRINT '(4x,"Turbulence, model 1, frequency is ",e10.3," s^-1")', Freq_turb_i_1_s1                          !@#$
                  END SELECT
               END DO

        END SELECT
        PRINT '(/2x,"Colliding portion (%): ",f8.5)', 100.0 * maxColfrac_spec(s)
     END DO

  END IF

  e_n_1_count = 0; e_n_2_count = 0; e_n_3_count = 0; i_n_1_count = 0; i_n_2_count = 0
  e_t_4_count = 0; i_t_3_count = 0                                                              !@#$
! 0.00054858 = m_e_kg / 1_amu_kg = 9.109534e-31 / 1.660565e-27 
  alpha_Vscl = SQRT(0.00054858_8 * T_neutral_eV / (T_e_eV * M_neutral_amu)) 
  
  DO s = 1, N_spec
     energy_factor_eV_s(s) = T_e_eV * Ms(s) / N_box_vel**2
  END DO
   
!  PRINT '(2x,"Finished :)")'

END SUBROUTINE INITIATE_MC_COLLISIONS 

!=========================================================================================================
! Note, this subroutine must be called before any first call of the FREQUENCY_OF_COLLISION function
!
SUBROUTINE CONFIG_READ_CRSECT_ARRAYS 

  USE ParallelOperationValues
  USE MCCollisions

  IMPLICIT NONE

  INTEGER ALLOC_ERR
  LOGICAL exists
  INTEGER j

! read cross sections of e-n elastic collisions from data file, if this kind of collisions is activated   
  IF (Colflag_kind_spec(1, 1).EQ.1) THEN       

     INQUIRE (FILE = 'ssc_crsect_en_elast.dat', EXIST = exists)
     
     IF (exists) THEN
     
        OPEN (9, FILE = 'ssc_crsect_en_elast.dat')

        PRINT '(2x,"Process ",i3," : e-n elastic collisions cross-sections data file is found. Reading the data file...")', &
                                                                                                              &Rank_of_process
        READ (9, '(2x,i4)') N_en_elast

        ALLOCATE(Energy_en_elast_eV(1:N_en_elast), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Energy_en_elast_eV !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF

        ALLOCATE(CrSect_en_elast_m2(1:N_en_elast), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE CrSect_en_elast_m2 !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        
        DO j = 1, N_en_elast
! was      READ (9, '(4x,f8.2,3x,e9.2)') Energy_en_elast_eV(j), CrSect_en_elast_m2(j)
           READ (9, '(4x,f9.3,2x,e10.3)') Energy_en_elast_eV(j), CrSect_en_elast_m2(j)
!           print '(4x,f8.2,3x,e9.2)', Energy_en_elast_eV(j), CrSect_en_elast_m2(j)
        END DO

        CLOSE (9, STATUS = 'KEEP')
     
     ELSE

        PRINT '(/2x,"Process ",i3," : Error! Cannot find file ssc_crsect_en_elast.dat")', Rank_of_process      
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
        
     END IF
  
  END IF                                        

! read cross sections of e-n excitation collisions from data file, if this kind of collisions is activated   
  IF (Colflag_kind_spec(2, 1).EQ.1) THEN           

     INQUIRE (FILE = 'ssc_crsect_en_excit.dat', EXIST = exists)
     
     IF (exists) THEN
     
        OPEN (9, FILE = 'ssc_crsect_en_excit.dat')

        PRINT '(2x,"Process ",i3," : e-n excitation collisions cross-sections data file is found. Reading the data file...")',&
                                                                                                               & Rank_of_process
        READ (9, '(2x,i4)') N_en_excit

        ALLOCATE(Energy_en_excit_eV(1:N_en_excit), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Energy_en_excit_eV !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF

        ALLOCATE(CrSect_en_excit_m2(1:N_en_excit), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE CrSect_en_excit_m2 !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        
        DO j = 1, N_en_excit
! was      READ (9, '(4x,f8.2,3x,e9.2)') Energy_en_excit_eV(j), CrSect_en_excit_m2(j)
           READ (9, '(4x,f8.2,3x,e10.3)') Energy_en_excit_eV(j), CrSect_en_excit_m2(j)
        END DO
        
        CLOSE (9, STATUS = 'KEEP')
        
        DO j = 1, N_en_excit                                 !
           IF (CrSect_en_excit_m2(j).GT.0.0) THEN            !
              Thresh_en_excit_eV = Energy_en_excit_eV(j)     ! finding the energy threshold value
              EXIT                                           !
           END IF                                            !
        END DO                                               !
     
     ELSE

        PRINT '(/2x,"Process ",i3," : Error! Cannot find file ssc_crsect_en_excit.dat !!!")', Rank_of_process      
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
        
     END IF
  
  END IF

! read cross sections of e-n ionization collisions from data file, if this kind of collisions is activated   
  IF (Colflag_kind_spec(3, 1).EQ.1) THEN           

     INQUIRE (FILE = 'ssc_crsect_en_ioniz.dat', EXIST = exists)
     
     IF (exists) THEN
     
        OPEN (9, FILE = 'ssc_crsect_en_ioniz.dat')

        PRINT '(2x,"Process ",i3," : e-n ionization collisions cross-sections data file is found. Reading the data file...")',&
                                                                                                              &  Rank_of_process
        READ (9, '(2x,i4)') N_en_ioniz

        ALLOCATE(Energy_en_ioniz_eV(1:N_en_ioniz), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Energy_en_ioniz_eV !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF

        ALLOCATE(CrSect_en_ioniz_m2(1:N_en_ioniz), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE CrSect_en_ioniz_m2 !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        
        DO j = 1, N_en_ioniz
! was      READ (9, '(4x,f7.1,4x,e10.3)') Energy_en_ioniz_eV(j), CrSect_en_ioniz_m2(j)
           READ (9, '(4x,f8.2,3x,e10.3)') Energy_en_ioniz_eV(j), CrSect_en_ioniz_m2(j)
        END DO
        
        CLOSE (9, STATUS = 'KEEP')
        
        DO j = 1, N_en_ioniz                                   !
           IF (CrSect_en_ioniz_m2(j).GT.0.0) THEN              ! 
              Thresh_en_ioniz_eV = Energy_en_ioniz_eV(j)       ! find the energy threshold value
              EXIT                                             !
           END IF                                              !
        END DO                                                 ! 
     
     ELSE

        PRINT '(/2x,"Process ",i3," : Error! Cannot find file ssc_crsect_en_ioniz.dat !!!")', Rank_of_process      
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
        
     END IF
  
  END IF

END SUBROUTINE CONFIG_READ_CRSECT_ARRAYS 

!=========================================================================================================
SUBROUTINE REMOVE_CRSECT_ARRAYS

  USE ParallelOperationValues
  USE MCCollisions
  IMPLICIT NONE
  INTEGER DEALLOC_ERR

  IF (ALLOCATED(Energy_en_elast_eV)) THEN
     DEALLOCATE(Energy_en_elast_eV, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Energy_en_elast_eV !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(Energy_en_excit_eV)) THEN
     DEALLOCATE(Energy_en_excit_eV, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Energy_en_excit_eV !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(Energy_en_ioniz_eV)) THEN
     DEALLOCATE(Energy_en_ioniz_eV, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Energy_en_ioniz_eV !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(CrSect_en_elast_m2)) THEN
     DEALLOCATE(CrSect_en_elast_m2, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE CrSect_en_elast_m2 !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(CrSect_en_excit_m2)) THEN
     DEALLOCATE(CrSect_en_excit_m2, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE CrSect_en_excit_m2 !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(CrSect_en_ioniz_m2)) THEN
     DEALLOCATE(CrSect_en_ioniz_m2, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE CrSect_en_ioniz_m2 !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

END SUBROUTINE REMOVE_CRSECT_ARRAYS

!=========================================================================================================
! Note, this subroutine is called in the INITIATE_MC_COLLISIONS, before calling CONFIGURE_COLLISION_ARRAYS 
!
SUBROUTINE CONFIGURE_COLLSPEC_ARRAYS 

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_spec
  IMPLICIT NONE
  INTEGER ALLOC_ERR

!  PRINT '(4x,"Allocating species arrays for electron-neutral and ion-neutral collisions ...")'

  IF (.NOT.ALLOCATED(Ncolkind_spec)) THEN
     ALLOCATE(Ncolkind_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Ncolkind_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (.NOT.ALLOCATED(deltaEnergy_eV_spec)) THEN
     ALLOCATE(deltaEnergy_eV_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE deltaEnergy_eV_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE CONFIGURE_COLLSPEC_ARRAYS 

!=========================================================================================================
! Note, this subroutine is called in the INITIATE_MC_COLLISIONS, before calling CONFIGURE_COLLISION_ARRAYS 
!
SUBROUTINE SETVALUES_COLLSPEC_ARRAYS 

  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_spec, T_e_eV

  IMPLICIT NONE
!  INTEGER ALLOC_ERR
  INTEGER s          ! species
!  INTEGER i          ! group
  INTEGER k          ! kind of collision

!  PRINT '(4x,"Allocating species arrays for electron-neutral and ion-neutral collisions ...")'

  DO s = 1, N_spec
     deltaEnergy_eV_spec(s) = maxEnergy_eV_spec(s) / (Nenergyval_spec(s) - 1)
  END DO

! determine the number of activated kinds of collisions for different species  
  Ncolkind_spec = 0                                   
  DO s = 1, N_spec                                   
     DO k = 1, maxNcolkind_spec(s)                    
        IF (Colflag_kind_spec(k, s).EQ.1) THEN       
           Ncolkind_spec(s) = Ncolkind_spec(s) + 1   
        END IF                                        
     END DO                                       
  END DO

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE SETVALUES_COLLSPEC_ARRAYS 

!===================================================================================================
! 
SUBROUTINE REMOVE_COLLSPEC_ARRAYS 

  USE ParallelOperationValues
  USE MCCollisions
  IMPLICIT NONE
  INTEGER DEALLOC_ERR

!  PRINT '(4x,"Allocating species arrays for electron-neutral and ion-neutral collisions ...")'

  IF (ALLOCATED(Ncolkind_spec)) THEN
     DEALLOCATE(Ncolkind_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Ncolkind_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(deltaEnergy_eV_spec)) THEN
     DEALLOCATE(deltaEnergy_eV_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE deltaEnergy_eV_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE REMOVE_COLLSPEC_ARRAYS 

!=======================================================================================================
! Note, this subroutine is called in the INITIATE_MC_COLLISIONS, after calling SETVALUES_COLLSPEC_ARRAYS  
!
SUBROUTINE CONFIGURE_COLLISION_ARRAYS 

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_spec
  IMPLICIT NONE
  INTEGER ALLOC_ERR
  INTEGER s           ! species

!  PRINT '(4x,"Allocating arrays for electron-neutral and ion-neutral collisions ...")'

  ALLOCATE(maxColfrac_spec(1:N_spec), STAT=ALLOC_ERR)
  ALLOCATE(Colkind_of_spec(1:N_spec), STAT=ALLOC_ERR)

  DO s = 1, N_spec
     NULLIFY(Colkind_of_spec(s)%activated)
     IF (.NOT.ASSOCIATED(Colkind_of_spec(s)%activated)) THEN
        ALLOCATE(Colkind_of_spec(s)%activated(1:Ncolkind_spec(s)), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Colkind_of_spec(s)%activated !!! s = ",i1)', Rank_of_process, s
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

  ALLOCATE(Colprob_of_spec(1:N_spec), STAT=ALLOC_ERR)

  DO s = 1, N_spec
     NULLIFY(Colprob_of_spec(s)%kind_energy)
     IF (.NOT.ASSOCIATED(Colprob_of_spec(s)%kind_energy)) THEN
        ALLOCATE(Colprob_of_spec(s)%kind_energy(1:Ncolkind_spec(s),1:Nenergyval_spec(s)), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Colprob_of_spec(s)%kind_energy !!! s = ",i1)', Rank_of_process, s
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

  ALLOCATE(Numbers_of_particles(0:(N_of_processes-1)), STAT=ALLOC_ERR)      !^%
  ALLOCATE(Numbers_of_collisions(0:(N_of_processes-1)), STAT=ALLOC_ERR)     !^%
  

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE CONFIGURE_COLLISION_ARRAYS 
 
!===================================================================================================
! Note, this subroutine is called in the INITIATE_PARAMETERS  
SUBROUTINE SETVALUES_COLLISION_ARRAYS 

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_spec, delta_t_s
  IMPLICIT NONE

  REAL(8) FREQUENCY_OF_COLLISION

  INTEGER s   ! species
  INTEGER k   ! index for activated collisions
  INTEGER j   ! index for all collisions, index for energy 
  REAL energy_eV, max_f !, delta_t

!  PRINT '(4x,"Setting the values of arrays for electron-neutral and ion-neutral collisions ...")'

  DO s = 1, N_spec    

! calculate the array of indexes of active collisions Colkind_of_spec(1:N_spec)%activated(1:Ncolkind_spec(s))
     k = 0
     DO j = 1, maxNcolkind_spec(s)              
        IF(Colflag_kind_spec(j, s).EQ.1) THEN
           k = k + 1              
           IF (k.GT.Ncolkind_spec(s)) THEN
              PRINT '(2x,"Process ",i3," : Error in determining Colkind_of_spec(",i2,")%activated")', Rank_of_process, s
              PRINT '(2x,"Current array index ",i2," exceeds the upper limit ",i2)', k, Ncolkind_spec(s)
              PRINT '(2x,"Program will be terminated now :(")'
              STOP
           END IF
           Colkind_of_spec(s)%activated(k) = j     
        END IF                
     END DO  
! if the species has no activated collisions, take the next species
     IF (Ncolkind_spec(s).EQ.0) CYCLE

! calculate the array of probabilities Colprob_of_spec(1:N_spec)%kind_energy(1:Ncolkind_spec(s),1:Nenergyval_spec(s))
     max_f = 0.0_8
     DO j = 1, Nenergyval_spec(s)
        energy_eV = (j-1) * deltaEnergy_eV_spec(s)
        DO k = 1, Ncolkind_spec(s)
      
           Colprob_of_spec(s)%kind_energy(k, j) = FREQUENCY_OF_COLLISION(energy_eV, Colkind_of_spec(s)%activated(k), s)
           IF (k.GT.1) THEN
              Colprob_of_spec(s)%kind_energy(k, j) = Colprob_of_spec(s)%kind_energy(k, j) + & 
                                                   & Colprob_of_spec(s)%kind_energy(k-1, j)
           END IF

        END DO
        IF (Colprob_of_spec(s)%kind_energy(Ncolkind_spec(s), j).GT.max_f) THEN
           max_f = Colprob_of_spec(s)%kind_energy(Ncolkind_spec(s), j)
        END IF
     END DO
! renormalize the probability (must be not greater than 1)
     DO j = 1, Nenergyval_spec(s)
        DO k = 1, Ncolkind_spec(s)
           Colprob_of_spec(s)%kind_energy(k, j) = Colprob_of_spec(s)%kind_energy(k, j) / max_f
        END DO
!!!!!!!!print '(2x,i4,2x,6(1x,e12.5))', j, Colprob_of_spec(s)%kind_energy(1:Ncolkind_spec(s), j)
     END DO

! calculating the array of fractions of particles, colliding each time step (including the NULL collisions), 
! for different groups and species maxColfrac_spec(1:N_of_groups,1:N_spec)
     maxColfrac_spec(s) = 1.0 - EXP(- max_f)

  END DO

!  PRINT '(2x,"Finished :)")'

  IF (Rank_of_process.EQ.0) THEN                                              ! server node saves the collision frequencies

     IF (Colflag_kind_spec(1, 1).EQ.1) THEN                                   ! electron-neutral, elastic
        OPEN (99, FILE = '_en_elast_coll_freqs.dat')
        DO j = 1, Nenergyval_spec(1)
           energy_eV = (j-1) * deltaEnergy_eV_spec(1)
           WRITE (99, '(2(2x,e12.5))') energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 1, 1) / delta_t_s
        END DO
        CLOSE (99, STATUS = 'KEEP')
     END IF

     IF (Colflag_kind_spec(2, 1).EQ.1) THEN                                   ! electron-neutral, excitation
        OPEN (99, FILE = '_en_excit_coll_freqs.dat')
        DO j = 1, Nenergyval_spec(1)
           energy_eV = (j-1) * deltaEnergy_eV_spec(1)
           WRITE (99, '(2(2x,e12.5))') energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 2, 1) / delta_t_s
        END DO
        CLOSE (99, STATUS = 'KEEP')
     END IF

     IF (Colflag_kind_spec(3, 1).EQ.1) THEN                                   ! electron-neutral, ionization
        OPEN (99, FILE = '_en_ioniz_coll_freqs.dat')
        DO j = 1, Nenergyval_spec(1)
           energy_eV = (j-1) * deltaEnergy_eV_spec(1)
           WRITE (99, '(2(2x,e12.5))') energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 3, 1) / delta_t_s
!print '(2(2x,e12.5))', energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 3, 1) / delta_t_s 
        END DO
        CLOSE (99, STATUS = 'KEEP')
     END IF

     IF (Colflag_kind_spec(4, 1).EQ.1) THEN                                   ! electron, turbulent
        OPEN (99, FILE = '_e_turb_coll_freqs.dat')
        DO j = 1, Nenergyval_spec(1)
           energy_eV = (j-1) * deltaEnergy_eV_spec(1)
           WRITE (99, '(2(2x,e12.5))') energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 4, 1) / delta_t_s
        END DO
        CLOSE (99, STATUS = 'KEEP')
     END IF

  END IF

END SUBROUTINE SETVALUES_COLLISION_ARRAYS 

!===================================================================================================
SUBROUTINE REMOVE_COLLISION_ARRAYS

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_spec
  IMPLICIT NONE
  INTEGER DEALLOC_ERR
  INTEGER s    ! species

!  PRINT '("Deleting arrays for electron-neutral and ion-neutral collisions ...")'

  IF (ALLOCATED(maxColfrac_spec)) THEN
     DEALLOCATE(maxColfrac_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE maxColfrac_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  DO s = 1, N_spec
     IF (ASSOCIATED(Colkind_of_spec(s)%activated)) THEN
        DEALLOCATE(Colkind_of_spec(s)%activated, STAT=DEALLOC_ERR)
        IF(DEALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Colkind_of_spec(s)%activated !!! s = ",i1)', Rank_of_process, s
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

  IF (ALLOCATED(Colkind_of_spec)) THEN
     DEALLOCATE(Colkind_of_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Colkind_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  DO s = 1, N_spec
     IF (ASSOCIATED(Colprob_of_spec(s)%kind_energy)) THEN
        DEALLOCATE(Colprob_of_spec(s)%kind_energy, STAT=DEALLOC_ERR)
        IF(DEALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Colprob_of_spec(s)%kind_energy !!! s = ",i1)', Rank_of_process, s
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

  IF (ALLOCATED(Colprob_of_spec)) THEN
     DEALLOCATE(Colprob_of_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Colprob_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF


  IF (ALLOCATED(Numbers_of_particles)) THEN
     DEALLOCATE(Numbers_of_particles, STAT=DEALLOC_ERR)
     IF (DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Numbers_of_particles !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(Numbers_of_collisions)) THEN
     DEALLOCATE(Numbers_of_collisions, STAT=DEALLOC_ERR)
     IF (DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Numbers_of_collisions !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  NULLIFY(Collided_particle)   

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE REMOVE_COLLISION_ARRAYS

!==========================================================
! 
SUBROUTINE PROCESS_COLLISIONS_WITH_ATOMS

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues

  USE rng_wrapper

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTERFACE
     RECURSIVE SUBROUTINE Node_Killer(node)
       USE MCCollisions
       TYPE (binary_tree), POINTER :: node
     END SUBROUTINE Node_Killer
  END INTERFACE

  INTEGER s             ! species 
  INTEGER random_j      ! number of selected particle
  REAL    real_random_j ! real value used to avoid errors if the number becomes more than the largest 32-bit integer 
  INTEGER j             ! counter of attempts to collide

  REAL    R_collided ! number of collided particles during one timestep, can be less than 1
  INTEGER I_collided ! integer part of R_collided
  REAL    F_collided ! fractional part of R_collided

! one should definitely take I_collided particles as the collided, and 
! 1 particle with the probability F_collided

  INTEGER start_attempt ! start index in the cycle, which takes the random numbers of particles
                        ! it is 1 if F_collided = 0 and it is 0 if F_collided > 0

  REAL(8)    random_n   ! to account for the nonuniform neutral density
  REAL(8)    random_r   ! to identify the type of collision
  INTEGER k          ! counter of kinds of collisions 
 

  REAL(8) energy_eV     ! energy of selected particle [eV]

  INTEGER indx_energy   ! 2-nd index of the array Colprob_of_spec(s)%kind_energy(k, indx_energy)  
  INTEGER indx_coll     ! index of selected collision type

  REAL(8) Vx_n, Vy_n, Vz_n  ! velociity components of the colliding neutral atom

  LOGICAL search_again
  LOGICAL already_stored
  INTEGER ALLOC_ERR
  LOGICAL Find_in_stored_list

  INTEGER n                           ! process number
  INTEGER ierr

  INTEGER N_of_collisions        ! number of collisions for a client process

  REAL(8) a0, a1

! function
  REAL(8) Shape_of_Natom  ! defines shape of the neutral atom density profile, 0<=Shape_of_Natom<=1

  DO s = 1, N_spec

     IF (Ncolkind_spec(s).EQ.0) CYCLE                  ! skip if no collision for the species

     N_of_collisions = 0

     IF (Rank_of_process.EQ.0) THEN                                                              ! server >>>

        Numbers_of_particles = 0
        N_part = 0
! Receive data from clients
        CALL MPI_GATHER(N_part(s), 1, MPI_INTEGER, Numbers_of_particles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        DO n = 2, N_of_processes - 1
           Numbers_of_particles(n) = Numbers_of_particles(n) + Numbers_of_particles(n - 1) 
        END DO
        N_part(s) = Numbers_of_particles(N_of_processes - 1)      ! ??? was (n-1)

        R_collided = maxColfrac_spec(s) * N_part(s)      
        I_collided = INT(R_collided)
        F_collided = R_collided - I_collided

        IF (F_collided.GT.0.000001) THEN
           start_attempt = 0                 ! start with statistical collisions
        ELSE
           start_attempt = 1                 ! skip the statistical collisions if the fractional part is too small (to avoid overflow)
        END IF

!============= Randomly take some particle "I_collided(+1)" times, define number of collisions for each process
        Numbers_of_collisions = 0
        DO j = start_attempt, I_collided

           IF (j.EQ.0) THEN                                                     ! statistical part
              real_random_j = REAL(well_random_number() * N_part(s) / F_collided)       ! the value of the right-hand part can exceed the maximal 32-bit integer number
              IF ((real_random_j.GT.N_part(s)).OR.(real_random_j.LT.0.0)) CYCLE ! that is why we use the real value and check the negative value here too
              random_j = INT(real_random_j)                                          ! initialize the integer value
              IF (random_j.GT.N_part(s)) CYCLE                                  ! double check (maybe that's silly...)
           ELSE                                                                 ! mandatory part
              random_j = INT(well_random_number() * N_part(s))                         ! note, we assume that N_part(s) is a trusted 32-bit integer value
           END IF                                                           
           IF (random_j.LT.1) random_j = 1                                      ! to be sure that we are within the array boundaries
           
           DO n = 1, N_of_processes - 1
              IF (random_j.LE.Numbers_of_particles(n)) THEN 
                 Numbers_of_collisions(n) = Numbers_of_collisions(n) + 1
                 EXIT
              END IF
           END DO

        END DO

! TRANSMIT DATA TO THE CLIENTS
        CALL MPI_SCATTER(Numbers_of_collisions, 1, MPI_INTEGER, N_of_collisions, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     ELSE                                                                                                      ! client >>>

! send data to the server
        CALL MPI_GATHER(N_part(s), 1, MPI_INTEGER, Numbers_of_particles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! Allocate binary tree to store the numbers of particles which have collided already
        NULLIFY(Collided_particle)
        IF (.NOT.ASSOCIATED(Collided_particle)) THEN
           ALLOCATE(Collided_particle, STAT=ALLOC_ERR)
           IF (ALLOC_ERR.NE.0) THEN
              PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Collided_particle !!!")', Rank_of_process
              PRINT  '(2x,"Program will be terminated now :(")'
              STOP
           END IF
        END IF
        NULLIFY(Collided_particle%Larger)
        NULLIFY(Collided_particle%Smaller)

! receive data from the server
        CALL MPI_SCATTER(Numbers_of_collisions, 1, MPI_INTEGER, N_of_collisions, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!============= Randomly take some (always different) particle "I_collided(+1)" times
        DO j = 1, N_of_collisions

!------------- Determine the kind of collision for the selected particle
           random_r = well_random_number()
           random_n = well_random_number()

           search_again = .TRUE.                              ! if I_collided >= 1 then
           DO WHILE (search_again)                            ! search will be repeated until a number will be successfully obtained
              random_j = INT(well_random_number() * N_part(s))
              random_j = MIN(MAX(random_j, 1), N_part(s))
              already_stored = Find_in_stored_list(random_j)   ! 
              IF (already_stored) CYCLE                                  ! skip particle if it has collided already  
              search_again = .FALSE.
           END DO

! #### NEW: first check whether collision is allowed by the neutral density profile #####
           IF (random_n.LE.Shape_of_Natom(species(s)%part(random_j)%X)) THEN      

              SELECT CASE (s)
              CASE (1)

! if we are here, we can proceed with processing the colliding particle as usual           
! for electrons the kind of collision is determined by the energy of electron
                 energy_eV = energy_factor_eV_s(s) * (species(s)%part(random_j)%VX**2 + &
                                                   &  species(s)%part(random_j)%VY**2 + & 
                                                   &  species(s)%part(random_j)%VZ**2)
           
                 indx_energy = energy_eV / deltaEnergy_eV_spec(s) + 1

                 a1 = (energy_eV - (indx_energy-1) * deltaEnergy_eV_spec(s)) / deltaEnergy_eV_spec(s)
                 a0 = 1.0_8 - a1

                 IF (indx_energy.GE.Nenergyval_spec(s)) THEN
                    indx_energy = Nenergyval_spec(s)-1
                    a0 = 0.0_8
                    a1 = 1.0_8
                 END IF
              
                 DO k = Ncolkind_spec(s), 1, -1 
                    IF (random_r.GT. &
                         & (a0 * Colprob_of_spec(s)%kind_energy(k, indx_energy) + &
                         &  a1 * Colprob_of_spec(s)%kind_energy(k, indx_energy+1)) &
                         & ) EXIT
                 END DO

                 indx_coll = k + 1        
                 CALL COLLIDE_ELECTRON(indx_coll, random_j, energy_eV)

              CASE (2)   
! for ions the kind of collision is determined by the energy of relative motion ion-neutral
! that's why we must take the neutral's velocity from Maxwell distribution:
                 CALL GetMaxwellVelocity(Vx_n) 
                 CALL GetMaxwellVelocity(Vy_n)       ! #### FUTURE WORK - CHECK THIS FOR CORRELATIONS ########
                 CALL GetMaxwellVelocity(Vz_n) 
! Use the factor above to obtain the dim-less velocity (V/V_te) of the produced neutral atom
                 Vx_n = Vx_n * alpha_Vscl
                 Vy_n = Vy_n * alpha_Vscl
                 Vz_n = Vz_n * alpha_Vscl
! Take the energy of ion at the frame where the neutral atom is at rest
                 energy_eV = energy_factor_eV_s(s) * ((species(s)%part(random_j)%VX - Vx_n)**2 + &
                                                   &  (species(s)%part(random_j)%VY - Vy_n)**2 + & 
                                                   &  (species(s)%part(random_j)%VZ - Vz_n)**2)
                 
                 indx_energy = energy_eV / deltaEnergy_eV_spec(s) + 1

                 a1 = (energy_eV - (indx_energy-1) * deltaEnergy_eV_spec(s)) / deltaEnergy_eV_spec(s)
                 a0 = 1.0_8 - a1

                 IF (indx_energy.GE.Nenergyval_spec(s)) THEN
                    indx_energy = Nenergyval_spec(s)-1
                    a0 = 0.0_8
                    a1 = 1.0_8
                 END IF

                 DO k = Ncolkind_spec(s), 1, -1 
                    IF (random_r.GT. &
                         & (a0 * Colprob_of_spec(s)%kind_energy(k, indx_energy) + &
                         &  a1 * Colprob_of_spec(s)%kind_energy(k, indx_energy+1)) &
                         & ) EXIT
                 END DO

                 indx_coll = k + 1
                 CALL COLLIDE_ION(indx_coll, random_j, Vx_n, Vy_n, Vz_n) 

              END SELECT    !------------- End of determination of the kind of collision
           END IF

        END DO           !============= End of random particle taking

        CALL Node_Killer(Collided_particle)

     END IF

  END DO                 !------------- End of the external cycle over species

!print *, 'exit  PROCESS_COLLISIONS_WITH_ATOMS' ------------------------------------------------------------------------------------

END SUBROUTINE PROCESS_COLLISIONS_WITH_ATOMS

!-----------------------------------------------------------------
! function's value equals .TRUE. if the particle is already stored,
! otherwise function's value is .FALSE. (i.e. particle does not collide yet)
LOGICAL FUNCTION Find_in_stored_list(number)

  USE MCCollisions
  IMPLICIT NONE

  INTEGER number
  TYPE (binary_tree), POINTER :: current

  Find_in_stored_list = .FALSE.

  current => Collided_particle

  DO 

     IF (number.GT.current%number) THEN
        IF (ASSOCIATED(current%Larger)) THEN
           current => current%Larger               ! 
           CYCLE                                   ! go to the next node, with larger "number"
        ELSE
           EXIT
        END IF
     END IF

     IF (number.LT.current%number) THEN
        IF (ASSOCIATED(current%Smaller)) THEN
           current => current%Smaller              ! 
           CYCLE                                   ! go to the next node, with smaller "number"
        ELSE
           EXIT
        END IF
     END IF

     Find_in_stored_list = .TRUE.                  ! number.EQ.current%number
     EXIT                                          ! if we are here, then we found the match
     
  END DO

END FUNCTION Find_in_stored_list

!-----------------------------------------------------------------
! subroutine adds number to the binary tree
! we assume that there are no nodes in the tree with the same value yet
SUBROUTINE Add_to_stored_list(number)

  USE ParallelOperationValues
  USE MCCollisions

  IMPLICIT NONE
  INTEGER number
  TYPE (binary_tree), POINTER :: current
  INTEGER ALLOC_ERR

  current => Collided_particle                  ! start from the head node of the binary tree

  DO                                            ! go through the allocated nodes to the end of the branch

     IF (number.GT.current%number) THEN         
        IF (ASSOCIATED(current%Larger)) THEN        
           current => current%Larger               ! 
           CYCLE                                   ! go to the next node, with larger "number"
        ELSE
           ALLOCATE(current%Larger, STAT=ALLOC_ERR)
           IF (ALLOC_ERR.NE.0) THEN
              PRINT '(/2x,"Process ",i3," : Add_to_stored_list: Error in ALLOCATE current%Larger !!!")', Rank_of_process
              PRINT  '(2x,"Program will be terminated now :(")'
              STOP
           END IF
           current => current%Larger
           EXIT
        END IF
     END IF

     IF (number.LT.current%number) THEN
        IF (ASSOCIATED(current%Smaller)) THEN        
           current => current%Smaller              ! 
           CYCLE                                   ! go to the next node, with smaller "number"
        ELSE
           ALLOCATE(current%Smaller, STAT=ALLOC_ERR)
           IF (ALLOC_ERR.NE.0) THEN
              PRINT '(/2x,"Process ",i3," : Add_to_stored_list: Error in ALLOCATE current%Smaller !!!")', Rank_of_process
              PRINT  '(2x,"Program will be terminated now :(")'
              STOP
           END IF
           current => current%Smaller
           EXIT
        END IF
     END IF
     
  END DO

  current%number = number                       ! store the number
  NULLIFY(current%Larger)
  NULLIFY(current%Smaller)

END SUBROUTINE Add_to_stored_list

!---------------------------------------------------
! this subroutine kills the nodes of the binary tree
RECURSIVE SUBROUTINE Node_Killer(node)

  USE ParallelOperationValues
  USE MCCollisions
  IMPLICIT NONE
  
  TYPE (binary_tree), POINTER :: node
  INTEGER DEALLOC_ERR

  IF (ASSOCIATED(node%Larger))  CALL Node_Killer(node%Larger)
  IF (ASSOCIATED(node%Smaller)) CALL Node_Killer(node%Smaller)

  DEALLOCATE(node, STAT=DEALLOC_ERR)
  IF (DEALLOC_ERR.NE.0) THEN
     PRINT '(/2x,"Process ",i3," : Error in Node_Killer !!!")', Rank_of_process
     PRINT  '(2x,"Program will be terminated now :(")'
     STOP
  END IF

  RETURN

END SUBROUTINE Node_Killer


!*******************************************************************************************
! This module is used by two subroutines, calculating the arbitrary velocity
! according to the isotropic maxwell distribution
MODULE MaxwellVelocity
  REAL(8) v(0:9000)          !(0:300)
  REAL(8) v_inj(0:4500)      !(0:150)
  REAL(8) :: U_max = 3.0_8
  INTEGER :: R_max     = 9000 !300
  INTEGER :: R_max_inj = 4500 !150
END MODULE MaxwellVelocity

!-------------------------------------------------------------------------------------------
! Prepares the tabulated values of integral of the maxwell distribution function
!  
SUBROUTINE PrepareMaxwellDistribIntegral

  USE ParallelOperationValues
  USE MaxwellVelocity
  USE CurrentProblemValues, ONLY : N_box_vel
  IMPLICIT NONE

  INTEGER i
  INTEGER N_pnts
  INTEGER count
  REAL(8) V_min, V_max
  REAL(8) F(0:180003)  !(0:30003)            ! to be sure that we overcome V_max
  REAL(8) temp
  REAL(8) dV

  LOGICAL check1, check2

  check1 = .FALSE.
  check2 = .FALSE.

! ------- for symmetrical maxwellian
  N_pnts = 180000  !30000
  V_min  = -U_max
  V_max  =  U_max
  dV = (V_max - V_min) / N_pnts

  F = 0.0_8
  DO i = 1, N_pnts + 3
     F(i) = F(i-1) + EXP( - (V_min + (REAL(i)-0.5_8) * dV)**2 )
  END DO

  temp = F(N_pnts)
  F = F * R_max / temp   ! normalize integral such that F(N_pnts) = R_max

  v(0) = V_min
  count = 0
  DO i = 1, N_pnts + 3
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        v(count) = V_min + i * dV
        IF (count.EQ.R_max) THEN
           check1 = .TRUE.
           EXIT
        END IF
     END IF
  END DO

  v = v * N_box_vel

!--------- for asymmetrical maxwellian * v (used for injection, v > 0)

  N_pnts = 90000   !15000
  V_min  = 0.0_8
  V_max  = U_max
  dV = (V_max - V_min) / N_pnts

  F = 0.0_8
  DO i = 1, N_pnts + 3
     temp = V_min + (REAL(i)-0.5_8) * dV
     F(i) = F(i-1) + EXP( - temp**2 ) * temp
  END DO

  temp = F(N_pnts)
  F(1:(N_pnts+3)) = F(1:(N_pnts+3)) * R_max_inj / temp   ! normalize integral such that F(N_pnts) = R_max_inj

  v_inj(0) = V_min
  count = 0
  DO i = 1, N_pnts + 3
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        v_inj(count) = V_min + i * dV
        IF (count.EQ.R_max_inj) THEN
           check2 = .TRUE.
           EXIT
        END IF
     END IF
  END DO

  v_inj = v_inj * N_box_vel

  IF (check1.AND.check2) THEN
     PRINT '(2x,"Process ",i3," : Integrals for producing maxwell distributions are successfully obtained ...")', &
                                                                                                  & Rank_of_process
  ELSE
     PRINT '(2x,"Process ",i3," : ERROR in PrepareMaxwellDistribIntegral !!!")', Rank_of_process
     PRINT '(2x,"The initialization in PrepareMaxwellDistribIntegral is not performed !!!")'
     PRINT '(2x,"The program will be terminated now :(")'
     STOP
  END IF

END SUBROUTINE PrepareMaxwellDistribIntegral

!-------------------------------------------------------------------------------------------
!  
SUBROUTINE GetMaxwellVelocity(U) 

  USE MaxwellVelocity

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) U

  REAL(8) R
  INTEGER indx
  
  R = R_max * well_random_number()

  indx = INT(R)

  IF (indx.LT.R_max) THEN
     U = v(indx) + (R - indx) * (v(indx+1) - v(indx))
  ELSE
     U = v(R_max)
  END IF
  RETURN
  
END SUBROUTINE GetMaxwellVelocity

!-------------------------------------------------------------------------------------------
!
!  
SUBROUTINE GetInjMaxwellVelocity(U) 

  USE MaxwellVelocity

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) U

  REAL(8) R
  INTEGER indx
  
  R = R_max_inj * well_random_number()

  indx = INT(R)

  IF (indx.LT.R_max_inj) THEN
     U = v_inj(indx) + (R - indx) * (v_inj(indx+1) - v_inj(indx))
  ELSE
     U = v_inj(R_max_inj)
  END IF
  RETURN
  
END SUBROUTINE GetInjMaxwellVelocity

!--------------------------------------------------
!
REAL(8) FUNCTION Shape_of_Natom(x)

  USE CurrentProblemValues, ONLY : N_cells, delta_x_m

  IMPLICIT NONE

  REAL(8) x
  
!  Shape_of_Natom = MIN( 400.0_8**(((x-DBLE(N_cells)) * delta_x_m) / 1.0_8), 1.0_8)

!  Shape_of_Natom = MIN( 10.0_8**(((x-DBLE(N_cells)) * delta_x_m) / 0.05_8), 1.0_8)  ! was 400.0_8 instead of 10.0_8

!  Shape_of_Natom = MIN( 4.0_8**(((x-DBLE(N_cells)) * delta_x_m) / 0.1_8), 1.0_8)  ! was 400.0_8 instead of 10.0_8


! here 400 and 1 are because the density decreases 400 times 
! if one goes 1 meter leftward from the right end (N_cells*delta_x_m) 
! also, shape_of_natom should not exceed 1

  Shape_of_Natom = 1.0_8  ! use this in case of uniform neutral density 

END FUNCTION Shape_of_Natom







