!============================================================================
! This subroutine calculates the mesh values of the streaming charge density
!.  
SUBROUTINE CALCULATE_STR_CHARGE_DENSITY 

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE

  INCLUDE 'mpif.h'

!  INTEGER s                             ! species
  INTEGER ierr, ALLOC_ERR, DEALLOC_ERR
  REAL(8), ALLOCATABLE :: rbufer(:)     ! for Q_strm_spec
  REAL(8), ALLOCATABLE :: rbufer2(:)   
  INTEGER ibufer(1:2)                   ! for Q_left and Q_right
  INTEGER ibufer2(1:2)  

!  LOGICAL flag
!  INTEGER stattus(MPI_STATUS_SIZE)
!  INTEGER n, source


  IF (.NOT.ALLOCATED(rbufer)) THEN
     ALLOCATE(rbufer(1:(2*N_nodes)), STAT = ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(2x,"Process ",i3," : CALCULATE_STR_CHARGE_DENSITY : Error in ALLOCATE rbufer !!!")', Rank_of_process
        PRINT '(2x,"The program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (.NOT.ALLOCATED(rbufer2)) THEN
     ALLOCATE(rbufer2(1:(2*N_nodes)), STAT = ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(2x,"Process ",i3," : CALCULATE_STR_CHARGE_DENSITY : Error in ALLOCATE rbufer2 !!!")', Rank_of_process
        PRINT '(2x,"The program will be terminated now :(")'
        STOP
     END IF
  END IF

  rbufer2 = 0.0_8
  ibufer2 = 0

  IF (Rank_of_process.GT.0) THEN                                                             ! client >>>

! send Q_left and Q_right to the server process
     ibufer(1) = Q_left
     ibufer(2) = Q_right
     CALL MPI_REDUCE(ibufer, ibufer2, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! send Q_strm_spec to the server
     rbufer(1:N_nodes) = Q_strm_spec(0:N_cells, 1)
     IF (N_spec.EQ.1) THEN
        CALL MPI_REDUCE(rbufer, rbufer2, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     ELSE
        rbufer((N_nodes+1):(2*N_nodes)) = Q_strm_spec(0:N_cells, 2)
        CALL MPI_REDUCE(rbufer, rbufer2, 2*N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  ELSE                                                                                       ! server >>>

! receive (accumulate) Q_left and Q_right from all client processes
     ibufer = 0
     CALL MPI_REDUCE(ibufer2, ibufer, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     Q_left  = ibufer(1)
     Q_right = ibufer(2)

     IF (BC_flag.EQ.1) THEN
! floating wall, direct accumulation
        full_Q_left = DBLE(Q_left)
        full_Q_right = DBLE(Q_right)
     ELSE
! constant given potential or external circuit, the accumulated value is the correction
        full_Q_left  = full_Q_left  + DBLE(Q_left)
        full_Q_right = full_Q_right + DBLE(Q_right)
     END IF

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! receive (accumulate) Q_strm_spec from all client processes
     rbufer = 0.0_8
     IF (N_spec.EQ.1) THEN
        CALL MPI_REDUCE(rbufer2, rbufer, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        Q_strm_spec(0:N_cells, 1) = rbufer(1:N_nodes)
        Q_stream(0:N_cells)   =     Qs(1) * Q_strm_spec(0:N_cells, 1) 
        Xi(0:N_cells)         =   K_Xi(1) * Q_strm_spec(0:N_cells, 1)
     ELSE
        CALL MPI_REDUCE(rbufer2, rbufer, 2*N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        Q_strm_spec(0:N_cells, 1) = rbufer(1:N_nodes)
        Q_strm_spec(0:N_cells, 2) = rbufer((N_nodes+1):(2*N_nodes))
        Q_stream(0:N_cells)   =     Qs(1) * Q_strm_spec(0:N_cells, 1) +   Qs(2) * Q_strm_spec(0:N_cells, 2)
        Xi(0:N_cells)         =   K_Xi(1) * Q_strm_spec(0:N_cells, 1) + K_Xi(2) * Q_strm_spec(0:N_cells, 2)
     END IF

  END IF

  IF (ALLOCATED(rbufer)) THEN
     DEALLOCATE(rbufer, STAT = DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(2x,"Process ",i3," : CALCULATE_STR_CHARGE_DENSITY : Error in DEALLOCATE rbufer !!!")', Rank_of_process
        PRINT *, 'The program will be terminated now :('
        STOP
     END IF
  END IF

  IF (ALLOCATED(rbufer2)) THEN
     DEALLOCATE(rbufer2, STAT = DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(2x,"Process ",i3," : CALCULATE_STR_CHARGE_DENSITY : Error in DEALLOCATE rbufer2 !!!")', Rank_of_process
        PRINT *, 'The program will be terminated now :('
        STOP
     END IF
  END IF
  
END SUBROUTINE CALCULATE_STR_CHARGE_DENSITY 
    
!=====================================================================================
! This subroutine calculates the mesh values of the potential F(1:N_nodes)
! at the integer nodes. Then it calculates the mesh values of the longitudinal electric field
! EZ(1:N_nodes)
! Besides, the energy of the longitudinal electric field is calculated too. IN FUTURE!!!
SUBROUTINE CALCULATE_STR_LONG_ELECTR_FIELD

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BeamInPlasma, ONLY : PeriodicBoundaryFlag
  IMPLICIT NONE

  INCLUDE 'mpif.h'

!  INTEGER igrid

  REAL(8) a_eq(0:N_cells,1:4)
  REAL(8) b_eq(0:N_cells)

  INTEGER i, j
  REAL(8) tmp
  INTEGER ierr

! function
  REAL(8) U_ext_vst

!  K_Q            = 1.0_8 / REAL(N_of_particles_cell * N_of_cells_debye)

  IF (Rank_of_process.EQ.0) THEN                                                          ! server >>>

     SELECT CASE (BC_flag)
     CASE (0)              ! potentials of walls (plasma boundaries) are given
        a_eq( 0, 1) =  0.0_8
        a_eq( 0, 2) =  1.0_8
        a_eq( 0, 3) =  0.0_8
        b_eq( 0)    = U_ext_vst(T_cntr * delta_t_s)
        
     CASE (1)              ! potentials of walls (plasma boundaries) are floating
        a_eq( 0, 1) =  0.0_8
        a_eq( 0, 2) =  1.0_8
        a_eq( 0, 3) = -1.0_8
        b_eq( 0)    = K_Q * (full_Q_left + Q_stream(0)) / (1.0_8 + Xi(0) + Xi(1))
 
     CASE (2)              ! there is an external contour with a battery and a resistor
        a_eq( 0, 1) =  0.0_8
        a_eq( 0, 2) =  1.0_8 + Xi(0) + Xi(1) + factor_SR * (1.0_8 + 2.0_8 * Xi(0))
        a_eq( 0, 3) = -1.0_8 - Xi(0) - Xi(1)
        b_eq( 0)    = factor_SR * U_ext * (1.0_8 + 2.0_8 * Xi(0)) + K_Q * (full_Q_left + Q_stream(0))

     END SELECT

     DO i = 1, N_cells - 1
        a_eq(i, 1) = 1.0_8 + Xi(i-1) + Xi(i)
        a_eq(i, 3) = 1.0_8 + Xi(i+1) + Xi(i)
        a_eq(i, 2) = - a_eq(i, 1) - a_eq(i, 3)
        b_eq(i)    = - K_Q * Q_stream(i)
     END DO

     a_eq(N_cells, 1) =  0.0_8
     a_eq(N_cells, 2) =  1.0_8
     a_eq(N_cells, 3) =  0.0_8
     b_eq(N_cells) = 0.0_8

     a_eq(0:N_cells, 4) = 0.0_8

! diagonalization . . .
     DO i = 1, N_cells !N_nodes
        IF (a_eq(i-1,2).NE.0.0_8) THEN
           tmp = a_eq(i, 1) / a_eq(i-1, 2)
           a_eq(i, 1) = 0.0_8
           a_eq(i, 2) = a_eq(i, 2) - tmp * a_eq(i-1, 3)
           b_eq(i)    = b_eq(i)    - tmp * b_eq(i-1) 
        ELSE
           IF (a_eq(i-1,3).EQ.0.0_8) THEN
              PRINT *, 'Error in Poisson system!', T_cntr
              PRINT '(2x,"Equation ",i5,": a(1) =",e12.5," a(2) =",e12.5," a(3) =",e12.5," b = ",e12.5)', &
                                                               & i, a_eq(i,1), a_eq(i,2), a_eq(i,3), b_eq(i)
              OPEN (41, FILE = 'error_a_b_eq.dat')
              DO j = -1, N_nodes
                 WRITE (41, '(2x,i5,4(2x,e12.5),6x,e12.5)') j, a_eq(j,1), a_eq(j,2), a_eq(j,3), a_eq(j,4), b_eq(j)
              END DO
              CLOSE (41, STATUS = 'KEEP')
              PRINT *, 'Program will be terminated now :('
              STOP
           END IF
           a_eq(i-1, 2) = a_eq(i, 1)
           a_eq(i, 1)   = 0.0_8
           
           tmp          = a_eq(i-1, 3)
           a_eq(i-1, 3) = a_eq(i, 2)
           a_eq(i, 2)   = tmp
           
           a_eq(i-1, 4) = a_eq(i, 3)
           a_eq(i, 3)   = 0.0_8
           
           tmp       = b_eq(i-1)
           b_eq(i-1) = b_eq(i)
           b_eq(i)   = tmp
        END IF
     END DO
  
     F(N_cells)   =  b_eq(N_cells) / a_eq(N_cells, 2)
     F(N_cells-1) = (b_eq(N_cells-1) - a_eq(N_cells-1, 3) * F(N_cells)) / a_eq(N_cells-1, 2)
     
     DO i = N_cells - 2, 0, -1
        F(i) = (b_eq(i) - a_eq(i, 3) * F(i+1) - a_eq(i, 4) * F(i+2)) / a_eq(i, 2)
     END DO

     EX(0)       = (-0.5_8 * K_Q * Q_stream(0)       + (F(0)         - F(1)      ) * (1.0_8 + Xi(0)         + Xi(1))       ) / (1.0_8 + 2.0_8 * Xi(0))
     EX(N_cells) = ( 0.5_8 * K_Q * Q_stream(N_cells) + (F(N_cells-1) - F(N_cells)) * (1.0_8 + Xi(N_cells-1) + Xi(N_cells)) ) / (1.0_8 + 2.0_8 * Xi(N_cells))

     DO i = 1, N_cells - 1
        EX(i) = 0.5_8 * (F(i - 1) - F(i + 1)) 
     END DO

     IF (PeriodicBoundaryFlag.EQ.1) THEN
        EX(0) = 0.5_8 * (F(N_cells-1) - F(1)) 
        EX(N_cells) = EX(0)
     END IF

  END IF
     
  CALL MPI_BCAST(EX, (N_cells+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)      ! server sends and clients receive array EX

  CALL MPI_BCAST(F, (N_cells+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)       ! server sends and clients receive array F

  DO i = 0, N_cells - 1      
     GradEX(i) = EX(i+1) - EX(i)
  END DO

END SUBROUTINE CALCULATE_STR_LONG_ELECTR_FIELD 

!-----------------------------------------------
!
REAL(8) FUNCTION U_ext_vst(t_s)

  USE CurrentProblemValues, ONLY : U_ext, U_scl_V

  IMPLICIT NONE

  REAL(8), PARAMETER :: U_rf_V = 1000.0_8
  REAL(8), PARAMETER :: f_rf_Hz = 1.356d7
  REAL(8), PARAMETER :: t_start_s = 1.0d-7

  REAL(8) t_s

  U_ext_vst = U_ext

!  IF (t_s.LT.t_start_s) RETURN

!  U_ext_vst = U_ext + (U_rf_V / U_scl_V) * SIN(6.283185307_8 * (t_s - t_start_s) * f_rf_Hz)

END FUNCTION U_ext_vst








