!==================================================================================
! This subroutine 
!
SUBROUTINE RESIZE_PARTICLE_ARRAYS

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  
  INTEGER s            ! species type

  INTEGER max_buf_len  ! length of bufer, necessary to store the values during the extension of arrays
  INTEGER N_max        ! Maximal number of particles in block after completion of exchange
  INTEGER ALLOC_ERR, DEALLOC_ERR

  REAL(8), ALLOCATABLE :: dbufer(:)  ! bufer for real(8) values
  INTEGER, ALLOCATABLE :: ibufer(:)  ! bufer for integer values

  TYPE(macroparticle), ALLOCATABLE :: bufer(:)

!  INTEGER count        ! counter (cycle) of injected particles

  INTEGER N_start, k, i

  TYPE(injected_particle), POINTER :: current

  max_buf_len = 0
! find the buffer length if it is necessary to increase the size of some blocks
  DO s = 1, N_spec
     N_max = N_part(s) + N_inject(s)     ! to be sure that we have space for injected particles
     IF (Length(s).GE.N_max) CYCLE
     IF (N_part(s).GT.max_buf_len) max_buf_len = N_part(s)
  END DO

! if necessary, allocate bufer
  ALLOCATE(bufer(1:max_buf_len), STAT=ALLOC_ERR)

! Increasing the array lengthes (if necessary)
  DO s = 1, N_spec
     N_max = N_part(s) + N_inject(s) 
 
     IF (Length(s).GE.N_max) CYCLE                                    ! do nothing if the present length is sufficient

     Length(s) = 1.1 * REAL(N_max) + 1                                ! calculate the new length                

     DO i = 1, N_part(s)
        bufer(i)%X       = species(s)%part(i)%X
        bufer(i)%Z       = species(s)%part(i)%Z
        bufer(i)%VX      = species(s)%part(i)%VX
        bufer(i)%VY      = species(s)%part(i)%VY
        bufer(i)%VZ      = species(s)%part(i)%VZ
        bufer(i)%AX      = species(s)%part(i)%AX
        bufer(i)%Tag     = species(s)%part(i)%Tag
        bufer(i)%prev_VX = species(s)%part(i)%prev_VX
     END DO
 
     DEALLOCATE(species(s)%part, STAT = DEALLOC_ERR)
     ALLOCATE(species(s)%part(1:Length(s)), STAT=ALLOC_ERR)

     DO i = 1, N_part(s)
        species(s)%part(i)%X       = bufer(i)%X   
        species(s)%part(i)%Z       = bufer(i)%Z
        species(s)%part(i)%VX      = bufer(i)%VX   
        species(s)%part(i)%VY      = bufer(i)%VY   
        species(s)%part(i)%VZ      = bufer(i)%VZ   
        species(s)%part(i)%AX      = bufer(i)%AX   
        species(s)%part(i)%Tag     = bufer(i)%Tag   
        species(s)%part(i)%prev_VX = bufer(i)%prev_VX   
     END DO
 
     DO i = N_part(s)+1, Length(s)
        species(s)%part(i)%X       = 0.0_8
        species(s)%part(i)%Z       = 0.0_8
        species(s)%part(i)%VX      = 0.0_8   
        species(s)%part(i)%VY      = 0.0_8   
        species(s)%part(i)%VZ      = 0.0_8   
        species(s)%part(i)%AX      = 0.0_8   
        species(s)%part(i)%Tag     = 0   
        species(s)%part(i)%prev_VX = 0.0_8   
     END DO
 
  END DO

  IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)

! end of arrays extension
!===========================                         
!print *, 'started recursive exchange', T_cntr
!print *, N_accept
!print *, Counter_add

  DO s = 1, N_spec

     NULLIFY(current)  

     SELECT CASE(s)
        CASE(1)                                 ! inject electrons into the first group
          current => Inj_electron
        CASE(2)                                 ! inject ions to the first group
          current => Inj_ion
     END SELECT

     N_start = N_part(s) + 1
     N_part(s) = N_part(s) + N_inject(s)

     IF (N_part(s).GT.Length(s)) THEN
        PRINT '(2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : While rearranging species ",i2," ...")', Rank_of_process, s
        PRINT '(2x,"Error! Number of particles is greater than the length!")'
        PRINT '(2x," N_part = ",i7," Length = ",i7)', N_part(s), Length(s)
        PRINT '(2x,"Program will be terminated now, sorry! :(")'
        STOP
     END IF

     DO k = N_start, N_part(s)

        IF(.NOT.ASSOCIATED(current)) THEN
           PRINT '(2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : While adding injected species ",i2," ...")', Rank_of_process, s
           PRINT '(2x,"Error! current is NOT associated!")'
           PRINT '(2x,"Program will be terminated now, sorry! :(")'
           STOP
        END IF
        species(s)%part(k)%X = current%X
        species(s)%part(k)%Z = 0.0_8
        species(s)%part(k)%VX = current%VX
        species(s)%part(k)%VY = current%VY
        species(s)%part(k)%VZ = current%VZ
        species(s)%part(k)%AX = 0.0_8  !current%AX
        species(s)%part(k)%Tag = current%Tag
        species(s)%part(k)%prev_VX = 0.0_8   

        IF (k.LT.N_part(s)) current => current%next
     END DO

  END DO

!print *, 'finished recursive exchange'              
!===========================

  max_buf_len = 0
! find the buffer length if it is necessary to decrease the size of some blocks
  DO s = 1, N_spec
     N_max = 0.8*REAL(Length(s))
     IF (N_part(s).GE.N_max) CYCLE
     IF (N_part(s).GT.max_buf_len) max_buf_len = N_part(s)   ! if for some species the particles occupy less than 80% of the assigned array
  END DO

! if necessary, allocate bufer
  ALLOCATE(bufer(1:max_buf_len), STAT=ALLOC_ERR)

! Decreasing the array lengthes, if necessary
  DO s = 1, N_spec
     N_max = 0.8*REAL(Length(s))
     IF (N_part(s).GE.N_max) CYCLE     ! if block is not too lengthy, take the next block
     IF (Length(s).LT.100) CYCLE     ! to avoid cycling 10-9-10 or 9-8-9, etc.

     Length(s) = 0.9*REAL(Length(s)) + 1   ! this is the new length

     DO i = 1, N_part(s)
        bufer(i)%X       = species(s)%part(i)%X
        bufer(i)%Z       = species(s)%part(i)%Z
        bufer(i)%VX      = species(s)%part(i)%VX
        bufer(i)%VY      = species(s)%part(i)%VY
        bufer(i)%VZ      = species(s)%part(i)%VZ
        bufer(i)%AX      = species(s)%part(i)%AX
        bufer(i)%Tag     = species(s)%part(i)%Tag
        bufer(i)%prev_VX = species(s)%part(i)%prev_VX
     END DO
  
     DEALLOCATE(species(s)%part, STAT = DEALLOC_ERR)
     ALLOCATE(species(s)%part(1:Length(s)), STAT=ALLOC_ERR)

     DO i = 1, N_part(s)
        species(s)%part(i)%X       = bufer(i)%X   
        species(s)%part(i)%Z       = bufer(i)%Z
        species(s)%part(i)%VX      = bufer(i)%VX   
        species(s)%part(i)%VY      = bufer(i)%VY   
        species(s)%part(i)%VZ      = bufer(i)%VZ   
        species(s)%part(i)%AX      = bufer(i)%AX   
        species(s)%part(i)%Tag     = bufer(i)%Tag   
        species(s)%part(i)%prev_VX = bufer(i)%prev_VX   
     END DO
 
     DO i = N_part(s)+1, Length(s)
        species(s)%part(i)%X       = 0.0_8
        species(s)%part(i)%Z       = 0.0_8
        species(s)%part(i)%VX      = 0.0_8   
        species(s)%part(i)%VY      = 0.0_8   
        species(s)%part(i)%VZ      = 0.0_8   
        species(s)%part(i)%AX      = 0.0_8   
        species(s)%part(i)%Tag     = 0   
        species(s)%part(i)%prev_VX = 0.0_8   
     END DO

  END DO

  IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)

!print *, 'finished REARRANGE_ACTIVE_BLOCK'

END SUBROUTINE  RESIZE_PARTICLE_ARRAYS
