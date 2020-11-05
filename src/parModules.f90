!========================================================================================================================
!
MODULE ParallelOperationValues

  INTEGER Rank_of_process
  INTEGER N_of_processes

END MODULE ParallelOperationValues

!========================================================================================================================
!
MODULE CurrentProblemValues

  INTEGER, PARAMETER :: N_spec = 2                   ! The number of species of particles used in simulations

  REAL(8), PARAMETER :: e_Cl     = 1.602189d-19      ! Charge of single electron [Cl]
  REAL(8), PARAMETER :: m_e_kg   = 9.109534d-31      ! Mass of single electron [kg]
  REAL(8), PARAMETER :: eps_0_Fm = 8.854188d-12      ! The dielectric constant [F/m]
  REAL(8), PARAMETER :: mu_0_Hm  = 1.256637d-6       ! The magnetic constant [H/m]
  REAL(8), PARAMETER :: amu_kg = 1.660565d-27        ! atomic mass unit [kg]

  REAL(8) R_ext_ohm       ! resistance of the resistor in the external circuit [ohm] 
  REAL(8) S_electrode_cm2 ! area of an electrode [cm^2]
  REAL(8) E_z_ext_Vm      ! The external accelerating Z-electric field [V/m]
  REAL(8) B_x_ext_Gs      ! The external radial X-magnetic field [Gauss]
  REAL(8) B_y_ext_Gs      ! The external radial Y-magnetic field [Gauss]
  REAL(8) U_ext_V         ! The externally maintained potential difference between the metal end plates [V]
  REAL(8) L_plasma_m      ! Plasma's length [m], aproximate, may be modified after mesh distribution
  REAL(8) N_plasma_m3     ! Plasma's density [m^-3], affects scaling
  REAL(8) M_i_amu         ! Mass of single ion [a.m.u.]
  REAL(8) T_e_eV          ! The electron temperature [eV], affects scaling
  REAL(8) T_i_eV          ! The ion temperature [eV]
  REAL(8) T_sim_ns        ! Duration of simulations [ns]

  LOGICAL ions_magnetized          ! flag shows whether to use the magnetic field in ion motion equations (true) or not (false)
  LOGICAL ions_sense_Ez            ! if the magnetic field is omitted for ions, they will be accelerated infinitely
                                   ! in the z-direction if the external electric field E_z_ext_Vm is applied,
                                   ! for some applications this may be not necessary, so this flag keeps (true) or removes (false) 
                                   ! the external electric field from the ion motion equation
                                   !
                                   ! WARNING, ions_sense_Ez is in effect only if ions_magnetized=false

  LOGICAL trace_z_shift            ! if it is true, shift of an electron in the Z-direction is calculated
                                   ! and once this shift exceeds a limit value, the particle is thermalized
                                   ! i.e. it is given a new velocity vector using distribution with parameters (temperatures, drifts)
                                   ! of the initial EVDF

  REAL(8) max_Z_shift              ! limiting value of the Z-shift which triggers thermalizetion of a particle [dim-less]

  LOGICAL trace_z_shift_ion        ! same as above but for ions
  REAL(8) max_Z_shift_ion          !

  REAL(8) VZi_0          ! ion Z-velocity at Z=0, use the ion acoustic speed for now

  REAL(8) Tx_e_eV         ! x-electron temperature for anisotropic electron distribution [eV]
  REAL(8) Tz_e_eV         ! z-electron temperature for anisotropic electron distribution [eV]
  REAL(8) N_distrib_m3    ! plasma density for electron distribution [m^-3], affects the number of particles, not the scaling
  LOGICAL use_ExB_for_electrons  ! these two flags command to use (true) or to not use (false) ExB drift
  LOGICAL use_ExB_for_ions       ! while setting the initial velocity distributions for electrons and ions

! we keep separate variables for electron and ion drift velocities to make sure that ion velocities may be set to zero
  REAL(8) Vx_e_drift        ! electron dim-less ExB-drift x-velocity, -Ez/By
  REAL(8) Vy_e_drift        ! electron dim-less ExB-drift y-velocity,  Ez/Bx
  REAL(8) Vx_i_drift        ! ion dim-less ExB-drift x-velocity, -Ez/By
  REAL(8) Vy_i_drift        ! ion dim-less ExB-drift y-velocity,  Ez/Bx

  INTEGER BC_flag                ! Flag, defines the boundary conditions for the potential
  INTEGER N_of_particles_cell    ! Number of macroparticles of one species per cell
  INTEGER N_of_cells_debye       ! Number of minimal cells per electron Debye length
  INTEGER N_max_vel              ! Factor defines the maximal expected electron velocity (in thermal velocities, used for calculation of timestep)
  INTEGER N_max_vel_distrib      ! Factor defines the maximal velocity for initial distribution (in thermal velocities)

  INTEGER N_box_vel              ! Number of velocity boxes per unit of initial thermal velocity
                                 ! Note, this value affects the scaling !

  REAL(8) W_plasma_s1   ! The electron plasma circular frequency [s^-1]
                        ! we assume that N_e = N_i and we will use only 
                        ! two species now - the electrons and the ions
  REAL(8) W_cycl_x_s1   ! The electron cyclotron frequency for Bx [s^-1]
  REAL(8) W_cycl_y_s1   ! The electron cyclotron frequency for By [s^-1]
  REAL(8) v_Te_ms       ! The electron thermal velocity [m/s]
  REAL(8) L_debye_m     ! The Debye length [m]

  REAL(8) E_scl_Vm      ! Scaling factor for E [V/m]
  REAL(8) U_scl_V       ! Scaling factor for the potential F [V]

  REAL(8) B_scl_Tl      ! Scaling factor for B [Tesla] 
  REAL(8) B_scl_Gs      ! Scaling factor for B [Gauss] - note, different units
  REAL(8) R_scl_nClm3   ! Scaling factor for charge density [nCl/m^3]
  REAL(8) N_scl_m3      ! Scaling factor for particles density [m^-3]
  REAL(8) V_scl_ms      ! Scaling factor for particle velocity [m/s]

  REAL(8) factor_SR     ! Dimensionless factor which appears when there is an external circuit

  REAL(8) U_ext         ! The externally maintained potential difference between the metal end plates [dim_less]

  INTEGER Q_left        ! Surface charge on the left  plasma-wall boundary [dim-less]
  INTEGER Q_right       ! Surface charge on the right plasma-wall boundary [dim-less]
  REAL(8) full_Q_left   ! Surface charge on the left  plasma-wall boundary which may account for external circuit effects [dim-less]
  REAL(8) full_Q_right  ! Surface charge on the right plasma-wall boundary which may account for external circuit effects [dim-less]

  INTEGER(8) N_of_part_initial   ! Initial number of macroparticles of one species in the whole system
                                 ! ###NEW### :: note that it is integer(8) now which allows to use much more particles

  REAL(8) delta_t_s     ! The time step [s] 
  REAL(8) delta_x_m     ! The cell size [m]

  INTEGER T_cntr        ! Counter of time steps
  INTEGER Start_T_cntr  ! Starting time step number
  INTEGER Max_T_cntr    ! Number of steps in time
  INTEGER N_nodes       ! Number of plasma nodes (node #0 is at z=0, and node #N_nodes-1=N_cells is at L_plasma_m)
  INTEGER N_cells       ! Number of cells in plasma (N_nodes = N_cells + 1)

  INTEGER SaveCheck_step            ! Time interval (in steps) for creating the checkpoints (full save)
  INTEGER Restore_from_checkpoint   ! Flag, controls the initialization
                                    ! if 0 then - ordinary initialization of particles, diagnostics and snapshots
                                    ! if 1 then - particles parameters (x, v, ax, tag) are being read from the checkpoint datafiles
                                    !             * assume that the diagnostic parameters for temporal output 
                                    !               (WriteOut_step, WriteStart_step, WriteAvg_step) 
                                    !               is not changed compared to the simulation, which produced the checkpoint
                                    !             * assume that the parameters of snapshot creation remains unchanged 
                                    !               for the groups, which precede or include the checkpoint; 
                                    !               the subsequent groups can be modified  
                                    !             * the upper time limit can be modified
                                    ! if 2 then -  particles parameters (x, v, ax, tag) are being read from the checkpoint datafiles
                                    !             * initialization of temporal diagnostics and snapshots occurs in regular manner,
                                    !               i.e. on the basis of data from files "ssc_initial.dat" and "ssc_snapshot.dat" only
!  CHARACTER(16) check_g_filename    ! name of file, where the general data will be saved at the checkpoint
!  CHARACTER(16) check_e_filename    ! name of file, where the parameters of electrons will be saved at the checkpoint
!  CHARACTER(16) check_i_filename    ! name of file, where the parameters of ions will be saved at the checkpoint

!  INTEGER I_random_seed   ! The integer number which will be used by the random number generator through the whole program
!  INTEGER I_random_seed_2 ! The integer number which will be used by the random number generator through the whole program in specifying type of collision
!  INTEGER I_random_seed_3 ! The integer number which will be used by the random number generator through the whole program in accounting for the nonuniform neutral density
!  INTEGER I_random_seed_4 ! The integer number which will be used by the random number generator through the whole program in specifying polar angle of scattering
!  INTEGER I_random_seed_5 ! The integer number which will be used by the random number generator through the whole program in accounting azimuthal angle of scattering

! SPECIES ARRAYS, HAVE FIXED DIMENSION (1:N_spec) ***** 
! particle parameters ------------------------

  REAL(8) Ms(1:N_spec)  ! normalized masses of species (by the electron mass)

  INTEGER Qs(1:N_spec)  ! normalized charges of species (by the electron charge ABS. value)

  REAL(8) QMs(1:N_spec) ! normalized charge-to-mass ratio of species (by that of electron ABS. value)

  REAL(8) VT(1:N_spec)  ! normalized thermal velocities of species (by that of the electrons)

! precalculated coefficients -----------------

  REAL(8) KVx         ! Coefficients for equations of motion (X-PRE-FIN-moving)  !
  REAL(8) K_Q         ! Coefficients used in Poisson's equation

  REAL(8) K_Xi(1:N_spec)   ! Coefficients used in interpolation of effective susceptibility 

  REAL(8)  K11(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (X-PRE-acceleration)
  REAL(8)  K12(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (X,Y-PRE-acceleration)
  REAL(8)  K13(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (X,Z-PRE-acceleration)
  REAL(8)  K22(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (Y-PRE-acceleration)
  REAL(8)  K23(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (Y,Z-PRE-acceleration)
  REAL(8)  K33(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (Z-PRE-acceleration)

  REAL(8)  A11(1:N_spec) ! Coefficients (before Ax^{n-1}) for equations of motion (X-PRE-acceleration)
  REAL(8)  A21(1:N_spec) ! Coefficients (before Ax^{n-1}) for equations of motion (Y-PRE-acceleration)
  REAL(8)  A31(1:N_spec) ! Coefficients (before Ax^{n-1}) for equations of motion (Z-PRE-acceleration)

  REAL(8)  A13(1:N_spec) ! Coefficients with Ez_external for equations of motion (X-PRE-acceleration)
  REAL(8)  A23(1:N_spec) ! Coefficients with Ez_external for equations of motion (Y-PRE-acceleration)
  REAL(8)  A33(1:N_spec) ! Coefficients with Ez_external for equations of motion (Z-PRE-acceleration)

  REAL(8) KvEx(1:N_spec) ! Coefficients before Ex^{n+1} for final VX correction
  REAL(8) KvEy(1:N_spec) ! Coefficients before Ex^{n+1} for final VY correction
  REAL(8) KvEz(1:N_spec) ! Coefficients before Ex^{n+1} for final VZ correction

  REAL(8) KxEx(1:N_spec) ! Coefficients before Ex^{n+1} for final X correction

! service data -------------------------------

  INTEGER N_part_client_spec(1:256,1:N_spec)     ! Number of particles in client processes, (1-st index - process, 2-nd index - species)
  INTEGER N_part(1:N_spec)                  ! Numbers of particles 
  INTEGER Length(1:N_spec)                  ! Lengthes of allocated arrays 
  INTEGER N_inject(1:N_spec)                ! Number (counter) of injected particles (due to ionization or SEE)

! ALLOCATABLE ARRAYS *****************************
! mesh values --------------------------------
  REAL(8), ALLOCATABLE :: EX(:)        ! The X-electric field

  REAL(8), ALLOCATABLE :: GradEX(:)    ! The X-gradient of X-electric field,                                   
                                       ! used for interpolation of mesh EX to the particle position

  REAL(8), ALLOCATABLE :: F(:)         ! The electric potential

  REAL(8), ALLOCATABLE :: Xi(:)        ! The effective dielectric susceptibility

  REAL(8), ALLOCATABLE :: Q_stream(:)  ! The interpolated streaming charge

  REAL(8), ALLOCATABLE :: Q_strm_spec(:,:)  ! Streaming charge in the nodes for all blocks of particles (is necessary because susceptibility depends on mass)

! allocatable array types --------------------

!  TYPE dalloc_array
!     REAL(8), POINTER :: part(:)
!  END TYPE dalloc_array

!  TYPE ialloc_array
!     INTEGER, POINTER :: part(:)
!  END TYPE ialloc_array

! particles dynamics ------------------------- 

!integer, allocatable :: collided(:)
!integer, allocatable :: selected(:)

!  TYPE(dalloc_array), ALLOCATABLE :: X_of_spec(:)       ! Array of pointers on arrays of X-coordinate of particles of blocks
!  TYPE(dalloc_array), ALLOCATABLE :: VX_of_spec(:)      ! Array of pointers on arrays of X-velocity of particles of blocks
!  TYPE(dalloc_array), ALLOCATABLE :: VY_of_spec(:)      ! Array of pointers on arrays of Y-velocity of particles of blocks
!  TYPE(dalloc_array), ALLOCATABLE :: VZ_of_spec(:)      ! Array of pointers on arrays of Z-velocity of particles of blocks
!  TYPE(dalloc_array), ALLOCATABLE :: AX_of_spec(:)      ! Array of pointers on arrays of smoothed X-acceleration of particles of blocks
!  TYPE(ialloc_array), ALLOCATABLE :: Tag_of_spec(:)     ! Array of pointers on arrays of tag of particles of blocks
!  TYPE(dalloc_array), ALLOCATABLE :: prev_VX_of_spec(:) ! Array of pointers on arrays of old X-velocity of particles of blocks

  TYPE macroparticle
     REAL(8) X
     REAL(8) Z
     REAL(8) VX
     REAL(8) VY
     REAL(8) VZ
     REAL(8) AX
     INTEGER Tag
     REAL(8) prev_VX
  END TYPE macroparticle

  TYPE macroparticle_array
     TYPE(macroparticle), ALLOCATABLE :: part(:)
  END TYPE macroparticle_array

  TYPE(macroparticle_array), ALLOCATABLE :: species(:)
!
!### should work like this: species(s)%part(k)%X

! LINKED LIST FOR INJECTED PARTICLES *****************

  TYPE injected_particle
     REAL(8) X
     REAL(8) VX
     REAL(8) VY
     REAL(8) VZ
!     REAL(8) AX
     INTEGER Tag
     TYPE (injected_particle), POINTER :: Next
  END TYPE injected_particle

  TYPE(injected_particle), POINTER :: Inj_electron
  TYPE(injected_particle), POINTER :: Current_electron
  TYPE(injected_particle), POINTER :: Inj_ion 
  TYPE(injected_particle), POINTER :: Current_ion 

  INTEGER, PARAMETER :: eTag_Coll_Neutral       =  2      ! tag to mark electron collided with neutral atom [elastic, excitation and ionization (both)]
  INTEGER, PARAMETER :: eTag_Emit_Left          =  1      ! tag to mark electrons reflected / backscattered / emitted from the left wall
  INTEGER, PARAMETER :: eTag_Emit_Right         = -1      ! tag to mark electrons reflected / backscattered / emitted from the right wall

!  LOGICAL Tags_must_be_nullified        ! flag, is true if it is necessary to clear all tags in CLEAR_TAGS, is modified in
                                        ! that subroutine and in snapshot procedures

!real(8) initial_electron_energy_eV           ! for testing purposes only, initial energy of electrons, [eV]
!real(8) initial_angle_of_incidence_deg       ! for testing purposes only, initial angle of incidence, [deg]

END MODULE CurrentProblemValues

!========================================================================================================================
! >>
!
! To add a NEW kind of collision for some species you must modify all lines containing !@#$ 
!
MODULE MCCollisions

  USE CurrentProblemValues, ONLY : N_spec

! allocatable array types
  TYPE ralloc_kind_energy
     REAL, POINTER :: kind_energy(:,:)
  END TYPE ralloc_kind_energy

  TYPE ialloc_activated
     INTEGER, POINTER :: activated(:) 
  END TYPE ialloc_activated

  INTEGER maxNcolkind_spec(1:N_spec)         ! maximal number of kinds of collisions for electrons,
                                        ! currently we have collisions for two plasma species only

!------- to read from file ---------
  REAL(8) M_neutral_amu                 ! The neutral component mass         [a.m.u.]
  REAL(8) N_neutral_m3                  ! The neutral component density      [m^-3]
  REAL(8) T_neutral_eV                  ! The neutral component temperature  [eV]

  REAL(8) Freq_turb_e_1_s1              ! Frequency of collisions, accounting the electron turbulence, model 1, [s^-1]
  REAL(8) Freq_turb_i_1_s1              ! Frequency of collisions, accounting the ion turbulence, model 1,      [s^-1]

  REAL(8) alpha_Vscl                    ! Scaling factor, when you take some velocity from normalized Maxwell distribution,
                                        ! and multiply the velocity by the above mentioned factor, 
                                        ! you obtain the dimensionless value in electron thermal velocities

  REAL(8) energy_factor_eV_s(1:N_spec)       ! Scaling factor for obtaining energy of particle in electron-Volts 
                                        ! from the squared particle dim-less velocity :  W_kin_eV(s) = energy_factor_eV_s(s) * v**2,
                                        ! here s - species index, v - dim-less absolute value of velocity (in units of v_Te / N_box_vel)

  INTEGER Colflag_kind_spec(1:4,1:N_spec)    ! For different species [electron(*,1), ion(*,2)],                !@#$
                                        ! the arrays of flags, indicating that some particular kind
                                        ! of collisions with neutrals is turned on
                                        ! 1%1 <=> e-n elastic,         model 1
                                        ! 2%1 <=> e-n excitation,      model 1
                                        ! 3%1 <=> e-n ionization,      model 1 
                                        ! 4%1 <=> turbulence,          model 1 
                                        ! 1%2 <=> i-n elastic,         model 1
                                        ! 2%2 <=> i-n charge exchange, model 1
                                        ! 3%2 <=> turbulence,          model 1
                                        ! 4%2 <=> yet empty

  INTEGER N_en_elast                             ! number of experimental values to be read 
  REAL, ALLOCATABLE :: Energy_en_elast_eV(:)     ! energies of known cross-secion values  for e-n elastic    collisions, [eV]
  REAL, ALLOCATABLE :: CrSect_en_elast_m2(:)     ! experimentally measured cross-sections for e-n elastic    collisions, [m^2]

  INTEGER N_en_excit                             ! number of experimental values to be read 
  REAL, ALLOCATABLE :: Energy_en_excit_eV(:)     ! energies of known cross-secion values  for e-n excitation collisions, [eV]
  REAL, ALLOCATABLE :: CrSect_en_excit_m2(:)     ! experimentally measured cross-sections for e-n excitation collisions, [m^2]
  REAL(8) Thresh_en_excit_eV                     ! energy threshold, [eV]

  INTEGER N_en_ioniz                             ! number of experimental values to be read 
  REAL, ALLOCATABLE :: Energy_en_ioniz_eV(:)     ! energies of known cross-secion values  for e-n ionization collisions, [eV]
  REAL, ALLOCATABLE :: CrSect_en_ioniz_m2(:)     ! experimentally measured cross-sections for e-n ionization collisions, [m^2]
  REAL(8) Thresh_en_ioniz_eV                     ! energy threshold, [eV]

  REAL(8) maxEnergy_eV_spec(1:N_spec)        ! maximal value of energy of colliding plasma species, [eV]

  INTEGER Nenergyval_spec(1:N_spec)          ! number of values of energy 
                                                     ! where the collisions probabilities are calculated 
!------- to calculate --------------

  INTEGER, ALLOCATABLE :: Ncolkind_spec(:)           ! Number of kinds of collisions for different plasma species
                                                     ! Ncolkind_spec(s) < maxNcolkind_spec(s)

  REAL(8), ALLOCATABLE :: deltaEnergy_eV_spec(:)     ! energy step for different colliding plasma species, dim-less

  REAL, ALLOCATABLE :: maxColfrac_spec(:)    ! The MAXIMAL RELATIVE number of particles of one species, !@
                                                   ! colliding at the current time level (including NULL-collisions)

  TYPE(ialloc_activated), ALLOCATABLE :: Colkind_of_spec(:)      ! Set of indexes determining those kinds of collisions,
                                                                 ! which are turned on in file "msc_partcols.dat" 

  TYPE(ralloc_kind_energy), ALLOCATABLE :: Colprob_of_spec(:)    ! For different species, the pointers to the 2-d arrays of
                                                                 ! probabilities of different kinds of activated collisions  
                                                                 ! versus the energy of the particle of that species

  INTEGER, ALLOCATABLE :: Numbers_of_particles(:)
  INTEGER, ALLOCATABLE :: Numbers_of_collisions(:)

!------- counters, for diagnostics -----------
  INTEGER e_n_1_count               ! e-n, elastic, model 1
  INTEGER e_n_2_count               ! e-n, excitation, model 1
  INTEGER e_n_3_count               ! e-n, ionization, model 1
  INTEGER e_t_4_count               ! e-turbulence, model 1               !@#$
  INTEGER i_n_1_count               ! i-n, elastic, model 1
  INTEGER i_n_2_count               ! i-n, charge exchange, model 1
  INTEGER i_t_3_count               ! i-turbulence, model 1               !@#$

! LINKED LIST, which store the numbers of particles participating in collisions

  TYPE binary_tree
     INTEGER number
     TYPE (binary_tree), POINTER :: Larger
     TYPE (binary_tree), POINTER :: Smaller     
  END TYPE binary_tree

  TYPE(binary_tree), POINTER :: Collided_particle

END MODULE MCCollisions

!========================================================================================================================
!
MODULE LangevinCollisions

  INTEGER Accounted_Langevin                         ! Flag, turns on/off the Langevin e-e collisions
  INTEGER T_skip_Lang                                ! Timesteps to skip between the application of collisions
  INTEGER added_cells                                ! number of cells added from each side to the cell-of-interest to
                                                     ! in order to calculate EVDF in the cell-of-interest

!  INTEGER N_of_left_node_Lang                        ! Left boundary for the region where teh EVDF will be calculated
!  INTEGER N_of_right_node_Lang                       ! Right -"-
  INTEGER W_max_Lang                                 ! Maximal velocity accounted for calculation of coefficients 
                                                     ! [in units of the scale thermal velocity]

  INTEGER N_box_w_Lang                               ! number of velocity boxes for calculation of EVDF   

  REAL(8) Factor_F_drag                              ! precalculated factors for the Langevin coefficients
  REAL(8) Factor_D_1_sqrt                            !
  REAL(8) Factor_D_3_sqrt                            !

  REAL(8), ALLOCATABLE :: Fd_w(:,:)                  ! tabulated distribution function over velocity module (in the electron flow frame)
  REAL(8), ALLOCATABLE :: F_drag(:,:)                ! tabulated drag force                                 (-"-)
  REAL(8), ALLOCATABLE :: D_1_sqrt(:,:)              ! tabulated square root of diffusion coefficient D1 in velocity space (-"-)
  REAL(8), ALLOCATABLE :: D_3_sqrt(:,:)              ! tabulated square root of diffusion coefficient D1 in velocity space (-"-)

  REAL(8), ALLOCATABLE ::  V_ion_threshold(:)        ! threshold value of velocity, below which Langevin coefficients for e-i 
  REAL(8), ALLOCATABLE ::  V_ion_threshold_2(:)      ! collisions do not change
  REAL(8), ALLOCATABLE ::  V_ion_threshold_sqrt(:)   ! is equal to the doubled initial ion thermal velocity

  REAL(8), ALLOCATABLE :: Vx_avg(:)                  !
  REAL(8), ALLOCATABLE :: Vy_avg(:)                  ! The averaged {flow} electron velocity
  REAL(8), ALLOCATABLE :: Vz_avg(:)                  !
  REAL(8), ALLOCATABLE :: energy_before(:)           ! total kinetic energy before collisions
  REAL(8), ALLOCATABLE :: energy_after(:)            ! total kinetic energy after collisions, before the corrections are applied
  REAL(8), ALLOCATABLE :: Log_coul(:)                ! Coulomb logarithm

  REAL(8), ALLOCATABLE :: Density_e_sqrt(:)          ! array of square root of electron density in the middles of cells, 
                                                     ! used only by clients, the server node uses (sends to clients) 
                                                     ! the square root of the modified streaming density
END MODULE LangevinCollisions

!========================================================================================================================
! >>
MODULE SEEmission

! the left plasma boundary can be the usual dielectric wall or the plasma source

  INTEGER PlasmaSourceFlag            ! controls the properties of the left plasma boundary
                                      ! turns on / off the refluxing at the left wall (reflection with thermalization)

  INTEGER AddInjectionFlag            ! turns on / off replenishing of particles, which were lost at the right wall, by the
                                      ! particles of the plasma source at the left plasma boundary

  INTEGER const_add_N_to_inject_by_proc       ! number of macroparticles to be injected by each process aat each timestep
                                              ! in order to ensure additional injection with constant flux (individual for each client process)

  INTEGER variable_add_N_to_inject    ! in order to achieve even load on all nodes, 
                                      ! total_N_to_inject = const_add_N_to_inject_by_proc * (N_of_processes-1) + variable_add_N_to_inject
                                      ! the latter number is distributed between different processes all the time:
                                      ! one process injects one pair

! we have three possible processes when an electron hits the wall, producing a secondary electron:
! 1 - elastic backscattering
! 2 - inelastic backscattering
! 3 - true secondary emission

  INTEGER Emitted_model(1:3)          ! for each possible process determines the way of processing
     
  INTEGER Elast_refl_type             ! type of reflection for elastically reflected electron: specular / random

! ELASTIC, MODEL 1
  REAL(8) setD_see_elastic            ! the coefficient (ratio of emitted to incident electrons)
  REAL(8) minE_see_elastic_eV         ! LOWER energy boundary, [eV] 
  REAL(8) minE_see_elastic            ! LOWER energy boundary, [dimensionless] 
  REAL(8) maxE_see_elastic_eV         ! UPPER energy boundary, [eV]
  REAL(8) maxE_see_elastic            ! UPPER energy boundary, [dimensionless]

! ELASTIC, MODEL 2
  REAL(8) E_elast_0_eV                ! threshold energy, [eV]
  REAL(8) E_elast_0                   ! threshold energy, [dimensionless]
  REAL(8) E_elast_max_eV              ! maximum emission energy, [eV]
  REAL(8) E_elast_max                 ! maximum emission energy, [dimensionless]
  REAL(8) maxD_elast                  ! maximum emission yield
  REAL(8) dE_elast_eV                 ! rate of decaying (like half-width) at energies > E_elast_max_eV, [eV]
  REAL(8) dE_elast                    ! rate of decaying (like half-width) at energies > E_elast_max_eV, [dimensionless]
  REAL(8) Frac_elast_highenergy       ! fraction of total SEE yield at high energies, must be  ZERO if the model is not used
     
! INELASTIC, MODEL 1
  REAL(8) setD_see_inelastic          ! the coefficient (ratio of emitted to incident electrons)
  REAL(8) minE_see_inelastic_eV       ! LOWER energy boundary, [eV] 
  REAL(8) minE_see_inelastic          ! LOWER energy boundary, [dimensionless] 
  REAL(8) maxE_see_inelastic_eV       ! UPPER energy boundary, [eV]
  REAL(8) maxE_see_inelastic          ! UPPER energy boundary, [dimensionless]

! INELASTIC, MODEL 2
  REAL(8) Frac_inelastic              ! fraction of total SEE yield, must be ZERO if the model is not used

! TRUE SECONDARY, MODEL 1
  REAL(8) setD_see_true               ! the coefficient (ratio of emitted to incident electrons)
  REAL(8) minE_see_true_eV            ! LOWER energy boundary, [eV] 
  REAL(8) minE_see_true               ! LOWER energy boundary, [dimensionless] 
  REAL(8) maxE_see_true_eV            ! UPPER energy boundary, [eV]
  REAL(8) maxE_see_true               ! UPPER energy boundary, [dimensionless]

! TRUE SECONDARY, MODEL 2
  REAL(8) E_see_0_eV                  ! threshold energy, [eV]
  REAL(8) E_see_0                     ! threshold energy, [dimensionless]
  REAL(8) E_see_max_eV                ! maximal emission energy, [eV]
  REAL(8) E_see_max                   ! maximal emission energy, [dimensionless]
  REAL(8) maxD_see_classic            ! maximal emission coefficient (normal to the surface)
  REAL(8) k_smooth                    ! smoothness factor (from 0 for very rough to 2 for polished surfaces)

  REAL(8) T_see_true_eV               ! Temperature of injected true secondary electrons, [eV]

  INTEGER electron_reflux_count       ! counter of refluxed electrons
  INTEGER ion_reflux_count            ! counter of refluxed ions
  REAL(8) electron_reflux_energy      ! stores the energy of refluxed electrons
  REAL(8) ion_reflux_energy           ! stores the energy of refluxed ions

  INTEGER electron_reemit_count       ! counter of refluxed electrons
  INTEGER ion_reemit_count            ! counter of refluxed ions
  REAL(8) electron_reemit_energy      ! stores the energy of refluxed electrons
  REAL(8) ion_reemit_energy           ! stores the energy of refluxed ions

  INTEGER N_of_lost_ions                ! counter of ions lost at the right wall (this number will be reemitted)

  INTEGER Ion_interac_model        ! Flag (0/1), defines interaction of ions with the walls

  INTEGER ion_left_reflect_count        ! counter of ions reflected from the wall
  REAL(8) ion_left_reflect_energy       ! store the energy of ions reflected from the wall

  INTEGER ion_right_reflect_count        ! counter of ions reflected from the wall
  REAL(8) ion_right_reflect_energy       ! store the energy of ions reflected from the wall

  INTEGER see_left_elastic_count        ! counter of particles, injected as elastically reflected                  tag = 1
  INTEGER see_left_inelastic_count      ! counter of particles, injected as inelastically backscattered            tag = 1
  INTEGER see_left_true_count           ! counter of particles, injected as true secondary                         tag = 1
  REAL(8) see_left_elastic_energy       ! stores the energy of particles, injected as elastically reflected
  REAL(8) see_left_inelastic_energy     ! stores the energy of particles, injected as inelastically backscattered
  REAL(8) see_left_true_energy          ! stores the energy of particles, injected as true secondary      

  INTEGER see_right_elastic_count        ! counter of particles, injected as elastically reflected                 tag = -1
  INTEGER see_right_inelastic_count      ! counter of particles, injected as inelastically backscattered           tag = -1
  INTEGER see_right_true_count           ! counter of particles, injected as true secondary                        tag = -1
  REAL(8) see_right_elastic_energy       ! stores the energy of particles, injected as elastically reflected
  REAL(8) see_right_inelastic_energy     ! stores the energy of particles, injected as inelastically backscattered
  REAL(8) see_right_true_energy          ! stores the energy of particles, injected as true secondary      

  INTEGER prie_left_from_right_count     ! counter of electrons, which come to the left wall as secondaries emitted from the right wall             tag = -1
  INTEGER prie_left_after_coll_count     ! counter of electrons, which come to the left wall after collisions with neutrals (all kinds)             tag = 2
  REAL(8) prie_left_from_right_energy    ! stores the energy of electrons, which come to the left wall as secondaries emitted from the right wall
  REAL(8) prie_left_after_coll_energy    ! stores the energy of electrons, which come to the left wall after collisions with neutrals (all kinds)

  INTEGER prie_right_from_left_count     ! counter of electrons, which come to the right wall as secondaries emitted from the left wall             tag = 1
  INTEGER prie_right_after_coll_count    ! counter of electrons, which come to the right wall after collisions with neutrals (all kinds)            tag = 2
  REAL(8) prie_right_from_left_energy    ! stores the energy of electrons, which come to the right wall as secondaries emitted from the left wall
  REAL(8) prie_right_after_coll_energy   ! stores the energy of electrons, which come to the right wall after collisions with neutrals (all kinds)

  INTEGER sece_left_from_right_count     ! counter of electrons, which are produced at the left wall by secondaries from the right wall             
  INTEGER sece_left_after_coll_count     ! counter of electrons, which are produced at the left wall by electrons after collisions with neutrals (all kinds) 
  INTEGER sece_right_from_left_count     ! counter of electrons, which are produced at the right wall by secondaries from the right wall             
  INTEGER sece_right_after_coll_count    ! counter of electrons, which are produced at the right wall by electrons after collisions with neutrals (all kinds) 

END MODULE SEEmission

!========================================================================================================================
!
MODULE IonInducedSEEmission

! the simplest case
  REAL(8) setD_ionsee_true  ! constant emission coefficient, must be below 1
  REAL(8) minE_ionsee_eV    ! threshold energy of ions [eV], no emission if ion energy is below the threshold
  REAL(8) T_ionsee_eV       ! temperature of emitted electrons [eV]

! related diagnostics
  INTEGER ionsee_left_count         ! counter of electron particles, injected at the left wall  tag = 1
  INTEGER ionsee_right_count        ! counter of electron particles, injected at the right wall  tag = 1

  REAL(8) ionsee_left_energy        ! stores the energy of particles, injected at the left wall      
  REAL(8) ionsee_right_energy       ! stores the energy of particles, injected at the right wall      

END MODULE IonInducedSEEmission

!========================================================================================================================
!
MODULE ElectronInjection

  INTEGER BeamInjectFlag_left                 ! flag, which turns on / off the constant electron emission at the LEFT wall
  INTEGER BeamInjectFlag_right                ! - " - at the RIGHT wall

  REAL(8) Delay_of_injection_ns_left          ! left wall, delay between the start of simulation and the start of emission , [ns] ([timesteps] if negative)
  REAL(8) Delay_of_injection_ns_right         ! right wall, - " -

  REAL(8) Beam_energy_eV_left                 ! left wall, initial energy of cold beam / temperature of warm source, [eV]
  REAL(8) Beam_energy_eV_right                ! right wall, - " -

  REAL(8) Beam_J_Am2_left                     ! left wall, current density of injected stream at the moment of injection only
  REAL(8) Beam_J_Am2_right                    ! right wall, - " -

  REAL(8) VX_e_beam_left                      ! left wall, dimensionless initial velocity of electron beam
  REAL(8) VX_e_beam_right                     ! right wall, - " -

  INTEGER N_to_inject_total_left              ! left wall, total number of macroparticles to be emitted at each timestep
  INTEGER N_to_inject_total_right             ! right wall, - " -

  INTEGER const_N_to_inject_by_proc_left      ! left wall, constant number of particles to be emitted by each process each timestep 
  INTEGER const_N_to_inject_by_proc_right     ! right wall, - " -

  INTEGER variable_N_to_inject_left           ! left wall, variable number of particles to be emitted by each timestep (distributed between all processes)
  INTEGER variable_N_to_inject_right          ! right wall, - " -

  INTEGER inject_every_this_many_timesteps_left  ! left wall, a particle must be emitted every this many timesteps (1=continuous, 2=every other, etc)
  INTEGER inject_every_this_many_timesteps_right ! right wall, - " -

  INTEGER Delay_of_injection_T_cntr_left      ! left wall, calculated delay between the start of simulation and the start of emission [timesteps]
  INTEGER Delay_of_injection_T_cntr_right     ! right wall, - " -
  
  REAL(8) Gap_betw_part_left                  ! left wall, distance between neighbour beam particles
  REAL(8) Gap_betw_part_right                 ! right wall, - " -

  INTEGER UseSmartTagsFlag                    ! flag, controls the modification of the particle's tag related with the change of particle's velocity
                                              ! normally this flag is 1 and the particles, which change the direction of propagation 
                                              ! (velocity sign) inside the plasma, obtain a zero tag, like the confined plasma background;
                                              ! if this flag is 0, the change of the particle's direction of propagation will not modify the particle's tag.

END MODULE ElectronInjection

!========================================================================================================================
!
MODULE BeamInPlasma

  INTEGER PeriodicBoundaryFlag
  INTEGER eBeamInPlasmaFlag          ! flag, which turns on / off the electron beam in plasma
  REAL(8) Start_ebeam_ns             ! delay between the start of simulation and the moment when the beam appears, [ns] ([timesteps] if negative)
  REAL(8) Energy_ebeam_eV            ! initial average energy of beam (energy of electron, which is at rest in the frame related with the beam), [eV]
  REAL(8) Te_ebeam_eV                ! temperature of the warm beam (for distribution in the frame related with the beam)
  REAL(8) Alfa_ebeam                 ! relative (to the scale plasma density) density of the beam 
  INTEGER N_ebeam_total              ! total number of beam macroparticles
  INTEGER N_ebeam_proc               ! number of beam macroparticles in each process (can be different for different processes)
  INTEGER Start_ebeam_T_cntr         ! calculated delay between the start of simulation and the appearance of the beam [timesteps]

  INTEGER, ALLOCATABLE :: N_ebeam_proc_array(:)  ! in this array the server process will store the number of beam particles in each client process
  REAL(8), ALLOCATABLE ::  X_ebeam_part(:)       ! used during the initial distribution, X-coordinates of beam particles
  REAL(8), ALLOCATABLE :: VX_ebeam_part(:)       ! used during the initial distribution, X-velocities of beam particles

END MODULE BeamInPlasma

!========================================================================================================================
!
MODULE Diagnostics

  USE CurrentProblemValues, ONLY : N_spec

  INTEGER WriteOut_step             ! Time interval (in steps, initially given in ???) for writing into the file
  INTEGER WriteStart_step           ! Time (in steps, initially given in ???) when the collecting of data for writing must be started
  INTEGER WriteAvg_step             ! Time interval (-"-) for averaging during writing into the file
  INTEGER TextOut_avg_prds_to_skip  ! Periods of averaging to be skipped between text outputs, 0 = each, 1 = skip one, 2 = skip two, etc.  

  INTEGER Start_diag_Tcntr          ! value of T_cntr when we start collecting data for a new diagnostic point
  INTEGER Finish_diag_Tcntr         ! value of T_cntr when we finish collecting data for that diagnostic point
  REAL(8) Averaging_factor          ! inversed number of timesteps used for averaging
  INTEGER N_of_saved_records        ! counter of records in each data file (versus time)
  INTEGER Save_check_Tcntr          ! value of T_cntr when we create a checkpoint, must always coincide with one of Finish_diag_Tcntr values

  INTEGER T_cntr_to_continue        ! value of T_cntr when a checkpoint was created, this checkpoint will be used to continue or start a simulation

  INTEGER N_of_probes                           ! number of probes for time dependencies

  INTEGER, ALLOCATABLE :: probe_node(:)         ! locations of probes where time dependencies will be saved (numbers of nodes)

  REAL(8) J_scl_Am2     ! Scaling factor for electric current density [nA/m^2]
                        ! This value is placed here because it includes Averaging_factor

! DIAGNOSTIC MESH-SPECIES ARRAYS ****************************

  REAL(8), ALLOCATABLE :: VX2_cell(:,:)         ! collects VX^2 for each cell (here cell is the space between two neighboring nodes) 
  REAL(8), ALLOCATABLE :: VY2_cell(:,:)         ! - " - VY^2  
  REAL(8), ALLOCATABLE :: VZ2_cell(:,:)         ! - " - VZ^2 

  REAL(8), ALLOCATABLE :: VX_cell(:,:)          ! - " - VX 
  REAL(8), ALLOCATABLE :: VY_cell(:,:)          ! - " - VY 
  REAL(8), ALLOCATABLE :: VZ_cell(:,:)          ! - " - VZ 

  INTEGER, ALLOCATABLE :: Npart_cell(:,:)       ! - " - number of particles 

  REAL(8), ALLOCATABLE :: QF_mesh(:)            ! collects potential energy in each node (usual distribution b/w nodes)
  REAL(8), ALLOCATABLE :: QVX_mesh(:)           ! - " - electrical current density in X-direction
  REAL(8), ALLOCATABLE :: QVY_mesh(:)           ! - " - electrical current density in Y-direction
  REAL(8), ALLOCATABLE :: QVZ_mesh(:)           ! - " - electrical current density in Z-direction

! energy values, dimensional
  REAL(8) Energy_full_eV                        ! full energy of system, [eV]

  REAL(8) Energy_pot_eV                         ! potential energy, [eV]                                           
  REAL(8) Energy_pot                            ! - " -                                                             IS USED FOR ACCUMULATION

  REAL(8) Energy_kin_eV(1:N_spec)                    ! kinetic energy, 1 = electrons, 2 = ions, [eV]

  INTEGER N_collect_part(1:N_spec)                   ! number of collected particles, 1 = electrons, 2 = ions

  REAL(8) Avg_kin_energy_eV(1:N_spec)                ! average kinetic energy of one particle (1=electron,2=ion), [eV]
  REAL(8) Avg_kin_energy_x_eV(1:N_spec)              ! - " -, X-motion,
  REAL(8) Avg_kin_energy_y_eV(1:N_spec)              ! - " -, Y-motion, 
  REAL(8) Avg_kin_energy_z_eV(1:N_spec)              ! - " -, Z-motion, 
  REAL(8) Avg_kin_energy_x(1:N_spec)                 ! - " -, X-motion,                                                 ARE USED FOR ACCUMULATION     
  REAL(8) Avg_kin_energy_y(1:N_spec)                 ! - " -, Y-motion,                                                 ARE USED FOR ACCUMULATION 
  REAL(8) Avg_kin_energy_z(1:N_spec)                 ! - " -, Z-motion,                                                 ARE USED FOR ACCUMULATION 

  REAL(8) Energy_wall_eV(1:N_spec)                   ! kinetic energy of (1=electrons,2=ions), which hit the walls, [eV]
  REAL(8) Energy_emit_eV(1:N_spec)                   ! kinetic energy of particles emitted/reflected from the walls (1=electrons,2=ions), [eV]
  REAL(8) Energy_coll_eV(1:N_spec)                   ! kinetic energy of (1=electrons,2=ions), lost in collisions of electrons with neutrals, [eV]

  INTEGER Rate_number_leftwall(1:N_spec)             ! number of (1=electrons,2=ions), which hit the  LEFT wall during the examined time interval
  INTEGER Rate_number_rightwall(1:N_spec)            ! number of (1=electrons,2=ions), which hit the RIGHT wall during the examined time interval
  REAL(8) Rate_number_leftwall_ns1(1:N_spec)         ! - " -, [1/ns]
  REAL(8) Rate_number_rightwall_ns1(1:N_spec)        ! - " -, [1/ns]

  INTEGER Rate_number_leftemit(1:N_spec)             ! number of (1=electrons,2=ions), which were emitted from the  LEFT wall during the examined time interval
  INTEGER Rate_number_rightemit(1:N_spec)            ! number of (1=electrons,2=ions), which were emitted from the RIGHT wall during the examined time interval
  REAL(8) Rate_number_leftemit_ns1(1:N_spec)         ! - " -, [1/ns] 
  REAL(8) Rate_number_rightemit_ns1(1:N_spec)        ! - " -, [1/ns] 

  REAL(8) Yield_left(1:N_spec)                       ! relative emission coefficient (1=electrons,2=ions), LEFT wall
  REAL(8) Yield_right(1:N_spec)                      ! relative emission coefficient (1=electrons,2=ions), RIGHT wall
  REAL(8) eYield_left_elast,  eYield_left_inelast,  eYield_left_true   ! relative coefficients for components of electron emission
  REAL(8) eYield_right_elast, eYield_right_inelast, eYield_right_true  !

  REAL(8) Rate_energy_leftwall(1:N_spec)                  ! kinetic energy of (1=electrons,2=ions), which hit the  LEFT wall 
  REAL(8) Rate_energy_rightwall(1:N_spec)                 ! kinetic energy of (1=electrons,2=ions), which hit the RIGHT wall
  REAL(8) Avg_energy_leftwall_eV(1:N_spec)           ! average kinetic energy of (1=electrons,2=ions), which hit the  LEFT wall, [eV] 
  REAL(8) Avg_energy_rightwall_eV(1:N_spec)          ! average kinetic energy of (1=electrons,2=ions), which hit the RIGHT wall, [eV]

  REAL(8) Rate_energy_leftemit(1:N_spec)                  ! kinetic energy of (1=electrons,2=ions), which were emitted from the  LEFT wall during the time interval
  REAL(8) Rate_energy_rightemit(1:N_spec)                 ! kinetic energy of (1=electrons,2=ions), which were emitted from the RIGHT wall during the time interval
  REAL(8) Avg_energy_leftemit_eV(1:N_spec)           ! average kinetic energy - " -, [eV]
  REAL(8) Avg_energy_rightemit_eV(1:N_spec)          ! average kinetic energy - " -, [eV]

  REAL(8) Energy_heat_eV                        ! energy, transferred into the plasma due to Joule heating (electric current), [eV]

  REAL(8) Rate_energy_full_eVns1                ! rate of change of full energy, [eV/ns]
  REAL(8) Rate_energy_pot_eVns1                 ! rate of change of potential energy, [eV/ns]
  REAL(8) Rate_energy_kin_eVns1(1:N_spec)            ! rate of cange of kinetic energy, 1 = electrons, 2 = ions, [eV/ns]

  REAL(8) Rate_energy_wall_eVns1(1:N_spec)           ! rate of change of (1=electrons,2=ions) energy due to collisions with the walls, [eV/ns]
  REAL(8) Rate_energy_emit_eVns1(1:N_spec)           ! rate of change of energy due to emission of particles from the wall (1=electrons,2=ions), [eV/ns]
  REAL(8) Rate_energy_coll_eVns1(1:N_spec)           ! rate of change of (1=electrons,2=ions) energy due to collisions with neutral particles, [eV/ns]

  REAL(8) Rate_energy_heat_eVns1                ! rate of change of energy due to Joule heating, [eV/ns]

! energy rates, dim-less
  REAL(8) Rate_energy_wall(1:N_spec)                 ! rate of change of energy of (1=electrons,2=ions) due to collisions with the walls
  REAL(8) Rate_energy_emit(1:N_spec)                 ! rate of change of energy of (1=electrons,2=ions) emitted/reflected from the walls 
  REAL(8) Rate_energy_coll(1:N_spec)                 ! rate of change of energy of (1=electrons,2=ions) due to collisions with neutral particles
  REAL(8) Rate_energy_heat                      ! rate of change of energy due to Joule heating

! flow values
  REAL(8) VY_e_avg_ms                           ! averaged (over all particles) Y-velocity for electrons, [m/s]           ARE USED FOR ACCUMULATION
  REAL(8) VZ_e_avg_ms                           ! averaged (over all particles) Z-velocity for electrons, [m/s]           ARE USED FOR ACCUMULATION
  REAL(8) JY_avg_Am2                            ! averaged (over cross section) Y-electric current density, [A/m^2] 
  REAL(8) JZ_avg_Am2                            ! averaged (over cross section) Z-electric current density, [A/m^2]

! density
  REAL(8) N_avg_m3(1:N_spec)                         ! averaged (over cross section) electron (1) and ion(2) density, [1/m^3]

  REAL(8) prev_Q_left   ! Surface charge on the left  plasma-wall boundary at previous diagnostics time-step [dim-less]
  REAL(8) prev_Q_right  ! Surface charge on the right plasma-wall boundary at previous diagnostics time-step [dim-less]

! dim-less quantities accumulated only during one timestep (the last timestep of the averaging period)
  REAL(8) VY_recent(1:N_spec)                        ! for Y-velocity, separate for electrons (1) and ions (2)
  REAL(8) VZ_recent(1:N_spec)                        ! for Z-velocity, separate for electrons (1) and ions (2)

! energy factors
  REAL(8) N_in_macro                            ! number of electrons/ions in one macroparticle
  REAL(8) Factor_energy_eV
  REAL(8) Factor_energy_macro_eV
  REAL(8) Factor_energy_pot_eV
  REAL(8) Factor_rate_ns1
  REAL(8) Factor_rate_macro_ns1
  REAL(8) Factor_rate_eVns1
  REAL(8) Factor_Joule_energy_eV

  REAL(8) Init_energy_full_eV

  INTEGER text_output_counter                   ! is used for skipping text output 

  REAL(8) f_factor(1:N_spec)                         ! factors for calculation of frequencies of collisions of electrons (1) and ions (2) with neutrals

  INTEGER, PARAMETER :: length_of_fbufer=600    ! length of the bufers below   
  REAL(8) fmid_bufer(1:length_of_fbufer)        ! bufer for calculation of time-averaged potential at the midplane
  REAL(8) fwall_bufer(1:length_of_fbufer)       ! bufer for calculation of time-averaged potential at the left wall (x=0)

  

END MODULE Diagnostics

!========================================================================================================================
!
MODULE Snapshots

  INTEGER N_of_all_snaps                        ! number of all snapshots
  INTEGER current_snap                          ! index of current snapshot (which must be created)
  INTEGER counter_of_profiles_to_save           ! number of profiles of the potential to be saved
  CHARACTER(4) snapnumber_txt                   ! prefix in the snapshot files, describes the number of the snapshot

  INTEGER Flag_pp_creation                      ! flag, turns on/off the creation of phase planes using the fraction of all particles
  INTEGER N_to_skip                             ! the fraction mention above is controlled by this parameter (number of particles to skip) for bulk electrons
  INTEGER N_to_skip_left                        ! -""- for electrons emitted at the left wall
  INTEGER N_to_skip_right                       ! -""- for electrons emitted at the right wall
  INTEGER N_to_skip_ion                         ! -""- for ions

  INTEGER N_of_all_vdf_locs                     ! number of locations for calculations of velocity distribution functions (VDF)

  REAL(8) Ve_x_max               ! Maximal x-velocity (normal) for electron velocity distribution (initially given in V_therm)
  REAL(8) Ve_yz_max              ! Maximal y,z-velocity (parallel) for electron velocity distribution (initially given in V_therm)

  INTEGER N_box_Vx_e             ! (Total number of velocity boxes - 1) for positive x-velocity (normal) for electrons
  INTEGER N_box_Vyz_e            ! (Total number of velocity boxes - 1) for positive y,z-velocity (parallel) for electrons
  INTEGER N_box_Vx_e_low         ! = -N_box_Vx_e                   !
  INTEGER N_box_Vx_e_top         ! =  N_box_Vx_e + 1               ! this will save three operations of addition and
  INTEGER N_box_Vyz_e_low        ! = -N_box_Vyz_e                  ! three operations of sign changing for each electron
  INTEGER N_box_Vyz_e_top        ! =  N_box_Vyz_e + 1              !

  INTEGER N_box_Vx_i             ! (Number of velocity boxes - 1) for positive x-velocity (normal) for ions
  INTEGER N_box_Vx_i_low         ! = -N_box_Vx_i                   ! this will save one addition and one sign changing for each ion
  INTEGER N_box_Vx_i_top         ! =  N_box_Vx_i + 1               ! 
  
  INTEGER, ALLOCATABLE :: Tcntr_snapshot(:)     ! timesteps when the snapshot files are written

  INTEGER, ALLOCATABLE :: Vdf_location_bnd(:)   ! Spatial boundaries, defining the locations for calculation of VDFs (numbers of nodes)

  LOGICAL Accumulate_wall_df                    ! flag, showing whether it is necessary (or not) to add particles to the wall distribution functions

! DIAGNOSTIC DISTRIBUTION FUNCTION ARRAYS *******************

  REAL(8), ALLOCATABLE :: evx_mid_of_box(:)     ! middles of   x-velocity boxes (in units of V_th_e), electrons
  REAL(8), ALLOCATABLE :: evyz_mid_of_box(:)    ! middles of y,z-velocity boxes (in units of V_th_e), electrons

  INTEGER, ALLOCATABLE :: ep_2vdf_lw(:,:)       ! 2-d (|v_norm|,|v_par|) distribution of incident (primary) electrons at the  left wall
  INTEGER, ALLOCATABLE :: ep_2vdf_rw(:,:)       ! 2-d (|v_norm|,|v_par|) distribution of incident (primary) electrons at the right wall
  INTEGER, ALLOCATABLE :: es_2vdf_lw(:,:)       ! 2-d (|v_norm|,|v_par|) distribution of emitted electrons at the  left wall
  INTEGER, ALLOCATABLE :: es_2vdf_rw(:,:)       ! 2-d (|v_norm|,|v_par|) distribution of emitted (secondary) electrons at the right wall

! common
!  INTEGER, ALLOCATABLE :: e_2vxvydf_loc(:,:,:)    ! 2-d (v_x,v_y) distribution of electrons at certain location (region) inside the plasma
!  INTEGER, ALLOCATABLE :: e_2vxvzdf_loc(:,:,:)    ! 2-d (v_x,v_z) distribution of electrons at certain location (region) inside the plasma
  INTEGER, ALLOCATABLE :: e_vxdf_loc(:,:)         ! distribution of electrons over v_x at certain location inside the plasma
  INTEGER, ALLOCATABLE :: e_vydf_loc(:,:)         ! distribution of electrons over v_y at certain location inside the plasma
  INTEGER, ALLOCATABLE :: e_vzdf_loc(:,:)         ! distribution of electrons over v_z at certain location inside the plasma
! emitted from the left wall, which do not mix with plasma
!  INTEGER, ALLOCATABLE :: ebl_2vxvydf_loc(:,:,:)  ! 2-d (v_x,v_y) distribution of electrons at certain location (region) inside the plasma
!  INTEGER, ALLOCATABLE :: ebl_2vxvzdf_loc(:,:,:)  ! 2-d (v_x,v_z) distribution of electrons at certain location (region) inside the plasma
  INTEGER, ALLOCATABLE :: ebl_vxdf_loc(:,:)       ! distribution of electrons over v_x at certain location inside the plasma
  INTEGER, ALLOCATABLE :: ebl_vydf_loc(:,:)       ! distribution of electrons over v_y at certain location inside the plasma
  INTEGER, ALLOCATABLE :: ebl_vzdf_loc(:,:)       ! distribution of electrons over v_z at certain location inside the plasma
! emitted from the right wall, which do not mix with plasma
!  INTEGER, ALLOCATABLE :: ebr_2vxvydf_loc(:,:,:)  ! 2-d (v_x,v_y) distribution of electrons at certain location (region) inside the plasma
!  INTEGER, ALLOCATABLE :: ebr_2vxvzdf_loc(:,:,:)  ! 2-d (v_x,v_z) distribution of electrons at certain location (region) inside the plasma
  INTEGER, ALLOCATABLE :: ebr_vxdf_loc(:,:)       ! distribution of electrons over v_x at certain location inside the plasma
  INTEGER, ALLOCATABLE :: ebr_vydf_loc(:,:)       ! distribution of electrons over v_y at certain location inside the plasma
  INTEGER, ALLOCATABLE :: ebr_vzdf_loc(:,:)       ! distribution of electrons over v_z at certain location inside the plasma

  REAL(8), ALLOCATABLE :: ivx_mid_of_box(:)     ! middles of x-velocity boxes (in units of V_th_e), ions
  INTEGER, ALLOCATABLE :: i_vxdf_loc(:,:)       ! distribution of ions over v_x at certain location inside the plasma
  INTEGER, ALLOCATABLE :: i_vydf_loc(:,:)       ! distribution of ions over v_y at certain location inside the plasma
  INTEGER, ALLOCATABLE :: i_vzdf_loc(:,:)       ! distribution of ions over v_z at certain location inside the plasma

END MODULE Snapshots

!----------------------------------
MODULE TestParticles

  INTEGER Leftmost_initial_test_part_node       ! left boundary of the initial distribution of test particles
  INTEGER N_of_initial_test_part_spatial_layers ! number of spatial layers in the initial distribution
  REAL(8) d_spatial                             ! difference between initial positions of neighbor layers

  REAL(8) Min_initial_test_part_energy_eV       ! minimal initial energy of test particles [eV]
  INTEGER N_of_initial_test_part_energy_layers  ! number of energy layers in the initial distribution
  REAL(8) d_energy_eV                           ! difference between initial energies of neighbor layers [eV] 

  INTEGER N_of_test_parts                       ! number of test particles
  INTEGER N_of_test_part_runs                   ! number of times when the test particles will be launched
  INTEGER current_test_run                      ! the number of the current test run

  INTEGER N_of_test_part_save_skip              ! number of timesteps to be skipped between saving test particles
  INTEGER test_part_skip_save_counter           ! counter of timesteps when saving the test particles was skipped

  TYPE TestParticle
     REAL(8) x
     REAL(8) vx
     REAL(8) ax
     INTEGER active              !  1 = active, 0 = non-active
  END TYPE TestParticle

! allocatable arrays
  INTEGER, ALLOCATABLE :: Test_part_run_limits(:,:)    ! initial and final time steps for each run 

  TYPE(TestParticle), ALLOCATABLE :: Test_part(:)            ! test particle parameters

END MODULE TestParticles

