! **************************************************************************************************
!> \brief Maximal Orbital Analysis (MOA)
!> \par History
!>      - module to calculate MOA occupations based on
!         Dupuis, et al., Journal of Computational Chemistry, 2019
!>      01.2020 created
!> \author Pavan Behara
! **************************************************************************************************


MODULE moa
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix,&
                                              write_fm_with_basis_info
   USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_init_random,&
                                              cp_fm_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_set_element,&
                                              cp_fm_to_fm,&
                                              cp_fm_write_formatted,&
                                              cp_fm_write_unformatted
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_deallocate_matrix, dbcsr_get_block_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_set, dbcsr_setname, dbcsr_type
   USE mathlib,                         ONLY: jacobi, diag, eigsrt
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_vals_val_set,&
                                              section_vals_get,&
                                              section_vals_val_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_duplicate
   USE kinds,                           ONLY: default_path_length, &
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: diag
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: nso
   USE particle_methods,                ONLY: get_particle_set
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              mo_set_p_type, &
                                              get_mo_set
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE scf_control_types,               ONLY: scf_control_type
   USE vab,                             ONLY: write_real_2dmatrix_to, &
                                              write_real_1dmatrix_to
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN =   'moa'

   PUBLIC :: maximal_orbital_analysis

CONTAINS
! **************************************************************************************************
!> \brief Perform Maximal Orbital Analysis for atoms and fragments
!>        {Check population analyses routine for better understanding of code
!>        structure}
!> \param qs_env ...quickstep environment post the scf calculation which
!         contains all the information needed to perform MOA and print to cubes
!         if needed
!> \param output_unit ...
!> \param print_level ...
!>      01.2020 created
!> \author Pavan Behara
! **************************************************************************************************
   SUBROUTINE maximal_orbital_analysis(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'maximal_orbital_analysis', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: handle, iatom, ibfn, ii, ispin, iw, icount, &
                                                            ikind, natom, nfrg, nkind, nmo, &
                                                            nsgf, nspin, z_val, ind_count, ifrg, lfrg, z_frg
      INTEGER, DIMENSION(:), POINTER                     :: len_frg, frg_indices, nomoa_in, &
                                                            nomoaa_in, nomoab_in, ind_frg, mo_list
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: nomoa, nomoaa, nomoab
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: real_atomblock, v_atomblock, &
                                                            tmp_block, v_mat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: diag_atomblock, diag_mat, popmoa, chgmoa, spnmoa, &
                                                            popfrg, chgfrg, spnfrg
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: fm_S, fm_T, fm_TI, fm_tmp, &
                                                            fm_O, fm_OI
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrixkp_s
      TYPE(dbcsr_type), POINTER                          :: sm_s
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: mo_print_sec
      LOGICAL                                            :: if_debug

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set)
      NULLIFY (qs_kind_set)
      NULLIFY (fmstruct)
      NULLIFY (fm_S)
      NULLIFY (fm_T)
      NULLIFY (fm_TI)
      NULLIFY (fm_tmp)
      NULLIFY (fm_O)
      NULLIFY (fm_OI)
      NULLIFY (matrixkp_s)
      NULLIFY (mos)
      NULLIFY (particle_set)
      NULLIFY (scf_control)
      NULLIFY (sm_s)
      NULLIFY (para_env)
      NULLIFY (blacs_env)
      NULLIFY (mo_print_sec)
      NULLIFY (real_atomblock)
      NULLIFY (v_atomblock)
      NULLIFY (tmp_block)
      NULLIFY (v_mat)
      NULLIFY (diag_atomblock)
      NULLIFY (diag_mat)
      NULLIFY (popmoa)
      NULLIFY (chgmoa)
      NULLIFY (spnmoa)
      NULLIFY (popfrg)
      NULLIFY (chgfrg)
      NULLIFY (spnfrg)


      if_debug = .FALSE.

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_s_kp=matrixkp_s, &
                      mos=mos, &
                      natom=natom, &
                      particle_set=particle_set, &
                      scf_control=scf_control, &
                      para_env=para_env, &
                      blacs_env=blacs_env)

      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(qs_kind_set))
      CPASSERT(ASSOCIATED(matrixkp_s))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(scf_control))

      IF (SIZE(matrixkp_s, 2) > 1) THEN

         CPWARN("Maximal orbital analysis not implemented for k-points.")

      ELSE
         logger => cp_get_default_logger()
         iw = cp_print_key_unit_nr(logger, qs_env%input,"DFT%SCF%PRINT%ITERATION_INFO", &
                                   extension=".scfLog")

         ! Write headline
         IF (iw > 0) THEN
            WRITE(iw, *) "-----------------------------------------------"
            WRITE(iw, *) "            MAXIMAL ORBITAL ANALYSIS           "
            WRITE(iw, *) "-----------------------------------------------"
         END IF

         sm_s => matrixkp_s(1, 1)%matrix ! Overlap matrix in sparse format
         nspin = SIZE(mos)               ! no. of spins
         nkind = SIZE(atomic_kind_set)   ! no. of atom kinds
         
         ! Input flags to MOA section
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%NFRG", i_val=nfrg)
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%LFRG", i_vals=len_frg)
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%IFRG", i_vals=frg_indices) 
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%NOMOA", i_vals=nomoa_in) 
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%NOMOAA", i_vals=nomoaa_in)
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%NOMOAB", i_vals=nomoab_in)

         ! Check if the MO_CUBES print section is active
         IF (cp_print_key_should_output(logger%iter_info, qs_env%input, &
                                     "PROPERTIES%MOA%MO_CUBES") /= 0) THEN
            mo_print_sec => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%MOA")
            CALL section_vals_val_get(mo_print_sec, "MO_CUBES%MO_LIST", &
                                               i_vals=mo_list)
         END IF

         ! NOMOA is not so relevant parameter in case of HF/DFT calculations.
         ! When doing CAS, etc., the orbital space has to be delineated between
         ! the virtual orbitals, active space and other orbitals. In that case
         ! to maintain the same energy unitary transformation is done only in
         ! respective categories (cannot rotate active space and other orbitals
         ! at the same time which results in a different energy)

         ! for RKS calculations
         CALL get_mo_set(mos(1)%mo_set, nmo=nmo)
         IF((SIZE(nomoa_in) .EQ. 1) .AND. (nomoa_in(1) .EQ. -100000)) THEN
                ALLOCATE(nomoa(nmo))
                nomoa(:) = 0
         ELSE IF (SIZE(nomoa_in) .LT. nmo) THEN
                IF(iw>0) WRITE(iw,*) "Entered the else-if part of nomoa"
                nomoa = [nomoa_in, [(0,ii=1,nmo-SIZE(nomoa_in))]]
         END IF
         IF(iw>0) WRITE(iw,*) "SIZE of NOMOA is ", SIZE(NOMOA), " and nmo is ", nmo
         IF(iw>0) WRITE(iw, *) "Size of Alpha orbs is", nmo

         ! alpha electrons
         IF((SIZE(nomoaa_in) .EQ. 1) .AND. (nomoaa_in(1) .EQ. -100000)) THEN
                ALLOCATE(nomoaa(nmo))
                nomoaa = nomoa
         ELSE IF (SIZE(nomoaa_in) .LT. nmo) THEN
                nomoaa = [nomoaa_in, [(0,ii=1,nmo-SIZE(nomoaa_in))]]
         END IF

         IF(nspin > 1) THEN
                 CALL get_mo_set(mos(2)%mo_set, nmo=nmo)
                 IF(iw>0) WRITE(iw, *) "Size of Beta orbs is", nmo
                 !beta electrons
                 IF((SIZE(nomoab_in) .EQ. 1) .AND. (nomoab_in(1) .EQ. -100000)) THEN
                        ALLOCATE(nomoab(nmo))
                        IF(nmo .LT. SIZE(nomoa)) THEN
                                nomoab = nomoa(1:nmo)
                        ELSE
                                nomoab = nomoa
                        END IF
                 ELSE IF (SIZE(nomoab_in) .LT. nmo) THEN
                        nomoab = [nomoab_in, [(0,ii=1,nmo-SIZE(nomoab_in))]]
                 END IF        
         END IF

         ! Printing out the input section for information
         IF(iw>0) WRITE(iw, *) "INPUT TAGS:"
         IF(iw>0) WRITE(iw,*) "NFRG = ", nfrg
         IF(iw>0) WRITE(iw,*) "LFRG = ", len_frg(:)
         IF(iw>0) WRITE(iw,*) "IFRG = ", frg_indices(:)
         IF(iw>0) WRITE(iw,*) "NOMOA = ", nomoa(:)
         IF(iw>0) WRITE(iw,*) "NOMOAA = ", nomoaa(:)
         IF(iw>0) WRITE(iw,*) "NOMOAB = ", nomoab(:)
         IF(iw>0) WRITE(iw,*) "IF MO_CUBE PRINT IS ACTIVE THEN"
         IF(iw>0) WRITE(iw,*) "MO_LIST = ", mo_list(:)
         IF(iw>0) WRITE(iw,*) "-------------------------"

         CALL get_mo_set(mos(1)%mo_set, nmo=nmo)
         IF(SIZE(nomoa) .NE. nmo) THEN
           IF(iw>0) WRITE(iw,*) "NOMOA must be of SIZE nmo which is ", nmo
         END IF
         IF(SIZE(nomoa) .NE. nmo) CPABORT("CHECK NOMOA INPUT: size must be same as no. of molecular orbitals")
         IF(SIZE(nomoaa) .NE. nmo) CPABORT("CHECK NOMOAA INPUT:")
         IF(nspin > 1) THEN
                CALL get_mo_set(mos(2)%mo_set, nmo=nmo)
                IF(SIZE(nomoab) .NE. nmo) CPABORT("CHECK NOMOAB INPUT:")
         END IF

         ! Get the total number of contracted spherical Gaussian basis functions
         CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

         ! Provide full size work matrices
         ! Overlap matrix S and P*P(t) are of the size (no. of basis fns) x
         ! (no. of basis fns) = (nsgf x nsgf)
         ! nsgf = no. of spherical gaussian functions

         CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                  para_env=para_env, &
                                  context=blacs_env, &
                                  nrow_global=nsgf, &
                                  ncol_global=nsgf)

         ! Creating empty full matrices fm_S with size nsgf by nsgf
         ! fm_ = full matrix
         ! fm_S  = Overlap matrix -S- in Sph. Harm. basis  
         ! S is block orthonormalized to get T 
         ! 
         !  fm_T  = block-orthonormalizing transformation -T- (= (EV)*(1/sqrt(diag))*(EV(t)))
         ! fm_TI = Inverse of transformation matrix T
         !
         ! fm_O  = orthonormalize blocked Orthonormalized -S-
         ! In other words: matrix T, which is block-orthonormal of S, is orthonormalized
         ! fm_OI = inverse of transformation matrix O

         CALL cp_fm_create(matrix=fm_S, &
                           matrix_struct=fmstruct, &
                           name="S MATRIX")
         CALL cp_fm_create(matrix=fm_T, &
                           matrix_struct=fmstruct, &
                           name="T MATRIX")
         CALL cp_fm_set_all(matrix=fm_T, alpha=0.0_dp)

         CALL cp_fm_create(matrix=fm_TI, &
                           matrix_struct=fmstruct, &
                           name="TI MATRIX")
         CALL cp_fm_set_all(matrix=fm_TI, alpha=0.0_dp)

         CALL cp_fm_create(matrix=fm_O, &
                           matrix_struct=fmstruct, &
                           name="O MATRIX")
         CALL cp_fm_set_all(matrix=fm_O, alpha=0.0_dp)

         CALL cp_fm_create(matrix=fm_OI, &
                           matrix_struct=fmstruct, &
                           name="OI MATRIX")
         CALL cp_fm_set_all(matrix=fm_OI, alpha=0.0_dp)

         CALL cp_fm_create(matrix=fm_tmp, &
                           matrix_struct=fmstruct, &
                           name="Temp MATRIX")
         CALL cp_fm_set_all(matrix=fm_tmp, alpha=0.0_dp)

         CALL cp_fm_struct_release(fmstruct=fmstruct)

         !Copy overlap matrix from sparse to full matrix format
         CALL copy_dbcsr_to_fm(sm_s, fm_S)

         ! Loop over all atoms 
         ! Essentially looping over particle_set
         icount = 0
         DO iatom=1,SIZE(particle_set)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
            CALL get_gto_basis_set(gto_basis_set=basis_set, nsgf=ibfn) ! no. of basis functions in iatom

            ALLOCATE(real_atomblock(ibfn, ibfn))
            ALLOCATE(diag_atomblock(ibfn))
            ALLOCATE(v_atomblock(ibfn,ibfn))
            ALLOCATE(tmp_block(ibfn,ibfn))

            CALL cp_fm_get_submatrix(fm_S, real_atomblock, icount+1, &
                                 icount+1, ibfn, ibfn)
            IF(iw>0) WRITE(iw,*) "No. of spherical basis fns is ", ibfn, &
               "for atom kind", ikind, "atom name", qs_kind_set(ikind)%name 
            CALL write_real_2dmatrix_to(real_atomblock, 124, "Atom block")
            ! Jacobi diagonalization of matrix real_atomblock
            ! Syntax diag(n, a, d, v): Diagonalize matrix a. The eigenvalues are returned in vector d
            ! and the eigenvectors are returned in matrix v.
            CALL diag(ibfn, real_atomblock, diag_atomblock, v_atomblock)
            CALL write_real_1dmatrix_to(diag_atomblock, 124, "daig_atomblock")
            CALL write_real_2dmatrix_to(v_atomblock, 124, "v_ block")

            ! Transformation matrix T is (EV)*(1/sqrt(diag))*(EV(t))
            ! Its inverse matrix TI is (EV)*(sqrt(diag))*(EV(t))
            ! where EV is the eigen vectors matrix v_atomblock
           
            DO ii=1,ibfn
               tmp_block(:,ii) = v_atomblock(:,ii)/sqrt(diag_atomblock(ii))
            END DO
            tmp_block = MATMUL(tmp_block, transpose(v_atomblock)) 
            IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "transformation matrix -T-")
            END IF
            CALL cp_fm_set_submatrix(fm_T, tmp_block, icount+1, &
                                 icount+1, ibfn, ibfn)

            DO ii=1,ibfn
               tmp_block(:,ii) = v_atomblock(:,ii)*sqrt(diag_atomblock(ii))
            END DO
            tmp_block = MATMUL(tmp_block, transpose(v_atomblock))
            IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "Inverse transformation matrix -TI-")
            END IF
            CALL cp_fm_set_submatrix(fm_TI, tmp_block, icount+1, &
                                 icount+1, ibfn, ibfn)

            DEALLOCATE(real_atomblock, diag_atomblock, v_atomblock, tmp_block)
            icount = icount + ibfn
         END DO
         
         IF(if_debug) THEN
                CALL cp_fm_write_formatted(fm_S, iw, "S in Sph. Harm. basis")
                CALL cp_fm_write_formatted(fm_T, iw, "fm_T")
                CALL cp_fm_write_formatted(fm_TI, iw, "fm_TI")
         END IF

         ! Doing T*S*T to get S in  Block-Orthonormal Sph. Harmonic Basis
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T, fm_S, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_T, 0.0_dp, fm_S)
         IF(if_debug) THEN
                CALL cp_fm_write_formatted(fm_S, iw, "S in Block-Orthonormal Sph. Harm. basis")
         END IF
         ! -----------

         ALLOCATE(tmp_block(nsgf, nsgf))
         ALLOCATE(diag_mat(nsgf))
         ALLOCATE(v_mat(nsgf, nsgf))
         CALL cp_fm_get_submatrix(fm_S, tmp_block, 1, 1, nsgf, nsgf)
         CALL diag(nsgf, tmp_block, diag_mat, v_mat)
         IF(if_debug) THEN
                CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of block-orthonormal S")
                CALL write_real_2dmatrix_to(v_mat, 124, "eigen vectors of block-orthonormal S")
         END IF

         ! Transformation matrix O is (EV)*(1/sqrt(diag))*(EV(t))
         ! Its inverse matrix OI is (EV)*(sqrt(diag))*(EV(t))
         ! where EV is the eigen vectors matrix v_atomblock

         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)/sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "transformation matrix -O-")
         END IF
         CALL cp_fm_set_submatrix(fm_O, tmp_block, 1, 1, nsgf, nsgf)

         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)*sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "Inverse transformation matrix -OI-")
         END IF
         CALL cp_fm_set_submatrix(fm_OI, tmp_block, 1, 1, nsgf, nsgf)
         IF(if_debug) THEN
                CALL cp_fm_write_formatted(fm_O, iw, "fm_O")
                CALL cp_fm_write_formatted(fm_OI, iw, "fm_OI")
         END IF

         ! Doing O*S*O to get S in  Orthonormalized Block-Orthonormal Sph. Harmonic Basis
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_O, fm_S, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_O, 0.0_dp, fm_S)
         IF(if_debug) THEN
                CALL cp_fm_write_formatted(fm_S, iw, "S in Orthonormalized Block-Orthonormal Sph. Harm. basis")
         END IF
         
         ! Ready to start MOA analysis, we have all the transformations we need
         ALLOCATE(popmoa(natom))
         ALLOCATE(chgmoa(natom))
         ALLOCATE(spnmoa(natom))
         popmoa(:) = 0.0_dp
         chgmoa(:) = 0.0_dp
         spnmoa(:) = 0.0_dp
         
         ! Calling orbmoa routine to do MOA analysis
         ! Inputs are S, V, T, SI, TI, OI
         ! S = Ortho. Block-Ortho. Sph. Harm. basis
         ! V = molecular orbital matrix
         ! T, TI, O, OI = transformations clearly listed above
         ! Outputs U, PP
         ! popmoa = population
         ! chgmoa = charge
         ! spnmoa = spin
         ! nmo = no .of molecular orbitals
         ! for nspin=2 popmoa, chgmoa, spnmoa will be reused to get proper chg,
         ! spn, pop by combining alpha and beta contributions
         
         
         DO ispin=1,nspin
            CALL get_mo_set(mos(ispin)%mo_set, nmo=nmo) 
            WRITE(144,*) "nmo in maximal orbital analysis routine is ", nmo
            CALL orbmoa(fm_S, mos(ispin)%mo_set%mo_coeff, fm_TI, fm_OI, &
                        popmoa, chgmoa, spnmoa, particle_set, qs_kind_set, nmo, nsgf,&
                        nspin, ispin, iw, para_env, blacs_env, if_debug)
         END DO
!
!        ----- PRINT -MOA- CHARGES ... -----
!
         IF(iw>0) THEN
            WRITE(iw,*) '  MOA for Atoms'
            WRITE(iw,*) '  -------------'
            WRITE(iw,*) '         ATOM            NUC       POP', &
                               '         CHRG          SPIN'
            WRITE(iw,*) '     ---------------------------',&
                               '-------------------------' 
            DO iatom=1,natom
               CALL get_atomic_kind(particle_set(iatom)%atomic_kind, element_symbol=element_symbol, z=z_val)
               WRITE(iw," (I9, A9, I9, 3(1X,F12.5)) ")  &
                    iatom, element_symbol, z_val, &
                    popmoa(iatom), chgmoa(iatom), spnmoa(iatom)
            ENDDO
         ENDIF
         
         IF(nfrg .GT. 0) THEN
         ! Call the routine that does MOA analysis by fragment
         ! nfrg, len_frg, frg_indices
            IF(iw>0) WRITE(iw,*) "-------"
            IF(iw>0) WRITE(iw,*) "Entered MOA fragment analysis section"
            ALLOCATE(popfrg(nfrg))
            ALLOCATE(chgfrg(nfrg))
            ALLOCATE(spnfrg(nfrg))
            popfrg(:) = 0.0_dp
            chgfrg(:) = 0.0_dp
            spnfrg(:) = 0.0_dp

            ind_count = 1
            DO ifrg=1,nfrg
               lfrg = len_frg(ifrg)
               IF(iw>0) WRITE(iw,*) "nfrg, lfrg ", nfrg, lfrg
               ALLOCATE(ind_frg(lfrg))
               ind_frg = frg_indices(ind_count:ind_count+lfrg-1)
               IF(lfrg .GT. 0) THEN
                  CALL moax_frg(qs_env, ifrg, lfrg, ind_frg, &
                         nomoa, nomoaa, nomoab, &
                         popfrg, chgfrg, spnfrg, z_frg, if_debug)
               ELSE
                  IF(iw>0) WRITE(iw,*) "ERROR: CHECK INPUT FOR length of fragement LFRG ARGUMENT"
               END IF
               ind_count = ind_count + lfrg 
               DEALLOCATE(ind_frg)
            END DO
         END IF
         IF(iw>0) THEN
            WRITE(iw,*) '  MOA for Fragments'
            WRITE(iw,*) '  -------------'
            WRITE(iw,*) ' FRAG   POP', &
                               '         CHRG          SPIN'
            WRITE(iw,*) '     ---------------------------',&
                               '-------------------------'
            DO ifrg=1,nfrg
               WRITE(iw," (I5, 3(1X,F12.5)) ")  &
                    ifrg, &
                    popfrg(ifrg), chgfrg(ifrg), spnfrg(ifrg)
            ENDDO
         ENDIF

         IF(iw>0) WRITE(iw,*) "..... END of MOA analysis for Fragments ....."
         IF(iw>0) WRITE(iw,*) " "
         IF(iw>0) WRITE(iw,*) "..... END of MOA analysis ....."


         DEALLOCATE(tmp_block, diag_mat, v_mat, popmoa, chgmoa, spnmoa, popfrg, chgfrg, spnfrg) 
         
         
      END IF
      CALL cp_fm_release(fm_S)
      CALL cp_fm_release(fm_T)
      CALL cp_fm_release(fm_TI)
      CALL cp_fm_release(fm_tmp)
      CALL cp_fm_release(fm_O)
      CALL cp_fm_release(fm_OI)
      CALL timestop(handle)
   END SUBROUTINE maximal_orbital_analysis

! **************************************************************************************************
!> \brief Transformations to do MOA for atoms
!>        depending on the number of spins this may be called twice
!> \param fm_S ... Overlap matrix in Ortho. Block-Ortho. Sph. Harm. basis
!> \param fm_V ... Molecular Orbitals at the end of SCF calculation in Sph. Harm. basis
!> \param fm_TI, fm_OI ... transformation matrices listed in 
!>        maximal_orbital_analysis routine
!> \param popmoa, chgmoa, spnmoa ... population, charge and spin from moa orbitals
!> \param particle_set, qs_kind_set ... objects in qs environment to get atom info
!> \param norb ... no. of molecular orbitals (or nmo)
!> \param nsgf ... no. of spherical gaussian functions
!> \param nspin ... no. of spins
!> \param ispin ... i^th spin being passed (required for spnmoa, chgmoa, popmoa updates)
!> \param iw ... write to file descriptor
!> \param para_env, blac_env ... quickstep env variables
!>      01.2020 created
!> \author Pavan Behara
! **************************************************************************************************
   SUBROUTINE orbmoa(fm_S, fm_V, fm_TI, fm_OI, &
                     popmoa, chgmoa, spnmoa, particle_set, &
                     qs_kind_set, norb, nsgf, nspin, ispin, iw, &
                     para_env, blacs_env, if_debug)
      TYPE(cp_fm_type), POINTER                           :: fm_S, fm_V, fm_TI, fm_OI
      REAL(KIND=dp), DIMENSION(:), POINTER, INTENT(INOUT) :: popmoa, chgmoa, spnmoa
      TYPE(particle_type), DIMENSION(:), POINTER          :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER           :: qs_kind_set
      INTEGER, INTENT(IN)                                 :: norb, nsgf, nspin, ispin, iw
      TYPE(cp_para_env_type), POINTER                     :: para_env
      TYPE(cp_blacs_env_type), POINTER                    :: blacs_env
      LOGICAL, INTENT(IN)                                 :: if_debug
      CHARACTER(LEN=*), PARAMETER :: routineN = 'orbmoa', &
         routineP = moduleN//':'//routineN
      INTEGER                                             :: handle, iatom, ibas, &
                                                             ibfn, icount, ikind, &
                                                             occ, z_atom
      TYPE(cp_fm_struct_type), POINTER                    :: fmstruct
      TYPE(cp_fm_type), POINTER                           :: fm_P, fm_PP, fm_V2, fm_tmp
      TYPE(gto_basis_set_type), POINTER                   :: basis_set
      REAL(KIND=dp)                                       :: pop_i
      REAL(KIND=dp), DIMENSION(:, :), POINTER             :: pp_atomblock, U_mat
      REAL(KIND=dp), DIMENSION(:), POINTER                :: diag_mat

      CALL timeset(routineN, handle)

      NULLIFY (fmstruct)
      NULLIFY (fm_P)
      NULLIFY (fm_PP)
      NULLIFY (fm_tmp)
      NULLIFY (fm_V2)
      NULLIFY (pp_atomblock)
      NULLIFY (U_mat)
      NULLIFY (diag_mat)

      IF(nspin .EQ. 1) THEN
        occ = 2
      ELSEIF(nspin .EQ. 2) THEN
        occ = 1
      END IF     

      CALL cp_fm_struct_create(fmstruct=fmstruct, &
                         para_env=para_env, &
                         context=blacs_env, &
                         nrow_global=nsgf, &
                         ncol_global=norb)

      ! Creating empty full matrices fm_tmp, fm_P with size nsgf by norb
      CALL cp_fm_create(matrix=fm_tmp, &
                        matrix_struct=fmstruct, &
                        name="Temp MATRIX")
      CALL cp_fm_set_all(matrix=fm_tmp, alpha=0.0_dp)

      CALL cp_fm_create(matrix=fm_P, &
                        matrix_struct=fmstruct, &
                        name="P = S*V MATRIX")
      CALL cp_fm_struct_release(fmstruct=fmstruct)


      CALL cp_fm_struct_create(fmstruct=fmstruct, &
                         para_env=para_env, &
                         context=blacs_env, &
                         nrow_global=nsgf, &
                         ncol_global=nsgf)

      ! Creating empty full matrices fm_PP with size nsgf by nsgf
      CALL cp_fm_create(matrix=fm_PP, &
                        matrix_struct=fmstruct, &
                        name="PP = P*P(t) MATRIX")
      CALL cp_fm_struct_release(fmstruct=fmstruct)

      CALL cp_fm_create(matrix=fm_V2, &
                        matrix_struct=fm_V%matrix_struct, &
                        name="Working V MATRIX")
      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_V, iw, "V in Sph. Harm. basis")
      END IF

      ! Multiplying orbitals with TI on the left  TI*V to get V in  Block-Orthonormal Sph. Harmonic Basis
      ! placing in a temp matrix since output should not overlap with the inputs
      ! in gemm routine (matrix_c should not overlap matrix_a, matrix_b where c = a*b)
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_TI, fm_V, 0.0_dp, fm_tmp)
      CALL cp_fm_to_fm(fm_tmp, fm_V2)
      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_V2, iw, "Orbitals -V- in Block-Orthonormal Sph. Harm. basis")
      END IF

      ! Multiplying orbitals with OI on the left  OI*V to get V in
      ! Orthonormalized Block-Orthonormal Sph. Harmonic Basis
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_OI, fm_V2, 0.0_dp, fm_tmp)
      CALL cp_fm_to_fm(fm_tmp, fm_V2)
      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_V2, iw, "Orbitals -V- in Orthonormalized Block-Orthonorm Sph. Harm. basis")
      END IF

      ! P = S*V
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_S, fm_V2, 0.0_dp, fm_P)
      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_P, iw, "Matrix P = S*V")
      END IF

      ! PP = P*P(t)
      CALL cp_gemm('N', 'T', nsgf, nsgf, norb, 1.0_dp, fm_P, fm_P, 0.0_dp, fm_PP)
      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_PP, iw, "Matrix PP = P*P(t)")
      END IF
      
      icount = 0

      DO iatom=1,SIZE(particle_set)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind, z=z_atom)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
         CALL get_gto_basis_set(gto_basis_set=basis_set, nsgf=ibfn) ! no. of basis functions in iatom
         
         ! ----- Get the atom block of - P*P(t) - -----
         ALLOCATE(pp_atomblock(ibfn, ibfn))
         CALL cp_fm_get_submatrix(fm_PP, pp_atomblock, icount+1, &
                              icount+1, ibfn, ibfn)
         IF(iw>0) WRITE(iw,*) "No. of spherical basis fns is ", ibfn, &
            "for atom kind", ikind, "atom name", qs_kind_set(ikind)%name
         IF(if_debug) THEN
                CALL write_real_2dmatrix_to(pp_atomblock, 124, "P*P(t) Atom block")
         END IF

         ! ----- Get Maximal Orbitals by diagonalization of - P*P(t) - -----
         ALLOCATE(U_mat(SIZE(pp_atomblock,1),SIZE(pp_atomblock,1)))
         ALLOCATE(diag_mat(MIN(SIZE(pp_atomblock,1),SIZE(pp_atomblock,2))))
   
         ! Doing diagonalization to get maximal overlap, MOA orbitals are in U_mat in block
         ! orthonormal Sph. Harm Basis on iatom
         !CALL jacobi(pp_atomblock, diag_mat, U_mat)
         ! Diagonalize matrix a
         CALL diag(ibfn, pp_atomblock, diag_mat, U_mat)

         ! Sort eigenvalues and eigenvector in ascending order
         CALL eigsrt(ibfn, diag_mat, U_mat)


         IF(if_debug) THEN
                CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of atomblock- P*P(t)")
                CALL write_real_2dmatrix_to(U_mat, 124, "U_mat of SVD of atomblock- P*P(t)")
         END IF

         ! Get MOA atomic population on iatom
         DO ibas=1,ibfn
            IF(ABS(diag_mat(ibas)) .LT. 1E-10) THEN
              diag_mat(ibas) = 0.0_dp
            END IF
            diag_mat(ibas) = SQRT(diag_mat(ibas))
            IF(iw>0) WRITE(iw," (A50, I3, F12.6, 2X, F12.6, I3, F12.6) ") &
                     "-moa- orb. -i- , max. ovlp., occ. number, occ.", &
                     ibas, diag_mat(ibas), diag_mat(ibas)**2, occ, diag_mat(ibas)**2*occ
      
         END DO         
         pop_i = 0.0_dp
         DO ibas=1,ibfn
            pop_i = pop_i + diag_mat(ibas)**2
         END DO
         pop_i = pop_i*occ
         
         IF(ispin .EQ. 1) THEN
           popmoa(iatom) = pop_i
         ELSE
           chgmoa(iatom) = pop_i
         END IF

         IF(ispin .EQ. 1) THEN
           spnmoa(iatom) =  0.0_dp
           popmoa(iatom) =  popmoa(iatom)
           chgmoa(iatom) = -popmoa(iatom) + REAL(z_atom, kind=dp)
         ELSE
           spnmoa(iatom) =  popmoa(iatom) - chgmoa(iatom)
           popmoa(iatom) =  popmoa(iatom) + chgmoa(iatom)
           chgmoa(iatom) = -popmoa(iatom) + REAL(z_atom, kind=dp)
         END IF
     

         icount = icount + ibfn
         DEALLOCATE(pp_atomblock, diag_mat, U_mat)
      END DO

      CALL cp_fm_release(fm_P)
      CALL cp_fm_release(fm_PP)
      CALL cp_fm_release(fm_tmp)
      CALL cp_fm_release(fm_V2)

      CALL timestop(handle)

   END SUBROUTINE orbmoa

! **************************************************************************************************
!> \brief Transformations to do MOA for fragments
!>        called nfrg number of times
!> \param qs_env ... quickstep environment containing all the information
!> \param ifrg ... i^th fragment
!> \param lfrg ... length of fragment
!> \param frg_indices ... indices of atoms in the fragment
!> \param nomoa, nomoaa, nomoab ... exclude this list of orbitals from moa
!> \param popfrg, chgfrg, spnfrg ... population, charge and spin for the fragment
!> \param z_frg ... nuclear charge on the fragment
!>      01.2020 created
!> \author Pavan Behara
! **************************************************************************************************
   SUBROUTINE moax_frg(qs_env, ifrg, lfrg, frg_indices, &
                       nomoa, nomoaa, nomoab, popfrg, chgfrg, spnfrg, z_frg, if_debug)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: ifrg, lfrg
      INTEGER, DIMENSION(:), POINTER                     :: frg_indices
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: nomoa, nomoaa, nomoab 
      REAL(KIND=dp), DIMENSION(:), POINTER               :: popfrg, chgfrg, spnfrg
      INTEGER, INTENT(OUT)                               :: z_frg
      LOGICAL, INTENT(IN)                                :: if_debug
      CHARACTER(LEN=*), PARAMETER :: routineN = 'moax_frg', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: handle, iatom, jatom, iw, natom, &
                                                            nkind, nspin, nsgf, &
                                                            z_val, z_frg1, z_frg2, &
                                                            ikind, jkind, ibfn, jbfn, bf_count, &
                                                            icount, jcount, ii, ispin, nmo
      LOGICAL                                            :: iat_yes, jat_yes
      LOGICAL, DIMENSION(:), POINTER                     :: atom_in_frg, bf_frg
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: real_atomblock, v_atomblock, &
                                                            tmp_block, v_mat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: diag_atomblock, diag_mat
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: fm_S, fm_S_working, fm_S_frg, fm_T, fm_TI, &
                                                            fm_T2, fm_T3, fm_T4, fm_TI2, &
                                                            fm_TI3, fm_TI4, fm_tmp
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrixkp_s
      TYPE(dbcsr_type), POINTER                          :: sm_s
      TYPE(gto_basis_set_type), POINTER                  :: basis_set, basis_set_i, basis_set_j
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: mo_print_sec

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set)
      NULLIFY (qs_kind_set)
      NULLIFY (fmstruct)
      NULLIFY (fm_S)
      NULLIFY (fm_S_frg)
      NULLIFY (fm_T2)
      NULLIFY (fm_TI2)
      NULLIFY (fm_T)
      NULLIFY (fm_TI)
      NULLIFY (fm_T3)
      NULLIFY (fm_TI3)
      NULLIFY (fm_T4)
      NULLIFY (fm_TI4)
      NULLIFY (fm_tmp)
      NULLIFY (matrixkp_s)
      NULLIFY (mos)
      NULLIFY (particle_set)
      NULLIFY (scf_control)
      NULLIFY (sm_s)
      NULLIFY (para_env)
      NULLIFY (blacs_env)
      NULLIFY (mo_print_sec)
      NULLIFY (real_atomblock)
      NULLIFY (v_atomblock)
      NULLIFY (tmp_block) 
      NULLIFY (v_mat)
      NULLIFY (diag_atomblock)
      NULLIFY (diag_mat)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_s_kp=matrixkp_s, &
                      mos=mos, &
                      natom=natom, &
                      particle_set=particle_set, &
                      scf_control=scf_control, &
                      para_env=para_env, &
                      blacs_env=blacs_env)

      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(qs_kind_set))
      CPASSERT(ASSOCIATED(matrixkp_s))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(scf_control))

      ! Get the total number of contracted spherical Gaussian basis functions
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

      IF (SIZE(matrixkp_s, 2) > 1) THEN

         CPWARN("Maximal orbital analysis not implemented for k-points.")

      ELSE
         logger => cp_get_default_logger()
         iw = cp_print_key_unit_nr(logger, qs_env%input,"DFT%SCF%PRINT%ITERATION_INFO", &
                                   extension=".scfLog")

         ! Write headline
         IF (iw > 0) THEN
            WRITE(iw, *) "-----------------------------------------------"
            WRITE(iw, *) "            MOA Fragment analysis           "
            WRITE(iw, *) "-----------------------------------------------"
         END IF


         sm_s => matrixkp_s(1, 1)%matrix ! Overlap matrix in sparse format
         nspin = SIZE(mos)
         nkind = SIZE(atomic_kind_set)
         
         ! Checking if atom belongs to the fragment or not
         ALLOCATE(atom_in_frg(natom))
         atom_in_frg(:) = .FALSE.
      
         DO iatom=1,lfrg
            atom_in_frg(frg_indices(iatom))=.TRUE.
         END DO
         IF(iw>0) WRITE(iw,*) "Atoms in fragement, ", ifrg," are: ",frg_indices(:)
          
         ! Nuclear charge of the fragment and rest of the system
         ! z_frg1 is for the defined fragment
         ! z_frg2 is for rest of the system
         z_frg1 = 0
         z_frg2 = 0
         DO iatom=1,natom
            IF (atom_in_frg(iatom)) THEN
               CALL get_atomic_kind(particle_set(iatom)%atomic_kind, element_symbol=element_symbol, z=z_val)
               z_frg1 = z_frg1+z_val
            ELSE 
               CALL get_atomic_kind(particle_set(iatom)%atomic_kind, element_symbol=element_symbol, z=z_val)
               z_frg2 = z_frg2+z_val
            END IF
         ENDDO
         z_frg = z_frg1

         IF(iw>0) WRITE(iw,*) "Total nuclear charge on fragement, ", ifrg," is: ", z_frg1
         IF(iw>0) WRITE(iw,*) "Total nuclear charge on rest of the system is: ", z_frg2

         ! Basis functions belong to fragment or not
         ! bf_frg has the truth value for all nsgf basis fns
         ALLOCATE(bf_frg(nsgf))
         bf_count = 1
         DO iatom=1,natom
            IF (atom_in_frg(iatom)) THEN
                CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
                CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
                CALL get_gto_basis_set(gto_basis_set=basis_set, nsgf=ibfn) ! no. of basis functions in iatom
                bf_frg(bf_count:bf_count+ibfn-1) = .TRUE.
                bf_count = bf_count+ibfn
            ELSE 
                CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
                CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
                CALL get_gto_basis_set(gto_basis_set=basis_set, nsgf=ibfn) ! no. of basis functions in iatom
                bf_frg(bf_count:bf_count+ibfn-1) = .FALSE.
                bf_count = bf_count+ibfn
            END IF
         END DO
         IF(iw>0) WRITE(iw,*) "Basis functions in the fragment or not, for frag. #", ifrg, bf_frg(:)


         ! Provide full size work matrices
         ! Overlap matrix S and P*P(t) are of the size (no. of basis fns) x
         ! (no. of basis fns)
         ! nsgf = no. of spherical gaussian functions

         CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                  para_env=para_env, &
                                  context=blacs_env, &
                                  nrow_global=nsgf, &
                                  ncol_global=nsgf)

         ! Creating empty full matrices fm_S with size nsgf by nsgf
         CALL cp_fm_create(matrix=fm_S, &
                           matrix_struct=fmstruct, &
                           name="S MATRIX")
         CALL cp_fm_create(matrix=fm_S_frg, &
                           matrix_struct=fmstruct, &
                           name="S-frg MATRIX")
         CALL cp_fm_create(matrix=fm_S_working, &
                           matrix_struct=fmstruct, &
                           name="S-working MATRIX")

         !     -T2, TI2-  to     Atom Ortho Sph. Harm. 
         !     -T3, TI3-  to     Frag Ortho Atom Ortho Sph. Harm. 
         !     -T4, TI4-  to     Molec Ortho Frag Ortho Atom Ortho Sph. Harm. 

         CALL cp_fm_create(matrix=fm_T2, &
                           matrix_struct=fmstruct, &
                           name="T2 MATRIX")
         CALL cp_fm_set_all(matrix=fm_T2, alpha=0.0_dp)

         CALL cp_fm_create(matrix=fm_TI2, &
                           matrix_struct=fmstruct, &
                           name="TI2 MATRIX")
         CALL cp_fm_set_all(matrix=fm_TI2, alpha=0.0_dp)

         CALL cp_fm_create(matrix=fm_T, &
                           matrix_struct=fmstruct, &
                           name="T MATRIX")
         CALL cp_fm_set_all(matrix=fm_T, alpha=0.0_dp)

         CALL cp_fm_create(matrix=fm_TI, &
                           matrix_struct=fmstruct, &
                           name="TI MATRIX")
         CALL cp_fm_set_all(matrix=fm_TI, alpha=0.0_dp)


         CALL cp_fm_create(matrix=fm_T3, &
                           matrix_struct=fmstruct, &
                           name="T3 MATRIX")
         CALL cp_fm_set_all(matrix=fm_T3, alpha=0.0_dp)

         CALL cp_fm_create(matrix=fm_TI3, &
                           matrix_struct=fmstruct, &
                           name="TI3 MATRIX")
         CALL cp_fm_set_all(matrix=fm_TI3, alpha=0.0_dp)

         CALL cp_fm_create(matrix=fm_T4, &
                           matrix_struct=fmstruct, &
                           name="T4 MATRIX")
         CALL cp_fm_set_all(matrix=fm_T4, alpha=0.0_dp)

         CALL cp_fm_create(matrix=fm_TI4, &
                           matrix_struct=fmstruct, &
                           name="TI4 MATRIX")
         CALL cp_fm_set_all(matrix=fm_TI4, alpha=0.0_dp)

         CALL cp_fm_create(matrix=fm_tmp, &
                           matrix_struct=fmstruct, &
                           name="Temp MATRIX")
         CALL cp_fm_set_all(matrix=fm_tmp, alpha=0.0_dp)

         CALL cp_fm_struct_release(fmstruct=fmstruct)

         !Copy overlap matrix from sparse to full matrix format
         CALL copy_dbcsr_to_fm(sm_s, fm_S)

         ! Loop over all atoms
         ! Essentially looping over particle_set
         icount = 0
         DO iatom=1,SIZE(particle_set)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
            CALL get_gto_basis_set(gto_basis_set=basis_set, nsgf=ibfn) ! no. of basis functions in iatom
            ALLOCATE(real_atomblock(ibfn, ibfn))
            ALLOCATE(diag_atomblock(ibfn))
            ALLOCATE(v_atomblock(ibfn,ibfn))
            ALLOCATE(tmp_block(ibfn,ibfn))
            CALL cp_fm_get_submatrix(fm_S, real_atomblock, icount+1, &
                                 icount+1, ibfn, ibfn)
            IF(iw>0) WRITE(iw,*) "No. of spherical basis fns is ", ibfn, &
               "for atom kind", ikind, "atom name", qs_kind_set(ikind)%name
            IF(if_debug) THEN
                CALL write_real_2dmatrix_to(real_atomblock, 124, "Atom block")
            END IF
            ! Jacobi diagonalization of matrix real_atomblock
            ! Syntax diag(n, a, d, v): Diagonalize matrix a. The eigenvalues are returned in vector d
            ! and the eigenvectors are returned in matrix v.
            CALL diag(ibfn, real_atomblock, diag_atomblock, v_atomblock)
            IF(if_debug) THEN
                CALL write_real_1dmatrix_to(diag_atomblock, 124, "diag_atomblock")
                CALL write_real_2dmatrix_to(v_atomblock, 124, "v_ block")
            END IF

            ! Transformation matrix T is (EV)*(1/sqrt(diag))*(EV(t))
            ! Its inverse matrix TI is (EV)*(sqrt(diag))*(EV(t))
            ! where EV is the eigen vectors matrix v_atomblock
            DO ii=1,ibfn
               tmp_block(:,ii) = v_atomblock(:,ii)/sqrt(diag_atomblock(ii))
            END DO
            tmp_block = MATMUL(tmp_block, transpose(v_atomblock))
            IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "transformation matrix -T-")
            END IF
            CALL cp_fm_set_submatrix(fm_T2, tmp_block, icount+1, &
                                 icount+1, ibfn, ibfn)

            DO ii=1,ibfn
               tmp_block(:,ii) = v_atomblock(:,ii)*sqrt(diag_atomblock(ii))
            END DO
            tmp_block = MATMUL(tmp_block, transpose(v_atomblock))
            IF(if_debug) THEN 
                CALL write_real_2dmatrix_to(tmp_block, 124, "Inverse transformation matrix -TI-")
            END IF
            CALL cp_fm_set_submatrix(fm_TI2, tmp_block, icount+1, &
                                 icount+1, ibfn, ibfn)

            DEALLOCATE(real_atomblock, diag_atomblock, v_atomblock, tmp_block)
            icount = icount + ibfn
         END DO
        
         IF(if_debug) THEN
                CALL cp_fm_write_formatted(fm_S, iw, "S in Atom Sph. Harm. basis")
                CALL cp_fm_write_formatted(fm_T2, iw, "fm_T -T2- to Atom Ortho Sph. Harm.")
                CALL cp_fm_write_formatted(fm_TI2, iw, "fm_TI -TI2- to Atom Ortho Sph. Harm.")
         END IF

         ! Doing T*S*T to get S in  Atom Ortho Sph. Harmonic Basis
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T2, fm_S, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_T2, 0.0_dp, fm_S)
         IF(if_debug) THEN
                CALL cp_fm_write_formatted(fm_S, iw, "S in Atom Ortho Sph. Harm. basis")
         END IF

         ALLOCATE(tmp_block(nsgf, nsgf))
         ALLOCATE(diag_mat(nsgf))
         ALLOCATE(v_mat(nsgf, nsgf))
         CALL cp_fm_get_submatrix(fm_S, tmp_block, 1, 1, nsgf, nsgf)
         CALL diag(nsgf, tmp_block, diag_mat, v_mat)
         IF(if_debug) THEN
                CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of S in Atom Ortho Sph. Harm.")
                CALL write_real_2dmatrix_to(v_mat, 124, "eigen vectors of S in Atom Ortho Sph. Harm.")
         END IF

         !!' -T- to Molec Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)/sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "-T- to Molec Ortho Atom Ortho Sph. Harm. basis")
         END IF
         CALL cp_fm_set_submatrix(fm_T, tmp_block, 1, 1, nsgf, nsgf)

         !!' -TI- (inv.) to Molec Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)*sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "-TI- Molec Ortho Atom Ortho Sph. Harm. basis")
         END IF
         CALL cp_fm_set_submatrix(fm_TI, tmp_block, 1, 1, nsgf, nsgf)


         ! Doing T*S*T for Frag. S block orthonromal to get Frag. S in  Orthonormalized Sph. Harmonic Basis
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T, fm_S, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_T, 0.0_dp, fm_S_working)
         IF(if_debug) THEN
                CALL cp_fm_write_formatted(fm_S_working, iw, "S in Molec Ortho Atom Ortho Sph. Harm. basis")
         END IF

         !!' eigen vectors and values of -S- in Molec Ortho Atom Ortho Sph. Harm. basis'
         CALL cp_fm_get_submatrix(fm_S_working, tmp_block, 1, 1, nsgf, nsgf)
         CALL diag(nsgf, tmp_block, diag_mat, v_mat)
         IF(if_debug) THEN
                CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of S in Molec Ortho Atom Ortho Sph. Harm. basis")
                CALL write_real_2dmatrix_to(v_mat, 124, "eigen vectors of S in Molec Ortho Atom Ortho Sph. Harm. basis")
         END IF


         ! ----- ... On the way to Frag Ortho Atom Ortho Sph. Harm. basis -----

         !! -S- in Atom Ortho Sph. Harm. basis
         CALL cp_fm_get_submatrix(fm_S, tmp_block, 1, 1, nsgf, nsgf)
         icount = 0
         DO iatom=1,natom
            iat_yes = atom_in_frg(iatom)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set_i)
            CALL get_gto_basis_set(gto_basis_set=basis_set_i, nsgf=ibfn) ! no. of basis functions in iatom
            jcount = 0
            DO jatom=1,natom
               jat_yes = atom_in_frg(jatom)
               CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
               CALL get_qs_kind(qs_kind_set(jkind), basis_set=basis_set_j)
               CALL get_gto_basis_set(gto_basis_set=basis_set_j, nsgf=jbfn) ! no. of basis functions in jatom
               IF(iw>0) WRITE(iw,*) "IAT, JAT = ", iatom, jatom, iat_yes, jat_yes
               ! Zero out (iatom, jatom) block if iat_yes or jat_yes is false or
               ! iatom.NE.jatom
               IF(iat_yes.AND.jat_yes) THEN
               
               ELSE
                 IF(iatom .NE. jatom) THEN
                   tmp_block(icount+1:icount+ibfn, jcount+1:jcount+jbfn) = 0.0_dp
                 END IF
               END IF
               jcount = jcount + jbfn
            END DO
            icount = icount + ibfn 
         END DO
         
         IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "S in Frag. blocked Atom Ortho Sph. Harm. basis")
         END IF
         CALL diag(nsgf, tmp_block, diag_mat, v_mat)
         IF(if_debug) THEN
                CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of S in Frag. blocked Atom Ortho Sph. Harm. basis")
                CALL write_real_2dmatrix_to(v_mat, 124, "eigen vectors of S in Frag. blocked Atom Ortho Sph. Harm. basis")
         END IF
         CALL cp_fm_set_submatrix(fm_S_frg, tmp_block, 1, 1, nsgf, nsgf)
      
         ! Transformation matrix T is (EV)*(1/sqrt(diag))*(EV(t))
         ! Its inverse matrix TI is (EV)*(sqrt(diag))*(EV(t))
         ! where EV is the eigen vectors matrix S Frag. block-orthonormal
       
         !!' -T- to Frag Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)/sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "-T3- to Frag Ortho Atom Ortho Sph. Harm. basis")
          END IF
         CALL cp_fm_set_submatrix(fm_T3, tmp_block, 1, 1, nsgf, nsgf)

         !!' -TI- (inv.) to Frag Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)*sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "-TI3- Frag Ortho Atom Ortho Sph. Harm. basis")
         END IF
         CALL cp_fm_set_submatrix(fm_TI3, tmp_block, 1, 1, nsgf, nsgf)


         ! Doing T*S*T for Frag. S block orthonromal to get Frag. S in  Orthonormalized Sph. Harmonic Basis
         ! Note: Here S is in Atom Ortho Sph. Harm. Basis
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T3, fm_S, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_T3, 0.0_dp, fm_S_frg)
         IF(if_debug) THEN
                CALL cp_fm_write_formatted(fm_S_frg, iw, "S in Frag Ortho Atom Ortho Sph. Harm. basis")
         END IF

         !!' eigen vectors and values of -S- in Frag Ortho Atom Ortho Sph. Harm. basis'
         CALL cp_fm_get_submatrix(fm_S_frg, tmp_block, 1, 1, nsgf, nsgf)
         CALL diag(nsgf, tmp_block, diag_mat, v_mat)
         IF(if_debug) THEN
                CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of S in Frag Ortho Atom Ortho Sph. Harm. basis")
                CALL write_real_2dmatrix_to(v_mat, 124, "eigen vectors of S in Frag Ortho Atom Ortho Sph. Harm. basis")
         END IF

         !!' -T- to Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)/sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "-T4- to Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis")
         END IF
         CALL cp_fm_set_submatrix(fm_T4, tmp_block, 1, 1, nsgf, nsgf)

         !!' -TI- (inv.) to Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)*sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         IF(if_debug) THEN
                CALL write_real_2dmatrix_to(tmp_block, 124, "-TI4- to Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis")
         END IF
         CALL cp_fm_set_submatrix(fm_TI4, tmp_block, 1, 1, nsgf, nsgf)

         ! Doing T*S*T for Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis 
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T4, fm_S_frg, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_T4, 0.0_dp, fm_S_frg)
         IF(if_debug) THEN
                CALL cp_fm_write_formatted(fm_S_frg, iw, "-S- in Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis")
         END IF

         popfrg(ifrg) = 0.0_dp
         chgfrg(ifrg) = 0.0_dp
         spnfrg(ifrg) = 0.0_dp
         
         ! Calling orbmoa_frg routine to do MOA analysis

         mo_print_sec => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%MOA")
         DO ispin=1,nspin
            IF(iw>0) THEN
                WRITE(iw,*) " ----- now ... the analysis -----"
                WRITE(iw,*) " -MOA- analysis for -fragment- =", ifrg
                WRITE(iw,*) " -------------------------------"
                WRITE(iw,*) "---------------------"
                IF(ispin.EQ.1) WRITE(iw,*) "alpha electrons"
                IF(ispin.EQ.2) WRITE(iw,*) "beta electrons"
                WRITE(iw,*) "---------------------"
            END IF
            CALL get_mo_set(mos(ispin)%mo_set, nmo=nmo)
            WRITE(144,*) "nmo is ", nmo
            !     -T2, TI2-  to     Atom Ortho Sph. Harm.
            !     -T3, TI3-  to     Frag Ortho Atom Ortho Sph. Harm.
            !     -T4, TI4-  to     Molec Ortho Frag Ortho Atom Ortho Sph. Harm.
            
            ! alpha electrons exclude orbs
            IF(ispin .EQ. 1) THEN
                IF(SIZE(nomoa) .EQ. SIZE(nomoaa)) THEN
                        nomoa = nomoaa
                ELSE
                        DEALLOCATE(nomoa)
                        ALLOCATE(nomoa(SIZE(nomoaa)))
                        nomoa = nomoaa
                END IF
            END IF
            
            ! beta electrons exclude orbs
            IF(ispin .EQ. 2) THEN
                IF(SIZE(nomoa) .EQ. SIZE(nomoab)) THEN
                        nomoa = nomoab
                ELSE
                        DEALLOCATE(nomoa)
                        ALLOCATE(nomoa(SIZE(nomoab)))
                        nomoa = nomoab
                END IF
            END IF

            CALL orbmoa_frg(fm_S_frg, mos(ispin)%mo_set%mo_coeff, fm_T2, fm_TI2, &
                            fm_T3, fm_TI3, fm_T4, fm_TI4, &
                            popfrg, chgfrg, spnfrg, bf_frg, z_frg1, nmo, nsgf,&
                            nspin, ispin, ifrg, iw, nomoa, para_env, blacs_env, &
                            mo_print_sec, qs_env, if_debug)
         END DO
         
         DEALLOCATE(atom_in_frg, bf_frg, tmp_block, diag_mat, v_mat)
      END IF 
      CALL cp_fm_release(fm_S)
      CALL cp_fm_release(fm_S_frg)
      CALL cp_fm_release(fm_T2)
      CALL cp_fm_release(fm_TI2)
      CALL cp_fm_release(fm_T)
      CALL cp_fm_release(fm_TI)
      CALL cp_fm_release(fm_T3)
      CALL cp_fm_release(fm_TI3)
      CALL cp_fm_release(fm_T4)
      CALL cp_fm_release(fm_TI4)
      CALL cp_fm_release(fm_tmp)
      CALL timestop(handle)

   END SUBROUTINE moax_frg

! **************************************************************************************************
!> \brief Transformations to do MOA for fragments
!> \param fm_S ... Overlap matrix in Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis
!> \param fm_V ... molecular orbitals in Sph. Harm. basis 
!>        (transformations are performed on this in this routine to match S as 
!>           Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis)
!> \param fm_TI2, TI3, TI4 ... transformation matrices from the calling routine
!> \param popfrg, chgfrg, spnfrg ... population, charge and spin for the fragment
!> \param bf_frg ... basis functions that belong to the fragment
!> \param z_frg ... nuclear charge on the fragment
!> \param norb ... no. of molecular orbitals (or nmo)
!> \param nsgf ... no. of spherical gaussian functions
!> \param nspin ... no. of spins
!> \param ispin ... i^th spin being passed (required for spnmoa, chgmoa, popmoa updates)
!> \param ifrg ... i^th fragment
!> \param iw ... write to file descriptor
!> \param nomoa ... list of molecular orbitals excluded from moa
!> \param para_env, blac_env ... quickstep env variables
!> \param mo_print_sec ... print section pointing to MO_CUBES in the input hierarchy
!> \param qs_env ... quickstep environment holding all the information 
!>      01.2020 created
!> \author Pavan Behara
! **************************************************************************************************
   SUBROUTINE orbmoa_frg(fm_S, fm_V, fm_T2, fm_TI2, fm_T3, fm_TI3, fm_T4, fm_TI4, &
                     popfrg, chgfrg, spnfrg, bf_frg, z_frg,  &
                     norb, nsgf, nspin, ispin, ifrg, iw, nomoa, &
                     para_env, blacs_env, mo_print_sec, qs_env, if_debug)
      TYPE(cp_fm_type), POINTER                           :: fm_S, fm_V, fm_T2, fm_TI2, fm_T3, fm_TI3, &
                                                             fm_T4, fm_TI4
      REAL(KIND=dp), DIMENSION(:), POINTER, INTENT(INOUT) :: popfrg, chgfrg, spnfrg
      LOGICAL, DIMENSION(:), POINTER                      :: bf_frg
      INTEGER                                             :: z_frg
      INTEGER, INTENT(IN)                                 :: norb, nsgf, nspin, ispin, ifrg, iw
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: nomoa
      TYPE(cp_para_env_type), POINTER                     :: para_env
      TYPE(cp_blacs_env_type), POINTER                    :: blacs_env
      TYPE(section_vals_type), POINTER                    :: mo_print_sec
      TYPE(qs_environment_type), POINTER                  :: qs_env
      LOGICAL, INTENT(IN)                                 :: if_debug
      CHARACTER(LEN=*), PARAMETER :: routineN = 'orbmoa_frg', &
         routineP = moduleN//':'//routineN
      INTEGER                                             :: handle, ibfn, &
                                                             jbfn, jorb, icount,  &
                                                             occ, imoa, nbf
      TYPE(cp_fm_struct_type), POINTER                    :: fmstruct
      TYPE(cp_fm_type), POINTER                           :: fm_P, fm_PP, fm_V2, fm_tmp, &
                                                             fm_rt1, fm_rt2, moa_coeff, fm_U_sph_harm_bas
      REAL(KIND=dp)                                       :: dum, pop_i, norm_of_diag
      REAL(KIND=dp), DIMENSION(:, :), POINTER             :: PP, U_mat, U_mat_truncated
      REAL(KIND=dp), DIMENSION(:), POINTER                :: diag_mat
      REAL(KIND=dp), EXTERNAL                             :: dnrm2
      TYPE(cp_logger_type), POINTER                       :: logger
      INTEGER, DIMENSION(:), POINTER                      :: mo_list

      CALL timeset(routineN, handle)

      NULLIFY (fmstruct)
      NULLIFY (fm_P)
      NULLIFY (fm_PP)
      NULLIFY (fm_tmp)
      NULLIFY (fm_rt1) !tmp array for reverse transforms
      NULLIFY (fm_rt2)
      NULLIFY (moa_coeff)
      NULLIFY (fm_U_sph_harm_bas)
      NULLIFY (PP)
      NULLIFY (U_mat)
      NULLIFY (U_mat_truncated)
      NULLIFY (diag_mat)

      IF(nspin .EQ. 1) THEN
        occ = 2
      ELSEIF(nspin .EQ. 2) THEN
        occ = 1
      END IF

      CALL cp_fm_struct_create(fmstruct=fmstruct, &
                         para_env=para_env, &
                         context=blacs_env, &
                         nrow_global=nsgf, &
                         ncol_global=norb)

      ! Creating empty full matrices fm_tmp, fm_P with size nsgf by norb
      CALL cp_fm_create(matrix=fm_tmp, &
                        matrix_struct=fmstruct, &
                        name="Temp MATRIX")
      CALL cp_fm_set_all(matrix=fm_tmp, alpha=0.0_dp)

      CALL cp_fm_create(matrix=fm_P, &
                        matrix_struct=fmstruct, &
                        name="P = S*V MATRIX")
      CALL cp_fm_struct_release(fmstruct=fmstruct)


      CALL cp_fm_struct_create(fmstruct=fmstruct, &
                         para_env=para_env, &
                         context=blacs_env, &
                         nrow_global=nsgf, &
                         ncol_global=nsgf)

      ! Creating empty full matrices fm_PP with size nsgf by nsgf
      CALL cp_fm_create(matrix=fm_PP, &
                        matrix_struct=fmstruct, &
                        name="PP = P*P(t) MATRIX")
      CALL cp_fm_create(matrix=fm_U_sph_harm_bas, &
                        matrix_struct=fmstruct, &
                        name="U_mat in Sph. Harm. basis")
      CALL cp_fm_create(matrix=fm_rt1, &
                        matrix_struct=fmstruct, &
                        name="TEMP MATRIX1 of size nsgf by nsgf")
      CALL cp_fm_create(matrix=fm_rt2, &
                        matrix_struct=fmstruct, &
                        name="TEMP MATRIX2 of size nsgf by nsgf")
      CALL cp_fm_struct_release(fmstruct=fmstruct)

      CALL cp_fm_create(matrix=fm_V2, &
                        matrix_struct=fm_V%matrix_struct, &
                        name="Working V MATRIX")
      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_V, iw, "V in Sph. Harm. basis")
      END IF
      ! Multiplying orbitals with TI2 on the left  TI2*V to get V in  Atom Ortho Sph. Harmonic Basis
      ! placing in a temp matrix since output should not overlap with the inputs
      ! in gemm routine (matrix_c should not overlap matrix_a, matrix_b where c = a*b)
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_TI2, fm_V, 0.0_dp, fm_tmp)
      CALL cp_fm_to_fm(fm_tmp, fm_V2)
      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_V2, iw, "Orbitals -V- in Atom Ortho Sph. Harm. basis")
      END IF

      ! Multiplying orbitals with TI3 on the left  TI3*V to get 
      ! V in Frag Ortho Atom Ortho Sph. Harmonic Basis
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_TI3, fm_V2, 0.0_dp, fm_tmp)
      CALL cp_fm_to_fm(fm_tmp, fm_V2)
      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_V2, iw, "Orbitals -V- in Frag Ortho Atom Ortho Sph. Harm. basis")
      END IF

      ! Multiplying orbitals with TI4 on the left  TI4*V to get V in
      ! Molec Ortho Frag Ortho Atom Ortho Sph. Harmonic Basis
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_TI4, fm_V2, 0.0_dp, fm_tmp)
      CALL cp_fm_to_fm(fm_tmp, fm_V2)
      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_V2, iw, "Orbitals -V- in Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis")

        CALL cp_fm_write_formatted(fm_S, iw, "-S- in Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis")
      END IF

      IF(iw>0) WRITE(iw,*) " ----- now ... the analysis -----"
      IF(iw>0) WRITE(iw,*) " -MOA- analysis for -fragment- =", ifrg 
      IF(iw>0) WRITE(iw,*) " -------------------------------"

      ! ----- calculate P ----
      ! P = S*V for all
      ! -S, V - are in Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_S, fm_V2, 0.0_dp, fm_P)
      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_P, iw, "Matrix P = S*V for all")
      END IF

      ! P_fragement
      DO ibfn = 1,nsgf
         DO jorb = 1,norb
            IF(.NOT.bf_frg(ibfn)) THEN
               CALL cp_fm_set_element(fm_P, ibfn, jorb, 0.0_dp)
            END IF
         END DO
      END DO
      IF(if_debug) THEN       
        CALL cp_fm_write_formatted(fm_P, iw, "P_fragment")
      END IF

      ! P_frag. for non-excluded molec. orbitals
      DO jorb = 1,norb
         WRITE(144,*) "jorb and nomoa(jorb) are", jorb, nomoa(jorb)
         IF(nomoa(jorb).NE.0) THEN
           DO ibfn = 1,nsgf
               CALL cp_fm_set_element(fm_P, ibfn, jorb, 0.0_dp)
           END DO
         END IF
      END DO

      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_P, iw, "P_fragment for non-excluded orbitals")
      END IF
      
      ! PP = P*P(t)
      CALL cp_gemm('N', 'T', nsgf, nsgf, norb, 1.0_dp, fm_P, fm_P, 0.0_dp, fm_PP)
      ALLOCATE(PP(nsgf, nsgf))
      ! Prepare -P*P(t)- 
      CALL cp_fm_get_submatrix(fm_PP, PP, 1, 1, nsgf, nsgf)
      DO ibfn = 1,nsgf
         DO jbfn = 1,ibfn
            IF(bf_frg(ibfn) .AND. bf_frg(jbfn)) THEN
                dum = PP(ibfn, jbfn)
            ELSE
              IF(jbfn.EQ.ibfn) THEN
                dum = 10.0_dp
              ELSE
                dum = 0.0_dp
              END IF
            END IF
            PP(ibfn, jbfn) = dum
         END DO
      END DO

      IF(if_debug) THEN
        CALL write_real_2dmatrix_to(PP, 124, "P*P(t) for fragment and non-excluded Molec Orbs")
      END IF  

      CALL cp_fm_set_submatrix(fm_PP, PP, 1, 1, nsgf, nsgf)

      IF(if_debug) THEN
        CALL cp_fm_write_formatted(fm_PP, iw, "Matrix PP = P*P(t)")
      END IF

      ! Get Maximal Orbitals -U-
      ALLOCATE(U_mat(SIZE(PP,1),SIZE(PP,1)))
      ALLOCATE(diag_mat(MIN(SIZE(PP,1),SIZE(PP,2))))
      ! Doing diagonalization to get maximal overlap, in -MOA- orbs -U- in Molec Ortho Frag Ortho ---
      ! --- Atom Ortho Sph. Harm. basis
      U_mat(:, :) = 0.0_dp
      diag_mat(:) = 0.0_dp
      IF(if_debug) THEN
        CALL write_real_1dmatrix_to(diag_mat, 124, "before jacobi diag_mat for fragment MOA orbs")
        CALL write_real_2dmatrix_to(U_mat, 124, "before jacobi U_mat for fragment MOA orbs")
      END IF

      !!CALL jacobi(PP, diag_mat, U_mat)
      ! Diagonalize matrix a
      CALL diag(nsgf, PP, diag_mat, U_mat)
      IF(if_debug) THEN
        CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of PP after diag before eigsrt for fragment MOA orbs")
        CALL write_real_2dmatrix_to(U_mat, 124, "U_mat of SVD of P*P(t) after diag before eigsrt for fragment MOA orbs")
      END IF

      ! Sort eigenvalues and eigenvector in ascending order
      CALL eigsrt(nsgf, diag_mat, U_mat)

      IF(if_debug) THEN
        CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of P*P(t) for fragment MOA orbs")
        CALL write_real_2dmatrix_to(U_mat, 124, "U_mat of SVD of P*P(t) for fragment MOA orbs")
      END IF
     
      ! Get MOA fragment population 
      icount = 0
      DO ibfn = 1,nsgf
         IF(diag_mat(ibfn) .GE. 2) THEN
           diag_mat(ibfn) = 0.0_dp
           DO jbfn = 1,nsgf
              U_mat(jbfn,ibfn) = 0.0_dp
           END DO
         ELSE 
           IF(ABS(diag_mat(ibfn)).LT.1E-10) THEN
             diag_mat(ibfn) = 0.0_dp
           END IF
           diag_mat(ibfn) = SQRT(diag_mat(ibfn))
           icount = icount + 1
           IF(iw>0) WRITE(iw," (A50, I3, F12.6, 2X, F12.6, I3, F12.6) ") &
                  "-moa- orb. -i- , max. ovlp., occ. number, occ.", &
                  icount, diag_mat(ibfn), diag_mat(ibfn)**2, occ, diag_mat(ibfn)**2*occ
         END IF
      END DO
      nbf = icount
     
      !IF(iw>0) CALL write_real_2dmatrix_to(U_mat, 124, "U_mat after zeroing (eig_val > 2) columns")
      IF(iw>0) WRITE(iw," (A94) ") "                                                                                      --------"
      norm_of_diag = dnrm2(SIZE(diag_mat), diag_mat, 1)
      IF(iw>0) WRITE(iw," (A82, F12.6) ") &
                  "-moa- population on fragment", norm_of_diag**2*occ

      ! population on fragment
      pop_i = 0.0_dp
      DO ibfn=1,nsgf
         pop_i = pop_i + diag_mat(ibfn)**2
      END DO
      pop_i = pop_i*occ

      IF(ispin .EQ. 1) THEN
        popfrg(ifrg) = pop_i
      ELSE
        chgfrg(ifrg) = pop_i
      END IF

      IF(ispin .EQ. 1) THEN
        spnfrg(ifrg) =  0.0_dp
        popfrg(ifrg) =  popfrg(ifrg)
        chgfrg(ifrg) = -popfrg(ifrg) + REAL(z_frg, kind=dp)
      ELSE
        spnfrg(ifrg) =  popfrg(ifrg) - chgfrg(ifrg)
        popfrg(ifrg) =  popfrg(ifrg) + chgfrg(ifrg)
        chgfrg(ifrg) = -popfrg(ifrg) + REAL(z_frg, kind=dp)
      END IF


      !--- -Trf- from Molec Ortho Frag Ortho Atom Ortho Sph. Harm. ---
      !            to             Frag Ortho Atom Ortho Sph. Harm.
      !            to                        Atom Ortho Sph. Harm.
      !            to                             Atom  Sph. Harm.
      ! U_mat right now is in Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis 
      ! transformations (T2*T3*T4)*U_mat to get in Sph. Harm. basis


      ALLOCATE(U_mat_truncated(nsgf, nbf))

      CALL cp_fm_struct_create(fmstruct=fmstruct, &
                         para_env=para_env, &
                         context=blacs_env, &
                         nrow_global=nsgf, &
                         ncol_global=nbf)
      CALL cp_fm_create(matrix=moa_coeff, &
                        matrix_struct=fmstruct, &
                        name="MOA coeff. of size nsgf by nbf")
      CALL cp_fm_struct_release(fmstruct=fmstruct)
      
      ! T3*T4
      CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T3, fm_T4, 0.0_dp, fm_rt1)
      ! T2*(T3*T4)
      CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T2, fm_rt1, 0.0_dp, fm_rt2)
      ! Taking U_mat into temp matrix
      CALL cp_fm_set_submatrix(fm_rt1, U_mat, 1, 1, nsgf, nsgf)
      ! (T2*T3*T4)*U
      CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_rt2, fm_rt1, 0.0_dp, fm_U_sph_harm_bas)
      

      IF(iw>0) WRITE(iw,*) "nbf value is ", nbf
      CALL cp_fm_get_submatrix(fm_U_sph_harm_bas, U_mat_truncated, 1, 1, nsgf, nbf)
      !IF(iw>0) CALL write_real_2dmatrix_to(U_mat_truncated, iw, "U_mat_truncated") 
      ! Now U_mat is in Sph. Harm. basis

      CALL cp_fm_set_submatrix(moa_coeff, U_mat_truncated, 1, 1, nsgf, nbf) 
      ! Printing to cubes
      logger => cp_get_default_logger()
      IF (cp_print_key_should_output(logger%iter_info, qs_env%input, &
                                     "PROPERTIES%MOA%MO_CUBES") /= 0) THEN
         CALL section_vals_val_get(mo_print_sec, "MO_CUBES%MO_LIST", &
                                               i_vals=mo_list)
      END IF
      DO imoa=1,SIZE(mo_list)
         CALL write_moa_cube(qs_env, logger, mo_print_sec, moa_coeff, mo_list(imoa), ispin, ifrg)
      END DO

      DEALLOCATE(PP, diag_mat, U_mat, U_mat_truncated)

      CALL cp_fm_release(fm_P)
      CALL cp_fm_release(fm_PP)
      CALL cp_fm_release(fm_tmp)
      CALL cp_fm_release(fm_rt1) !tmp array for reverse transforms
      CALL cp_fm_release(fm_rt2)
      CALL cp_fm_release(moa_coeff)
      CALL cp_fm_release(fm_U_sph_harm_bas)

      CALL timestop(handle)
   END SUBROUTINE orbmoa_frg


! **************************************************************************************************
!> \brief take mo_coeff and save given electronic state to cube files
!         repurposed save_mo_cube() subroutine in et_coupling_proj.F
!> \param qs_env ... QuickStep environment containing all system data
!> \param logger ... output logger
!> \param input ...  block of print section in moa input
!> \param mo_coeff ... electronic states data
!> \param im ... i^th vector in the moa orbitals
!> \param is spin ID
! **************************************************************************************************
   SUBROUTINE write_moa_cube(qs_env, logger, input, mo_coeff, im, is, ifrg)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      INTEGER                                            :: im, is, ifrg

      CHARACTER(len=*), PARAMETER :: routineN = 'save_mo_cube', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename
      CHARACTER(LEN=default_string_length)               :: title
      INTEGER                                            :: unit_nr
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: wf_g, wf_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys

      ! Initialization
      NULLIFY (particles)
      NULLIFY (subsys)

      NULLIFY (pw_env)
      NULLIFY (pw_pools)
      NULLIFY (auxbas_pw_pool)

      NULLIFY (atomic_kind_set)
      NULLIFY (cell)
      NULLIFY (dft_control)
      NULLIFY (particle_set)
      NULLIFY (qs_kind_set)

      ! Name of the cube file
      WRITE (filename, '(A4,I5.5,A6,I1.1,A5,I5.5)') 'MOA_', im, '_spin_', is, "ifrg_", ifrg
      ! Open the file
      unit_nr = cp_print_key_unit_nr(logger, input, 'MO_CUBES', extension='.cube', &
                                     middle_name=TRIM(filename), file_position='REWIND', log_filename=.FALSE.)
      ! Title of the file
      WRITE (title, *) 'MOA WAVEFUNCTION ', im, ' spin ', is

      ! List of all atoms
      CALL get_qs_env(qs_env, subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)

      ! Grids for wavefunction
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_g%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      ! Calculate the grid values
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, dft_control=dft_control, particle_set=particle_set)
      CALL calculate_wavefunction(mo_coeff, im, wf_r, wf_g, atomic_kind_set, &
                                  qs_kind_set, cell, dft_control, particle_set, pw_env)
      CALL cp_pw_to_cube(wf_r%pw, unit_nr, title, particles=particles, &
                         stride=section_get_ivals(input, 'MO_CUBES%STRIDE'))

      ! Close file
      CALL cp_print_key_finished_output(unit_nr, logger, input, 'MO_CUBES')

      ! Clean memory
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g%pw)

   END SUBROUTINE write_moa_cube

END MODULE moa
