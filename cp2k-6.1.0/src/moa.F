! **************************************************************************************************
!> \brief Maximal Orbital Analysis
!> \par History
!>      - module to calculate MOA occupations based on
!         Dupuis, et al., Journal of Computational Chemistry, 2019
!
!> \author  Pavan Behara
! **************************************************************************************************


MODULE moa
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
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
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_deallocate_matrix, dbcsr_get_block_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_set, dbcsr_setname, dbcsr_type
   USE mathlib,                         ONLY: diag
   USE input_section_types,             ONLY: section_vals_val_set,&
                                              section_vals_get,&
                                              section_vals_val_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_duplicate
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: diag
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: nso
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              mo_set_p_type, &
                                              get_mo_set
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
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
!> \brief Perform Maximal Orbital Analysis
!>        Check population analyses routine for code guidance
!> \param qs_env ...
!> \param output_unit ...
!> \param print_level ...
! **************************************************************************************************
   SUBROUTINE maximal_orbital_analysis(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'maximal_orbital_analysis', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: handle, iatom, ibfn, ii, ispin, iw, icount, ikind, &
                                                            natom, nfrg, nkind, nmo, &
                                                            nsgf, nspin, z_val, ind_count, ifrg, lfrg, z_frg
      INTEGER, DIMENSION(:), POINTER                     :: len_frg, frg_indices, nomoa, &
                                                            nomoaa, nomoab, ind_frg
      LOGICAL                                            :: ifcub, ifuno, ifpun
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

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set)
      NULLIFY (qs_kind_set)
      NULLIFY (fmstruct)
      NULLIFY (fm_S)
      NULLIFY (fm_T)
      NULLIFY (fm_TI)
      NULLIFY (fm_tmp)
      NULLIFY (matrixkp_s)
      NULLIFY (mos)
      NULLIFY (particle_set)
      NULLIFY (scf_control)
      NULLIFY (sm_s)
      NULLIFY (para_env)
      NULLIFY (blacs_env)
      

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
         nspin = SIZE(mos)
         nkind = SIZE(atomic_kind_set)
         

         ! Input flags to MOA section
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%NFRG", i_val=nfrg)
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%LFRG", i_vals=len_frg)
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%IFRG", i_vals=frg_indices) 
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%NOMOA", i_vals=nomoa) 
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%NOMOAA", i_vals=nomoaa)
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%NOMOAB", i_vals=nomoab)
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%IFUNO", l_val=ifuno)
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%IFPUN", l_val=ifpun)
         CALL section_vals_val_get(qs_env%input,"PROPERTIES%MOA%IFCUB", l_val=ifcub)

         ! Printing out the input tags for information 
         IF(iw>0) WRITE(iw, *) "INPUT TAGS:"
         IF(iw>0) WRITE(iw,*) "NFRG = ", nfrg
         IF(iw>0) WRITE(iw,*) "LFRG = ", len_frg(:)
         IF(iw>0) WRITE(iw,*) "IFRG = ", frg_indices(:)
         IF(iw>0) WRITE(iw,*) "NOMOA = ", nomoa(:)
         IF(iw>0) WRITE(iw,*) "NOMOAA = ", nomoaa(:)
         IF(iw>0) WRITE(iw,*) "NOMOAB = ", nomoab(:)
         IF(iw>0) WRITE(iw,*) "IFUNO = ", ifuno
         IF(iw>0) WRITE(iw,*) "IFPUN = ", ifpun
         IF(iw>0) WRITE(iw,*) "IFCUB = ", ifcub
         IF(iw>0) WRITE(iw,*) "-------------------------"


         ! Get the total number of contracted spherical Gaussian basis functions
         CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

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
            !CALL write_real_2dmatrix_to(tmp_block, 124, "transformation matrix -T-")
            CALL cp_fm_set_submatrix(fm_T, tmp_block, icount+1, &
                                 icount+1, ibfn, ibfn)

            DO ii=1,ibfn
               tmp_block(:,ii) = v_atomblock(:,ii)*sqrt(diag_atomblock(ii))
            END DO
            tmp_block = MATMUL(tmp_block, transpose(v_atomblock))
            !CALL write_real_2dmatrix_to(tmp_block, 124, "Inverse transformation matrix -TI-")
            CALL cp_fm_set_submatrix(fm_TI, tmp_block, icount+1, &
                                 icount+1, ibfn, ibfn)

            DEALLOCATE(real_atomblock, diag_atomblock, v_atomblock, tmp_block)
            icount = icount + ibfn
         END DO
         !CALL cp_fm_write_formatted(fm_S, iw, "S in Sph. Harm. basis")
         !CALL cp_fm_write_formatted(fm_T, iw, "fm_T")
         !CALL cp_fm_write_formatted(fm_TI, iw, "fm_TI")

         ! Doing T*S*T to get S in  Orthonormalized Sph. Harmonic Basis
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T, fm_S, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_T, 0.0_dp, fm_S)
         !CALL cp_fm_write_formatted(fm_S, iw, "S in Block-Orthonormal Sph. Harm. basis")
         
         ALLOCATE(tmp_block(nsgf, nsgf))
         ALLOCATE(diag_mat(nsgf))
         ALLOCATE(v_mat(nsgf, nsgf))
         CALL cp_fm_get_submatrix(fm_S, tmp_block, 1, 1, nsgf, nsgf)
         CALL diag(nsgf, tmp_block, diag_mat, v_mat)
         CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of block-orthonormal S")
         CALL write_real_2dmatrix_to(v_mat, 124, "eigen vectors of block-orthonormal S")

         ! Transformation matrix O is (EV)*(1/sqrt(diag))*(EV(t))
         ! Its inverse matrix OI is (EV)*(sqrt(diag))*(EV(t))
         ! where EV is the eigen vectors matrix v_atomblock

         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)/sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         CALL write_real_2dmatrix_to(tmp_block, 124, "transformation matrix -O-")
         CALL cp_fm_set_submatrix(fm_O, tmp_block, 1, 1, nsgf, nsgf)

         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)*sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         CALL write_real_2dmatrix_to(tmp_block, 124, "Inverse transformation matrix -OI-")
         CALL cp_fm_set_submatrix(fm_OI, tmp_block, 1, 1, nsgf, nsgf)
         !CALL cp_fm_write_formatted(fm_O, iw, "fm_O")
         !CALL cp_fm_write_formatted(fm_OI, iw, "fm_OI")
         ! Doing T*S*T to get S in  Orthonormalized Sph. Harmonic Basis
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_O, fm_S, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_O, 0.0_dp, fm_S)
         !CALL cp_fm_write_formatted(fm_S, iw, "S in Orthonormalized Block-Orthonormal Sph. Harm. basis")
         
         ! Ready to start MOA analysis, we have all the transformations we need
         ALLOCATE(popmoa(natom))
         ALLOCATE(chgmoa(natom))
         ALLOCATE(spnmoa(natom))
         popmoa(:) = 0.0_dp
         chgmoa(:) = 0.0_dp
         spnmoa(:) = 0.0_dp
         
         ! Calling orbmoa routine to do MOA analysis
         ! Inputs are S, V, T, SI, TI, OI
         ! Outputs U, PP
         ! for nspin=2 popmoa, chgmoa, spnmoa will be reused to get proper chg,
         ! spn, pop by combining alpha and beta contributions

         
         DO ispin=1,nspin
            CALL get_mo_set(mos(ispin)%mo_set, nmo=nmo) 
            CALL orbmoa(fm_S, mos(ispin)%mo_set%mo_coeff, fm_T, fm_TI, fm_O, fm_OI, &
                        popmoa, chgmoa, spnmoa, particle_set, qs_kind_set, nmo, nsgf,&
                        nspin, ispin, iw, para_env, blacs_env)
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
                  CALL moax_frg(qs_env, nfrg, ifrg, lfrg, ind_frg, &
                         ifpun, ifcub, ifuno, nomoa, nomoaa, nomoab, &
                         popfrg, chgfrg, spnfrg, z_frg)
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

      CALL timestop(handle)
   END SUBROUTINE maximal_orbital_analysis

   SUBROUTINE orbmoa(fm_S, fm_V, fm_T, fm_TI, fm_O, fm_OI, &
                     popmoa, chgmoa, spnmoa, particle_set, &
                     qs_kind_set, norb, nsgf, nspin, ispin, iw, &
                     para_env, blacs_env)
      TYPE(cp_fm_type), POINTER                           :: fm_S, fm_V, fm_T, fm_TI, fm_O, fm_OI
      REAL(KIND=dp), DIMENSION(:), POINTER, INTENT(INOUT) :: popmoa, chgmoa, spnmoa
      TYPE(particle_type), DIMENSION(:), POINTER          :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER           :: qs_kind_set
      INTEGER, INTENT(IN)                                 :: norb, nsgf, nspin, ispin, iw
      TYPE(cp_para_env_type), POINTER                     :: para_env
      TYPE(cp_blacs_env_type), POINTER                    :: blacs_env
      CHARACTER(LEN=*), PARAMETER :: routineN = 'orbmoa', &
         routineP = moduleN//':'//routineN
      INTEGER                                             :: err_info, handle, iatom, ibas, &
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

      !CALL cp_fm_write_formatted(fm_V, iw, "V in Sph. Harm. basis")

      ! Multiplying orbitals with TI on the left  TI*V to get V in  Block-Orthonormal Sph. Harmonic Basis
      ! placing in a temp matrix since output should not overlap with the inputs
      ! in gemm routine (matrix_c should not overlap matrix_a, matrix_b where c = a*b)
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_TI, fm_V, 0.0_dp, fm_tmp)
      CALL cp_fm_to_fm(fm_tmp, fm_V2)
      !CALL cp_fm_write_formatted(fm_V2, iw, "Orbitals -V- in Block-Orthonormal Sph. Harm. basis")

      ! Multiplying orbitals with OI on the left  OI*V to get V in
      ! Orthonormalized Block-Orthonormal Sph. Harmonic Basis
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_OI, fm_V2, 0.0_dp, fm_tmp)
      CALL cp_fm_to_fm(fm_tmp, fm_V2)
      !CALL cp_fm_write_formatted(fm_V2, iw, "Orbitals -V- in Orthonormalized Block-Orthonorm Sph. Harm. basis")

      ! P = S*V
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_S, fm_V2, 0.0_dp, fm_P)
      !CALL cp_fm_write_formatted(fm_P, iw, "Matrix P = S*V")

      ! PP = P*P(t)
      CALL cp_gemm('N', 'T', nsgf, nsgf, norb, 1.0_dp, fm_P, fm_P, 0.0_dp, fm_PP)
      !CALL cp_fm_write_formatted(fm_PP, iw, "Matrix PP = P*P(t)")
      
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
         CALL write_real_2dmatrix_to(pp_atomblock, 124, "P*P(t) Atom block")

         ! ----- Get Maximal Orbitals by diagonalization of - P*P(t) - -----
         ALLOCATE(U_mat(SIZE(pp_atomblock,1),SIZE(pp_atomblock,1)))
         ALLOCATE(diag_mat(MIN(SIZE(pp_atomblock,1),SIZE(pp_atomblock,2))))
   
         ! Doing diagonalization to get maximal overlap, MOA orbitals are in U_mat in block
         ! orthonormal Sph. Harm Basis on iatom
         CALL diag(ibfn, pp_atomblock, diag_mat, U_mat)
         CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of atomblock- P*P(t)")
         CALL write_real_2dmatrix_to(U_mat, 124, "U_mat of SVD of atomblock- P*P(t)")

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
      
      CALL timestop(handle)

   END SUBROUTINE orbmoa
!
!
   SUBROUTINE moax_frg(qs_env, nfrg, ifrg, lfrg, frg_indices, ifpun, &
                       ifcub, ifuno, nomoa, nomoaa, nomoab, popfrg, chgfrg, spnfrg, z_frg)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: nfrg, ifrg, lfrg
      INTEGER, DIMENSION(:), POINTER                     :: frg_indices
      LOGICAL                                            :: ifcub, ifuno, ifpun
      INTEGER, DIMENSION(:), POINTER                     :: nomoa, nomoaa, nomoab 
      REAL(KIND=dp), DIMENSION(:), POINTER               :: popfrg, chgfrg, spnfrg
      INTEGER, INTENT(OUT)                               :: z_frg
      CHARACTER(LEN=*), PARAMETER :: routineN = 'moax_frg', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: handle, iatom, jatom, iw, natom, &
                                                            nkind, nspin, nsgf, ind_count, &
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
            CALL write_real_2dmatrix_to(real_atomblock, 124, "Atom block")
            ! Jacobi diagonalization of matrix real_atomblock
            ! Syntax diag(n, a, d, v): Diagonalize matrix a. The eigenvalues are returned in vector d
            ! and the eigenvectors are returned in matrix v.
            CALL diag(ibfn, real_atomblock, diag_atomblock, v_atomblock)
            CALL write_real_1dmatrix_to(diag_atomblock, 124, "diag_atomblock")
            CALL write_real_2dmatrix_to(v_atomblock, 124, "v_ block")


            ! Transformation matrix T is (EV)*(1/sqrt(diag))*(EV(t))
            ! Its inverse matrix TI is (EV)*(sqrt(diag))*(EV(t))
            ! where EV is the eigen vectors matrix v_atomblock

            DO ii=1,ibfn
               tmp_block(:,ii) = v_atomblock(:,ii)/sqrt(diag_atomblock(ii))
            END DO
            tmp_block = MATMUL(tmp_block, transpose(v_atomblock))
            !CALL write_real_2dmatrix_to(tmp_block, 124, "transformation matrix -T-")
            CALL cp_fm_set_submatrix(fm_T2, tmp_block, icount+1, &
                                 icount+1, ibfn, ibfn)

            DO ii=1,ibfn
               tmp_block(:,ii) = v_atomblock(:,ii)*sqrt(diag_atomblock(ii))
            END DO
            tmp_block = MATMUL(tmp_block, transpose(v_atomblock))
            !CALL write_real_2dmatrix_to(tmp_block, 124, "Inverse transformation matrix -TI-")
            CALL cp_fm_set_submatrix(fm_TI2, tmp_block, icount+1, &
                                 icount+1, ibfn, ibfn)

            DEALLOCATE(real_atomblock, diag_atomblock, v_atomblock, tmp_block)
            icount = icount + ibfn
         END DO
         !CALL cp_fm_write_formatted(fm_S, iw, "S in Atom Sph. Harm. basis")
         !CALL cp_fm_write_formatted(fm_T2, iw, "fm_T -T2- to Atom Ortho Sph. Harm.")
         !CALL cp_fm_write_formatted(fm_TI2, iw, "fm_TI -TI2- to Atom Ortho Sph. Harm.")

         ! Doing T*S*T to get S in  Atom Ortho Sph. Harmonic Basis
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T2, fm_S, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_T2, 0.0_dp, fm_S)
         !CALL cp_fm_write_formatted(fm_S, iw, "S in Atom Ortho Sph. Harm. basis")

         ALLOCATE(tmp_block(nsgf, nsgf))
         ALLOCATE(diag_mat(nsgf))
         ALLOCATE(v_mat(nsgf, nsgf))
         CALL cp_fm_get_submatrix(fm_S, tmp_block, 1, 1, nsgf, nsgf)
         CALL diag(nsgf, tmp_block, diag_mat, v_mat)
         CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of S in Atom Ortho Sph. Harm.")
         CALL write_real_2dmatrix_to(v_mat, 124, "eigen vectors of S in Atom Ortho Sph. Harm.")

         !!' -T- to Molec Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)/sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         CALL write_real_2dmatrix_to(tmp_block, 124, "-T- to Molec Ortho Atom Ortho Sph. Harm. basis")
         CALL cp_fm_set_submatrix(fm_T, tmp_block, 1, 1, nsgf, nsgf)

         !!' -TI- (inv.) to Molec Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)*sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         CALL write_real_2dmatrix_to(tmp_block, 124, "-TI- Molec Ortho Atom Ortho Sph. Harm. basis")
         CALL cp_fm_set_submatrix(fm_TI, tmp_block, 1, 1, nsgf, nsgf)


         ! Doing T*S*T for Frag. S block orthonromal to get Frag. S in  Orthonormalized Sph. Harmonic Basis
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T, fm_S, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_T, 0.0_dp, fm_S_working)
         !CALL cp_fm_write_formatted(fm_S_working, iw, "S in Molec Ortho Atom Ortho Sph. Harm. basis")

         !!' eigen vectors and values of -S- in Molec Ortho Atom Ortho Sph. Harm. basis'
         CALL cp_fm_get_submatrix(fm_S_working, tmp_block, 1, 1, nsgf, nsgf)
         CALL diag(nsgf, tmp_block, diag_mat, v_mat)
         CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of S in Molec Ortho Atom Ortho Sph. Harm. basis")
         CALL write_real_2dmatrix_to(v_mat, 124, "eigen vectors of S in Molec Ortho Atom Ortho Sph. Harm. basis")


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

         CALL write_real_2dmatrix_to(tmp_block, 124, "S in Frag. blocked Atom Ortho Sph. Harm. basis")
         CALL diag(nsgf, tmp_block, diag_mat, v_mat)
         CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of S in Frag. blocked Atom Ortho Sph. Harm. basis")
         CALL write_real_2dmatrix_to(v_mat, 124, "eigen vectors of S in Frag. blocked Atom Ortho Sph. Harm. basis")
         CALL cp_fm_set_submatrix(fm_S_frg, tmp_block, 1, 1, nsgf, nsgf)
      
         ! Transformation matrix T is (EV)*(1/sqrt(diag))*(EV(t))
         ! Its inverse matrix TI is (EV)*(sqrt(diag))*(EV(t))
         ! where EV is the eigen vectors matrix S Frag. block-orthonormal
       
         !!' -T- to Frag Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)/sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         CALL write_real_2dmatrix_to(tmp_block, 124, "-T3- to Frag Ortho Atom Ortho Sph. Harm. basis")
         CALL cp_fm_set_submatrix(fm_T3, tmp_block, 1, 1, nsgf, nsgf)

         !!' -TI- (inv.) to Frag Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)*sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         CALL write_real_2dmatrix_to(tmp_block, 124, "-TI3- Frag Ortho Atom Ortho Sph. Harm. basis")
         CALL cp_fm_set_submatrix(fm_TI3, tmp_block, 1, 1, nsgf, nsgf)


         ! Doing T*S*T for Frag. S block orthonromal to get Frag. S in  Orthonormalized Sph. Harmonic Basis
         ! Note: Here S is in Atom Ortho Sph. Harm. Basis
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T3, fm_S, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_T3, 0.0_dp, fm_S_frg)
         !CALL cp_fm_write_formatted(fm_S_frg, iw, "S in Frag Ortho Atom Ortho Sph. Harm. basis")

         !!' eigen vectors and values of -S- in Frag Ortho Atom Ortho Sph. Harm. basis'
         CALL cp_fm_get_submatrix(fm_S_frg, tmp_block, 1, 1, nsgf, nsgf)
         CALL diag(nsgf, tmp_block, diag_mat, v_mat)
         CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of S in Frag Ortho Atom Ortho Sph. Harm. basis")
         CALL write_real_2dmatrix_to(v_mat, 124, "eigen vectors of S in Frag Ortho Atom Ortho Sph. Harm. basis")
          

         !!' -T- to Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)/sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         CALL write_real_2dmatrix_to(tmp_block, 124, "-T4- to Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis")
         CALL cp_fm_set_submatrix(fm_T4, tmp_block, 1, 1, nsgf, nsgf)

         !!' -TI- (inv.) to Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis'
         DO ii=1,nsgf
            tmp_block(:,ii) = v_mat(:,ii)*sqrt(diag_mat(ii))
         END DO
         tmp_block = MATMUL(tmp_block, transpose(v_mat))
         CALL write_real_2dmatrix_to(tmp_block, 124, "-TI4- to Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis")
         CALL cp_fm_set_submatrix(fm_TI4, tmp_block, 1, 1, nsgf, nsgf)

         ! Doing T*S*T for Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis 
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_T4, fm_S_frg, 0.0_dp, fm_tmp)
         CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, fm_tmp, fm_T4, 0.0_dp, fm_S_frg)
         !CALL cp_fm_write_formatted(fm_S_frg, iw, "-S- in Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis")

         popfrg(ifrg) = 0.0_dp
         chgfrg(ifrg) = 0.0_dp
         spnfrg(ifrg) = 0.0_dp
         
         ! Calling orbmoa_frg routine to do MOA analysis
         ! Inputs are S, V, T, SI, TI, OI
         ! Outputs U, PP
         ! for nspin=2 popmoa, chgmoa, spnmoa will be reused to get proper chg,
         ! spn, pop by combining alpha and beta contributions


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
            CALL orbmoa_frg(fm_S_frg, mos(ispin)%mo_set%mo_coeff, fm_T2, fm_TI2, fm_T3, fm_TI3, fm_T4, fm_TI4, &
                        popfrg, chgfrg, spnfrg, bf_frg, z_frg1, particle_set, qs_kind_set, nmo, nsgf,&
                        nspin, ispin, ifrg, iw, nomoa, para_env, blacs_env)
         END DO
         
         DEALLOCATE(atom_in_frg, bf_frg, tmp_block, diag_mat, v_mat)
      END IF 
      CALL timestop(handle)

   END SUBROUTINE moax_frg

   SUBROUTINE orbmoa_frg(fm_S, fm_V, fm_T2, fm_TI2, fm_T3, fm_TI3, fm_T4, fm_TI4, &
                     popfrg, chgfrg, spnfrg, bf_frg, z_frg, particle_set, &
                     qs_kind_set, norb, nsgf, nspin, ispin, ifrg, iw, nomoa, &
                     para_env, blacs_env)
      TYPE(cp_fm_type), POINTER                           :: fm_S, fm_V, fm_T2, fm_TI2, fm_T3, fm_TI3, fm_T4, fm_TI4
      REAL(KIND=dp), DIMENSION(:), POINTER, INTENT(INOUT) :: popfrg, chgfrg, spnfrg
      LOGICAL, DIMENSION(:), POINTER                      :: bf_frg
      INTEGER                                             :: z_frg
      TYPE(particle_type), DIMENSION(:), POINTER          :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER           :: qs_kind_set
      INTEGER, INTENT(IN)                                 :: norb, nsgf, nspin, ispin, ifrg, iw
      INTEGER, DIMENSION(:), POINTER                      :: nomoa
      TYPE(cp_para_env_type), POINTER                     :: para_env
      TYPE(cp_blacs_env_type), POINTER                    :: blacs_env
      CHARACTER(LEN=*), PARAMETER :: routineN = 'orbmoa_frg', &
         routineP = moduleN//':'//routineN
      INTEGER                                             :: err_info, handle, iatom, ibfn, &
                                                             jbfn, iorb, jorb, icount, ikind, &
                                                             occ, z_atom
      TYPE(cp_fm_struct_type), POINTER                    :: fmstruct
      TYPE(cp_fm_type), POINTER                           :: fm_P, fm_PP, fm_V2, fm_tmp
      TYPE(gto_basis_set_type), POINTER                   :: basis_set
      REAL(KIND=dp)                                       :: dum, pop_i
      REAL(KIND=dp), DIMENSION(:, :), POINTER             :: PP, U_mat
      REAL(KIND=dp), DIMENSION(:), POINTER                :: diag_mat

      CALL timeset(routineN, handle)

      NULLIFY (fmstruct)
      NULLIFY (fm_P)
      NULLIFY (fm_PP)
      NULLIFY (fm_tmp)

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

      !CALL cp_fm_write_formatted(fm_V, iw, "V in Sph. Harm. basis")
      ! Multiplying orbitals with TI2 on the left  TI2*V to get V in  Block-Orthonormal Sph. Harmonic Basis
      ! placing in a temp matrix since output should not overlap with the inputs
      ! in gemm routine (matrix_c should not overlap matrix_a, matrix_b where c = a*b)
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_TI2, fm_V, 0.0_dp, fm_tmp)
      CALL cp_fm_to_fm(fm_tmp, fm_V2)
      !CALL cp_fm_write_formatted(fm_V2, iw, "Orbitals -V- in Atom Ortho Sph. Harm. basis")

      ! Multiplying orbitals with TI3 on the left  TI3*V to get V in
      ! Orthonormalized Block-Orthonormal Sph. Harmonic Basis
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_TI3, fm_V2, 0.0_dp, fm_tmp)
      CALL cp_fm_to_fm(fm_tmp, fm_V2)
      !CALL cp_fm_write_formatted(fm_V2, iw, "Orbitals -V- in Frag Ortho Atom Ortho Sph. Harm. basis")

      ! Multiplying orbitals with TI4 on the left  TI4*V to get V in
      ! Orthonormalized Block-Orthonormal Sph. Harmonic Basis
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_TI4, fm_V2, 0.0_dp, fm_tmp)
      CALL cp_fm_to_fm(fm_tmp, fm_V2)
      !CALL cp_fm_write_formatted(fm_V2, iw, "Orbitals -V- in Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis")

      !CALL cp_fm_write_formatted(fm_S, iw, "-S- in Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis")

      IF(iw>0) WRITE(iw,*) " ----- now ... the analysis -----"
      IF(iw>0) WRITE(iw,*) " -MOA- analysis for -fragment- =", ifrg 
      IF(iw>0) WRITE(iw,*) " -------------------------------"

      ! ----- calculate P ----
      ! P = S*V for all
      ! -S, V - are in Molec Ortho Frag Ortho Atom Ortho Sph. Harm. basis
      CALL cp_gemm('N', 'N', nsgf, norb, nsgf, 1.0_dp, fm_S, fm_V2, 0.0_dp, fm_P)
      !CALL cp_fm_write_formatted(fm_P, iw, "Matrix P = S*V for all")

      ! P_fragement
      DO ibfn = 1,nsgf
         DO jorb = 1,norb
            IF(.NOT.bf_frg(ibfn)) THEN
               CALL cp_fm_set_element(fm_P, ibfn, jorb, 0.0_dp)
            END IF
         END DO
      END DO       
      !CALL cp_fm_write_formatted(fm_P, iw, "P_fragment")

      ! P_frag. for non-excluded molec. orbitals
      DO jorb = 1,norb
         IF(nomoa(jorb).NE.0) THEN
           DO ibfn = 1,nsgf
               CALL cp_fm_set_element(fm_P, ibfn, jorb, 0.0_dp)
           END DO
         END IF
      END DO
      !CALL cp_fm_write_formatted(fm_P, iw, "P_fragment for non-excluded orbitals")
      
      ! PP = P*P(t)
      CALL cp_gemm('N', 'T', nsgf, nsgf, norb, 1.0_dp, fm_P, fm_P, 0.0_dp, fm_PP)
      ALLOCATE(PP(nsgf, nsgf))
      ! Prepare -P*P(t)- 
      CALL cp_fm_get_submatrix(fm_PP, PP, 1, 1, nsgf, nsgf)
      DO ibfn = 1,nsgf
         DO jbfn = 1,nsgf
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
      CALL write_real_2dmatrix_to(PP, 124, "P*P(t) for fragment and non-excluded Molec Orbs")
      CALL cp_fm_set_submatrix(fm_PP, PP, 1, 1, nsgf, nsgf)
      !CALL cp_fm_write_formatted(fm_PP, iw, "Matrix PP = P*P(t)")

      ! Get Maximal Orbitals -U-
      ALLOCATE(U_mat(SIZE(PP,1),SIZE(PP,1)))
      ALLOCATE(diag_mat(MIN(SIZE(PP,1),SIZE(PP,2))))

      ! Doing diagonalization to get maximal overlap, MOA orbitals are in U_mat in block
      ! orthonormal Sph. Harm Basis on iatom
      CALL diag(nsgf, PP, diag_mat, U_mat)
      CALL write_real_1dmatrix_to(diag_mat, 124, "eigen values of P*P(t) for fragment MOA orbs")
      CALL write_real_2dmatrix_to(U_mat, 124, "U_mat of SVD of P*P(t) for fragment MOA orbs")

      ! Get MOA fragment population 
      icount = 0
      DO ibfn = 1,nsgf
         IF(diag_mat(ibfn) .GE. 2) THEN
           diag_mat(ibfn) = 0.0_dp
           DO jbfn = 1,nsgf
              U_mat(ibfn,jbfn) = 0.0_dp
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
      DEALLOCATE(PP, diag_mat, U_mat)

      CALL timestop(handle)
   END SUBROUTINE orbmoa_frg
   
!
END MODULE moa
