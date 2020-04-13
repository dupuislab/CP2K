&FORCE_EVAL
  METHOD QS
  &DFT  
     UKS
     BASIS_SET_FILE_NAME  ../../../cp2k/data/BASIS_MOLOPT
     POTENTIAL_FILE_NAME  /projects/academic/mdupuis2/pavan/POTENTIALS/GTH_POTENTIALS_5-12-10
     WFN_RESTART_FILE_NAME ./p_1-RESTART.wfn
     CHARGE -1
     MULTIPLICITY 2
     PLUS_U_METHOD MULLIKEN
    &MGRID
      CUTOFF 650
      NGRIDS   10
      REL_CUTOFF   60
    &END MGRID
    &QS
       EPS_DEFAULT 1.0E-12
       METHOD GPW
       MAP_CONSISTENT
       EXTRAPOLATION ASPC
       EPS_PGF_ORB 1E-12
    &END QS
     &SCF
       MAX_SCF    500
       EPS_SCF    ${EPSSCF}
       EPS_DIIS 0.5
       SCF_GUESS  RESTART
       &MIXING  T
            METHOD BROYDEN_MIXING
       &END MIXING
       &DIAGONALIZATION
               ALGORITHM STANDARD
       &END
       &OUTER_SCF
         MAX_SCF 500
         EPS_SCF ${EPSSCF}
       &END
       &PRINT
         #&PROGRAM_RUN_INFO
         #  FILENAME ./p_state1.out
         #&END
         &RESTART
           FILENAME ./p_state1
           &EACH
               QS_SCF 100
           &END
           ADD_LAST NUMERIC
         &END
       &END

     &END SCF

     &XC
        &XC_FUNCTIONAL PBE
        &END XC_FUNCTIONAL
     &END XC
   &PRINT
    &MULLIKEN ON
    &END
    &LOWDIN OFF
    &END
    &XRAY_DIFFRACTION_SPECTRUM OFF
    &END
    &PDOS OFF
      FILENAME PDOS-TMAA
      NLUMO -1
    &END PDOS
    &E_DENSITY_CUBE OFF
    &END
    &MO_CUBES                 ! Controls which MOs are written to cube-files.
      WRITE_CUBE FALSE
    &END MO_CUBES
    &TOT_DENSITY_CUBE OFF
    &END TOT_DENSITY_CUBE
   &END PRINT
  &END DFT

  &SUBSYS
    &CELL
   A   ${CELL_A}        0.0000000000        0.0000000000
   B   0.0000000000       ${CELL_B}      0.0000000000
   C   0.0000000000        0.0000000000      ${CELL_C}
    &END CELL
    &TOPOLOGY
      CONNECTIVITY OFF
      COORD_FILE_NAME ${XYZ}
      COORDINATE      xyz
    &END TOPOLOGY
   
    &KIND Ti
      ELEMENT Ti
      BASIS_SET SZV-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q12
      &DFT_PLUS_U
         L 2
         U_MINUS_J 0.367
      &END DFT_PLUS_U
    &END KIND
    &KIND O
      BASIS_SET SZV-MOLOPT-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
    
  &END SUBSYS
&END FORCE_EVAL
