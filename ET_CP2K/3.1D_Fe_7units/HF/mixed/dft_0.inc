&FORCE_EVAL
  METHOD QS
  &DFT  
     UKS
     BASIS_SET_FILE_NAME  /projects/academic/mdupuis2/pavan/BASIS_SETS/GTH_BASIS_SETS_5-12-10
     POTENTIAL_FILE_NAME  /projects/academic/mdupuis2/pavan/POTENTIALS/GTH_POTENTIALS_5-12-10
     WFN_RESTART_FILE_NAME ./p_0-RESTART.wfn
     CHARGE 6
     MULTIPLICITY 35
    &MGRID
      CUTOFF 500
      NGRIDS   5
      REL_CUTOFF   60
    &END MGRID
    &QS
       EPS_DEFAULT 1.0E-12
       METHOD GPW
       MAP_CONSISTENT
       EXTRAPOLATION ASPC
       EPS_PGF_ORB 1E-10
    &END QS
     &SCF
       MAX_SCF    500
       EPS_SCF    ${EPSSCF}
       SCF_GUESS  RESTART
       &OT
       &END
       &OUTER_SCF
         MAX_SCF 500
         EPS_SCF ${EPSSCF}
       &END
       &PRINT
         #&PROGRAM_RUN_INFO
         #  FILENAME ./p_state0.out
         #&END
         &RESTART
           FILENAME ./p_state0
           &EACH
                JUST_ENERGY 1
           &END
           ADD_LAST NUMERIC
         &END
       &END
     &END SCF

     &XC
        &XC_FUNCTIONAL NONE
        &END XC_FUNCTIONAL
        &HF
        &INTERACTION_POTENTIAL
                POTENTIAL_TYPE TRUNCATED
                CUTOFF_RADIUS [angstrom]  10.6711
                T_C_G_DATA ./t_c_g.dat
        &END
        &MEMORY
            MAX_MEMORY 5000
        &END
        &END
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
   &KIND Fe
      ELEMENT Fe
      BASIS_SET mod-DZV
      POTENTIAL GTH-BLYP-q16
   &END KIND
   &KIND O
      BASIS_SET mod-DZV
      POTENTIAL GTH-BLYP-q6
   &END KIND

   &KIND H
      BASIS_SET mod-DZV
      POTENTIAL GTH-BLYP-q1
   &END KIND
  &END SUBSYS
&END FORCE_EVAL
