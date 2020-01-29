&FORCE_EVAL
  METHOD QS
  STRESS_TENSOR NONE
  &DFT  
     UKS
     BASIS_SET_FILE_NAME  /projects/academic/mdupuis2/pavan/BASIS_SETS/GTH_BASIS_SETS_5-12-10
     POTENTIAL_FILE_NAME  /projects/academic/mdupuis2/pavan/POTENTIALS/GTH_POTENTIALS_5-12-10
     WFN_RESTART_FILE_NAME ./p_0-RESTART.wfn
     CHARGE 6
     MULTIPLICITY 35
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
       EPS_PGF_ORB 1E-10
    &END QS
     &SCF
       MAX_SCF    1000 
       EPS_SCF    ${EPSSCF}
       SCF_GUESS  RESTART
       &OT
       &END
       &OUTER_SCF
         MAX_SCF 500
         EPS_SCF ${EPSSCF}
       &END
       &PRINT
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
        &XC_FUNCTIONAL BLYP
        &END XC_FUNCTIONAL
     &END XC
   &PRINT
    &MULLIKEN ON
    &END
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
      &DFT_PLUS_U
         L 2
         U_MINUS_J 0.18375  #5 eV, conv factor = 27.2107
      &END DFT_PLUS_U
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
