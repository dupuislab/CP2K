&FORCE_EVAL
  METHOD QS
  &DFT  
     UKS
     BASIS_SET_FILE_NAME  /projects/academic/mdupuis2/pavan/BASIS_SETS/GTH_BASIS_SETS_5-12-10
     POTENTIAL_FILE_NAME  /projects/academic/mdupuis2/pavan/POTENTIALS/GTH_POTENTIALS_5-12-10
     PLUS_U_METHOD MULLIKEN
     WFN_RESTART_FILE_NAME ./p_0-RESTART.wfn
     CHARGE -1
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
    &END QS
     &SCF
       MAX_SCF    100 
       EPS_SCF    ${EPSSCF}
       SCF_GUESS  RESTART
       &OT
       &END
       &OUTER_SCF
         MAX_SCF 50
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
     A       14.7248573303         0.0000000000         0.0000000000
     B       -7.3624286652        12.7521005152         0.0000000000
     C        0.0000000000         0.0000000000        13.2674207687
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
         U_MINUS_J 0.2205  #6 eV, conv factor = 27.2107
      &END DFT_PLUS_U
   &END KIND
   &KIND O
      BASIS_SET mod-DZV
      POTENTIAL GTH-BLYP-q6
   &END KIND
  &END SUBSYS
&END FORCE_EVAL
