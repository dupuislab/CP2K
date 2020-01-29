&FORCE_EVAL
  METHOD QS
  &DFT
     UKS
     BASIS_SET_FILE_NAME  ./HE_onlyS
     POTENTIAL_FILE_NAME  ./POTENTIAL
     WFN_RESTART_FILE_NAME ./p_0-RESTART.wfn
     CHARGE 1
     MULTIPLICITY  2
    &QS
       EPS_DEFAULT 1.0E-12
       METHOD GAPW
       MAP_CONSISTENT
       EXTRAPOLATION ASPC
    &END QS
     &SCF
      EPS_SCF 1.0E-7
      SCF_GUESS ATOMIC
      MAX_SCF 60
      &DIAGONALIZATION
              ALGORITHM STANDARD
      &END
      &OUTER_SCF
        EPS_SCF 1.0E-6
        MAX_SCF 6
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
    &POISSON
      PERIODIC NONE
      PSOLVER  MT
    &END POISSON
    &XC                        ! Parametes needed to compute the electronic exchange potential
      &XC_FUNCTIONAL NONE      ! No xc functional
      &END XC_FUNCTIONAL
      &HF                      ! Hartree Fock exchange. In this case is 100% (no fraction specified).
        &INTERACTION_POTENTIAL
                POTENTIAL_TYPE COULOMB
        &END
        &SCREENING             ! Screening of the electronic repulsion up to the given threshold.              
          EPS_SCHWARZ 1E-12  ! Threshold specification
        &END SCREENING
      &END HF
    &END XC
   &PRINT
    &MULLIKEN ON
    &END
   &END PRINT
  &END DFT

  &SUBSYS
    &CELL
            PERIODIC NONE
            ABC ${CELL_A}  ${CELL_B}  ${CELL_C}
    &END CELL
    &TOPOLOGY
      CONNECTIVITY OFF
      COORD_FILE_NAME ./state0.xyz
      COORDINATE      xyz
      &CENTER_COORDINATES
      &END
    &END TOPOLOGY
      &KIND He_plus
      ELEMENT He
      &BS
        &ALPHA
           L 0
           N 1
           NEL 0
        &END
        &BETA
           L 0
           N 1
           NEL -2
        &END
      &END
      BASIS_SET DZVP-S-ONLY
      POTENTIAL ALL
      LEBEDEV_GRID 200
      RADIAL_GRID 100
   &END KIND

   &KIND He
      ELEMENT He
      BASIS_SET DZVP-S-ONLY
      POTENTIAL ALL
      LEBEDEV_GRID 200
      RADIAL_GRID 100
    &END KIND
  &END SUBSYS
&END FORCE_EVAL

