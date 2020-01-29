&FORCE_EVAL
  METHOD QS
  &DFT
     UKS
     BASIS_SET_FILE_NAME  /projects/academic/mdupuis2/pavan/BASIS_SETS/GTH_BASIS_SETS_5-12-10
     POTENTIAL_FILE_NAME  /projects/academic/mdupuis2/pavan/POTENTIALS/GTH_POTENTIALS_5-12-10
     WFN_RESTART_FILE_NAME ./p_0-RESTART.wfn
     CHARGE 1
     MULTIPLICITY 2
     &POISSON
        PERIODIC NONE
        POISSON_SOLVER WAVELET
     &END
    &QS
       EPS_DEFAULT 1.0E-12
       METHOD GAPW
       MAP_CONSISTENT
       EXTRAPOLATION ASPC
    &END QS
     &SCF
      EPS_SCF 1.0E-6
      SCF_GUESS RESTART
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
           BACKUP_COPIES 0
           FILENAME ./p_state0
           &EACH
               JUST_ENERGY 1
           &END
           ADD_LAST NUMERIC
         &END
       &END

     &END SCF
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
            ABC ${CELL_A}  ${CELL_A}  ${CELL_A}
    &END CELL
    &TOPOLOGY
     COORD_FILE_FORMAT XYZ
     COORD_FILE_NAME ./state0.xyz
     &CENTER_COORDINATES ON
     &END CENTER_COORDINATES
    &END TOPOLOGY
    &KIND Zn_plus
      ELEMENT Zn 
      &BS
        &ALPHA
           L 2
           N 3
           NEL 0
        &END
        &BETA
           L 2
           N 3
           NEL -2
        &END
      &END
    BASIS_SET DZVP-MOLOPT-SR-GTH-q12
    POTENTIAL GTH-PBE-q12
      LEBEDEV_GRID 500
      RADIAL_GRID 2500
    &END KIND

  &KIND Zn
    BASIS_SET DZVP-MOLOPT-SR-GTH-q12
    POTENTIAL GTH-PBE-q12
      LEBEDEV_GRID 500
      RADIAL_GRID 2500
  &END KIND

  &END SUBSYS
&END FORCE_EVAL

