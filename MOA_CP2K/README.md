# MOA-CP2K

Maximal Orbital Analysis (MOA) for better interpretation of electronic interactions based on Dupuis et al.'s, work is available.

Testcases for few molecular systems are in this folder. More test cases will be added soon. Input section is as follows:<br />

> 
>>     &PROPERTIES
            &MOA
                 NFRG 2             ! No. of fragments
                 LFRG 1 2           ! Length of fragments 
                 IFRG 1 2 3         ! Atoms in fragments by indices 
                 NOMOA 1 1 1 0 0    ! Exclude MOs if index is non-zero
                 NOMOAA 1 2 3 0 0   ! Exclude alpha MOs if index is non-zero
                 NOMOAB 1 1 1 0 0   ! Exclude beta MOs if index is non-zero
                 &MO_CUBES          ! build wfn from MOA orbitals and print individual cubes                             
                         STRIDE 2 2 2                  ! The stride (X,Y,Z) used to write the cube file
                         MO_LIST 1 2 3               ! List of MOs to print as cube
                 &END
         &END MOA
    &END PROPERTIES


