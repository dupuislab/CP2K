# ET-CP2K 
### VAB Calculation

Code addition to a cloned verison of CP2K-6.1 to calculate electron coupling matrix element VAB for molecular and periodic systems. <sup>[1]</sup> Calculation method is based on Farazdel & Dupuis work. <sup>[2]<sup> 

Testcases are in the folder ET_CP2K. Input section needs to be specified as follows:<br />
Hartree Fock options are a replication of Guidon et al.'s HF section from CP2K manual https://manual.cp2k.org/cp2k-6_1-branch/CP2K_INPUT/FORCE_EVAL/DFT/XC/HF.html

>
>>     &FORCE_EVAL
        ...
        &MIXED
          &VAB
          DO_VAB .TRUE.
          &HF !
              FRACTION 1.0  !MUST BE 1.0 FOR VAB_CALCULATION
      	     !Other optional sections for setting up HF calculation as in
        	   !Guidon, et al.’s papers on HF implementation in CP2K  111,134,139
              HF_INFO
              INTERACTION_POTENTIAL
              LOAD_BALANCE
              MEMORY
              PERIODIC
              SCREENING 
            &END
          &END
        &END
        ...
    &END

# MOA-CP2K

Maximal Orbital Analysis (MOA) for better interpretation of electronic interactions based on Dupuis et al.'s, work <sup> [3] </sup> is available.

Testcases are in the folder MOA_CP2K. Input section is as follows:<br />

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

#### Installation instructions
To download the code
 
`git clone https://github.com/dupuislab/CP2K.git`

To build the code with the default Makefile and an arch file use 

`make -j8 ARCH=local VERSION=popt`

where arch file local.popt exists in the arch folder. 

For more detailed instructions on the dependencies and how to modify Makefiles to suit your compute nodes please visit the CP2K website

`https://www.cp2k.org/howto`

References:

1. PK Behara, M Dupuis, Electron transfer in extended systems: characterization by periodic density functional theory including the electronic coupling, Phys. Chem. Chem. Phys., 2020, Advance Article.(https://dx.doi.org/10.1039/C9CP05133C)
2.    Abbas Farazdel, Michel Dupuis, Enrico Clementi, and Ari Aviram Journal of the American Chemical Society 1990 112 (11), 4206-4214 DOI: 10.1021/ja00167a016
3.    Michel Dupuis, Meghana Nallapu, Journal of Computational Chemistry 2019, 40, 39–50

