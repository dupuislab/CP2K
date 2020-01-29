### VAB Calculation

Refer to the paper for a detailed summary on the calculation method and results for molecular and periodic systems. <sup>[1]</sup>  

Testcases in this folder show how to setup and run VAB calculation. All the parameters are set to optimal as per our knowledge. If further fine tuning is needed for any system please do so and report back the improvements if you can.
The procedure followed is <br />
1. Initial and final charged states are obtained by elongating the bonds around the required site by around 0.15 A and optimizing the geometry within a specific level of theory (HF or DFT or others).<br />
2. Geometries are interpolated between these two states and electronic coupling is calcuated at the halfway point.<br />
3. Wherever possible wavefunctions of already optimized initial and final states are used as initial guess.<br />

 Input section needs to be specified as follows:<br />
Hartree Fock options are a replication of Guidon et al.'s HF section from CP2K manual https://manual.cp2k.org/cp2k-6_1-branch/CP2K_INPUT/FORCE_EVAL/DFT/XC/HF.html

>
>>     &FORCE_EVAL
        ...
        &MIXED
          &VAB
          DO_VAB .TRUE.
          &HF !
              FRACTION 1.0  !Must be 1.0
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

References:

1. PK Behara, M Dupuis, Electron transfer in extended systems: characterization by periodic density functional theory including the electronic coupling, Phys. Chem. Chem. Phys., 2020, Advance Article.(https://dx.doi.org/10.1039/C9CP05133C)
