# Identifying porous cages in the Cambridge Structural Database using topological data analysis

This repository contains the scripts and datasets used and presented in "Identifying porous cages subsets in the Cambridge Structural Database using topological data analysis" by Aurelia Li, Rocio Bueno-Perez and David Fairen-Jimenez, Adsorption and Advanced Materials Group (aam.ceb.cam.ac.uk), led by David Fairen-Jimenez from the Department of Chemical Engineering and Biotechnology, University of Cambridge.

## List of datasets
- The __Latest datasets__ folder contains:
    - **_MOC_Mar2022.gcd_**: list of MOCs identified in the CSD with updates up to March 2022.
    - **_OC_Mar2022.gcd_**: list of OCs identified in the CSD with updates up to March 2022.
- The __Previous datasets__ folder contains the first lists of structures obtained during the writing and peer-review process of the paper:
    - **_MOC_Nov2019.gcd_**: list of MOCs identified in the CSD with updates up to November 2019.
    - **_OC_Nov2019.gcd_**: list of OCs identified in the CSD with updates up to November 2019.
    - **_MOC_update_Nov2019-Mar2022.gcd_**: list of MOCs identified in the CSD between November 2019 and March 2022.
    - **_OC_update_Nov2019-Mar2022.gcd_**: list of MOCs identified in the CSD between November 2019 and March 2022.
    - **_CARBOXYLATE_Mar2022.gcd_**: list of carboxylate-based structures in the CSD with updates up to March 2022. 

_Disclaimer_: these structures were programmatically found in the CSD, and were not all visually checked. Please let us know if any structure shouldn't be included and if other structures should be.

## Scripts:
- __TDA_structural_data_preparation.py__: given a list of refcodes (a GCD list), this script 
    - removes any guest molecules and solvents, 
    - performs quick checks to eliminate structures unlikely to be cages,
    - writes the resulting structure to a CIF, 
    - writes the corresponding coordinates to a CSV, 
    - returns the reduced list of potential cage candidates.
- __TDA_persistence_landscapes.py__: given the list of potential cages and the .csv coordinates obtained previously, this script calculates the persistence landscapes and outputs the result as an NPY (numpy) file.
- __TDA_classification.py__: given the NPY files obtained previously, this script classifies the different structures using hierarchical clustering. A random forest seciton is also included. 
- _Bonus_: TDA_noise_removal.py. This srcipt contains functions that were used to identify and discard noisy structures. 
