This repository contains files associated with our work entitled "Accurate Modeling of Acid Dissociation: Collective Variables, Hybrid Exchange-Correlation Functionals, and Long-Range Machine Learning Potentials".

## Repository Structure
```
pKaCalculations/
├── 0-data/
│   ├── training-data.xyz          # MLIP training dataset
│   ├── test-data.xyz              # MLIP test dataset
│   ├── vasp2mace.py               # Script to extract data from VASP output
│   └── DFT/
│       ├── INCAR		           # Example DFT calculation files
│       ├── POSCAR
│       └── OUTCAR                 
│
├── 1-mliap_training/
│   ├── a-mace/                    # Representative files for training using MACE
│   │   ├── config.yml
│   │   ├── run_train.py
│   │   └── MACE_models/
│   └── b-maceles/                 # Representative files for training using MACELES
│       ├── config.yml
│       ├── run_train.py
│       └── MACE_models/
│
├── 2-umbrella_sampling/
│   ├── a-water/		            # Representative files for ASE-MACELES simulation
│   ├── b-formic_acid/
│   ├── c-acetic_acid/
│   ├── d-glycine_O/
│   ├── e-glycine_N/
│   ├── f-aspartic_O1/
│   ├── g-aspartic_O2/
│   ├── h-aspartic_N/
│   ├── i-sample_lammps_mace/      # Representative files for LAMMPS-MACE simulation
│   └── umbrella_integration.sh    # Script to calculate free energy from umbrella sampling simulations
│
├── 3-collective_variables/	
│   ├── pKa.cpp			            # Custom PLUMED source files, implementing CV for acid dissociations cases
│   ├── pKaBase.cpp
│   └── pKaBase.h                  
│
└── README.md
```
