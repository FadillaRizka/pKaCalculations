# pKaCalculations
This repository contains files associated with our work entitled "Accurate Modeling of Acid Dissociation: Collective Variables, Hybrid Exchange-Correlation Functionals, and Long-Range Machine Learning Potentials".  

## Repository Structure

```text
pKaCalculations/
├── 0-data/
│   ├── training-data.xyz          # MLIP training dataset
│   ├── test-data.xyz              # MLIP test dataset
│   ├── vasp2mace.py               # Script for extracting training data from VASP outputs
│   └── DFT/                       # Representative DFT calculation files
│       ├── INCAR
│       ├── POSCAR
│       └── OUTCAR
│
├── 1-mliap_training/
│   ├── a-mace/                    # Representative files for MACE training
│   │   ├── config.yml
│   │   ├── run_train.py
│   │   └── MACE_models/
│   └── b-maceles/                 # Representative files for MACE-LES training
│       ├── config.yml
│       ├── run_train.py
│       └── MACE_models/
│
├── 2-umbrella_sampling/
│   ├── a-water/                   # Representative files for ASE–MACE-LES simulations
│   ├── b-formic_acid/
│   ├── c-acetic_acid/
│   ├── d-glycine_O/
│   ├── e-glycine_N/
│   ├── f-aspartic_O1/
│   ├── g-aspartic_O2/
│   ├── h-aspartic_N/
│   ├── i-sample_lammps_mace/      # Representative files for LAMMPS–MACE simulations
│   └── umbrella_integration.sh    # Script for umbrella integration free-energy calculations
│
├── 3-collective_variables/
│   ├── pKa.cpp                    # Custom PLUMED source code for acid dissociation CVs
│   ├── pKaBase.cpp
│   └── pKaBase.h
│
└── README.md
```
