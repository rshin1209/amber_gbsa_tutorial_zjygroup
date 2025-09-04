# GBSA Preparation and Submission Tool 
## ZJYGroup At Vanderbilt University

This repository provides a Python script (`gbsa_prep.py`) to **prepare and run MM-GBSA or QM/MM-GBSA calculations** with [AmberTools](https://ambermd.org/AmberTools.php).  
It automates:

- Stripping **complex/receptor/ligand** prmtops with `cpptraj`
- Preparing stripped/autoimaged trajectory (`md.nc`)
- Generating **MMPBSA input** files (MM or QM/MM)
- Creating a **SLURM job script**
- Optionally submitting the job with `sbatch`

---

## Features

- Reads a JSON configuration file with all parameters
- Supports both **MM** and **QM/MM** GBSA workflows
- Workspace created automatically as `<prefix>_gbsa`
- Safety checks for required executables (`cpptraj`, `mpirun`, `sbatch`)
- Includes a `--dry-run` mode to prepare inputs without submission

---

## Requirements

- Python ≥ 3.6  
- Amber/AmberTools in your `PATH`  
  - `cpptraj`  
  - `MMPBSA.py.MPI`
- SLURM (`sbatch`) environment for job submission

Make sure Amber is sourced before running:

```bash
source /home/shaoq1/bin/amber_env/amber-accre.sh
```

# File Directory
## Before running:
```bash
./
├── gbsa_prep.py
├── input.json
└── <prefix>/
    ├── md.nc
    └── solvated_complex.prmtop
```

## After running:
```bash
./
├── gbsa_prep.py
├── input.json
├── <prefix>/
│   ├── md.nc
│   └── solvated_protein.prmtop
└── <prefix>_gbsa/
    ├── MMgbsa/
    │   ├── submit.job
    │   └── MMgbsa.in
    ├── complex.prmtop
    ├── complex_prmtop.in
    ├── ligand.prmtop
    ├── ligand_prmtop.in
    ├── md.nc
    ├── receptor.prmtop
    ├── receptor_prmtop.in
    └── strip_traj.in
```

## 1. Prepare JSON input
Example input.json:
```bash
{
  "directory": "//nobackup/yang_lab/shinw3/lasso/FC_wt",
  "complex_residues": "1-723",
  "receptor_residues": "18-723",
  "ligand_residues": "1-17",
  "level_of_theory": "MM",
  "startframe": 5000,
  "endframe": 50000,
  "interval": 100,
  "igb": 2,
  "saltcon": 0.150,
  "qm_residues": "1-17,108-120",
  "qmcharge_com": 0,
  "qmcharge_rec": -1,
  "qmcharge_lig": 1,
  "submit_job": "True"
}
```
## 2. Run the script
```bash
source /home/shaoq1/bin/amber_env/amber-accre.sh
python gbsa_prep.py -i input.json
```
## 3. Job Submission Options

- If `"submit_job": "True"` in your JSON file, the script will **automatically submit the job** with `sbatch`.  
- If `"submit_job": "False"`, the script will only prepare the input files. You will then need to submit manually:

```bash
cd ./<prefix>_gbsa/MMgbsa/
sbatch submit.job
