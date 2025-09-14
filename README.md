# GBSA Preparation and Submission Tool 
## ZJYang Lab At Vanderbilt University

This repository provides a Python script (`gbsa_prep.py`) to **prepare and run MM-GBSA or QM/MM-GBSA calculations** with [AmberTools](https://ambermd.org/AmberTools.php) on the **ACCRE cluster at Vanderbilt University**.
It automates:

- Stripping **complex/receptor/ligand** prmtops with `cpptraj`
- Preparing stripped/autoimaged trajectory (`md.nc`)
- Generating **GBSA input** files (MM or semiempirical-QM/MM)
- Creating a **SLURM job script**
- Optionally submitting the job with `sbatch`

---

## Features

- Reads a JSON configuration file with all parameters
- Supports both **MM** and **semiempirical-QM/MM** GBSA workflows
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
module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 ambertools/25.0
```

## File Directory
### Before running:
```bash
./
├── gbsa_prep.py
├── input.json
└── <prefix>/
    ├── md.nc
    └── solvated_complex.prmtop
```

### After running:
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
  "directory": "./FC_wt",
  "directory_description": "File system path to the working directory where .nc and .prmtop files for this simulation are stored",

  "complex_residues": "1-723",
  "complex_residues_description": "Residue range belonging to complex",

  "receptor_residues": "18-723",
  "receptor_residues_description": "Residue range belonging to receptor",

  "ligand_residues": "1-17",
  "ligand_residues_description": "Residue range belonging to ligand",

  "level_of_theory": "MM",
  "level_of_theory_description": "Method used for calculation: MM or semiempirical (AM1, PM6, DFTB3, etc.)",

  "startframe": 5000,
  "startframe_description": "Starting frame index from trajectory",

  "endframe": 50000,
  "endframe_description": "Last frame index from trajectory",

  "interval": 100,
  "interval_description": "Frame sampling interval",

  "igb": 2,
  "igb_description": "Generalized Born solvent model parameter (igb = 2 recommended)",

  "saltcon": 0.150,
  "saltcon_description": "Salt concentration (mM)",

  "qm_residues": "1-17,108-120",
  "qm_residues_description": "Residues included in QM region",

  "qmcharge_com": 0,
  "qmcharge_com_description": "Total charge of the complex. qmcharge_com = qmcharge_rec + qmcharge_lig",

  "qmcharge_rec": -1,
  "qmcharge_rec_description": "Charge of the receptor residues in QM region",

  "qmcharge_lig": 1,
  "qmcharge_lig_description": "Charge of the ligand residues in QM region",

  "submit_job": "True",
  "submit_job_description": "Automatically submit job after GBSA preparation (True or False)"
}
```
## 2. Run the script
```bash
module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 ambertools/25.0
python gbsa_prep.py -i input.json
```
## 3. Job Submission Options

- If `"submit_job": "True"` in your JSON file, the script will **automatically submit the job** with `sbatch`.  
- If `"submit_job": "False"`, the script will only prepare the input files. You will then need to submit manually:

```bash
cd ./<prefix>_gbsa/MMgbsa/
sbatch submit.job
```
## Contact
For any questions, please contact wook.shin@vanderbilt.edu.
