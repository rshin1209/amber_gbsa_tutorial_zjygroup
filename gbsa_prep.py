#!/usr/bin/env python3
"""
Prepare and run (QMMM-)MMGBSA from a JSON config.

- Creates a stripped COMPLEX topology (complex.prmtop) via cpptraj
- Strips and autoimages the trajectory to a local md.nc
- Writes an MMPBSA input (MM or QM/MM)
- Creates a SLURM submission script and optionally submits it

IMPORTANT CHANGE vs your old version:
- We DO NOT create receptor.prmtop / ligand.prmtop anymore.
- We let MMPBSA.py build receptor/ligand from complex.prmtop using -m OR -n.

Mask semantics (MMPBSA):
  -m / --receptor-mask : atoms/residues to STRIP FROM COMPLEX to create RECEPTOR
  -n / --ligand-mask   : atoms/residues to STRIP FROM COMPLEX to create LIGAND
You cannot use both -m and -n at the same time.

This script will, by default, use:
  -m ":<ligand_residues>"
(i.e., "strip ligand from complex to make receptor")
If you want the opposite mode, set in JSON:
  "mmpbsa_mask_mode": "n"
Then it will use:
  -n ":<receptor_residues>"

Usage:
    python prepare_gbsa.py -i input.json
"""

import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, Any, Optional


# ---------------------------
# Utilities
# ---------------------------

def run(cmd: str, cwd: Optional[Path] = None):
    """Run a shell command, raising on failure and echoing the command."""
    cwd_str = str(cwd) if cwd else os.getcwd()
    print(f"[run] {cwd_str}$ {cmd}")
    subprocess.run(shlex.split(cmd), cwd=str(cwd) if cwd else None, check=True)

def ensure_cmd(name: str):
    """Check that a command exists in PATH."""
    from shutil import which
    if which(name) is None:
        raise RuntimeError(f"Required command '{name}' not found in PATH.")

def coerce_bool(x) -> bool:
    if isinstance(x, bool):
        return x
    if isinstance(x, (int, float)):
        return bool(x)
    if isinstance(x, str):
        return x.strip().lower() in {"1", "true", "yes", "y", "t"}
    return False

def require_keys(d: Dict[str, Any], keys):
    missing = [k for k in keys if k not in d]
    if missing:
        raise KeyError(f"Missing required config key(s): {', '.join(missing)}")

def pick_single_file(sim_dir: Path, pattern: str) -> Path:
    """Pick a single file matching pattern; error if none or multiple."""
    hits = sorted(sim_dir.glob(pattern))
    if len(hits) == 0:
        raise FileNotFoundError(f"No files found matching {sim_dir / pattern}")
    if len(hits) > 1:
        # You can relax this if you truly have multiple prmtops, but most workflows don't.
        msg = "\n".join(str(p) for p in hits[:10])
        raise RuntimeError(
            f"Multiple files found matching {sim_dir / pattern}. "
            f"Please keep exactly one or change the script.\nFound:\n{msg}"
        )
    return hits[0]


# ---------------------------
# Writers
# ---------------------------

def write_cpptraj_prmtop_in(out_path: Path, source_prmtop: Path, residue_mask: str, out_prmtop: str):
    """
    Write a cpptraj input that loads the source parm, strips to residue_mask,
    removes box info, and writes the stripped topology.
    NOTE: residue_mask should be residue ranges only, e.g. "1-530,1062"
    """
    text = f"""parm {source_prmtop}
parmstrip !(:{residue_mask})
parmbox nobox
parmwrite out {out_prmtop}
run
quit
"""
    out_path.write_text(text)

def write_cpptraj_striptraj_in(out_path: Path, source_prmtop: Path,
                               source_nc: Path, residue_mask: str,
                               start: int, end: int, interval: int,
                               out_nc: str):
    """
    Write a cpptraj input that loads parm, reads trajectory frames with range/interval,
    autoimages, strips to residue_mask, and writes a boxless trajectory.
    """
    text = f"""parm {source_prmtop}
trajin {source_nc} {start} {end} {interval}
autoimage
strip !(:{residue_mask})
trajout {out_nc} nobox
run
quit
"""
    out_path.write_text(text)

def write_mmgbsa_input(out_path: Path, data: Dict[str, Any]):
    """
    Write MMPBSA input file for MM or QM/MM based on level_of_theory.
    """
    level = str(data["level_of_theory"]).upper()

    header = (
        "Input file for running MM-GBSA\n"
        if level == "MM" else
        "Input file for running QMMM-GBSA\n"
    )

    general = """&general
   verbose=1,
/
"""

    if level == "MM":
        gb = f"""&gb
   igb={data["igb"]},
   saltcon={data["saltcon"]},
/
"""
    else:
        gb = f"""&gb
   igb={data["igb"]},
   saltcon={data["saltcon"]},
   ifqnt=1,
   qm_theory="{data["level_of_theory"]}",
   qm_residues="{data["qm_residues"]}",
   qmcharge_com={data["qmcharge_com"]},
   qmcharge_rec={data["qmcharge_rec"]},
   qmcharge_lig={data["qmcharge_lig"]},
/
"""
    out_path.write_text(header + general + gb)

def write_slurm_job(out_path: Path, job_name: str, nprocs: int,
                    level: str,
                    mask_mode: str,
                    receptor_residues: str,
                    ligand_residues: str):
    """
    Write a SLURM script that runs MMPBSA.py.MPI using -m OR -n (no -rp/-lp).
    mask_mode:
      - "m" => use -m ":<ligand_residues>"  (strip ligand to build receptor)
      - "n" => use -n ":<receptor_residues>" (strip receptor to build ligand)
    """
    mask_mode = str(mask_mode).strip().lower()
    if mask_mode not in {"m", "n"}:
        raise ValueError('mmpbsa_mask_mode must be "m" or "n".')

    if mask_mode == "m":
        # receptor is made by stripping ligand from complex
        mask_arg = f'-m ":{ligand_residues}"'
    else:
        # ligand is made by stripping receptor from complex
        mask_arg = f'-n ":{receptor_residues}"'

    text = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name={job_name}
#SBATCH --partition=production
#SBATCH --ntasks={nprocs}
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --account=yang_lab

set -euo pipefail

source /home/shaoq1/bin/amber_env/amber-accre.sh

echo "[$(date)] Running MMPBSA.py.MPI ({level}, mask_mode={mask_mode})..."
mpirun -np {nprocs} "$AMBERHOME/bin/MMPBSA.py.MPI" -O \\
  -i ./*.in \\
  -cp ../complex.prmtop \\
  -y ../md.nc \\
  {mask_arg} \\
  > progress.log 2>&1

echo "[$(date)] Done."
"""
    out_path.write_text(text)
    out_path.chmod(0o755)


# ---------------------------
# Pipeline
# ---------------------------

def main():
    ap = argparse.ArgumentParser(description="Prepare and (optionally) run (QMMM-)MMGBSA from JSON.")
    ap.add_argument("-i", "--input", required=True, help="Path to JSON config.")
    ap.add_argument("--procs", type=int, default=16, help="MPI ranks for MMPBSA.py.MPI")
    ap.add_argument("--dry-run", action="store_true", help="Prepare everything but do not submit.")
    args = ap.parse_args()

    # Load config
    with open(args.input, "r") as f:
        data = json.load(f)

    # Basic validation
    require_keys(data, [
        "directory",
        "complex_residues",
        "receptor_residues",
        "ligand_residues",
        "level_of_theory",
        "startframe",
        "endframe",
        "interval",
        "igb",
        "saltcon",
        "submit_job",
    ])

    level = str(data["level_of_theory"]).upper()
    if level != "MM":
        require_keys(data, ["qm_residues", "qmcharge_com", "qmcharge_rec", "qmcharge_lig"])

    # External tool checks (fail early with a helpful error)
    ensure_cmd("cpptraj")
    # MMPBSA.py.MPI path is resolved inside the SLURM job via $AMBERHOME

    # Resolve simulation dir
    sim_dir = Path(data["directory"]).expanduser().resolve()
    if not sim_dir.exists():
        raise FileNotFoundError(f"directory not found: {sim_dir}")

    # Pick a single prmtop; avoids glob ambiguity inside cpptraj
    prmtop_path = pick_single_file(sim_dir, "*.prmtop")

    # Trajectory
    source_nc = sim_dir / "md.nc"
    if not source_nc.exists():
        raise FileNotFoundError(f"Trajectory not found: {source_nc}")

    # Output workspace: <sim_dir_name>_gbsa in CWD
    prefix = sim_dir.name + "_gbsa"
    work = Path.cwd() / prefix
    work.mkdir(parents=True, exist_ok=True)

    # Mask mode selection:
    #   default "m" => use -m ":ligand_residues"
    #   set "mmpbsa_mask_mode": "n" in JSON if you want -n ":receptor_residues"
    mask_mode = str(data.get("mmpbsa_mask_mode", "m")).strip().lower()

    # 1) Create stripped complex.prmtop (this is the COMPLEX MMPBSA will strip from)
    write_cpptraj_prmtop_in(
        work / "complex_prmtop.in",
        prmtop_path,
        data["complex_residues"],
        "./complex.prmtop",
    )
    run("cpptraj -i complex_prmtop.in", cwd=work)

    # 2) Prepare stripped/autoimaged trajectory matching complex_residues
    write_cpptraj_striptraj_in(
        work / "strip_traj.in",
        prmtop_path,
        source_nc,
        data["complex_residues"],
        int(data["startframe"]),
        int(data["endframe"]),
        int(data["interval"]),
        "./md.nc",
    )
    run("cpptraj -i strip_traj.in", cwd=work)

    # 3) Prepare MMPBSA input + job dir
    mmgbsa_dirname = f"{level}gbsa"   # e.g., MMgbsa or PM6gbsa
    mmgbsa_dir = work / mmgbsa_dirname
    mmgbsa_dir.mkdir(parents=True, exist_ok=True)

    mmgbsa_in = mmgbsa_dir / f"{level}gbsa.in"
    write_mmgbsa_input(mmgbsa_in, data)

    # 4) SLURM script (inside mmgbsa_dir): uses -m OR -n (no -rp/-lp)
    job_script = mmgbsa_dir / "submit.job"
    write_slurm_job(
        job_script,
        job_name=prefix,
        nprocs=args.procs,
        level=level,
        mask_mode=mask_mode,
        receptor_residues=data["receptor_residues"],
        ligand_residues=data["ligand_residues"],
    )

    # 5) Optionally submit
    submit = coerce_bool(data["submit_job"])
    if submit and not args.dry_run:
        ensure_cmd("sbatch")
        run("sbatch submit.job", cwd=mmgbsa_dir)
        print(f"\nSubmitted: {job_script}")
    else:
        print("\nPreparation complete.")
        print(f"To run later:\n  (cd {mmgbsa_dir} && sbatch submit.job)")
        if mask_mode == "m":
            print(f"Using MMPBSA strip mode: -m \":{data['ligand_residues']}\"  (strip ligand to build receptor)")
        else:
            print(f"Using MMPBSA strip mode: -n \":{data['receptor_residues']}\" (strip receptor to build ligand)")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
