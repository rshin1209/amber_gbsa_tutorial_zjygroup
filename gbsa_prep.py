"""
Prepare and run (QMMM-)MMGBSA from a JSON config.

- Creates stripped topologies (complex/receptor/ligand) via cpptraj
- Strips and autoimages the trajectory to a local md.nc
- Writes an MMPBSA input (MM or QM/MM)
- Creates a SLURM submission script and optionally submits it

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
from typing import Dict, Any

# ---------------------------
# Utilities
# ---------------------------

def run(cmd, cwd=None):
    """Run a shell command, raising on failure and echoing the command."""
    print(f"[run] {cwd or os.getcwd()}$ {cmd}")
    subprocess.run(shlex.split(cmd), cwd=cwd, check=True)

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

# ---------------------------
# Writers
# ---------------------------

def write_cpptraj_prmtop_in(out_path: Path, source_prmtop_glob: str, residue_mask: str, out_prmtop: str):
    """
    Write a cpptraj input that loads the source parm, strips to residue_mask,
    removes box info, and writes the stripped topology.
    """
    text = f"""parm {source_prmtop_glob}
parmstrip !(:{residue_mask})
parmbox nobox
parmwrite out {out_prmtop}
run
quit
"""
    out_path.write_text(text)

def write_cpptraj_striptraj_in(out_path: Path, source_prmtop_glob: str,
                               source_nc: str, residue_mask: str,
                               start: int, end: int, interval: int,
                               out_nc: str):
    """
    Write a cpptraj input that loads parm, reads trajectory frames with range/interval,
    autoimages, strips to residue_mask, and writes a boxless trajectory.
    """
    # cpptraj trajin syntax: trajin <file> [start [stop [offset]]]
    text = f"""parm {source_prmtop_glob}
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
    Uses the same shape as the user's original template.
    """
    level = str(data["level_of_theory"]).upper()

    header = (
        "Input file for running MM-GBSA\n"
        if level == "MM" else
        "Input file for running QMMM-GBSA\n"
    )

    # &general block
    general = (
        f"""&general
   verbose=1,
/
"""
    )

    # &gb block (with or without QM params)
    if level == "MM":
        gb = f"""&gb
   igb={data["igb"]},
   saltcon={data["saltcon"]},
/
"""
    else:
        # Keep the user's original keys/shape for QM flags
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

def write_slurm_job(out_path: Path, job_name: str, amber_env_script: str, nprocs: int):
    """
    Write a SLURM script that runs MMPBSA.py.MPI with the prepared inputs.
    Assumes this script lives inside the mmgbsa work dir and that inputs are one dir up.
    """
    text = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name={job_name}
#SBATCH --partition=production
#SBATCH --ntasks={nprocs}
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --account=yang_lab_csb

set -euo pipefail

source {amber_env_script}

echo "[$(date)] Running MMPBSA.py.MPI..."
mpirun -np {nprocs} "$AMBERHOME/bin/MMPBSA.py.MPI" -O \\
  -i ./*.in \\
  -cp ../complex.prmtop \\
  -rp ../receptor.prmtop \\
  -lp ../ligand.prmtop \\
  -y ../md.nc \\
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
    ap.add_argument("--amber-env", default="/home/shaoq1/bin/amber_env/amber-accre.sh",
                    help="Path to shell script that sets AMBER environment.")
    ap.add_argument("--procs", type=int, default=8, help="MPI ranks for MMPBSA.py.MPI")
    ap.add_argument("--dry-run", action="store_true", help="Prepare everything but do not submit.")
    args = ap.parse_args()

    # Load config
    with open(args.input, "r") as f:
        data = json.load(f)

    # Basic validation
    require_keys(data, [
        "directory", "complex_residues", "receptor_residues", "ligand_residues",
        "level_of_theory", "startframe", "endframe", "interval",
        "igb", "saltcon", "submit_job"
    ])
    level = str(data["level_of_theory"]).upper()
    if level != "MM":
        require_keys(data, ["qm_residues", "qmcharge_com", "qmcharge_rec", "qmcharge_lig"])

    # External tool checks (fail early with a helpful error)
    ensure_cmd("cpptraj")
    # MMPBSA.py.MPI path is resolved inside the SLURM job via $AMBERHOME

    # Resolve paths
    sim_dir = Path(data["directory"]).expanduser().resolve()
    if not sim_dir.exists():
        raise FileNotFoundError(f"directory not found: {sim_dir}")

    # Identify a prmtop to load in cpptraj (use glob, as per your original approach)
    # We'll pass the glob string directly to cpptraj, matching your pipeline.
    prmtop_glob = str(sim_dir / "*.prmtop")

    # Output workspace: <prefix> alongside this script / CWD
    prefix = sim_dir.name + "_gbsa"
    work = Path.cwd() / prefix
    work.mkdir(parents=True, exist_ok=True)

    # 1) Create stripped prmtops
    write_cpptraj_prmtop_in(work / "complex_prmtop.in", prmtop_glob, data["complex_residues"], "./complex.prmtop")
    write_cpptraj_prmtop_in(work / "receptor_prmtop.in", prmtop_glob, data["receptor_residues"], "./receptor.prmtop")
    write_cpptraj_prmtop_in(work / "ligand_prmtop.in",  prmtop_glob, data["ligand_residues"],  "./ligand.prmtop")

    run("cpptraj -i complex_prmtop.in", cwd=work)
    run("cpptraj -i receptor_prmtop.in", cwd=work)
    run("cpptraj -i ligand_prmtop.in",  cwd=work)

    # 2) Prepare stripped/autoimaged trajectory
    source_nc = str(sim_dir / "md.nc")
    if not Path(source_nc).exists():
        raise FileNotFoundError(f"Trajectory not found: {source_nc}")

    write_cpptraj_striptraj_in(
        work / "strip_traj.in",
        prmtop_glob,
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

    # 4) SLURM script (inside mmgbsa_dir)
    job_script = mmgbsa_dir / "submit.job"
    write_slurm_job(job_script, job_name=prefix, amber_env_script=args.amber_env, nprocs=args.procs)

    # 5) Optionally submit
    submit = coerce_bool(data["submit_job"])
    if submit and not args.dry_run:
        ensure_cmd("sbatch")
        run("sbatch submit.job", cwd=mmgbsa_dir)
        print(f"\nSubmitted: {job_script}")
    else:
        print("\nPreparation complete.")
        print(f"To run later:\n  (cd {mmgbsa_dir} && sbatch submit.job)")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)


