from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from rinrus import paths


def _read_sdf_coords(path: Path) -> list[tuple[float, float, float]]:
    lines = path.read_text().splitlines()
    if len(lines) < 4:
        raise ValueError(f"Invalid SDF file: {path}")
    counts = lines[3]
    natoms = None
    try:
        natoms = int(counts[0:3])
    except ValueError:
        tokens = counts.split()
        if tokens and tokens[0].isdigit():
            natoms = int(tokens[0])
    if natoms is None or natoms <= 0:
        raise ValueError(f"Unable to read atom count from SDF: {path}")
    start = 4
    coords: list[tuple[float, float, float]] = []
    for i in range(natoms):
        line = lines[start + i]
        tokens = line.split()
        if len(tokens) < 3:
            raise ValueError(f"Invalid atom line in SDF: {path}")
        coords.append((float(tokens[0]), float(tokens[1]), float(tokens[2])))
    return coords


def _read_pdb_coords(path: Path) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    for line in path.read_text().splitlines():
        record = line[:6]
        if record not in ("ATOM  ", "HETATM"):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        coords.append((x, y, z))
    return coords


def _read_ligand_coords(path: Path) -> list[tuple[float, float, float]]:
    suffix = path.suffix.lower()
    if suffix in (".sdf", ".mol"):
        return _read_sdf_coords(path)
    if suffix == ".pdb":
        return _read_pdb_coords(path)
    raise ValueError(f"Unsupported ligand format: {path}")


def _read_pdb_atoms(path: Path, include_hetatm: bool) -> list[tuple[str, int, str, float, float, float]]:
    atoms: list[tuple[str, int, str, float, float, float]] = []
    for line in path.read_text().splitlines():
        record = line[:6]
        if record == "HETATM" and not include_hetatm:
            continue
        if record not in ("ATOM  ", "HETATM"):
            continue
        resn = line[17:20].strip()
        if resn in ("HOH", "WAT", "H2O"):
            continue
        chain = line[21].strip()
        try:
            resi = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        atoms.append((chain, resi, resn, x, y, z))
    return atoms


def _format_seed(residues: list[tuple[str, int]]) -> str:
    parts = []
    for chain, resi in residues:
        parts.append(f"{chain}:{resi}" if chain else f":{resi}")
    return ",".join(parts)


def _merge_seed_strings(existing: str, extra: str) -> str:
    def _split(seed: str) -> list[str]:
        return [s for s in seed.replace(" ", "").split(",") if s]

    seen = set()
    merged: list[str] = []
    for item in _split(existing) + _split(extra):
        if item not in seen:
            seen.add(item)
            merged.append(item)
    return ",".join(merged)


def _select_seed_from_ligand(
    pdb_path: Path, ligand_path: Path, cutoff: float | None, include_hetatm: bool
) -> str:
    ligand_coords = _read_ligand_coords(ligand_path)
    if not ligand_coords:
        raise ValueError(f"No coordinates found in ligand file: {ligand_path}")
    pdb_atoms = _read_pdb_atoms(pdb_path, include_hetatm=include_hetatm)
    if not pdb_atoms:
        raise ValueError(f"No atoms found in PDB file: {pdb_path}")

    if cutoff is None:
        cutoff = 4.0
    cutoff_sq = cutoff * cutoff
    selected = set()
    for chain, resi, _resn, x, y, z in pdb_atoms:
        for lx, ly, lz in ligand_coords:
            dx = x - lx
            dy = y - ly
            dz = z - lz
            if (dx * dx + dy * dy + dz * dz) <= cutoff_sq:
                selected.add((chain, resi))
                break
    if not selected:
        raise ValueError("No residues found within cutoff of ligand")
    residues = sorted(selected, key=lambda item: (item[0], item[1]))
    return _format_seed(residues)


def _append_env_path(env: dict, key: str, value: Path) -> None:
    value_str = str(value)
    current = env.get(key, "")
    entries = [e for e in current.split(os.pathsep) if e]
    if value_str not in entries:
        entries.insert(0, value_str)
    env[key] = os.pathsep.join(entries)


def _runtime_env(bin_path: Path, lib3_path: Path, home_path: Path) -> dict:
    env = os.environ.copy()
    _append_env_path(env, "PATH", bin_path)
    _append_env_path(env, "PYTHONPATH", lib3_path)
    env.setdefault("RINRUS_HOME", str(home_path))
    return env


def _write_driver_input(options: dict[str, str]) -> Path:
    lines = []
    for key, value in options.items():
        if value is None:
            continue
        lines.append(f"{key}: {value}")
    content = "\n".join(lines) + "\n"

    tmp = tempfile.NamedTemporaryFile(mode="w", prefix="rinrus_", suffix=".inp", delete=False)
    try:
        tmp.write(content)
        return Path(tmp.name)
    finally:
        tmp.close()


def _run_driver(args: argparse.Namespace) -> int:
    home_path = paths.rinrus_home()
    bin_path = paths.bin_dir()
    lib3_path = paths.lib3_dir()

    if not (bin_path / "RINRUS_driver.py").is_file():
        raise FileNotFoundError(f"RINRUS driver not found under {bin_path}")
    if not lib3_path.is_dir():
        raise FileNotFoundError(f"RINRUS lib3 not found under {lib3_path}")

    env = _runtime_env(bin_path, lib3_path, home_path)

    if args.driver_input:
        if args.ligand:
            raise ValueError("--ligand cannot be used with --input")
        driver_input = Path(args.driver_input)
        if not driver_input.is_file():
            raise FileNotFoundError(f"Driver input file not found: {driver_input}")
        cmd = [sys.executable, str(bin_path / "RINRUS_driver.py"), "-i", str(driver_input)]
        result = subprocess.run(cmd, env=env)
        return result.returncode

    options: dict[str, str] = {}
    for item in args.set or []:
        if "=" not in item:
            raise ValueError(f"Invalid --set entry (expected key=value): {item}")
        key, value = item.split("=", 1)
        options[key.strip().lower()] = value.strip()

    if args.pdb:
        options["pdb"] = args.pdb
    if args.seed:
        options["seed"] = args.seed
    if args.rin_program:
        options["rin_program"] = args.rin_program
    if args.model:
        options["model"] = args.model
    if args.qm_input_format:
        options["qm_input_format"] = args.qm_input_format
    if args.seed_charge is not None:
        options["seed_charge"] = str(args.seed_charge)
    if args.multiplicity is not None:
        options["multiplicity"] = str(args.multiplicity)
    if args.ph is not None:
        options["model_prot_ph"] = str(args.ph)
    if args.path_to_scripts:
        options["path_to_scripts"] = args.path_to_scripts

    if args.auto is not None:
        if not args.ligand:
            raise ValueError("--auto requires --ligand")
        options.setdefault("rin_program", "distance")
        options.setdefault("model", "max")
        options.setdefault("dist_type", "closest")
        if args.ligand_cutoff is None:
            args.ligand_cutoff = args.auto

    if args.ligand:
        if "pdb" not in options:
            raise ValueError("--ligand requires --pdb")
        ligand_seed = _select_seed_from_ligand(
            Path(options["pdb"]),
            Path(args.ligand),
            args.ligand_cutoff,
            args.ligand_include_hetatm,
        )
        if "seed" in options:
            options["seed"] = _merge_seed_strings(options["seed"], ligand_seed)
        else:
            options["seed"] = ligand_seed

    missing = [k for k in ("pdb", "seed", "rin_program", "model") if k not in options]
    if missing:
        raise ValueError(f"Missing required options: {', '.join(missing)}")

    driver_input = _write_driver_input(options)
    try:
        cmd = [sys.executable, str(bin_path / "RINRUS_driver.py"), "-i", str(driver_input)]
        result = subprocess.run(cmd, env=env)
        if result.returncode == 0 and args.output:
            output_path = Path(args.output)
            if output_path.exists() and output_path.is_dir():
                target_dir = output_path
            else:
                target_dir = output_path.parent if output_path.suffix else output_path
            if not target_dir.exists():
                raise FileNotFoundError(f"Output directory does not exist: {target_dir}")
            candidates = []
            for item in Path.cwd().iterdir():
                if not item.is_file():
                    continue
                match = re.match(r"res_(\d+)_h\.pdb$", item.name)
                if match:
                    candidates.append((int(match.group(1)), item))
            if not candidates:
                raise FileNotFoundError("No res_*_h.pdb files found to copy")
            candidates.sort(key=lambda entry: entry[0])
            _, source = candidates[-1]
            if output_path.exists() and output_path.is_dir():
                dest = output_path / source.name
            elif output_path.suffix:
                dest = output_path
            else:
                dest = output_path / source.name
            shutil.copyfile(source, dest)
        return result.returncode
    finally:
        try:
            driver_input.unlink()
        except OSError:
            pass


def _run_script(args: argparse.Namespace) -> int:
    home_path = paths.rinrus_home()
    bin_path = paths.bin_dir()
    lib3_path = paths.lib3_dir()
    env = _runtime_env(bin_path, lib3_path, home_path)

    script = args.script
    script_path = bin_path / script
    if not script_path.is_file() and not script.endswith(".py"):
        script_path = bin_path / f"{script}.py"
    if not script_path.is_file():
        raise FileNotFoundError(f"Script not found under {bin_path}: {script}")

    cmd = [sys.executable, str(script_path)]
    if args.args:
        cmd.extend(args.args)
    result = subprocess.run(cmd, env=env)
    return result.returncode


def main() -> None:
    parser = argparse.ArgumentParser(prog="rinrus", description="RINRUS command line interface")
    subparsers = parser.add_subparsers(dest="command", required=True)

    driver = subparsers.add_parser("driver", help="Run the RINRUS driver")
    driver.add_argument("-i", "--input", dest="driver_input", help="Driver input file")
    driver.add_argument("--pdb", help="Input PDB filename")
    driver.add_argument("--seed", help="Seed fragment(s), e.g. A:300,A:301")
    driver.add_argument("--ligand", help="Ligand file (SDF/MOL/PDB) used to auto-select seed residues")
    driver.add_argument("--ligand-cutoff", type=float, default=None, help="Ligand contact cutoff in Angstrom")
    driver.add_argument(
        "--ligand-include-hetatm",
        action="store_true",
        help="Include HETATM residues as potential seeds",
    )
    driver.add_argument("--auto", type=float, help="Shortcut for ligand-driven run (sets cutoff/rin_program/model/dist_type)")
    driver.add_argument("--rin-program", dest="rin_program", help="probe, arpeggio, distance, or manual")
    driver.add_argument("--model", help="all, max, maximal, or a number")
    driver.add_argument("--qm-input-format", dest="qm_input_format", help="gaussian, orca, qchem, gau-xtb, psi4-fsapt")
    driver.add_argument("--seed-charge", dest="seed_charge", type=int, help="Seed charge")
    driver.add_argument("--multiplicity", type=int, help="Spin multiplicity")
    driver.add_argument("--ph", type=float, help="Use pdb2pqr to protonate models at this pH")
    driver.add_argument("--output", help="Copy final res_<N>_h.pdb to this path")
    driver.add_argument("--path-to-scripts", dest="path_to_scripts", help="Override bin directory path")
    driver.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Additional driver options (repeatable)",
    )

    run_script = subparsers.add_parser("run", help="Run a bundled script from bin/")
    run_script.add_argument("script", help="Script name (e.g. probe2rins.py)")
    run_script.add_argument("args", nargs=argparse.REMAINDER)

    args = parser.parse_args()
    try:
        if args.command == "driver":
            code = _run_driver(args)
        else:
            code = _run_script(args)
    except (FileNotFoundError, ValueError) as exc:
        print(f"rinrus: {exc}", file=sys.stderr)
        sys.exit(2)
    sys.exit(code)
