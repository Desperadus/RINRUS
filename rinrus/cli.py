from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from rinrus import paths


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
    if args.path_to_scripts:
        options["path_to_scripts"] = args.path_to_scripts

    missing = [k for k in ("pdb", "seed", "rin_program", "model") if k not in options]
    if missing:
        raise ValueError(f"Missing required options: {', '.join(missing)}")

    driver_input = _write_driver_input(options)
    try:
        cmd = [sys.executable, str(bin_path / "RINRUS_driver.py"), "-i", str(driver_input)]
        result = subprocess.run(cmd, env=env)
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
    driver.add_argument("--rin-program", dest="rin_program", help="probe, arpeggio, distance, or manual")
    driver.add_argument("--model", help="all, max, maximal, or a number")
    driver.add_argument("--qm-input-format", dest="qm_input_format", help="gaussian, orca, qchem, gau-xtb, psi4-fsapt")
    driver.add_argument("--seed-charge", dest="seed_charge", type=int, help="Seed charge")
    driver.add_argument("--multiplicity", type=int, help="Spin multiplicity")
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
