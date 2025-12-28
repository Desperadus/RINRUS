from __future__ import annotations

import os
import sys
from pathlib import Path


def _expand_path(path: str | Path) -> Path:
    return Path(path).expanduser().resolve()


def rinrus_home() -> Path:
    env_home = os.environ.get("RINRUS_HOME")
    if env_home:
        return _expand_path(env_home)

    here = Path(__file__).resolve().parent
    for parent in [here] + list(here.parents):
        if (parent / "bin" / "RINRUS_driver.py").is_file():
            return parent

    share_root = Path(sys.prefix) / "share" / "rinrus"
    if (share_root / "bin" / "RINRUS_driver.py").is_file():
        return share_root

    raise FileNotFoundError(
        "Unable to locate RINRUS home. Set RINRUS_HOME to the repo root."
    )


def bin_dir() -> Path:
    env_bin = os.environ.get("RINRUS_BIN")
    if env_bin:
        return _expand_path(env_bin)
    return rinrus_home() / "bin"


def lib3_dir() -> Path:
    env_lib3 = os.environ.get("RINRUS_LIB3")
    if env_lib3:
        return _expand_path(env_lib3)
    return rinrus_home() / "lib3"


def template_dir() -> Path:
    env_template = os.environ.get("RINRUS_TEMPLATE_DIR")
    if env_template:
        return _expand_path(env_template)
    return rinrus_home() / "template_files"
