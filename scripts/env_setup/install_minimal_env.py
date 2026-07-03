#!/usr/bin/env python3
"""Create or update the cross-platform scMechanism-Agent conda environment."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def build_conda_command(conda: str, env_name: str, env_file: str, create: bool) -> list[str]:
    if create:
        return [conda, "env", "create", "-n", env_name, "-f", env_file]
    return [conda, "env", "update", "-n", env_name, "-f", env_file, "--prune"]


def run(command: list[str]) -> int:
    print("+ " + " ".join(command))
    return subprocess.call(command)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-name", default="scmechanism-agent")
    parser.add_argument("--env-file", default="environment.yml")
    parser.add_argument("--conda", default=shutil.which("conda") or "conda")
    parser.add_argument("--create", action="store_true", help="Use conda env create instead of env update")
    parser.add_argument("--profile", choices=["minimal", "extended"], default="minimal")
    parser.add_argument("--skip-r-packages", action="store_true", help="Only create/update the conda environment")
    parser.add_argument("--skip-check", action="store_true", help="Do not run the minimal environment check after install")
    args = parser.parse_args()

    env_file = Path(args.env_file)
    if not env_file.exists():
        print(f"Environment file not found: {env_file}", file=sys.stderr)
        return 2

    rc = run(build_conda_command(args.conda, args.env_name, str(env_file), args.create))
    if rc != 0:
        return rc

    if not args.skip_r_packages:
        r_installer = Path(__file__).with_name("install_r_packages.R")
        rc = run(
            [
                args.conda,
                "run",
                "-n",
                args.env_name,
                "Rscript",
                str(r_installer),
                "--profile",
                args.profile,
            ]
        )
        if rc != 0:
            return rc

    if args.skip_check:
        return 0

    checker = Path(__file__).with_name("check_environment.py")
    return run([args.conda, "run", "-n", args.env_name, "python", str(checker), "--profile", args.profile])


if __name__ == "__main__":
    sys.exit(main())
