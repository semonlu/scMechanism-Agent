#!/usr/bin/env python3
"""Self-tests for the cross-platform environment setup helpers."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def load_module(rel_path: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class EnvSetupToolTests(unittest.TestCase):
    def test_minimal_profile_contains_core_seurat_runtime_packages(self):
        check_env = load_module("scripts/env_setup/check_environment.py", "check_environment_for_test")

        profile = check_env.get_profile("minimal")

        self.assertIn("scanpy", profile.python_imports)
        self.assertIn("Seurat", profile.r_packages)
        self.assertIn("SingleR", profile.r_packages)
        self.assertIn("celldex", profile.r_packages)
        self.assertNotIn("CellChat", profile.r_packages)

    def test_extended_profile_adds_optional_course_modules(self):
        check_env = load_module("scripts/env_setup/check_environment.py", "check_environment_extended_for_test")

        profile = check_env.get_profile("extended")

        self.assertIn("CellChat", profile.r_packages)
        self.assertIn("monocle3", profile.r_packages)
        self.assertIn("clusterProfiler", profile.r_packages)
        self.assertIn("org.Mm.eg.db", profile.r_packages)

    def test_json_status_marks_missing_imports_as_failed(self):
        check_env = load_module("scripts/env_setup/check_environment.py", "check_environment_status_for_test")

        report = check_env.build_report(
            profile=check_env.EnvironmentProfile(
                name="unit",
                python_imports=["definitely_missing_python_package_for_scmechanism"],
                r_packages=[],
                required_tools=[],
            ),
            json_ready=True,
        )

        self.assertEqual(report["status"], "FAILED")
        self.assertFalse(report["python"][0]["available"])

    def test_install_command_uses_cross_platform_conda_update(self):
        installer = load_module("scripts/env_setup/install_minimal_env.py", "install_minimal_env_for_test")

        command = installer.build_conda_command(
            conda="conda",
            env_name="scmechanism-agent",
            env_file="environment.yml",
            create=False,
        )

        self.assertEqual(
            command,
            ["conda", "env", "update", "-n", "scmechanism-agent", "-f", "environment.yml", "--prune"],
        )


if __name__ == "__main__":
    unittest.main()
