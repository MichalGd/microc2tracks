"""Small shared-launcher/access regressions; no environment installs or real data."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from check_shared_access import check_access
from workflow_common import VERSION, sha256


def bash_path(path):
    value = Path(path).resolve().as_posix()
    return "/" + value[0].lower() + value[2:] if os.name == "nt" else value


class SharedInstallationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.release = self.root / "read only release"
        self.scripts = self.release / "scripts"
        self.scripts.mkdir(parents=True)
        shutil.copy2(ROOT / "scripts/shared_launcher.sh", self.scripts / "shared_launcher.sh")
        self.bash = Path("C:/Program Files/Git/bin/bash.exe") if os.name == "nt" else Path("/bin/bash")
        self.prefix = self.root / "shared environment"
        (self.prefix / "bin").mkdir(parents=True)
        shim = self.prefix / "bin/python"
        shim.write_text("#!/usr/bin/env bash\nexec " + shlex.quote(Path(sys.executable).as_posix()) + ' "$@"\n', newline="\n")
        shim.chmod(0o755)
        self.config = self.root / "site.conf"
        self.config.write_text("ENV_PREFIX=" + shlex.quote(bash_path(self.prefix)) + "\nEXTRA_TOOL_PATH=''\n", newline="\n")
        self.cwd = self.root / "user project"
        self.cwd.mkdir()
        self.env = {**os.environ, "MICROC2TRACKS_LAUNCHER_CONFIG": str(self.config),
                    "PYTHONPATH": "unrelated-private-packages", "PYTHON": "must-not-be-used"}

    def launch(self, *args):
        return subprocess.run([str(self.bash), str(self.scripts / "shared_launcher.sh"), *args],
            cwd=self.cwd, env=self.env, capture_output=True, text=True)

    def test_launcher_preserves_arguments_cwd_and_pins_environment_without_activation(self):
        code = ('import json, os, sys\nprint(json.dumps({"args":sys.argv[1:], "cwd":os.getcwd(), '
                '"prefix":os.environ["CONDA_PREFIX"], "no_user_site":sys.flags.no_user_site, '
                '"no_bytecode":sys.dont_write_bytecode, "pythonpath":os.environ.get("PYTHONPATH"), '
                '"path":os.environ["PATH"]}))\n')
        (self.scripts / "pipeline.py").write_text(code)
        (self.scripts / "check_shared_access.py").write_text(code)
        for command in ("run", "preflight", "merge", "check-access"):
            result = self.launch(command, "-c", "config with spaces.conf")
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(result.stdout)
            if command == "check-access":
                self.assertEqual(payload["args"][0], "--env-prefix")
                self.assertEqual(payload["args"][-2:], ["-c", "config with spaces.conf"])
            else:
                self.assertEqual(payload["args"], [command, "-c", "config with spaces.conf"])
            self.assertEqual(Path(payload["cwd"]), self.cwd)
            self.assertTrue(payload["no_bytecode"])
            self.assertTrue(payload["no_user_site"])
            self.assertIsNone(payload["pythonpath"])
            if os.name == "nt":
                self.assertEqual(Path(payload["path"].split(os.pathsep)[0]), self.prefix / "bin")
            else:
                self.assertTrue(payload["path"].startswith(str(self.prefix) + "/bin:"))
        result = self.launch("--version")
        self.assertEqual(json.loads(result.stdout)["args"], ["run", "--version"])

    def test_inaccessible_shared_python_fails_without_using_fallback(self):
        self.config.write_text("ENV_PREFIX='/definitely_missing_microc2tracks_environment'\n")
        result = self.launch("--version")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("shared Python is inaccessible", result.stderr)

    def test_actual_version_runs_without_writes_in_release(self):
        for path in (ROOT / "scripts").glob("*.py"):
            shutil.copy2(path, self.scripts / path.name)
        shutil.copytree(ROOT / "config", self.release / "config")
        before = {str(p.relative_to(self.release)): sha256(p) for p in self.release.rglob("*") if p.is_file()}
        # On POSIX this exercises a genuinely non-writable installed source tree.
        if os.name != "nt":
            for path in self.release.rglob("*"):
                path.chmod(0o555 if path.is_dir() else 0o444)
            self.release.chmod(0o555)
            self.addCleanup(self.restore_permissions)
        result = self.launch("--version")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "microc2tracks " + VERSION)
        after = {str(p.relative_to(self.release)): sha256(p) for p in self.release.rglob("*") if p.is_file()}
        self.assertEqual(before, after)
        self.assertFalse(list(self.release.rglob("__pycache__")))

    def restore_permissions(self):
        self.release.chmod(0o755)
        for path in self.release.rglob("*"):
            path.chmod(0o755 if path.is_dir() else 0o644)

    def test_reference_and_jar_access_checks_do_not_modify_inputs(self):
        import test_integration_features as fixtures
        fixture = fixtures.IntegrationFeatures("test_digest_is_complete_and_cache_reuses_verified_content")
        fixture.setUp()
        self.addCleanup(fixture.tearDown)
        jar = fixture.root / "existing.jar"
        jar.write_bytes(b"readability fixture, not an actual JAR")
        before = {str(p): sha256(p) for p in fixture.root.iterdir() if p.is_file()}
        with patch("check_shared_access.shutil.which", return_value=str(self.prefix / "bin/python")):
            _, warnings, errors = check_access(self.prefix, registry=fixture.registry, jar=jar)
        self.assertEqual(errors, [])
        self.assertFalse(any("not checked" in warning for warning in warnings))
        self.assertEqual(before, {str(p): sha256(p) for p in fixture.root.iterdir() if p.is_file()})
        Path(str(fixture.fasta) + ".pac").unlink()
        with patch("check_shared_access.shutil.which", return_value=str(self.prefix / "bin/python")):
            self.assertTrue(check_access(self.prefix, registry=fixture.registry, jar=jar)[2])
            self.assertEqual(check_access(self.prefix, registry=fixture.registry, jar=jar, refilter=True)[2], [])

    def test_access_failures_and_unchecked_resources_are_not_reported_as_verified(self):
        with patch("check_shared_access.shutil.which", return_value=None):
            checks, warnings, errors = check_access(self.prefix)
        self.assertTrue(any("pairtools" in row for row in errors))
        self.assertTrue(any("JAR was not checked" in row for row in warnings))
        self.assertTrue(any("Reference paths were not checked" in row for row in warnings))
        self.assertFalse(any("reference " in row for row in checks))
        with patch("check_shared_access.shutil.which", return_value=str(self.prefix / "bin/python")), \
             patch("check_shared_access.os.access", return_value=False):
            errors = check_access(self.prefix)[2]
        self.assertTrue(any("shared Python" in row for row in errors))
        self.assertTrue(any("release" in row for row in errors))


if __name__ == "__main__":
    unittest.main()
