#!/usr/bin/env python3
"""Check that the pipeline compiles under Nextflow's strict syntax (run before every commit).

Two checks, both must pass:

1. `nextflow inspect` (what `nf-core pipelines download` runs) compiles the whole pipeline
   with the profiles singularity,test,test_full; once with the minimum Nextflow version of
   the manifest (nextflowVersion '!>=X') and once with the Nextflow on PATH (or --versions).
   `nextflow lint` alone is not enough: it reported 'no errors' for a file with a stray
   closing brace that `nextflow run` / `inspect` rejected ('Script compilation failed').
2. `nextflow lint` of the whole repository: no errors, and no warnings except the accepted
   ones (ACCEPTED_WARNINGS: the `when:` section, which nf-core module lint requires).

Exit status 0 if both pass, 1 otherwise; the failing output is printed.

Examples
  scripts/check_syntax.py
  scripts/check_syntax.py --versions 25.10.4 26.10.2 --profiles singularity,test,test_full
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

PIPELINE_DIR = Path(__file__).resolve().parent.parent
# lint warnings that are accepted: nf-core module lint fails a module without the
# `when: task.ext.when == null || task.ext.when` block that strict syntax deprecates
ACCEPTED_WARNINGS = ["The `when` section is deprecated"]


def minimum_version():
    """Minimum Nextflow version of the manifest (nextflowVersion = '!>=X')."""
    text = (PIPELINE_DIR / "nextflow.config").read_text()
    match = re.search(r"nextflowVersion\s*=\s*'!?>=\s*([\d.]+)'", text)
    if not match:
        sys.exit("ERROR: nextflow.config: no nextflowVersion = '!>=X' in the manifest")
    return match.group(1)


def inspect(version, profiles):
    """(ok, output) of `nextflow inspect` with a Nextflow version (None: the one on PATH)."""
    env = dict(os.environ)
    if version:
        env["NXF_VER"] = version
    # run in a temporary folder: inspect writes .nextflow/ and a log
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [
            "nextflow",
            "inspect",
            "-format",
            "json",
            "-profile",
            profiles,
            str(PIPELINE_DIR / "main.nf"),
        ]
        result = subprocess.run(cmd, cwd=tmp, env=env, capture_output=True, text=True)
        log = Path(tmp, ".nextflow.log")
        details = log.read_text() if result.returncode and log.exists() else ""
    output = result.stderr + result.stdout[-2000:] + details[-4000:]
    return result.returncode == 0, output


def lint():
    """(ok, problems) of `nextflow lint` of the repository: errors and unaccepted warnings."""
    result = subprocess.run(
        ["nextflow", "lint", "-o", "concise", "."], cwd=PIPELINE_DIR, capture_output=True, text=True
    )
    lines = (result.stdout + result.stderr).splitlines()
    # findings start with 'Error' or 'Warn'; files of the nf-test work folders are skipped
    findings = [ln for ln in lines if ln.startswith(("Error", "Warn")) and "/.nf-test/" not in ln]
    problems = [ln for ln in findings if not any(a in ln for a in ACCEPTED_WARNINGS)]
    accepted = len(findings) - len(problems)
    return result.returncode == 0 and not problems, problems, accepted


def main():
    """Run inspect for every version and lint; report and set the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--versions", nargs="*", help="Nextflow versions (default: manifest minimum and PATH)"
    )
    parser.add_argument(
        "--profiles", default="singularity,test,test_full", help="profiles for inspect"
    )
    args = parser.parse_args()
    versions = args.versions or [minimum_version(), None]
    failed = False
    for version in versions:
        ok, output = inspect(version, args.profiles)
        label = version or "on PATH"
        print(
            f"{'OK  ' if ok else 'FAIL'} nextflow inspect -profile {args.profiles} (Nextflow {label})"
        )
        if not ok:
            failed = True
            print(output, file=sys.stderr)
    ok, problems, accepted = lint()
    print(
        f"{'OK  ' if ok else 'FAIL'} nextflow lint: {len(problems)} problems, {accepted} accepted warnings"
    )
    if not ok:
        failed = True
        print("\n".join(problems), file=sys.stderr)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
