#!/usr/bin/env python3
"""Build an offline bundle of the workflow for an HPC without internet access.

Run this on a machine with internet access, git, Nextflow and (for --images)
Singularity/Apptainer. The bundle folder contains:

  workflow-nxf-somatic-variantcomp-<tag>.tar.gz   the pipeline code at a git tag
  containers.txt                                  all container images the pipeline uses
  singularity/                                    (--images) the images, named like
                                                  Nextflow's Singularity cache
  nxf_plugins/                                    Nextflow plugins (nf-schema)
  README.md                                       how to run the bundle offline

The container list is taken from `nextflow inspect` for the test profiles that
together use every process of the pipeline.

Example:
  scripts/make_offline_bundle.py --tag v1.0.0 --outdir /data/bundles/variantcomp --images
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

PIPELINE_DIR = Path(__file__).resolve().parent.parent
NAME = "workflow-nxf-somatic-variantcomp"
INSPECT_PROFILES = ["test", "test_seqc2"]


def run(cmd, **kwargs):
    print("$ " + " ".join(str(c) for c in cmd), file=sys.stderr)
    return subprocess.run([str(c) for c in cmd], check=True, **kwargs)


def singularity_cache_name(url):
    """File name Nextflow uses for an image in NXF_SINGULARITY_CACHEDIR."""
    name = url.split("://", 1)[1] if "://" in url else url
    name = name.replace(":", "-").replace("/", "-")
    return name if name.endswith(".img") else name + ".img"


def list_containers(workdir):
    images = set()
    for profile in INSPECT_PROFILES:
        out = run(
            [
                "nextflow",
                "inspect",
                PIPELINE_DIR,
                "-profile",
                f"{profile},singularity",
                "-format",
                "json",
            ],
            capture_output=True,
            text=True,
            cwd=workdir,
        )
        # log lines may surround the JSON document
        data, _end = json.JSONDecoder().raw_decode(out.stdout[out.stdout.index("{") :])
        images |= {p["container"] for p in data["processes"] if p.get("container")}
    return sorted(images)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n".join(__doc__.split("\n")[1:]),
    )
    parser.add_argument("--tag", required=True, help="git tag (or commit) to bundle")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--images", action="store_true", help="pull the Singularity images")
    args = parser.parse_args()

    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    work = out / ".work"
    work.mkdir(exist_ok=True)

    # 1. pipeline code at the tag
    archive = out / f"{NAME}-{args.tag}.tar.gz"
    run(
        [
            "git",
            "-C",
            PIPELINE_DIR,
            "archive",
            "--format=tar.gz",
            f"--prefix={NAME}/",
            "-o",
            archive,
            args.tag,
        ]
    )

    # 2. containers (also downloads the nf-schema plugin into NXF_HOME)
    plugins = out / "nxf_plugins"
    env_home = work / "nxf_home"
    os.environ["NXF_HOME"] = str(env_home)
    images = list_containers(work)
    (out / "containers.txt").write_text("\n".join(images) + "\n")
    if (env_home / "plugins").exists():
        shutil.copytree(env_home / "plugins", plugins, dirs_exist_ok=True)

    # 3. images
    if args.images:
        engine = shutil.which("singularity") or shutil.which("apptainer")
        if not engine:
            sys.exit("--images needs singularity or apptainer on PATH")
        cache = out / "singularity"
        cache.mkdir(exist_ok=True)
        for url in images:
            target = cache / singularity_cache_name(url)
            if target.exists():
                continue
            source = url if "://" in url else f"docker://{url}"
            run([engine, "pull", "--name", target, source])

    (out / "README.md").write_text(f"""# {NAME} offline bundle ({args.tag})

Created with `scripts/make_offline_bundle.py` from the workflow repository.

## Install on the HPC

```bash
tar -xzf {archive.name}            # -> {NAME}/
export NXF_OFFLINE=true
export NXF_PLUGINS_DIR=$PWD/nxf_plugins
export NXF_SINGULARITY_CACHEDIR=$PWD/singularity
```

## Run

```bash
nextflow run {NAME} \\
    -profile singularity \\
    -c <your-hpc.config> \\
    -params-file params.yaml \\
    --outdir results
```

An example HPC config (PBS Pro, Singularity, resource limits) is in
`{NAME}/docs/hpc.md`. Container images: `containers.txt`.
""")
    shutil.rmtree(work / "work", ignore_errors=True)
    print(f"bundle: {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
