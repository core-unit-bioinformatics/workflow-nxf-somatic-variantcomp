#!/usr/bin/env python3
"""Download, prepare and run the SEQC2 HCC1395/HCC1395BL WES + WGS test dataset.

Subcommands:

  prepare  download the public SEQC2 files (with md5 verification where SEQC2
           provides checksums), subset them for the chosen mode, and write a
           ready-to-run dataset folder: reference, target BED, unfiltered
           MuTect2 VCFs, optional BAM slices, the SEQC2 truth set,
           samplesheet.tsv, checksums.md5 and README.md (provenance).
  run      run the workflow on a prepared dataset folder.
  all      prepare + run.

Modes (see seqc2_testdata.yaml):

  small    chr22 only, 3 pseudo-samples, BAM slices (two 1-Mb windows) for one.
           A few tens of MB, runs on a laptop / Raspberry Pi with -profile conda.
           This is the dataset committed to the test-data repository.
  full     genome-wide, 12 pseudo-samples, VCFs only by default (the full BAMs
           are 10-160 GB each; add --with-bams to download them). Meant for an HPC.

Pseudo-samples are SEQC2 center/replicate pairs (e.g. IL_1 = WGS_IL_1 + WES_IL_1).
All pairs are sequenced from the same DNA batch, so each pair behaves like one
tumor/normal "patient" sequenced with WGS and WES.

Required tools on PATH: curl, samtools, bcftools, tabix, bgzip, bedtools
(all in envs/variantcomp-dev.yaml).

Examples:
  seqc2_testdata.py prepare --mode small
  seqc2_testdata.py run --mode small            # -profile test_seqc2,conda
  seqc2_testdata.py all --mode full --outdir /scratch/seqc2_full --profile singularity
"""

import argparse
import datetime
import gzip
import hashlib
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_DIR = SCRIPT_DIR.parent.parent
DEFAULT_CONFIG = SCRIPT_DIR / "seqc2_testdata.yaml"
DEFAULT_TESTDATA = PIPELINE_DIR.parent / "workflow-nxf-somatic-variantcomp-testdata" / "seqc2"
REQUIRED_TOOLS = ["curl", "samtools", "bcftools", "tabix", "bgzip", "bedtools"]
DEFAULT_PROFILES = {"small": "test_seqc2,conda", "full": "test_full,singularity"}

log = logging.getLogger("seqc2_testdata")


# --------------------------------------------------------------------------- helpers


def run(cmd, stdout=None, cwd=None):
    """Run a command (list), log it and fail loudly."""
    log.info("$ %s", " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True, stdout=stdout, cwd=cwd)


def md5sum(path, chunk=1 << 20):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url, dest, md5=None):
    """Download url to dest (resumable); verify md5 if given. Skips finished files."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        if md5 is None or md5sum(dest) == md5:
            log.info("cached: %s", dest.name)
            return dest
        log.warning("md5 mismatch for cached %s, downloading again", dest.name)
        dest.unlink()
    part = dest.with_name(dest.name + ".part")
    run(["curl", "-fsSL", "--retry", "5", "--retry-delay", "10", "-C", "-", "-o", part, url])
    if md5 is not None:
        observed = md5sum(part)
        if observed != md5:
            part.unlink()
            raise RuntimeError(f"md5 mismatch for {url}: expected {md5}, got {observed}")
        log.info("md5 ok: %s", dest.name)
    part.rename(dest)
    return dest


def fetch_md5_table(url, cache_file):
    """Read a SEQC2 md5 list; returns {basename: md5}.

    The list is cached; if the server is unreachable, the cached copy is used.
    Without either, files are downloaded/used without md5 verification.
    """
    cache_file = Path(cache_file)
    out = subprocess.run(["curl", "-fsSL", "--retry", "3", url], capture_output=True)
    if out.returncode == 0:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        cache_file.write_bytes(out.stdout)
        text = out.stdout.decode()
    elif cache_file.exists():
        log.warning("could not fetch %s, using cached copy %s", url, cache_file)
        text = cache_file.read_text()
    else:
        log.warning("could not fetch %s: files will NOT be md5-verified", url)
        return {}
    table = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2:
            table[Path(parts[1]).name] = parts[0]
    return table


def index_vcf(vcf):
    run(["tabix", "-f", "-p", "vcf", vcf])


def subset_vcf(src, dest, contigs):
    """Copy (contigs=None) or subset a VCF to contigs; output bgzipped + indexed.

    Streams with -t (targets) because some SEQC2 VCFs are plain gzip, not BGZF,
    and therefore cannot be indexed.
    """
    cmd = ["bcftools", "view", "--no-version", "-Oz", "-o", dest]
    if contigs:
        cmd += ["-t", ",".join(contigs)]
    run([*cmd, src])
    index_vcf(dest)


def check_tools():
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        sys.exit(
            f"Missing tools on PATH: {', '.join(missing)}. "
            "Activate the dev env: conda activate variantcomp-dev"
        )


def git_describe(path):
    try:
        out = subprocess.run(
            ["git", "-C", str(path), "describe", "--tags", "--always", "--dirty"],
            check=True,
            capture_output=True,
        )
        return out.stdout.decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


# --------------------------------------------------------------------------- prepare


class Dataset:
    def __init__(self, cfg, mode, outdir, cache, pairs=None, with_bams=False):
        self.cfg = cfg
        self.mode = mode
        self.mcfg = cfg["modes"][mode]
        self.base = cfg["base_url"].rstrip("/")
        self.out = Path(outdir).resolve()
        self.cache = Path(cache).resolve() if cache else self.out / ".downloads"
        self.pairs = pairs or self.mcfg["pairs"]
        self.contigs = self.mcfg.get("contigs")
        self.bam_pairs = self.mcfg.get("bam_pairs") or []
        if with_bams:
            self.bam_pairs = list(self.pairs)
        unknown = sorted(set(self.pairs) - set(cfg["available_pairs"]))
        if unknown:
            raise SystemExit(f"Unknown SEQC2 pairs: {', '.join(unknown)}")
        self.sources = []  # (relative output, url, md5 or n/a)

    def url(self, rel):
        return f"{self.base}/{rel}"

    def fetch(self, rel, md5=None):
        dest = download(self.url(rel), self.cache / rel, md5)
        self.sources.append((rel, self.url(rel), md5 or "n/a"))
        return dest

    # reference ---------------------------------------------------------------

    def prepare_reference(self):
        ref_dir = self.out / "reference"
        ref_dir.mkdir(parents=True, exist_ok=True)
        rel = self.cfg["reference"]["fasta"]
        name = self.cfg["reference"]["name"]
        if self.contigs:
            fasta = ref_dir / f"{name}_{'_'.join(self.contigs)}.fa"
            if not Path(str(fasta) + ".gz.fai").exists():
                with fasta.open("w") as fh:
                    for contig in self.contigs:
                        run(["samtools", "faidx", self.url(rel), contig], stdout=fh)
                run(["bgzip", "-f", "-@", "2", fasta])
                run(["samtools", "faidx", f"{fasta}.gz"])
            self.sources.append(
                (f"reference/{fasta.name}.gz", self.url(rel), "n/a (contig subset)")
            )
            # samtools caches the remote index in the working directory
            for leftover in Path.cwd().glob(Path(rel).name + ".fai"):
                leftover.unlink()
            return Path(f"{fasta}.gz")
        fasta = self.fetch(rel)
        self.fetch(rel + ".fai")
        dest = ref_dir / Path(rel).name
        for suffix in ("", ".fai"):
            target = Path(str(dest) + suffix)
            if not target.exists():
                os.link(Path(str(fasta) + suffix), target)
        return dest

    # regions -----------------------------------------------------------------

    def prepare_wes_bed(self):
        src = self.fetch(self.cfg["wes_bed"])
        bed_dir = self.out / "bed"
        bed_dir.mkdir(parents=True, exist_ok=True)
        tmp = bed_dir / "wes.unsorted.bed"
        with gzip.open(src, "rt") as fin, tmp.open("w") as fout:
            for line in fin:
                if line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if self.contigs and fields[0] not in self.contigs:
                    continue
                fout.write("\t".join(fields[:3]) + "\n")
        sorted_tmp = bed_dir / "wes.sorted.bed"
        with sorted_tmp.open("w") as fh:
            run(["bedtools", "sort", "-i", tmp], stdout=fh)
        bed = bed_dir / "S07604624_Covered_human_all_v6_plus_UTR_hg38.bed"
        with bed.open("w") as fh:
            run(["bedtools", "merge", "-i", sorted_tmp], stdout=fh)
        tmp.unlink()
        sorted_tmp.unlink()
        return bed

    # VCFs ----------------------------------------------------------------------

    def prepare_vcfs(self):
        md5 = fetch_md5_table(self.url(self.cfg["vcf_md5"]), self.cache / self.cfg["vcf_md5"])
        vcf_dir = self.out / "vcf"
        vcf_dir.mkdir(parents=True, exist_ok=True)
        result = {}
        for pair in self.pairs:
            for tech, tcfg in self.cfg["technologies"].items():
                rel = self.cfg["vcf_pattern"].format(tech=tcfg["folder"], pair=pair)
                src = self.fetch(rel, md5.get(Path(rel).name))
                dest = vcf_dir / f"{pair}_{tech}.mutect2.vcf.gz"
                if not Path(str(dest) + ".tbi").exists():
                    subset_vcf(src, dest, self.contigs)
                result[(pair, tech)] = dest
        return result

    def prepare_truth(self):
        truth_dir = self.out / "truth"
        truth_dir.mkdir(parents=True, exist_ok=True)
        files = {}
        for key in ("snv", "indel"):
            src = self.fetch(self.cfg["truth"][key])
            dest = truth_dir / Path(self.cfg["truth"][key]).name
            if not Path(str(dest) + ".tbi").exists():
                subset_vcf(src, dest, self.contigs)
            files[key] = dest
        src = self.fetch(self.cfg["truth"]["regions"])
        dest = truth_dir / Path(self.cfg["truth"]["regions"]).name
        with open(src) as fin, dest.open("w") as fout:
            for line in fin:
                if not self.contigs or line.split("\t", 1)[0] in self.contigs:
                    fout.write(line)
        files["regions"] = dest
        return files

    # BAMs ----------------------------------------------------------------------

    def bam_regions(self, wes_bed):
        """WES targets +/- padding inside the BAM windows (small mode)."""
        windows = self.mcfg.get("bam_windows") or []
        if not windows:
            return None
        pad = int(self.mcfg.get("bam_padding", 0))
        regions = self.out / "bam" / "bam_regions.bed"
        regions.parent.mkdir(parents=True, exist_ok=True)
        with open(wes_bed) as fh:
            targets = [line.split("\t")[:3] for line in fh if line.strip()]
        with regions.open("w") as fout:
            for window in windows:
                chrom, span = window.split(":")
                start, end = (int(x.replace(",", "")) for x in span.split("-"))
                last = None
                for c, s, e in targets:
                    s, e = max(int(s) - pad, start - 1), min(int(e) + pad, end)
                    if c != chrom or e <= s:
                        continue
                    if last and s <= last[1]:
                        last[1] = max(last[1], e)
                        continue
                    if last:
                        fout.write(f"{chrom}\t{last[0]}\t{last[1]}\n")
                    last = [s, e]
                if last:
                    fout.write(f"{chrom}\t{last[0]}\t{last[1]}\n")
        return regions

    def prepare_bams(self, wes_bed):
        if not self.bam_pairs:
            return {}
        regions = self.bam_regions(wes_bed)
        bam_dir = self.out / "bam"
        bam_dir.mkdir(parents=True, exist_ok=True)
        result = {}
        for pair in self.bam_pairs:
            center, rep = pair.split("_")
            for tech, tcfg in self.cfg["technologies"].items():
                for role in ("T", "N"):
                    keys = dict(tech=tcfg["folder"], center=center, role=role, rep=rep)
                    bam_rel = self.cfg["bam_pattern"].format(**keys)
                    bai_rel = self.cfg["bai_pattern"].format(**keys)
                    dest = bam_dir / f"{pair}_{tech}_{role}.bam"
                    if regions is None:
                        bam = self.fetch(bam_rel)
                        bai = self.fetch(bai_rel)
                        if not dest.exists():
                            os.link(bam, dest)
                            os.link(bai, str(dest) + ".bai")
                    else:
                        fraction = self.mcfg.get("bam_subsample")
                        subsample = (
                            ["--subsample", str(fraction), "--subsample-seed", "42"]
                            if fraction
                            else []
                        )
                        if not Path(str(dest) + ".bai").exists():
                            run(
                                ["samtools", "view", "-b", "-M", "-L", regions, *subsample]
                                + ["-o", dest, "-X", self.url(bam_rel), self.url(bai_rel)]
                            )
                            run(["samtools", "index", dest])
                        note = f"slice of {regions.name}"
                        if fraction:
                            note += f", subsampled to {fraction} (seed 42)"
                        self.sources.append(
                            (f"bam/{dest.name}", self.url(bam_rel), f"n/a ({note})")
                        )
                    result[(pair, tech, role)] = dest
        # samtools caches remote indices in the working directory
        for leftover in Path.cwd().glob("*.bwa.dedup.bai"):
            leftover.unlink()
        return result

    # outputs -------------------------------------------------------------------

    def write_samplesheet(self, vcfs, bams, wes_bed):
        header = [
            "sample",
            "technology",
            "library_type",
            "tumor_id",
            "normal_id",
            "caller",
            "vcf",
            "vcf_index",
            "tumor_bam",
            "tumor_bai",
            "normal_bam",
            "normal_bai",
            "target_bed",
        ]
        rows = []
        for pair in self.pairs:
            for tech, tcfg in self.cfg["technologies"].items():
                vcf = vcfs[(pair, tech)].relative_to(self.out)
                tbam = bams.get((pair, tech, "T"))
                nbam = bams.get((pair, tech, "N"))

                def rel(p, suffix=""):
                    return str(Path(str(p) + suffix).relative_to(self.out)) if p else "n/a"

                rows.append(
                    [
                        f"SEQC2_{pair}",
                        tech,
                        tcfg["library_type"],
                        self.cfg["tumor_id"],
                        self.cfg["normal_id"],
                        self.cfg["caller"],
                        str(vcf),
                        f"{vcf}.tbi",
                        rel(tbam),
                        rel(tbam, ".bai"),
                        rel(nbam),
                        rel(nbam, ".bai"),
                        str(wes_bed.relative_to(self.out)) if tech == "wes" else "n/a",
                    ]
                )
        path = self.out / "samplesheet.tsv"
        with path.open("w") as fh:
            fh.write("\t".join(header) + "\n")
            fh.writelines("\t".join(r) + "\n" for r in rows)
        return path

    def write_params(self, fasta, truth):
        params = {
            "input": "samplesheet.tsv",
            "truth_technology": "wgs",
            "fasta": str(fasta.relative_to(self.out)),
            "known_variants_snv": str(truth["snv"].relative_to(self.out)),
            "known_variants_indel": str(truth["indel"].relative_to(self.out)),
            "known_variants_regions": str(truth["regions"].relative_to(self.out)),
        }
        path = self.out / "params.yaml"
        with path.open("w") as fh:
            fh.write("# Pipeline parameters, paths relative to this folder.\n")
            fh.write("# 'seqc2_testdata.py run' resolves them; not all keys may be used by\n")
            fh.write("# the current pipeline version (known_variants_* are for a later feature).\n")
            yaml.safe_dump(params, fh, sort_keys=False)
        return path

    def write_provenance(self, args):
        files = sorted(
            p
            for p in self.out.rglob("*")
            if p.is_file()
            and ".downloads" not in p.parts
            and p.name not in ("checksums.md5", "README.md")
        )
        with (self.out / "checksums.md5").open("w") as fh:
            for p in files:
                fh.write(f"{md5sum(p)}  {p.relative_to(self.out)}\n")
        today = datetime.date.today().isoformat()
        lines = [
            f"# SEQC2 HCC1395/HCC1395BL WES + WGS test dataset ({self.mode} mode)",
            "",
            "Created by `tests/seqc2/seqc2_testdata.py` from",
            "core-unit-bioinformatics/workflow-nxf-somatic-variantcomp. Do not edit by hand;",
            "regenerate instead.",
            "",
            f"- Created: {today}",
            f"- Workflow repo version: `{git_describe(PIPELINE_DIR)}`",
            f"- Command: `seqc2_testdata.py {' '.join(args)}`",
            f"- Mode: `{self.mode}`; pairs: {', '.join(self.pairs)}; contigs: "
            f"{', '.join(self.contigs) if self.contigs else 'all'}; "
            f"BAMs: {', '.join(self.bam_pairs) if self.bam_pairs else 'none'}",
            f"- Reference: {self.cfg['reference']['name']}",
            "",
            "## Content",
            "",
            "| Path | Content |",
            "|---|---|",
            "| `samplesheet.tsv` | Pipeline samplesheet (paths relative to this folder) |",
            "| `params.yaml` | Pipeline parameters (paths relative to this folder) |",
            "| `reference/` | Reference FASTA (+ .fai/.gzi) |",
            "| `bed/` | WES target BED (Agilent SureSelect v6 + UTR, hg38 liftover), merged |",
            "| `vcf/` | Unfiltered MuTect2 (GATK 3.8) tumor/normal VCFs per pseudo-sample |",
            "| `bam/` | BAM slices (WES targets +/- padding in the BAM windows, subsampled) |",
            "| `truth/` | SEQC2 somatic reference call set v1.2.1 + high-confidence regions |",
            "| `checksums.md5` | md5 of every file in this folder |",
            "",
            "## Sources",
            "",
            "SEQC2 Somatic Mutation Working Group, SRA SRP162370. Please cite Fang et al. 2021,",
            "*Nat Biotechnol* 39:1151 (doi:10.1038/s41587-021-00993-6) and Xiao et al. 2021,",
            "*Nat Biotechnol* 39:1141 (doi:10.1038/s41587-021-00994-5).",
            "",
            "| File | URL | md5 (SEQC2) |",
            "|---|---|---|",
        ]
        lines += [f"| `{rel}` | {url} | {md5} |" for rel, url, md5 in self.sources]
        (self.out / "README.md").write_text("\n".join(lines) + "\n")

    def prepare(self, args):
        check_tools()
        self.out.mkdir(parents=True, exist_ok=True)
        log.info("Preparing SEQC2 %s dataset in %s", self.mode, self.out)
        fasta = self.prepare_reference()
        wes_bed = self.prepare_wes_bed()
        vcfs = self.prepare_vcfs()
        truth = self.prepare_truth()
        bams = self.prepare_bams(wes_bed)
        self.write_samplesheet(vcfs, bams, wes_bed)
        self.write_params(fasta, truth)
        self.write_provenance(args)
        log.info("Done: %s", self.out)


# --------------------------------------------------------------------------- run


def run_pipeline(dataset_dir, rundir, profile, extra_config, resume, nextflow_args):
    dataset_dir = Path(dataset_dir).resolve()
    rundir = Path(rundir).resolve()
    rundir.mkdir(parents=True, exist_ok=True)
    with (dataset_dir / "params.yaml").open() as fh:
        params = yaml.safe_load(fh)
    resolved = {}
    for key, value in params.items():
        candidate = dataset_dir / str(value)
        resolved[key] = str(candidate) if candidate.exists() else value
    resolved["outdir"] = str(rundir / "results")
    params_file = rundir / "params.yaml"
    with params_file.open("w") as fh:
        yaml.safe_dump(resolved, fh, sort_keys=False)
    cmd = ["nextflow", "run", PIPELINE_DIR, "-profile", profile, "-params-file", params_file]
    for config in extra_config or []:
        cmd += ["-c", Path(config).resolve()]
    if resume:
        cmd.append("-resume")
    cmd += nextflow_args or []
    run(cmd, cwd=rundir)


# --------------------------------------------------------------------------- CLI


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n".join(__doc__.split("\n")[1:]),
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def common(p):
        p.add_argument("--mode", choices=["small", "full"], default="small")
        p.add_argument(
            "--config",
            type=Path,
            default=DEFAULT_CONFIG,
            help="dataset configuration (default: %(default)s)",
        )
        p.add_argument(
            "--outdir", type=Path, help=f"dataset folder (default: {DEFAULT_TESTDATA}/<mode>)"
        )

    def prep(p):
        p.add_argument("--cache", type=Path, help="download cache (default: <outdir>/.downloads)")
        p.add_argument("--pairs", help="comma-separated SEQC2 pairs (default: from config)")
        p.add_argument("--with-bams", action="store_true", help="add (full) BAMs for all pairs")

    def runopts(p):
        p.add_argument("--rundir", type=Path, help="run folder (default: ./seqc2_<mode>_run)")
        p.add_argument(
            "--profile",
            help="Nextflow profile(s) (default: test_seqc2,conda for small; "
            "test_full,singularity for full)",
        )
        p.add_argument(
            "-c", "--nextflow-config", action="append", help="additional Nextflow config file(s)"
        )
        p.add_argument("--no-resume", action="store_true")
        p.add_argument(
            "nextflow_args",
            nargs=argparse.REMAINDER,
            help="further arguments passed to nextflow (after --)",
        )

    p_prep = sub.add_parser("prepare", help="download and prepare the dataset")
    common(p_prep)
    prep(p_prep)
    p_run = sub.add_parser("run", help="run the workflow on a prepared dataset")
    common(p_run)
    runopts(p_run)
    p_all = sub.add_parser("all", help="prepare + run")
    common(p_all)
    prep(p_all)
    runopts(p_all)

    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    with args.config.open() as fh:
        cfg = yaml.safe_load(fh)
    outdir = args.outdir or DEFAULT_TESTDATA / args.mode

    if args.command in ("prepare", "all"):
        pairs = args.pairs.split(",") if args.pairs else None
        Dataset(cfg, args.mode, outdir, args.cache, pairs, args.with_bams).prepare(argv)
    if args.command in ("run", "all"):
        extra = [a for a in (args.nextflow_args or []) if a != "--"]
        run_pipeline(
            outdir,
            args.rundir or Path(f"seqc2_{args.mode}_run"),
            args.profile or DEFAULT_PROFILES[args.mode],
            args.nextflow_config,
            not args.no_resume,
            extra,
        )


if __name__ == "__main__":
    main()
