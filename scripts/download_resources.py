#!/usr/bin/env python3
"""Download the annotation resources of the pipeline before an offline run.

The pipeline never downloads anything (it must run offline on an HPC). Run this script
once on a machine with internet access, copy the output folder to the HPC and pass the
written config to the pipeline (`-c <outdir>/resources.config`).

Resources (each can be switched off)
  vep_cache/   Ensembl VEP cache (indexed, offline) of the species and build, e.g.
               vep_cache/homo_sapiens/116_GRCh38/ -> --vep_cache, --vep_cache_version,
               --vep_genome (and --merged / --refseq in --vep_cache_args for those caches)
               Source: https://ftp.ensembl.org/pub/release-<v>/variation/indexed_vep_cache/
               Checked with the Ensembl CHECKSUMS file (BSD `sum`), if `sum` is on PATH.
               GRCh38 cache v116: about 27 GB download, about 30 GB unpacked.
  giab/<build>/ GIAB genome stratifications v3.5 of the build, in the release layout
               (LowComplexity/, Mappability/, ...) -> --giab_dir, --giab_files
               Source: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/
               genome-stratifications/v3.5/<build>@all/; md5-checked.
  gencode/     (--gtf) GENCODE basic gene annotation, sorted, bgzipped + tabix index,
               for VEP without a cache -> --vep_gtf (needs bgzip and tabix on PATH)

OUTPUT (in --outdir)
  resources.config        Nextflow config with a params block: the absolute paths of the
                          downloaded resources (vep_cache or vep_gtf, vep_cache_version,
                          vep_genome, vep_species, giab_dir, giab_files); use it with
                          `-c resources.config` next to your own `-params-file`
  resources.params.yaml   the same parameters as YAML (for -params-file, or to copy into
                          your own params file)
  README.md               what was downloaded, when, from where, with which checksum

Downloads are resumable (an interrupted run continues where it stopped) and finished
files are not downloaded again.

TLS certificates: the first existing file of --ca-bundle, $CURL_CA_BUNDLE, $SSL_CERT_FILE,
$REQUESTS_CA_BUNDLE, the Python package certifi and the usual system paths is used. This
avoids errors with conda builds whose built-in CA path was not relocated.

Examples
  scripts/download_resources.py --outdir /data/variantcomp_resources --genome GRCh38
  scripts/download_resources.py --outdir res --genome GRCh37 --no-vep --gtf
  scripts/download_resources.py --outdir res --dry-run
"""

import argparse
import datetime
import gzip
import hashlib
import logging
import os
import shutil
import ssl
import subprocess
import sys
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path

ENSEMBL_FTP = "https://ftp.ensembl.org/pub/release-{version}/variation/indexed_vep_cache"
GIAB_FTP = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "genome-stratifications/v3.5/{build}@all"
)
GENCODE_FTP = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release}"
# GIAB stratifications used by default (pipeline --giab_files): name -> path in the release
GIAB_FILES = {
    "homopolymer": "LowComplexity/{build}_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz",
    "tandem_repeat": "LowComplexity/{build}_AllTandemRepeats.bed.gz",
    "segmental_duplication": "SegmentalDuplications/{build}_segdups.bed.gz",
    "low_mappability": "Mappability/{build}_lowmappabilityall.bed.gz",
    "all_difficult": "Union/{build}_alldifficultregions.bed.gz",
}
# CA bundles of common Linux distributions (Debian/Ubuntu, RHEL/Fedora, openSUSE, Alpine)
SYSTEM_CA_BUNDLES = [
    "/etc/ssl/certs/ca-certificates.crt",
    "/etc/pki/tls/certs/ca-bundle.crt",
    "/etc/ssl/ca-bundle.pem",
    "/etc/ssl/cert.pem",
]

log = logging.getLogger("download_resources")
SSL_CONTEXT = None  # set in main()


# --------------------------------------------------------------------------- downloads


def find_ca_bundle(explicit=None):
    """Path of the first existing CA bundle: --ca-bundle, env vars, certifi, system paths."""
    candidates = [explicit] if explicit else []
    # environment variables only count if they point to a file
    candidates += [os.environ.get(v) for v in ("CURL_CA_BUNDLE", "SSL_CERT_FILE")]
    candidates.append(os.environ.get("REQUESTS_CA_BUNDLE"))
    try:
        import certifi

        candidates.append(certifi.where())
    except ImportError:
        pass
    candidates += SYSTEM_CA_BUNDLES
    for path in candidates:
        if path and Path(path).is_file():
            return str(path)
    if explicit:
        sys.exit(f"ERROR: --ca-bundle {explicit} does not exist")
    return None


def ssl_context(explicit=None):
    """SSL context with the CA bundle found by find_ca_bundle (Python defaults if none)."""
    bundle = find_ca_bundle(explicit)
    if bundle is None:
        log.warning("no CA bundle found (set --ca-bundle or install certifi); using defaults")
        return ssl.create_default_context()
    log.info("CA bundle: %s", bundle)
    return ssl.create_default_context(cafile=bundle)


def http_get(url, part, retries=5, delay=10):
    """Download url into the file part; resumes a partial file (HTTP Range); retries."""
    for attempt in range(1, retries + 1):
        # bytes already downloaded by an earlier, interrupted attempt
        offset = part.stat().st_size if part.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, context=SSL_CONTEXT, timeout=300) as resp:
                # 206: the server continues at offset; 200: it sends the whole file
                mode = "ab" if offset and resp.status == 206 else "wb"
                with open(part, mode) as fh:
                    shutil.copyfileobj(resp, fh, 1 << 20)
            return
        except urllib.error.HTTPError as err:
            # 416: nothing left to download (the partial file is complete)
            if err.code == 416 and offset:
                return
            # other client errors (404, 403, ...) won't go away by retrying
            if 400 <= err.code < 500:
                sys.exit(f"ERROR: download of {url} failed: HTTP {err.code}")
            reason = f"HTTP {err.code}"
        except (urllib.error.URLError, TimeoutError, ConnectionError, ssl.SSLError) as err:
            reason = str(getattr(err, "reason", err))
        log.warning("download of %s failed (%s); attempt %d/%d", url, reason, attempt, retries)
        if attempt < retries:
            time.sleep(delay)
    sys.exit(f"ERROR: download of {url} failed after {retries} attempts")


def http_text(url):
    """Content of a small text file at url (stops on error: checksums are required)."""
    try:
        with urllib.request.urlopen(url, context=SSL_CONTEXT, timeout=120) as resp:
            return resp.read().decode()
    except (urllib.error.URLError, TimeoutError, ConnectionError, ssl.SSLError) as err:
        sys.exit(f"ERROR: could not fetch {url}: {getattr(err, 'reason', err)}")


def http_size(url):
    """Size in bytes of the file at url (HEAD request), or None if unknown."""
    try:
        request = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(request, context=SSL_CONTEXT, timeout=60) as resp:
            length = resp.headers.get("Content-Length")
            return int(length) if length else None
    except (urllib.error.URLError, TimeoutError, ConnectionError, ssl.SSLError):
        return None


def md5sum(path, chunk=1 << 20):
    """md5 hex digest of a file."""
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url, dest, check=None):
    """Download url to dest (resumable) unless dest exists; check(path) verifies the file.

    check returns an error message or None; a failed check deletes the file and stops.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        log.info("exists: %s", dest)
        return dest
    part = dest.with_name(dest.name + ".part")
    log.info("download: %s -> %s", url, dest)
    http_get(url, part)
    problem = check(part) if check else None
    if problem:
        part.unlink()
        sys.exit(f"ERROR: {url}: {problem}; the file was deleted, run the script again")
    part.rename(dest)
    return dest


# --------------------------------------------------------------------------- resources


def bsd_sum(path):
    """BSD checksum and 1-kB block count of a file (`sum`, as in Ensembl CHECKSUMS)."""
    if shutil.which("sum") is None:
        return None
    out = subprocess.run(["sum", str(path)], check=True, capture_output=True, text=True)
    checksum, blocks = out.stdout.split()[:2]
    return int(checksum), int(blocks)


def safe_extract(archive, dest):
    """Extract a tar.gz below dest; members with absolute paths or '..' are refused."""
    dest = Path(dest).resolve()
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            target = (dest / member.name).resolve()
            if dest not in target.parents and target != dest:
                sys.exit(f"ERROR: {archive}: member {member.name} would be written outside {dest}")
        # the 'data' filter (Python >= 3.12, backported) also drops device files and links
        # out of the tree
        if hasattr(tarfile, "data_filter"):
            tar.extractall(dest, filter="data")
        else:
            tar.extractall(dest)


def vep_cache(args, out, sources):
    """Download and unpack the indexed VEP cache; return the params for the pipeline."""
    species = (
        args.species
        + {"ensembl": "", "merged": "_merged", "refseq": "_refseq"}[args.vep_cache_type]
    )
    name = f"{species}_vep_{args.vep_cache_version}_{args.genome}.tar.gz"
    base = ENSEMBL_FTP.format(version=args.vep_cache_version)
    url = f"{base}/{name}"
    cache_dir = out / "vep_cache"
    unpacked = cache_dir / species / f"{args.vep_cache_version}_{args.genome}"
    if args.dry_run:
        size = http_size(url)
        log.info("would download %s (%s GB)", url, f"{size / 1e9:.1f}" if size else "?")
    elif unpacked.is_dir():
        log.info("exists: %s", unpacked)
    else:
        # expected BSD checksum and block count from the Ensembl CHECKSUMS file
        expected = None
        for line in http_text(f"{base}/CHECKSUMS").splitlines():
            fields = line.split()
            if len(fields) == 3 and fields[2] == name:
                expected = (int(fields[0]), int(fields[1]))
        if expected is None:
            sys.exit(f"ERROR: {name} is not listed in {base}/CHECKSUMS; check --vep-cache-version")

        def check(path):
            observed = bsd_sum(path)
            if observed is None:
                log.warning("`sum` not on PATH: %s is not checksum-verified", name)
                return None
            return None if observed == expected else f"checksum {observed} != {expected}"

        archive = download(url, out / "downloads" / name, check)
        log.info("unpacking %s into %s", archive.name, cache_dir)
        safe_extract(archive, cache_dir)
        if not unpacked.is_dir():
            sys.exit(f"ERROR: {archive}: expected folder {unpacked} after unpacking")
        if not args.keep_archives:
            archive.unlink()
    sources.append(("vep_cache/", url, "Ensembl CHECKSUMS (BSD sum)"))
    params = {
        "vep_cache": str(cache_dir),
        "vep_cache_version": args.vep_cache_version,
        "vep_genome": args.genome,
        "vep_species": args.species,
    }
    # merged / RefSeq caches need the matching VEP flag
    if args.vep_cache_type != "ensembl":
        params["vep_cache_args"] = (
            f"--check_existing --var_synonyms --clin_sig_allele 1 --{args.vep_cache_type}"
        )
    return params


def giab(args, out, sources):
    """Download the GIAB stratifications (md5-checked); return the params for the pipeline."""
    base = GIAB_FTP.format(build=args.genome)
    files = {name: rel.format(build=args.genome) for name, rel in GIAB_FILES.items()}
    if args.giab_names:
        unknown = sorted(set(args.giab_names) - set(files))
        if unknown:
            sys.exit(f"ERROR: unknown GIAB name(s) {', '.join(unknown)}; known: {', '.join(files)}")
        files = {name: files[name] for name in args.giab_names}
    # local folder named after the build (no '@' in file names)
    giab_dir = out / "giab" / args.genome
    md5 = {}
    if not args.dry_run:
        for line in http_text(f"{base}/{args.genome}-genome-stratifications-md5s.txt").splitlines():
            fields = line.split()
            if len(fields) == 2:
                md5[fields[1].lstrip("./")] = fields[0]
    for rel in files.values():
        url = f"{base}/{rel}"
        if args.dry_run:
            log.info("would download %s", url)
            continue
        if rel not in md5:
            sys.exit(f"ERROR: {rel} is not listed in the GIAB md5 file of {args.genome}")

        def check(path, expected=md5[rel]):
            observed = md5sum(path)
            return None if observed == expected else f"md5 {observed} != {expected}"

        download(url, giab_dir / rel, check)
        sources.append((f"giab/{args.genome}/{rel}", url, md5[rel]))
    return {
        "giab_dir": str(giab_dir),
        "giab_files": ",".join(f"{name}={rel}" for name, rel in files.items()),
    }


def gencode_gtf(args, out, sources):
    """GENCODE basic GTF, sorted, bgzipped + tabix (VEP --gtf); return the params."""
    for tool in ("bgzip", "tabix"):
        if shutil.which(tool) is None:
            sys.exit(f"ERROR: --gtf needs {tool} on PATH (e.g. conda activate variantcomp-dev)")
    base = GENCODE_FTP.format(release=args.gencode_release)
    if args.genome == "GRCh38":
        name = f"gencode.v{args.gencode_release}.basic.annotation.gtf.gz"
        url = f"{base}/{name}"
    else:
        name = f"gencode.v{args.gencode_release}lift37.basic.annotation.gtf.gz"
        url = f"{base}/GRCh37_mapping/{name}"
    gtf_dir = out / "gencode"
    final = gtf_dir / name.replace(".gtf.gz", ".sorted.gtf.gz")
    if args.dry_run:
        log.info("would download %s", url)
    elif not Path(str(final) + ".tbi").exists():
        raw = download(url, out / "downloads" / name)
        # VEP --gtf needs a position-sorted, bgzipped, tabix-indexed GTF without headers
        records = []
        with gzip.open(raw, "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    fields = line.split("\t", 5)
                    records.append((fields[0], int(fields[3]), int(fields[4]), line))
        records.sort(key=lambda r: (r[0], r[1], r[2]))
        plain = final.with_suffix("")
        gtf_dir.mkdir(parents=True, exist_ok=True)
        with open(plain, "w") as fh:
            fh.writelines(r[3] for r in records)
        subprocess.run(["bgzip", "-f", str(plain)], check=True)
        subprocess.run(["tabix", "-f", "-p", "gff", str(final)], check=True)
        if not args.keep_archives:
            raw.unlink()
    sources.append((f"gencode/{final.name}", url, "n/a (not published)"))
    return {"vep_gtf": str(final)}


def write_outputs(args, out, params, sources):
    """resources.config / resources.params.yaml (pipeline parameters), README.md."""
    # the same parameters as YAML (-params-file) and as Nextflow config (-c)
    yaml_lines = [
        "# Annotation resources for workflow-nxf-somatic-variantcomp, written by",
        "# scripts/download_resources.py. Use with: -params-file resources.params.yaml",
    ]
    config_lines = [
        "// Annotation resources for workflow-nxf-somatic-variantcomp, written by",
        "// scripts/download_resources.py. Use with: -c resources.config",
        "params {",
    ]
    for key, value in params.items():
        if isinstance(value, str):
            yaml_lines.append(f"{key}: {value!r}")
            config_lines.append(f"    {key} = '{value}'")
        else:
            text = str(value).lower() if isinstance(value, bool) else str(value)
            yaml_lines.append(f"{key}: {text}")
            config_lines.append(f"    {key} = {text}")
    config_lines.append("}")
    (out / "resources.params.yaml").write_text("\n".join(yaml_lines) + "\n")
    (out / "resources.config").write_text("\n".join(config_lines) + "\n")
    readme = [
        "# Annotation resources",
        "",
        f"- Created: {datetime.date.today().isoformat()}",
        f"- Command: `download_resources.py {' '.join(sys.argv[1:])}`",
        f"- Genome build: {args.genome}",
        "",
        "| Path | URL | Checksum |",
        "|---|---|---|",
    ]
    readme += [f"| `{rel}` | {url} | {checksum} |" for rel, url, checksum in sources]
    readme += [
        "",
        "Cite: Ensembl VEP (McLaren et al. 2016, Genome Biol 17:122); GIAB genome",
        "stratifications (Dwarshuis et al. 2024, Nat Commun 15:9029); GENCODE (Frankish et",
        "al. 2023, Nucleic Acids Res 51:D942).",
    ]
    (out / "README.md").write_text("\n".join(readme) + "\n")
    log.info("wrote resources.config, resources.params.yaml and README.md in %s", out)


def main():
    """Parse the command line and download the selected resources."""
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n".join(__doc__.split("\n")[1:]),
    )
    parser.add_argument("--outdir", type=Path, required=True, help="resource folder")
    parser.add_argument("--genome", choices=["GRCh38", "GRCh37"], default="GRCh38")
    parser.add_argument("--species", default="homo_sapiens", help="VEP species (cache name)")
    parser.add_argument("--vep-cache-version", type=int, default=116, help="Ensembl release")
    parser.add_argument(
        "--vep-cache-type", choices=["ensembl", "merged", "refseq"], default="ensembl"
    )
    parser.add_argument("--no-vep", action="store_true", help="don't download the VEP cache")
    parser.add_argument("--no-giab", action="store_true", help="don't download GIAB files")
    parser.add_argument(
        "--giab-names",
        type=lambda s: s.split(","),
        help=f"comma-separated subset of the GIAB files ({', '.join(GIAB_FILES)})",
    )
    parser.add_argument("--gtf", action="store_true", help="GENCODE GTF for --vep_gtf")
    parser.add_argument("--gencode-release", type=int, default=46)
    parser.add_argument("--keep-archives", action="store_true", help="keep downloaded archives")
    parser.add_argument("--ca-bundle", help="CA certificate bundle for HTTPS")
    parser.add_argument("--dry-run", action="store_true", help="only list what would be done")
    args = parser.parse_args()
    if args.no_vep and args.no_giab and not args.gtf:
        parser.error("nothing to do: --no-vep and --no-giab without --gtf")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    global SSL_CONTEXT
    SSL_CONTEXT = ssl_context(args.ca_bundle)
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    params, sources = {}, []
    if not args.no_vep:
        params.update(vep_cache(args, out, sources))
    if args.gtf:
        gtf_params = gencode_gtf(args, out, sources)
        # the pipeline takes either a VEP cache or a GTF: the GTF is only a parameter
        # when no cache was downloaded
        if args.no_vep:
            params.update(gtf_params)
            params["vep_genome"] = args.genome
    if not args.no_giab:
        params.update(giab(args, out, sources))
    else:
        params["giab_stratifications"] = False
    if args.dry_run:
        log.info("dry run: nothing written")
        return
    write_outputs(args, out, params, sources)


if __name__ == "__main__":
    main()
