"""Unit tests for tests/seqc2/seqc2_testdata.py (offline parts only)."""

import importlib.util
from pathlib import Path

import yaml

SCRIPT = Path(__file__).resolve().parents[1] / "seqc2" / "seqc2_testdata.py"
spec = importlib.util.spec_from_file_location("seqc2_testdata", SCRIPT)
seqc2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(seqc2)


def load_config():
    with seqc2.DEFAULT_CONFIG.open() as fh:
        return yaml.safe_load(fh)


def test_config_modes_use_available_pairs():
    cfg = load_config()
    for mode in cfg["modes"].values():
        assert set(mode["pairs"]) <= set(cfg["available_pairs"])
        assert set(mode.get("bam_pairs") or []) <= set(mode["pairs"])


def test_bam_regions_pad_clip_and_merge(tmp_path):
    cfg = load_config()
    cfg["modes"]["small"]["bam_window"] = "chr22:1000-2000"
    cfg["modes"]["small"]["bam_padding"] = 50
    wes_bed = tmp_path / "wes.bed"
    wes_bed.write_text(
        "chr21\t1100\t1200\n"  # other contig: dropped
        "chr22\t900\t1010\n"  # clipped to window start
        "chr22\t1100\t1200\n"  # padded, overlaps the next one after padding
        "chr22\t1280\t1300\n"
        "chr22\t1990\t2100\n"  # clipped to window end
        "chr22\t3000\t3100\n"  # outside window: dropped
    )
    ds = seqc2.Dataset(cfg, "small", tmp_path / "out", None)
    regions = ds.bam_regions(wes_bed).read_text().splitlines()
    assert regions == ["chr22\t999\t1350", "chr22\t1940\t2000"]


def test_unknown_pairs_are_rejected(tmp_path):
    cfg = load_config()
    try:
        seqc2.Dataset(cfg, "small", tmp_path, None, pairs=["XX_9"])
    except SystemExit as err:
        assert "XX_9" in str(err)
    else:
        raise AssertionError("unknown pair accepted")
