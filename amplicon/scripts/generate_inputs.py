#!/usr/bin/env python3
"""
Generate QIIME2 manifest and metadata files from input/16S and input/ITS directories.
Metadata values are looked up from lib/metadata.tsv; samples not found there receive
placeholder values.

Outputs:
  manifests/16S_manifest.tsv
  manifests/ITS_manifest.tsv
  metadata/16S_metadata.tsv
  metadata/ITS_metadata.tsv
"""

import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_ROOT   = PROJECT_ROOT / "input"
LIB_METADATA = PROJECT_ROOT / "lib" / "metadata.tsv"
MANIFEST_DIR = PROJECT_ROOT / "manifests"
METADATA_DIR = PROJECT_ROOT / "metadata"

MANIFEST_DIR.mkdir(exist_ok=True)
METADATA_DIR.mkdir(exist_ok=True)

METADATA_COLUMNS = ["sample-id", "Country", "host", "material", "description"]
QIIME2_TYPES = {
    "sample-id":   "sample-id",
    "Country":     "categorical",
    "host":        "categorical",
    "material":    "categorical",
    "description": "categorical",
}


def load_lib_metadata(path: Path) -> dict[str, dict]:
    """Read lib/metadata.tsv and return a dict keyed by sample-id."""
    if not path.exists():
        print(f"[WARNING] lib metadata not found: {path}", file=sys.stderr)
        return {}
    lookup: dict[str, dict] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            row = dict(zip(header, fields))
            sid = row.get("sample-id", "").strip()
            if sid:
                lookup[sid] = row
    return lookup


def parse_samples(amplicon: str) -> list[dict]:
    """Scan input/{amplicon}/ and pair R1/R2 files; return sample records."""
    indir = INPUT_ROOT / amplicon
    if not indir.exists():
        print(f"[WARNING] Directory not found: {indir}", file=sys.stderr)
        return []

    samples = []
    for r1 in sorted(indir.glob("*_R1_001.fastq.gz")):
        r2 = Path(str(r1).replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))
        if not r2.exists():
            print(f"[WARNING] Missing R2 for {r1.name}", file=sys.stderr)
            continue
        stem = r1.name.replace("_R1_001.fastq.gz", "")
        sample_id = re.sub(rf"_{amplicon}_S\d+$", "", stem)
        samples.append({
            "sample_id": sample_id,
            "forward":   str(r1.resolve()),
            "reverse":   str(r2.resolve()),
        })
    return samples


def write_manifest(samples: list[dict], path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for s in samples:
            fh.write(f"{s['sample_id']}\t{s['forward']}\t{s['reverse']}\n")
    print(f"Wrote manifest: {path}  ({len(samples)} samples)")


def write_metadata(samples: list[dict], path: Path,
                   lib_meta: dict[str, dict]) -> None:
    rows = []
    for s in samples:
        sid = s["sample_id"]
        if sid in lib_meta:
            src = lib_meta[sid]
            row = {
                "sample-id":   sid,
                "Country":     src.get("Country",     "NA").strip(),
                "host":        src.get("host",        "NA").strip(),
                "material":    src.get("material",    "NA").strip(),
                "description": src.get("description", sid).strip(),
            }
        else:
            print(f"[WARNING] {sid} not in lib/metadata.tsv — using placeholders",
                  file=sys.stderr)
            row = {
                "sample-id":   sid,
                "Country":     "NA",
                "host":        "NA",
                "material":    "NA",
                "description": sid,
            }
        rows.append(row)

    with open(path, "w") as fh:
        fh.write("\t".join(METADATA_COLUMNS) + "\n")
        fh.write("\t".join(QIIME2_TYPES[c] for c in METADATA_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(c, "NA") for c in METADATA_COLUMNS) + "\n")
    print(f"Wrote metadata: {path}  ({len(rows)} samples)")


def main():
    lib_meta = load_lib_metadata(LIB_METADATA)
    if not lib_meta:
        print("[ERROR] No entries loaded from lib/metadata.tsv", file=sys.stderr)
        sys.exit(1)

    for amplicon in ("16S", "ITS"):
        samples = parse_samples(amplicon)
        if not samples:
            print(f"[ERROR] No samples found for {amplicon}", file=sys.stderr)
            continue
        write_manifest(samples, MANIFEST_DIR / f"{amplicon}_manifest.tsv")
        write_metadata(samples, METADATA_DIR / f"{amplicon}_metadata.tsv", lib_meta)


if __name__ == "__main__":
    main()
