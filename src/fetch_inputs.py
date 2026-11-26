"""Fetch the versioned arXiv source and extract the two inputs for this code."""

import argparse
import gzip
from hashlib import sha256
from io import BytesIO
from pathlib import Path
import shutil
import tarfile
import tempfile
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "data" / "upstream"
ARCHIVE = UPSTREAM / "2604.20838v4-source.tar"
URL = "https://arxiv.org/src/2604.20838v4"
ARCHIVE_SHA256 = "d638e1fd8a4ff1b0f34c3ab8bff16fea3d3b79696523ffcca00cdcc822c677b7"
PREFIX = "anc/distance_certificate"
INPUTS = ("data/affine.csv", "results/affine_certificates.json")


def read_members(source_bytes, names):
    """Read named files without unpacking the rest of the upstream archive."""
    files = {}
    with tarfile.open(fileobj=BytesIO(source_bytes), mode="r:") as source:
        for name in names:
            member = source.getmember(f"{PREFIX}/{name}")
            if not member.isfile():
                raise ValueError(f"Expected a regular file: {name}")
            with source.extractfile(member) as input_file:
                files[name] = input_file.read()
    return files


def install_inputs(archive):
    archive_bytes = archive.read_bytes()
    source_bytes = archive_bytes
    # A saved copy of the arXiv response may itself have been gzip-compressed.
    for _ in range(2):
        if not source_bytes.startswith(b"\x1f\x8b"):
            break
        source_bytes = gzip.decompress(source_bytes)
    if sha256(source_bytes).hexdigest() != ARCHIVE_SHA256:
        raise ValueError("Unexpected SHA-256 for the uncompressed arXiv v4 source")
    for name, contents in read_members(source_bytes, INPUTS).items():
        (UPSTREAM / Path(name).name).write_bytes(contents)
    if archive.resolve() != ARCHIVE.resolve() or source_bytes != archive_bytes:
        ARCHIVE.write_bytes(source_bytes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", nargs="?", type=Path, help="Local arXiv source archive")
    args = parser.parse_args()
    UPSTREAM.mkdir(parents=True, exist_ok=True)
    if args.archive is not None or ARCHIVE.exists():
        install_inputs(args.archive or ARCHIVE)
    else:
        with tempfile.TemporaryDirectory() as scratch:
            downloaded = Path(scratch) / ARCHIVE.name
            with urlopen(URL, timeout=60) as response, downloaded.open("wb") as output:
                shutil.copyfileobj(response, output)
            install_inputs(downloaded)
    print(f"Verified both arXiv:2604.20838v4 inputs in {UPSTREAM.relative_to(ROOT)}.")


if __name__ == "__main__":
    main()
