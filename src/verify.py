"""Certify increasing even weights on the fixed P=32 code. Run from the project root.

python3 -B src/verify.py --through 26 --workers 6
Only a fully completed pair of sectors advances the lower bound. All scratch
matrices, pattern catalogues, compiled kernels and per-job inputs are temporary.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import json
from pathlib import Path
import subprocess
import tempfile
import time

from fetch_inputs import INPUTS, ROOT, UPSTREAM
from prepare import construct, syndrome, verify_and_write

INPUT_FILES = tuple(UPSTREAM / Path(name).name for name in INPUTS)


def run_weight(weight, directory, pool, workers, control):
    catalogue = directory / "catalogue"
    subprocess.run([directory / "patterns", str(weight), catalogue], check=True)
    with catalogue.open() as f:
        count = int(next(f))
        chunks = []
        chunk_size = max(256, (count + 4 * workers - 1) // (4 * workers))
        for start in range(0, count, chunk_size):
            size = min(count - start, chunk_size)
            p = directory / f"chunk{start}"
            with p.open("w") as out:
                print(size, file=out)
                for _ in range(size):
                    out.write(next(f))
            chunks.append((start, p, size))
        if f.read():
            raise ValueError("Unexpected trailing pattern data")

    def search(job):
        sector, (start, path, size) = job
        args = [directory / "embed", directory / f"{sector}.matrix", path, str(weight)]
        if not control:
            args.append("--minimal")
        record = json.loads(subprocess.check_output(args))
        if not record["found"] and record["completed_patterns"] != size:
            raise ValueError("Partial shard")
        record["start"] = start
        return sector, record

    sectors = {
        s: {"completed_patterns": 0, "embeddings_tested": 0, "witnesses": []}
        for s in "XZ"
    }
    jobs = [(s, chunk) for chunk in chunks for s in "XZ"]
    for sector, rec in pool.map(search, jobs):
        summary = sectors[sector]
        for field in ("completed_patterns", "embeddings_tested"):
            summary[field] += rec[field]
        if rec["found"]:
            summary["witnesses"].append(rec["support"])
    for _, path, _ in chunks:
        path.unlink()
    return {
        "weight": weight,
        "catalogue_size": count,
        "minimality_pruning": not control,
        "sectors": sectors,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--through", type=int, required=True, choices=range(2, 31, 2))
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--output", type=Path, help="Path for the new JSONL run record")
    args = p.parse_args()
    if args.workers < 1:
        p.error("--workers must be positive")
    stem = f"affine_p32_{datetime.now():%Y%m%d_%H%M%S_%f}"
    output = args.output or ROOT / "data" / "runs" / (stem + ".jsonl")
    output.parent.mkdir(parents=True, exist_ok=True)
    if not all(path.is_file() for path in INPUT_FILES):
        p.error("run python3 -B src/fetch_inputs.py first")
    lower, upper = 2, 32
    with tempfile.TemporaryDirectory() as tmp, output.open("x", buffering=1) as out:
        directory = Path(tmp)
        tick = time.monotonic()
        hashes = verify_and_write(UPSTREAM, directory)
        for name in ("patterns", "embed"):
            subprocess.run(
                [
                    "c++",
                    "-O3",
                    "-std=c++17",
                    str(ROOT / "src" / f"{name}.cpp"),
                    "-o",
                    str(directory / name),
                ],
                check=True,
            )
        print(json.dumps({"matrix_sha256": hashes}), file=out)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for weight in range(2, args.through + 1, 2):
                r = run_weight(
                    weight, directory, pool, args.workers, control=weight <= 22
                )
                found = any(s["witnesses"] for s in r["sectors"].values())
                if found:
                    # Freshly verify each claimed quantum witness outside the C++ kernel.
                    rows, columns, _ = construct(UPSTREAM / "affine.csv")
                    for sector, v in r["sectors"].items():
                        basis = {}
                        for qs in rows[sector]:
                            b = sum(1 << q for q in qs)
                            while b and b.bit_length() in basis:
                                b ^= basis[b.bit_length()]
                            if b:
                                basis[b.bit_length()] = b
                        for qs in v["witnesses"]:
                            other = "Z" if sector == "X" else "X"
                            if (
                                len(set(qs)) != weight
                                or len(qs) != weight
                                or syndrome(columns[other], qs)
                            ):
                                raise ValueError("Invalid new support")
                            b = sum(1 << q for q in qs)
                            while b and b.bit_length() in basis:
                                b ^= basis[b.bit_length()]
                            if not b:
                                raise ValueError("New witness is a stabilizer")
                    upper = weight
                else:
                    lower = weight + 2
                r.update(
                    lower=lower,
                    upper=upper,
                    elapsed_seconds=round(time.monotonic() - tick, 3),
                )
                if lower > upper:
                    raise ValueError("Inconsistent bounds")
                print(json.dumps(r), file=out)
                print(
                    weight,
                    r["catalogue_size"],
                    lower,
                    upper,
                    r["elapsed_seconds"],
                    flush=True,
                )
                if lower == upper:
                    break
    print(f"Recorded interval [{lower}, {upper}] in {output}.")


if __name__ == "__main__":
    main()
