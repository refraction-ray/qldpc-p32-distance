"""Cross-check low-weight catalogues/counters and test known logical patterns at weight 32."""

import argparse
import hashlib
from collections import defaultdict
import json
from pathlib import Path
import subprocess
import tempfile

from fetch_inputs import ARCHIVE, ARCHIVE_SHA256, ROOT, UPSTREAM, read_members
from prepare import construct, verify_and_write


def positive_pattern(support, columns):
    """Pair incidences at each check, then extract the connected component(s)."""
    rows = defaultdict(list)
    for i, q in enumerate(support):
        for c, r in enumerate(columns[q]):
            rows[(c, r)].append(i)
    adj = [[-1, -1, -1] for _ in support]
    for (c, _), indices in rows.items():
        if len(indices) % 2:
            raise ValueError("Odd syndrome in control support")
        for a, b in zip(indices[::2], indices[1::2]):
            adj[a][c] = b
            adj[b][c] = a
    pending = set(range(len(adj)))
    components = []
    while pending:
        indices = [min(pending)]
        pending.remove(indices[0])
        for a in indices:
            for b in adj[a]:
                if b in pending:
                    indices.append(b)
                    pending.remove(b)
        ids = {v: i for i, v in enumerate(indices)}
        components.append([[ids[b] for b in adj[a]] for a in indices])
    return components


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--published",
        action="store_true",
        help="Compare through weight 22 using the downloaded arXiv source",
    )
    args = p.parse_args()
    with tempfile.TemporaryDirectory() as scratch:
        d = Path(scratch)
        for exe in ["patterns", "embed"]:
            subprocess.run(
                ["c++", "-O3", "-std=c++17", ROOT / f"src/{exe}.cpp", "-o", d / exe],
                check=True,
            )
        ref = UPSTREAM
        _, cols, _ = construct(ref / "affine.csv")
        verify_and_write(ref, d)
        cert = json.loads((ref / "affine_certificates.json").read_text())
        for s, other in (("X", "Z"), ("Z", "X")):
            components = positive_pattern(cert["sides"][s]["support"], cols[other])
            for minimal in (False, True):
                matched = False
                for graph in components:
                    w = len(graph)
                    line = [w] + [b for neighbours in graph for b in neighbours]
                    (d / "positive").write_text("1\n" + " ".join(map(str, line)) + "\n")
                    cmd = [d / "embed", d / f"{s}.matrix", d / "positive", str(w)]
                    if minimal:
                        cmd.append("--minimal")
                    witness = json.loads(subprocess.check_output(cmd))
                    if witness["found"]:
                        matched = True
                        print(
                            s,
                            "positive control found at weight",
                            w,
                            "minimal pruning",
                            minimal,
                            flush=True,
                        )
                if not matched:
                    raise ValueError("Positive control failed")
        if not args.published:
            return
        archive_bytes = ARCHIVE.read_bytes()
        if hashlib.sha256(archive_bytes).hexdigest() != ARCHIVE_SHA256:
            raise ValueError("The downloaded arXiv archive has changed")
        src = d / "published"
        names = ["src/colored_exact.cpp"]
        names.extend(
            f"results/affine_{sector}/{name}"
            for sector in "XZ"
            for name in ("matrix.txt", "filter.txt", "matrix.txt.filtered.results.jsonl")
        )
        for name, data in read_members(archive_bytes, names).items():
            path = src / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        subprocess.run(
            [
                "c++",
                "-O3",
                "-std=c++17",
                src / "src/colored_exact.cpp",
                "-o",
                d / "author",
            ],
            check=True,
        )
        subprocess.run(
            [
                "c++",
                "-O3",
                "-std=c++17",
                ROOT / "src/check_catalogues.cpp",
                "-o",
                d / "control",
            ],
            check=True,
        )
        for w in range(6, 23, 2):
            patterns = d / f"patterns{w}"
            subprocess.run([d / "patterns", str(w), patterns], check=True)
            brute = d / "unpruned_catalogue"
            subprocess.run([d / "control", str(w), brute], check=True)
            if brute.read_bytes() != patterns.read_bytes():
                raise ValueError("Independent complete-matching catalogue differs")
            with patterns.open() as file:
                count = int(file.readline())
            print("catalogue", w, count, flush=True)
            for s in "XZ":
                mat = d / f"author_{s}"
                srcdir = src / "results" / f"affine_{s}"
                mat.write_bytes((srcdir / "matrix.txt").read_bytes())
                ours = json.loads(
                    subprocess.check_output(
                        [d / "embed", d / f"{s}.matrix", patterns, str(w)]
                    )
                )
                if ours["found"]:
                    raise ValueError("Unexpected logical at control weight")
                if w <= 16:
                    published = [
                        json.loads(line)
                        for line in (srcdir / "matrix.txt.filtered.results.jsonl")
                        .read_text()
                        .splitlines()
                    ]
                    base = next(row for row in published if row["weight"] == w)
                    if base["embeddings"] != ours["embeddings_tested"]:
                        raise ValueError(
                            "Mismatch with published unrestricted catalogue"
                        )
                if w >= 18:
                    subprocess.run(
                        [
                            d / "author",
                            mat,
                            patterns,
                            str(w),
                            srcdir / "filter.txt",
                            "1800",
                        ],
                        check=True,
                        stdout=subprocess.DEVNULL,
                    )
                    raw = json.loads(
                        Path(str(mat) + ".filtered.results.jsonl").read_text()
                    )
                    if (
                        raw["found"]
                        or raw["patterns_completed"] != count
                        or raw["embeddings"] != ours["embeddings_tested"]
                    ):
                        raise ValueError("Mismatch with independent embedding kernel")
                print(s, w, count, ours["embeddings_tested"], flush=True)


if __name__ == "__main__":
    main()
