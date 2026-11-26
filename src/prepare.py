"""Rebuild the two full CSS matrices from Okada–Kasai's published lift coefficients."""

import csv
import hashlib
import json
from pathlib import Path

LIFT = 32
N = 512 * LIFT
M = 192 * LIFT


def matrix_digest(rows):
    data = json.dumps(rows, separators=(",", ":")).encode()
    return hashlib.sha256(data).hexdigest()


def syndrome(columns, support):
    odd = set()
    for q in support:
        odd.symmetric_difference_update(columns[q])
    return odd


def construct(coefficients):
    rows = {s: [[] for _ in range(M)] for s in "XZ"}
    base = {s: [[] for _ in range(192)] for s in "XZ"}
    with coefficients.open() as f:
        for record in csv.DictReader(f):
            side = {"HX": "X", "HZ": "Z"}[record["side"]]
            r, q, shift = (
                int(record[k]) for k in ("base_row", "base_col", "shift_mod_32")
            )
            if not (0 <= r < 192 and 0 <= q < 512 and 0 <= shift < LIFT):
                raise ValueError("Invalid lift coefficient")
            base[side][r].append(q)
            for t in range(LIFT):
                rows[side][r * LIFT + t].append(q * LIFT + (t + shift) % LIFT)
    columns = {}
    for side, matrix in rows.items():
        cols = [[] for _ in range(N)]
        for r, qs in enumerate(matrix):
            qs.sort()
            if len(set(qs)) != 8 or len(qs) != 8:
                raise ValueError("Unexpected row weight")
            for q in qs:
                cols[q].append(r)
        if any([r // (M // 3) for r in col] != [0, 1, 2] for col in cols):
            raise ValueError("Each qubit must see one check per colour")
        columns[side] = cols
    if any(syndrome(columns["Z"], row) for row in rows["X"]):
        raise ValueError("CSS commutation failed")
    return rows, columns, base


def coordinates(lines):
    """For each deleted colour, label the eight components of the other two."""
    labels = [[] for _ in range(512)]
    for colour in range(3):
        parent = list(range(512))

        def root(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for r, line in enumerate(lines):
            if r // 64 != colour:
                for q in line:
                    parent[root(q)] = root(line[0])
        names = {}
        for q in range(512):
            v = root(q)
            labels[q].append(names.setdefault(v, len(names)))
        if len(names) != 8:
            raise ValueError("The base is not the expected Cartesian product")
    coords = list(map(tuple, labels))
    if len(set(coords)) != 512:
        raise ValueError("Base coordinates are not one-to-one")
    for r, line in enumerate(lines):
        for c in range(3):
            expected = 8 if c == r // 64 else 1
            if len({coords[q][c] for q in line}) != expected:
                raise ValueError("A check fails the coordinate-line property")
    return coords


def verify_and_write(reference, directory):
    rows, columns, base = construct(reference / "affine.csv")
    cert = json.loads((reference / "affine_certificates.json").read_text())
    hashes = {s: matrix_digest(rows[s]) for s in "XZ"}
    if hashes != cert["matrix_sha256"]:
        raise ValueError("Matrix hash differs from the author's certificate")
    for logical, other in (("X", "Z"), ("Z", "X")):
        witness = cert["sides"][logical]["support"]
        separator = cert["sides"][logical]["separating_support"]
        if (
            len(witness) != 32
            or len(set(witness)) != 32
            or len(set(separator)) != len(separator)
            or any(not 0 <= q < N for q in witness + separator)
            or syndrome(columns[other], witness)
            or syndrome(columns[logical], separator)
            or len(set(witness) & set(separator)) % 2 != 1
        ):
            raise ValueError("Published weight-32 logical certificate failed")
        xyz = coordinates(base[other])
        with (directory / f"{logical}.matrix").open("w") as out:
            print(N, M, file=out)
            for block in (xyz, columns[other], rows[logical]):
                for row in block:
                    print(*row, file=out)
    return hashes


if __name__ == "__main__":
    import sys

    print(verify_and_write(Path(sys.argv[1]), Path(sys.argv[2])))
