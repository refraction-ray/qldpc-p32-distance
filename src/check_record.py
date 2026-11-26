"""Check the recorded search coverage and optionally compare a fresh run."""

import argparse
import json
from pathlib import Path
import tempfile

from fetch_inputs import ROOT, UPSTREAM
from prepare import verify_and_write
from verify import INPUT_FILES

RECORDED = ROOT / "data" / "recorded"


def read_records(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def layer_signature(layer):
    weight = layer["weight"]
    count = layer["catalogue_size"]
    sectors = layer["sectors"]
    require(set(sectors) == {"X", "Z"}, f"Missing sector at weight {weight}")
    require(
        all(
            record["completed_patterns"] == count and not record["witnesses"]
            for record in sectors.values()
        ),
        f"Incomplete exclusion at weight {weight}",
    )
    require(
        (layer["lower"], layer["upper"]) == (weight + 2, 32),
        f"Inconsistent bounds after weight {weight}",
    )
    return (
        count,
        layer["minimality_pruning"],
        tuple(sectors[sector]["embeddings_tested"] for sector in "XZ"),
    )


def check_archive():
    layers = read_records(RECORDED / "through28.jsonl")
    header, *shards = read_records(RECORDED / "weight30.jsonl")
    require([r["weight"] for r in layers] == list(range(2, 30, 2)), "Missing layer")
    expected = {layer["weight"]: layer_signature(layer) for layer in layers}
    count, shard_count = header["catalogue_size"], header["shards"]
    require(count >= shard_count > 0, "Invalid weight-30 header")
    records = {(r["sector"], r["shard"]): r for r in shards}
    require(
        len(shards) == len(records) == 2 * shard_count,
        "Duplicate or missing weight-30 shard",
    )
    embeddings = {"X": 0, "Z": 0}
    for sector in "XZ":
        for shard in range(shard_count):
            row = records[(sector, shard)]
            size = (count + shard_count - 1 - shard) // shard_count
            require(
                not row["found"]
                and row["completed_patterns"] == row["shard_size"] == size,
                f"Incomplete weight-30 exclusion in {sector} shard {shard}",
            )
            embeddings[sector] += row["embeddings_tested"]
    expected[30] = (
        count,
        True,
        tuple(embeddings[sector] for sector in "XZ"),
    )
    return expected


def compare_run(path, hashes, expected):
    header, *layers = read_records(path)
    require(header["matrix_sha256"] == hashes, "New run uses different matrices")
    weights = [layer["weight"] for layer in layers]
    require(
        bool(weights)
        and weights == list(range(2, weights[-1] + 1, 2))
        and weights[-1] <= 30,
        "New run is missing layers",
    )
    for layer in layers:
        weight = layer["weight"]
        require(
            layer_signature(layer) == expected[weight],
            f"New and archived records disagree at weight {weight}",
        )
    print(f"Fresh run matches the archived search through weight {weights[-1]}.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recomputed", type=Path, help="Compare this new run as well")
    args = parser.parse_args()
    if not all(path.is_file() for path in INPUT_FILES):
        parser.error("run python3 -B src/fetch_inputs.py first")
    with tempfile.TemporaryDirectory() as scratch:
        hashes = verify_and_write(UPSTREAM, Path(scratch))
    expected = check_archive()
    print("Both weight-32 witnesses verified and archive coverage checked through weight 30.")
    print("Recorded result: [[16384,4142,32]], d_X = d_Z = 32.")
    if args.recomputed:
        compare_run(args.recomputed, hashes, expected)


if __name__ == "__main__":
    main()
