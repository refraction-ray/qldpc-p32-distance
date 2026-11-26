# Exact distance of the Okada–Kasai P=32 code

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Input: arXiv:2604.20838v4](https://img.shields.io/badge/Input-arXiv%3A2604.20838v4-b31b1b.svg)](https://arxiv.org/abs/2604.20838v4)

This repository certifies **`[[16384,4142,32]]` with `d_X = d_Z = 32`** for the fixed affine-coset CSS matrices in [Okada and Kasai, arXiv:2604.20838v4](https://arxiv.org/html/2604.20838v4). The authors excluded weights through 16 and supplied the two weight-32 logical operators. The computation here excludes the intervening even weights on both sides. It applies to the specified matrices.

## How we found the reduction

1. **Start with the authors' pairing graphs.** Each qubit is incident to one check of each of three colours. A zero-syndrome support meets each check an even number of times, so its incidences can be paired within each check. This gives one perfect matching per colour. The graph of a minimum logical support must be connected. Otherwise each component would separately have zero syndrome, and at least one component would be a smaller logical support.
2. **Read coordinates from the published coefficients.** We asked GPT Astra 6.2 to try excluding weight 18. It deleted one colour at a time from each 512-position base graph and labelled the components left by the other two colours. Each deletion gave eight components. The three labels together distinguish all 512 positions. The input checker verifies that a check of colour `c` varies coordinate `c` across eight values and fixes the other two.
3. **Turn the coordinates into an exclusion rule.** A path avoiding colour `c` leaves coordinate `c` fixed. An edge of colour `c` must change it, including in the lifted matrix, so its endpoints cannot be joined by such a path. Applying this rule before embedding cuts the authors' weight-16 catalogue from **25,375 to 26 graphs**. Both sides still give the authors' respective embedding counts of 69,816 and 69,704. At weight 18 only 14 graphs survive, and neither side has an embedding.
4. **Enumerate the remaining graphs and embeddings.** Two matchings form alternating cycles. The code enumerates their possible sizes and adds the third matching, keeping one representative per colour-preserving isomorphism class. The common cyclic shift of the lift fixes one vertex to sheet zero. The embedder checks all remaining assignments in the original matrices. In the higher layers it also prunes any partial support meeting a weight-8 same-type stabilizer in at least five positions, since multiplying by that stabilizer would strictly reduce the weight of any completion.

## Search results

The recorded runs excluded every weight in the table in both the X and Z sectors. Each sector checked the same graph catalogue. The shard counts describe how the original runs divided the work.

| Weight | Candidate graphs per sector | Shards per sector |
|---:|---:|---:|
| 18 | 14 | 1 |
| 20 | 173 | 1 |
| 22 | 467 | 2 |
| 24 | 5,011 | 20 |
| 26 | 27,429 | 32 |
| 28 | 291,592 | 32 |
| 30 | 2,714,836 | 80 |

The odd column weight forces every pure CSS kernel support to have even weight. Any nontrivial mixed Pauli has a nontrivial pure component of no greater weight. The authors' weight-32 witnesses are verified by zero syndromes and odd overlap with separating logical operators.

## Reproduce

Use Python 3.9+ and a C++17 compiler available as `c++`. The Python code uses only the standard library. Run from the repository root.

```sh
python3 -B src/fetch_inputs.py
python3 -B src/check_record.py
python3 -B src/verify.py --through 18 --workers 4
```

The fetch command verifies the arXiv v4 archive and also accepts a local tar or tar.gz file. The record check validates the reconstructed matrices and both weight-32 witnesses, then checks coverage in the recorded search. The last command reruns the search through weight 18 and reports `[20, 32]`.

```sh
python3 -B src/check_controls.py --published
python3 -B src/verify.py --through 30 --workers 10 --output data/runs/full.jsonl
python3 -B src/check_record.py --recomputed data/runs/full.jsonl
```

The control command checks both weight-32 positive cases and cross-checks low weights through 22. Use a fresh output path for each run. The full command regenerates the graph catalogues and writes the counts per weight. The original weight-30 run took about an hour including catalogue generation. Allow 2 GB of temporary disk space. The comparison checks the new counts against the recorded search. The independent controls cover weights through 22. The recorded weight-24 through weight-30 exclusions come from this project's enumeration.

`src/` contains the search and checks. `data/recorded/` contains counts from the original runs, including the completion count for each weight-30 shard. Downloaded inputs and new run records go to `data/upstream/` and `data/runs/`, respectively.

## Citation and license

Cite Koki Okada and Kenta Kasai, *High-Girth Regular Quantum LDPC Codes from Affine-Coset Structures*, for the input matrices, pairing-graph approach, and weight-32 witnesses. For the distance computation, cite this repository at the revision used. The code in `src/` is available under the [Apache License 2.0](LICENSE). The upstream archive is obtained directly from [arXiv](https://arxiv.org/src/2604.20838v4).
