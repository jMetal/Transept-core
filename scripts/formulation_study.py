#!/usr/bin/env python3
"""Compare formulations (sets of objectives) by the alignments they return.

Usage:
    formulation_study.py <output-root> <dataset> --data-root <GPCRdb> --output <study.csv>
                         [--jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar]

For every run under <output-root>/ejecuciones/<dataset>/, the script evaluates its alignments
against the GPCRdb reference with EvaluationMain (once; the reports are kept next to the output
file) and applies two rules that do not look at the reference:

- R1: the alignment with the best lDDT (ties: best SOPwTP, if optimised, then the first row);
- R2: the knee, the distinct non-dominated solution closest to the ideal point after scaling
  every objective to [0, 1] over the run (not defined with a single objective).

It writes one line per run with the reference-based metrics of R1 and R2, the median and the best
pair F1 of all the returned alignments (the best is an oracle, an upper bound that no rule
without the reference can guarantee), and the objectives of the run. Run identifiers of the form
<configuration>_seed_<seed> let compare_runs.py compare configurations, e.g.:

    compare_runs.py study.csv --baseline nsgaii --indicator r1_pair_f1

The hypervolume cannot compare formulations: each one measures a different objective space.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
import sys
from pathlib import Path

LDDT = "lDDT"
SOPWTP = "SOPwTP"
METRICS = ("pair_f1", "tm_pair_f1", "tm_gap_rate", "coverage")


def read_points(run: Path) -> list[tuple[float, ...]]:
    """Objective values of FUN.csv, one tuple per row (higher is better for all)."""
    lines = (run / "FUN.csv").read_text(encoding="utf-8").splitlines()
    return [tuple(float(value) for value in line.split(",")) for line in lines if line.strip()]


def objective_names(run: Path) -> list[str]:
    """Objectives of the run, from metadata.json; older runs used SOPwTP and lDDT."""
    metadata = json.loads((run / "metadata.json").read_text(encoding="utf-8"))
    return list(metadata.get("objectives") or [SOPWTP, LDDT])


def best_lddt(points: list[tuple[float, ...]], names: list[str]) -> int:
    """Rule R1: row with the best lDDT, ties broken by SOPwTP (if any), then by row."""
    lddt = names.index(LDDT)
    sop = names.index(SOPWTP) if SOPWTP in names else None
    return max(
        range(len(points)),
        key=lambda row: (points[row][lddt], points[row][sop] if sop is not None else 0.0, -row),
    )


def knee(points: list[tuple[float, ...]]) -> int | None:
    """Rule R2: distinct non-dominated row closest to the ideal point after scaling."""
    if not points or len(points[0]) < 2:
        return None
    front: dict[tuple[float, ...], int] = {}
    for row, point in enumerate(points):
        dominated = any(
            all(o >= p for o, p in zip(other, point)) and other != point for other in points
        )
        if not dominated:
            front.setdefault(point, row)
    lows = [min(point[i] for point in front) for i in range(len(points[0]))]
    highs = [max(point[i] for point in front) for i in range(len(points[0]))]

    def distance(point: tuple[float, ...]) -> float:
        return math.sqrt(
            sum(((high - value) / (high - low)) ** 2
                for value, low, high in zip(point, lows, highs) if high > low)
        )

    return front[min(front, key=lambda point: (distance(point), front[point]))]


def evaluate(run: Path, report: Path, jar: str, data_root: Path, dataset: str) -> dict[str, dict]:
    """Reference-based metrics of every alignment of a run, keyed by file name."""
    if not report.is_file():
        subprocess.run(
            [
                "java", "-cp", jar, "org.tmm2structalign.evaluation.EvaluationMain",
                "--reference", str(data_root / "reference_alignments" / f"{dataset}.fasta"),
                "--topology",
                str(data_root / "sequences" / "tmregions"
                    / f"{dataset}_predicted_topologies.3line"),
                "--msa-root", str(run),
                "--output", str(report),
            ],
            check=True,
            capture_output=True,
        )
    metrics = {}
    with report.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            row["coverage"] = int(row["test_pairs"]) / int(row["reference_pairs"])
            metrics[row["msa_file"]] = row
    return metrics


def summarize(run: Path, metrics: dict[str, dict]) -> dict[str, object]:
    points = read_points(run)
    names = objective_names(run)
    line: dict[str, object] = {"run": run.name, "objectives": " ".join(names)}
    for rule, row in (("r1", best_lddt(points, names)), ("r2", knee(points))):
        line[f"{rule}_row"] = "" if row is None else row
        for metric in METRICS:
            value = None if row is None else metrics[f"MSASol{row}.fasta"][metric]
            line[f"{rule}_{metric}"] = "" if value is None else float(value)
    f1 = [float(metrics[f"MSASol{row}.fasta"]["pair_f1"]) for row in range(len(points))]
    line["median_pair_f1"] = statistics.median(f1)
    line["best_pair_f1"] = max(f1)
    line["solutions"] = len(points)
    return line


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_root", type=Path)
    parser.add_argument("dataset")
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--jar", default="target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar")
    arguments = parser.parse_args()

    runs = sorted(
        path for path in (arguments.output_root / "ejecuciones" / arguments.dataset).iterdir()
        if (path / "metadata.json").is_file()
    )
    if not runs:
        sys.exit(f"No finished runs of {arguments.dataset} under {arguments.output_root}")
    reports = arguments.output.parent / f"evaluations_{arguments.dataset}"
    reports.mkdir(parents=True, exist_ok=True)
    lines = []
    for run in runs:
        metrics = evaluate(
            run, reports / f"{run.name}.csv", arguments.jar, arguments.data_root, arguments.dataset)
        lines.append({"dataset": arguments.dataset, **summarize(run, metrics)})
    with arguments.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(lines[0]))
        writer.writeheader()
        writer.writerows(lines)
    print(f"{len(lines)} runs of {arguments.dataset} -> {arguments.output}")


if __name__ == "__main__":
    main()
