#!/usr/bin/env python3
"""Compare configurations from the per-run indicators of pareto_analysis.py or formulation_study.py.

Usage:
    compare_runs.py <indicators_by_run.csv> --baseline <configuration> [--indicator hypervolume]

Run identifiers must have the form <configuration>_seed_<seed>, e.g. nsgaii_seed_1. Runs of two
configurations with the same seed are paired. For every configuration the script prints the
median and interquartile range of the indicator, how many seeds it beats the baseline on, and the
two-sided p-value of the exact Wilcoxon signed-rank test against the baseline (no SciPy needed).
"""
from __future__ import annotations

import argparse
import csv
import itertools
import statistics
import sys
from collections import defaultdict

LOWER_IS_BETTER = {"igd_plus"}


def lower_is_better(indicator: str) -> bool:
    """IGD+ and TM-gap rates (also prefixed, e.g. r1_tm_gap_rate) are better when lower."""
    return indicator in LOWER_IS_BETTER or indicator.endswith("tm_gap_rate")


def wilcoxon_exact(first: list[float], second: list[float]) -> float:
    """Two-sided exact p-value of the Wilcoxon signed-rank test; zero differences are dropped."""
    differences = [a - b for a, b in zip(first, second) if a != b]
    if not differences:
        return 1.0
    order = sorted(range(len(differences)), key=lambda i: abs(differences[i]))
    ranks = [0.0] * len(differences)
    start = 0
    while start < len(order):
        end = start
        while end + 1 < len(order) and abs(differences[order[end + 1]]) == abs(differences[order[start]]):
            end += 1
        for position in range(start, end + 1):
            ranks[order[position]] = (start + end) / 2 + 1
        start = end + 1
    positive = sum(rank for rank, difference in zip(ranks, differences) if difference > 0)
    total = sum(ranks)
    observed = min(positive, total - positive)
    extreme = 0
    for signs in itertools.product((False, True), repeat=len(ranks)):
        statistic = sum(rank for rank, sign in zip(ranks, signs) if sign)
        if min(statistic, total - statistic) <= observed + 1e-12:
            extreme += 1
    return extreme / 2 ** len(ranks)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "indicators",
        help="indicators_by_run.csv of pareto_analysis.py or the output of formulation_study.py",
    )
    parser.add_argument("--baseline", required=True, help="Configuration to compare against")
    parser.add_argument("--indicator", default="hypervolume")
    args = parser.parse_args()

    values: dict[str, dict[str, float]] = defaultdict(dict)
    with open(args.indicators, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            configuration, separator, seed = row["run"].rpartition("_seed_")
            if not separator:
                sys.exit(f"Run '{row['run']}' does not follow <configuration>_seed_<seed>")
            if row[args.indicator] != "":  # undefined for this run, e.g. R2 with one objective
                values[configuration][seed] = float(row[args.indicator])
    if args.baseline not in values:
        sys.exit(f"No runs of the baseline '{args.baseline}'; found: {', '.join(sorted(values))}")

    lower = lower_is_better(args.indicator)
    baseline = values[args.baseline]
    print(f"{args.indicator} ({'lower' if lower else 'higher'} is better), "
          f"baseline {args.baseline}")
    print(f"{'configuration':20s} {'runs':>4s} {'median':>8s} {'IQR':>8s} {'wins':>7s} {'p':>7s}")
    for configuration in sorted(values, key=lambda name: (name != args.baseline, name)):
        runs = values[configuration]
        series = list(runs.values())
        quartiles = statistics.quantiles(series, n=4) if len(series) > 1 else [series[0]] * 3
        wins = p_value = ""
        if configuration != args.baseline:
            seeds = sorted(set(runs) & set(baseline))
            mine = [runs[seed] for seed in seeds]
            theirs = [baseline[seed] for seed in seeds]
            better = sum((m < t) if lower else (m > t) for m, t in zip(mine, theirs))
            wins = f"{better}/{len(seeds)}"
            p_value = f"{wilcoxon_exact(mine, theirs):.3f}"
        print(f"{configuration:20s} {len(series):4d} {statistics.median(series):8.3f} "
              f"{quartiles[2] - quartiles[0]:8.3f} {wins:>7s} {p_value:>7s}")


if __name__ == "__main__":
    main()
