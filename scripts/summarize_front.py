#!/usr/bin/env python3
"""Usage: summarize_front.py <run-directory> [<evaluation.csv>]

Lists the distinct non-dominated solutions of a run, from the best SOPwTP to the best lDDT,
marks the two extremes and the knee (closest to the ideal point after normalising both
objectives), and adds the reference-based metrics if an evaluation file is given.
"""
import csv
import math
import pathlib
import sys

run = pathlib.Path(sys.argv[1])
points = [tuple(map(float, line.split(","))) for line in open(run / "FUN.csv") if line.strip()]
metrics = {}
if len(sys.argv) > 2:
    metrics = {row["msa_file"]: row for row in csv.DictReader(open(sys.argv[2]))}

front = {}
for index, (sop, lddt) in enumerate(points):
    dominated = any(s >= sop and l >= lddt and (s, l) != (sop, lddt) for s, l in points)
    if not dominated:
        front.setdefault((sop, lddt), index)  # keep the first of repeated solutions

sops = [s for s, _ in front]
lddts = [l for _, l in front]
def distance(sop, lddt):
    return math.hypot((max(sops) - sop) / (max(sops) - min(sops)),
                      (max(lddts) - lddt) / (max(lddts) - min(lddts)))
knee = min(front, key=lambda point: distance(*point))

print(f"{len(points)} solutions, {len(front)} distinct non-dominated")
print("file            SOPwTP    lDDT    pair-F1  TM-pair-F1  TM-gap-rate  note")
for sop, lddt in sorted(front, reverse=True):
    name = f"MSASol{front[(sop, lddt)]}.fasta"
    note = ("best SOPwTP" if sop == max(sops) else "best lDDT" if lddt == max(lddts)
            else "knee" if (sop, lddt) == knee else "")
    row = metrics.get(name)
    extra = (f"{float(row['pair_f1']):8.4f}  {float(row['tm_pair_f1']):10.4f}"
             f"  {float(row['tm_gap_rate']):11.4f}" if row else " " * 33)
    print(f"{name:14s} {sop:8.0f}  {lddt:.4f}  {extra}  {note}")
