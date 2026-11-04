# Metrics

## Objectives

The two objectives that the search maximises are defined in the
[formulation tutorial](../tutorials/T2-problem-formulation.md): **SOPwTP**, a topology-aware sum
of pairs, and **lDDT**, a structural consistency score computed from AlphaFold2 distance matrices.
They do not use any reference alignment.

## Comparison with a reference alignment

`org.tmm2structalign.evaluation` compares an alignment with the GPCRdb reference alignment of the
dataset ([how to run it](cli.md#comparison-with-the-reference-alignments)).

### Pair precision, recall and F1

Every pair of residues aligned in the same column is recorded as (protein A, position A, protein B,
position B), with ungapped positions. Then

```text
Pair precision = correct test pairs / test pairs
Pair recall    = correct test pairs / reference pairs
Pair F1        = 2 * precision * recall / (precision + recall)
```

The comparison does not depend on the order of the sequences, the length of the alignment or where
the gaps are.

### TM-pair precision, recall and F1

The same definitions restricted to pairs whose two residues are transmembrane according to the
dataset's three-line topology file.

### TM-gap rate

For every transmembrane residue and every other protein, a TM-gap event occurs when the other
protein has a gap in that column:

```text
TM-gap rate = TM-residue-to-gap events / (TM residues * (N - 1))
```

It ranges from 0 to 1 regardless of the number of proteins; lower is better.

### Input checks

Identifiers are matched case-insensitively, after normalising UniProt `sp|ID|...` and `tr|ID|...`
headers, chain suffixes and `.pdb`, `.ent`, `.cif` or `.mmcif` extensions. Before scoring, the
evaluator checks that the topology, the reference and every alignment contain the same proteins
with the same ungapped sequences.

## BAliBASE SP and TC

A Java reimplementation of the BAliBASE `bali_score` program existed in an earlier version and was
removed in commit `23bb95e`; it is only available in the git history.
