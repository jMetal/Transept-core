# Operators

The crossovers and mutations of Transept core, with their names on the command line and in the
configurations, and their parameters. None of them ever separates two transmembrane residues that
were together: a gap can only appear inside a helix if a seed alignment already had one there.

The [formulation tutorial](../tutorials/T2-problem-formulation.md#step-8-crossover) explains how
each one works, with examples. The tutorial *Choosing the operators* of the Transept interface
(`docs/tutorials/operators.md` in the Transept repository) says when each one helps, from a pilot
study.

## Crossovers

| `--crossover` | In a configuration | Parameters | What it does |
| --- | --- | --- | --- |
| `TM_SAFE_SINGLE_POINT` (default) | `tmSafeSinglePoint` | `cutCandidates` (1–64, default 1) | Cuts the first parent at a column that splits no helix, and continues with the rest of every protein as aligned in the second parent, padding with gaps to keep that part aligned. With `cutCandidates` > 1, draws that many cuts and uses the one that adds the fewest gaps |
| `TM_SAFE_UNIFORM` | `tmSafeUniform` | none | Exchanges, at random, the blocks between the columns where both parents have the same residues on the left in every protein. Adds no gap |

Both use `--crossover-probability` (`crossoverProbability`, default 0.8). The crossover that fills
the initial population is always the single-point one, with `cutCandidates`. PAES and the random
searches do not use the crossover.

## Mutations

| `--mutation` | In a configuration | Parameters | What it does |
| --- | --- | --- | --- |
| `TM_SAFE_GAP_SHIFT` (default) | `tmSafeGapShift` | `maxShift` (0–200, default 0 = no limit) | Moves a run of gaps of one protein to another position that splits no helix |
| `MERGE_GAP_RUNS` | `mergeGapRuns` | none | Joins two consecutive gap runs of one protein |
| `RANDOM_CHOICE` | `randomChoice` | `maxShift`, `blockShift`, `realignWindow` | One of the other four mutations, chosen uniformly at random each time. Not evaluated yet |
| `TM_SAFE_BLOCK_SHIFT` | `tmSafeBlockShift` | `blockShift` (1–20, default 1) | Moves a block of gaps shared by several proteins up to `blockShift` columns, keeping those proteins aligned with each other (SAGA) |
| `TM_SAFE_WINDOW_REALIGN` | `tmSafeWindowRealign` | `realignWindow` (10–500, default 50) | Realigns the residues of one protein against the others inside a window of columns, by dynamic programming with the substitution matrices and gap penalties of SOPwTP. It approximates SOPwTP (the gaps of the other proteins are simplified), so it usually raises it but can lower it |

All use `--mutation-probability` (`mutationProbability`, default 0.2; 1.0 for PAES).

## Where to set them

- **Command line**, for every algorithm: `--crossover`, `--mutation`, `--crossover-probability`,
  `--mutation-probability`, `--max-shift` and `--cut-candidates` ([cli](cli.md)). `blockShift` and
  `realignWindow` keep their defaults there.
- **Configuration**, for NSGA-II, SPEA2 and the genetic algorithm: every parameter of the tables,
  e.g. `--configuration "--crossover tmSafeUniform --mutation tmSafeWindowRealign --realignWindow 100"`
  ([configuration](configuration.md)). A parameter that belongs to another operator, such as
  `maxShift` with `tmSafeBlockShift`, has no effect.
- A parameter can be given in one place only: `--mutation` together with `--mutation` in the
  configuration is an error.
