# Configuring NSGA-II

NSGA-II is built with [Evolver](https://github.com/jMetal/Evolver) from a **parameter space**:
a YAML file that lists the components of the algorithm and the values each one can take. A
**configuration** chooses one value for every parameter. The parameter space is
`src/main/resources/parameterSpaces/MSANSGAII.yaml` and the default configuration, the one used
when nothing is given, is `src/main/resources/defaultConfigurations/MSANSGAIIDefault.txt`. Both are
inside the JAR.

## Giving a configuration

A configuration is a list of `--parameter value` pairs, passed with `--configuration` (or, from
a file, `--configuration-file`). Only the parameters you give change; the rest keep their default
value:

```bash
java -jar target/TM-MSAStructAligner-1.0-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results --run-id random-selection --seed 1 \
  --configuration "--selection random --mutation mergeGapRuns --mutationProbability 0.4"
```

`metadata.json` records the complete configuration that was used (`configuration`), so a run can
be repeated with `--configuration "<that string>"`. Unknown parameter names and values outside
the space are errors.

## The parameter space

| Parameter | Values | Default | Meaning |
| --- | --- | --- | --- |
| `algorithmResult` | `population`, `externalArchive` | `population` | What the run returns: the final population, or an archive of the non-dominated solutions found during the search |
| `populationSizeWithArchive` | integer, 2–1000 | (none) | With `externalArchive`: the population of the search; the archive has `--population-size` solutions |
| `archiveType` | `crowdingDistanceArchive`, `unboundedArchive` | (none) | With `externalArchive`: bounded by crowding distance, or unbounded |
| `createInitialSolutions` | `seededDistinct` | `seededDistinct` | The seed alignments followed by distinct crossover children ([T2, Step 7](../tutorials/T2-problem-formulation.md#step-7-the-initial-population)) |
| `offspringPopulationSize` | integer, 1–1000 | the population size | Children per generation |
| `variation` | `crossoverAndMutationVariation` | | Crossover, then mutation |
| `crossover` | `tmSafeSinglePoint` | `tmSafeSinglePoint` | Single-point crossover that never cuts a transmembrane helix |
| `crossoverProbability` | 0.0–1.0 | 0.8 | Probability of crossing two parents |
| `mutation` | `tmSafeGapShift`, `mergeGapRuns`, `randomChoice` | `tmSafeGapShift` | The mutations of `--mutation` ([cli](cli.md)) |
| `mutationProbability` | 0.0–1.0 | 0.2 | Probability of mutating a child |
| `selection` | `tournament`, `random` | `tournament` | How parents are chosen |
| `selectionTournamentSize` | integer, 2–10 | 2 | With `tournament`: solutions compared per choice (by rank, then crowding distance) |

The population size is not part of the space: it is `--population-size`. The replacement is
NSGA-II's (non-dominated sorting and crowding distance). The default configuration is the
published one; with it, the results are identical to those of earlier versions with the same
seed.

The operator options of the command line (`--mutation`, `--crossover-probability`,
`--mutation-probability`) set the same parameters; give each parameter in one place only.

## Other algorithms

The other algorithms of `--algorithm` are still built directly with jMetal and are configured with
the operator options only. They will move to parameter spaces one by one, with the same
approach.
