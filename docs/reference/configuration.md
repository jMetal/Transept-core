# Configuring the algorithms

NSGA-II, SPEA2 and the genetic algorithm are built with [Evolver](https://github.com/jMetal/Evolver)
from a **parameter space**:
a YAML file that lists the components of the algorithm and the values each one can take. A
**configuration** chooses one value for every parameter. The spaces are in
`src/main/resources/parameterSpaces/` and the default configurations, used when nothing is given, in
`src/main/resources/defaultConfigurations/`; both are inside the JAR.

| `--algorithm` | Parameter space | Default configuration |
| --- | --- | --- |
| `NSGA_II` | `MSANSGAII.yaml` | `MSANSGAIIDefault.txt` (the published one) |
| `SPEA2` | `MSASPEA2.yaml` | `MSASPEA2Default.txt` |
| `GENETIC_ALGORITHM` | `MSAGA.yaml` | `MSAGADefault.txt` |

## Giving a configuration

A configuration is a list of `--parameter value` pairs, passed with `--configuration` (or, from
a file, `--configuration-file`). Only the parameters you give change; the rest keep their default
value:

```bash
java -jar target/transept-core-0.1.0-jar-with-dependencies.jar run \
  --data-root ../Dataset_Structural_GPCRs/GPCRdb --dataset classF \
  --seeds clustalw,kalign,mafft,muscle5,tcoffee,tmaligner \
  --output-root results --run-id random-selection --seed 1 \
  --configuration "--selection random --mutation mergeGapRuns --mutationProbability 0.4"
```

`metadata.json` records the complete configuration that was used (`configuration`), so a run can
be repeated with `--configuration "<that string>"`. Unknown parameter names and values outside
the space are errors.

## The parameter space of NSGA-II

| Parameter | Values | Default | Meaning |
| --- | --- | --- | --- |
| `algorithmResult` | `population`, `externalArchive` | `population` | What the run returns: the final population, or an archive of the non-dominated solutions found during the search |
| `populationSizeWithArchive` | integer, 2–1000 | (none) | With `externalArchive`: the population of the search; the archive has `--population-size` solutions |
| `archiveType` | `crowdingDistanceArchive`, `unboundedArchive` | (none) | With `externalArchive`: bounded by crowding distance, or unbounded |
| `createInitialSolutions` | `seededDistinct` | `seededDistinct` | The seed alignments followed by distinct crossover children ([T2, Step 7](../tutorials/T2-problem-formulation.md#step-7-the-initial-population)) |
| `offspringPopulationSize` | integer, 1–1000 | the population size | Children per generation |
| `variation` | `crossoverAndMutationVariation` | | Crossover, then mutation |
| `crossover` | `tmSafeSinglePoint`, `tmSafeUniform` | `tmSafeSinglePoint` | The crossovers of `--crossover` ([cli](cli.md)) |
| `crossoverProbability` | 0.0–1.0 | 0.8 | Probability of crossing two parents |
| `cutCandidates` | integer, 1–64 | 1 | Valid cuts the single-point crossover draws, of which it uses the one that adds the fewest gaps; 1 is a cut at random. It also applies to the crossover that fills the initial population, so it counts with `tmSafeUniform` too |
| `mutation` | `tmSafeGapShift`, `mergeGapRuns`, `randomChoice`, `tmSafeBlockShift`, `tmSafeWindowRealign` | `tmSafeGapShift` | The mutations of `--mutation` ([cli](cli.md)) |
| `mutationProbability` | 0.0–1.0 | 0.2 | Probability of mutating a child |
| `maxShift` | integer, 0–200 | 0 | With `tmSafeGapShift` or `randomChoice`: largest number of columns a gap run can move; 0 is no limit |
| `blockShift` | integer, 1–20 | 1 | With `tmSafeBlockShift` or `randomChoice`: largest number of columns a block of gaps can move |
| `realignWindow` | integer, 10–500 | 50 | With `tmSafeWindowRealign` or `randomChoice`: number of columns of the window that is realigned |
| `selection` | `tournament`, `random` | `tournament` | How parents are chosen |
| `selectionTournamentSize` | integer, 2–10 | 2 | With `tournament`: solutions compared per choice (by rank, then crowding distance) |

The population size is not part of the space: it is `--population-size`. The replacement is
NSGA-II's (non-dominated sorting and crowding distance). The default configuration is the
published one; with it, the results are identical to those of earlier versions with the same
seed.

The operator options of the command line (`--crossover`, `--mutation`, `--crossover-probability`,
`--mutation-probability`, `--max-shift`, `--cut-candidates`) set the same parameters; give each parameter in one place only.

## SPEA2

`MSASPEA2.yaml` is NSGA-II's space plus three parameters of Evolver's RDEMOEA, a generational
algorithm with a configurable ranking and density estimator:

| Parameter | Values | Default |
| --- | --- | --- |
| `ranking` | `strengthRanking`, `dominanceRanking` | `strengthRanking` |
| `densityEstimator` | `knn` (with `knnNeighborhoodSize` 1–5 and `knnNormalizeObjectives`), `crowdingDistance` | `knn`, 1 neighbour, normalised objectives |
| `replacement` | `rankingAndDensityEstimator` with `removalPolicy` `sequential` or `oneShot` | `sequential` |

The default is SPEA2: strength ranking, the distance to the nearest neighbour as density and the
worst solutions removed one by one. Unlike the original SPEA2, the population plays the role of
its archive. With `dominanceRanking`, `crowdingDistance` and `oneShot` it behaves as NSGA-II.

## The genetic algorithm

`MSAGA.yaml` has the initial population, offspring size, variation and selection of NSGA-II (the
tournament compares the single objective) and `gaReplacement`: `muPlusLambda` (default, parents and
children compete) or `muCommaLambda` (only the children survive; it needs more offspring than the
population size). It needs `--objectives LDDT`.

## Other algorithms

The other algorithms of `--algorithm` are still built directly with jMetal and are configured with
the operator options only. SMS-EMOA, MOEA/D and PAES will move to parameter spaces with the same
approach, so that all the sequential base algorithms can be configured automatically later. The
asynchronous algorithms (`ASYNC_NSGA_II`, `ASYNC_GENETIC_ALGORITHM`) are meant for production
runs on many cores and use the operator options; automatic configuration runs the base
algorithms sequentially and parallelises over configurations instead.
