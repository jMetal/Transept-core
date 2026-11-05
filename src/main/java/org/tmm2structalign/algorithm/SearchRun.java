package org.tmm2structalign.algorithm;

import java.util.List;
import java.util.Map;
import java.util.function.IntConsumer;
import org.tmm2structalign.alignment.MsaSolution;
import org.uma.jmetal.component.algorithm.EvolutionaryAlgorithm;
import org.uma.jmetal.component.algorithm.RandomSearchAlgorithm;
import org.uma.jmetal.util.observable.Observable;
import org.uma.jmetal.util.observer.Observer;

/**
 * Uniform view of a jMetal algorithm run, whatever its class: evolutionary algorithms and
 * {@link RandomSearchAlgorithm} expose their progress and results differently.
 */
public interface SearchRun {

  /** Runs the search to completion. */
  void run();

  /**
   * Returns the current result: the population for population-based algorithms, the archive for
   * PAES and the non-dominated solutions found so far for random search.
   */
  List<MsaSolution> front();

  /** Returns the number of evaluations performed so far. */
  int evaluations();

  /** Returns the wall-clock time of the last run, in milliseconds. */
  long runtimeMillis();

  /** Returns the name reported by jMetal. */
  String name();

  /**
   * Returns the configuration of the algorithm ({@code --parameter value ...}) if it was built
   * from a parameter space, or an empty string.
   */
  default String configuration() {
    return "";
  }

  /**
   * Registers a listener called with the number of evaluations every time the algorithm reports
   * progress.
   *
   * @param listener progress listener
   */
  void onProgress(IntConsumer listener);

  /** Adapts an evolutionary algorithm. */
  static SearchRun of(EvolutionaryAlgorithm<MsaSolution> algorithm) {
    return new SearchRun() {
      @Override
      public void run() {
        algorithm.run();
      }

      @Override
      public List<MsaSolution> front() {
        return algorithm.result();
      }

      @Override
      public int evaluations() {
        return algorithm.numberOfEvaluations();
      }

      @Override
      public long runtimeMillis() {
        return algorithm.totalComputingTime();
      }

      @Override
      public String name() {
        return algorithm.name();
      }

      @Override
      public void onProgress(IntConsumer listener) {
        algorithm.observable().register(progressObserver(listener));
      }
    };
  }

  /** Adapts a random search. */
  static SearchRun of(RandomSearchAlgorithm<MsaSolution> algorithm) {
    return new SearchRun() {
      @Override
      public void run() {
        algorithm.run();
      }

      @Override
      public List<MsaSolution> front() {
        return algorithm.result();
      }

      @Override
      public int evaluations() {
        return algorithm.evaluations();
      }

      @Override
      public long runtimeMillis() {
        return algorithm.totalComputingTime();
      }

      @Override
      public String name() {
        return algorithm.name();
      }

      @Override
      public void onProgress(IntConsumer listener) {
        algorithm.getObservable().register(progressObserver(listener));
      }
    };
  }

  /**
   * Attaches the configuration that built an algorithm to its run.
   *
   * @param run the run
   * @param configuration configuration string of the algorithm
   * @return the same run, reporting the configuration
   */
  static SearchRun withConfiguration(SearchRun run, String configuration) {
    return new SearchRun() {
      @Override
      public void run() {
        run.run();
      }

      @Override
      public List<MsaSolution> front() {
        return run.front();
      }

      @Override
      public int evaluations() {
        return run.evaluations();
      }

      @Override
      public long runtimeMillis() {
        return run.runtimeMillis();
      }

      @Override
      public String name() {
        return run.name();
      }

      @Override
      public String configuration() {
        return configuration;
      }

      @Override
      public void onProgress(IntConsumer listener) {
        run.onProgress(listener);
      }
    };
  }

  private static Observer<Map<String, Object>> progressObserver(IntConsumer listener) {
    return (Observable<Map<String, Object>> observable, Map<String, Object> data) -> {
      if (data.get("EVALUATIONS") instanceof Integer evaluations) {
        listener.accept(evaluations);
      }
    };
  }
}
