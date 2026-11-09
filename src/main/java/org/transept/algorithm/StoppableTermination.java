package org.transept.algorithm;

import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import org.uma.jmetal.component.catalogue.common.termination.Termination;

/**
 * A termination condition that is also met when the user asks the run to stop.
 *
 * <p>A stop requested before the algorithm has results waits until it has some: at least two
 * populations evaluated (the initial one and a generation), or a first solution in a random
 * search. Otherwise SPEA2, the asynchronous algorithms or the random searches could return
 * nothing.
 *
 * @param termination the termination of the algorithm, by evaluations or by time
 * @param stop the request to stop early
 */
public record StoppableTermination(Termination termination, StopRequest stop)
    implements Termination {

  /** Validates the arguments. */
  public StoppableTermination {
    if (termination == null || stop == null) {
      throw new IllegalArgumentException("termination and stop are required");
    }
  }

  @Override
  public boolean isMet(Map<String, Object> algorithmStatusData) {
    return termination.isMet(algorithmStatusData)
        || (hasResults(algorithmStatusData) && stop.isRequested());
  }

  /** Whether stopping now would leave results, judging by what the algorithm reports. */
  static boolean hasResults(Map<String, Object> status) {
    boolean ready = true;
    if (status.get("BEST_SOLUTIONS") instanceof List<?> best) {
      ready = !best.isEmpty();
    } else if (status.get("POPULATION") instanceof List<?> population
        && status.get("EVALUATIONS") instanceof Number evaluations) {
      ready = !population.isEmpty() && evaluations.longValue() >= 2L * population.size();
    }
    return ready;
  }

  @Override
  public OptionalDouble progress(Map<String, Object> algorithmStatusData) {
    return termination.progress(algorithmStatusData);
  }
}
