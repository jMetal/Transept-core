package org.transept.cli;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import org.transept.algorithm.AlgorithmType;
import org.transept.run.Json;
import picocli.CommandLine.Command;

/**
 * Lists the algorithms of {@code run --algorithm} and what they accept, as JSON, so that other
 * programs (the studio) do not keep their own copy of these facts.
 */
@Command(
    name = "algorithms",
    mixinStandardHelpOptions = true,
    description = "List the algorithms and the objectives, configuration and threads they accept.")
public final class AlgorithmsCommand implements Callable<Integer> {

  @Override
  public Integer call() {
    List<Map<String, Object>> algorithms = new ArrayList<>();
    for (AlgorithmType type : AlgorithmType.values()) {
      Map<String, Object> algorithm = new LinkedHashMap<>();
      algorithm.put("name", type.name());
      algorithm.put("minObjectives", type.minObjectives());
      algorithm.put("maxObjectives", type.maxObjectives());
      algorithm.put("parameterSpace", type.parameterSpace());
      algorithm.put("steadyState", type.isSteadyState());
      algorithm.put("asynchronous", type.isAsynchronous());
      algorithm.put("usesCrossover", type.usesCrossover());
      algorithms.add(algorithm);
    }
    System.out.println(Json.write(algorithms));
    return 0;
  }
}
