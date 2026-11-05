package org.transept.evolver;

import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import org.uma.evolver.parameter.ParameterManagement;
import org.uma.evolver.parameter.ParameterSpace;
import org.uma.evolver.parameter.yaml.YAMLParameterSpace;

/**
 * Configurations of the algorithms built from a parameter space: strings of the form
 * {@code --parameter value ...}, as Evolver writes and reads them.
 *
 * <p>A configuration is resolved in layers: the default configuration of the algorithm
 * ({@code defaultConfigurations/<space>Default.txt}), then the values that the caller derives from
 * other options, then the configuration given by the user. Later layers replace earlier values,
 * and names that the parameter space does not define are rejected, because Evolver would ignore
 * them silently.
 */
public final class MsaConfigurations {

  private MsaConfigurations() {}

  /**
   * Loads a parameter space of the MSA algorithms.
   *
   * @param resource parameter space resource, e.g. {@link MsaNsga2#PARAMETER_SPACE}
   * @return the parameter space
   */
  public static ParameterSpace parameterSpace(String resource) {
    return new YAMLParameterSpace(resource, new MsaParameterFactory());
  }

  /**
   * Reads the default configuration of a parameter space.
   *
   * @param resource parameter space resource, e.g. {@code MSANSGAII.yaml}
   * @return the default values, in their order
   * @throws IllegalArgumentException if there is no default configuration
   */
  public static Map<String, String> defaults(String resource) {
    String path = "defaultConfigurations/" + resource.replace(".yaml", "Default.txt");
    try (InputStream input = MsaConfigurations.class.getClassLoader().getResourceAsStream(path)) {
      if (input == null) {
        throw new IllegalArgumentException("No default configuration: " + path);
      }
      return parse(new String(input.readAllBytes(), StandardCharsets.UTF_8));
    } catch (IOException exception) {
      throw new UncheckedIOException(exception);
    }
  }

  /**
   * Parses a configuration.
   *
   * @param text configuration such as {@code --selection tournament --selectionTournamentSize 2}
   * @return the values by parameter name, in their order
   * @throws IllegalArgumentException if a token is not a {@code --name value} pair or a name is
   *     repeated
   */
  public static Map<String, String> parse(String text) {
    Map<String, String> values = new LinkedHashMap<>();
    String trimmed = text == null ? "" : text.strip();
    if (trimmed.isEmpty()) {
      return values;
    }
    String[] tokens = trimmed.split("\\s+");
    if (tokens.length % 2 != 0) {
      throw new IllegalArgumentException(
          "A configuration is a list of --name value pairs: " + text);
    }
    for (int index = 0; index < tokens.length; index += 2) {
      String name = tokens[index];
      if (!name.startsWith("--") || name.length() == 2 || tokens[index + 1].startsWith("--")) {
        throw new IllegalArgumentException("Expected --name value, found: " + name + " "
            + tokens[index + 1]);
      }
      if (values.put(name.substring(2), tokens[index + 1]) != null) {
        throw new IllegalArgumentException("Parameter given twice: " + name);
      }
    }
    return values;
  }

  /**
   * Formats a configuration.
   *
   * @param values values by parameter name
   * @return the configuration string
   */
  public static String format(Map<String, String> values) {
    return values.entrySet().stream()
        .map(entry -> "--" + entry.getKey() + " " + entry.getValue())
        .collect(Collectors.joining(" "));
  }

  /**
   * Checks that every name is a parameter of the space.
   *
   * @param values values by parameter name
   * @param space parameter space
   * @throws IllegalArgumentException if a name is unknown
   */
  public static void requireKnownNames(Map<String, String> values, ParameterSpace space) {
    Set<String> known =
        ParameterManagement.parameterFlattening(space.topLevelParameters()).stream()
            .map(parameter -> parameter.name())
            .collect(Collectors.toSet());
    for (String name : values.keySet()) {
      if (!known.contains(name)) {
        throw new IllegalArgumentException(
            "Unknown parameter --" + name + "; the parameter space defines " + known);
      }
    }
  }

  /**
   * Converts a configuration into the arguments that Evolver parses.
   *
   * @param values values by parameter name
   * @return the arguments, {@code --name value} pairs
   */
  public static String[] arguments(Map<String, String> values) {
    return format(values).split(" ");
  }
}
