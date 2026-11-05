package org.tmm2structalign.evolver;

import java.util.List;

/** Checks shared by the parameters of the parameter spaces. */
final class MsaParameters {

  private MsaParameters() {}

  static void requireKnown(List<String> values, List<String> known, String parameter) {
    for (String value : values) {
      if (!known.contains(value)) {
        throw new IllegalArgumentException(
            "Unknown value of " + parameter + ": " + value + "; valid: " + known);
      }
    }
  }
}
