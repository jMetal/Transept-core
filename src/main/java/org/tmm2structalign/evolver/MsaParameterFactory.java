package org.tmm2structalign.evolver;

import java.util.List;
import org.tmm2structalign.alignment.MsaSolution;
import org.uma.evolver.parameter.catalogue.ExternalArchiveParameter;
import org.uma.evolver.parameter.catalogue.selectionparameter.SelectionParameter;
import org.uma.evolver.parameter.factory.ParameterFactory;
import org.uma.evolver.parameter.type.CategoricalParameter;

/** Maps the names of a parameter space in YAML to the parameters of the MSA algorithms. */
public class MsaParameterFactory implements ParameterFactory<MsaSolution> {

  @Override
  public CategoricalParameter createParameter(String parameterName, List<String> values) {
    if (values == null || values.isEmpty()) {
      throw new IllegalArgumentException("Parameter " + parameterName + " has no values");
    }
    return switch (parameterName) {
      case "archiveType" -> new ExternalArchiveParameter<MsaSolution>("archiveType", values);
      case "createInitialSolutions" -> new CreateInitialSolutionsMsaParameter(values);
      case "crossover" -> new MsaCrossoverParameter(values);
      case "mutation" -> new MsaMutationParameter(values);
      case "selection" -> new SelectionParameter<MsaSolution>(values);
      case "variation" -> new MsaVariationParameter(values);
      default -> new CategoricalParameter(parameterName, values);
    };
  }
}
