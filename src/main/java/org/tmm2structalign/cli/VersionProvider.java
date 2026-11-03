package org.tmm2structalign.cli;

import java.util.Map;
import org.tmm2structalign.run.BuildInfo;
import picocli.CommandLine.IVersionProvider;

/** Reports the git commit of the build as the version. */
public final class VersionProvider implements IVersionProvider {

  @Override
  public String[] getVersion() {
    Map<String, Object> info = BuildInfo.describe();
    return new String[] {
      "TM-M2StructAlign, commit " + info.get("gitCommit") + " (" + info.get("gitBranch") + ")",
      "Java " + info.get("javaVersion")
    };
  }
}
