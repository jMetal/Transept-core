package org.tmm2structalign.run;

import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Properties;

/** Version information of the running build, read from {@code git.properties} if present. */
public final class BuildInfo {

  private static final String RESOURCE = "/git.properties";
  private static final String UNKNOWN = "unknown";

  private BuildInfo() {}

  /** Returns the git commit, branch and dirty flag of the build, or "unknown" values. */
  public static Map<String, Object> describe() {
    Properties properties = new Properties();
    try (InputStream input = BuildInfo.class.getResourceAsStream(RESOURCE)) {
      if (input != null) {
        properties.load(input);
      }
    } catch (IOException exception) {
      properties.clear();
    }
    Map<String, Object> info = new LinkedHashMap<>();
    info.put("gitCommit", properties.getProperty("git.commit.id", UNKNOWN));
    info.put("gitBranch", properties.getProperty("git.branch", UNKNOWN));
    info.put("gitDirty", properties.getProperty("git.dirty", UNKNOWN));
    info.put("javaVersion", System.getProperty("java.version"));
    return info;
  }
}
