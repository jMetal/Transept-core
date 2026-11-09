package org.transept.algorithm;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.LongSupplier;

/**
 * A request to stop a run early and keep its results, made by creating a file.
 *
 * <p>A user who does not want to wait for the whole budget creates the file (for example, Transept
 * writes {@code STOP} in the run directory). The run then ends at its next check as if the budget
 * were spent, and writes its results as usual. The file is looked for at most every {@link
 * #CHECK_INTERVAL_MILLIS} milliseconds, because some algorithms check their termination after every
 * evaluation.
 */
public final class StopRequest {

  /** Name of the file that requests the stop, inside the run directory. */
  public static final String FILE_NAME = "STOP";

  /** Shortest time between two looks for the file. */
  public static final long CHECK_INTERVAL_MILLIS = 200;

  /** A request that is never made, for runs that cannot be stopped this way. */
  public static final StopRequest NEVER = new StopRequest(null, System::nanoTime);

  private final Path file;
  private final LongSupplier nanoClock;
  private final AtomicBoolean requested = new AtomicBoolean();
  private volatile long lastCheck = Long.MIN_VALUE;

  StopRequest(Path file, LongSupplier nanoClock) {
    this.file = file;
    this.nanoClock = nanoClock;
  }

  /**
   * Returns a request made by creating a file.
   *
   * @param file the file whose existence requests the stop
   * @return the request
   */
  public static StopRequest byFile(Path file) {
    return new StopRequest(file, System::nanoTime);
  }

  /**
   * Tells whether the stop has been requested, looking for the file if enough time has passed.
   * Once requested, it stays requested.
   *
   * @return true if the run must stop
   */
  public boolean isRequested() {
    if (file != null && !requested.get()) {
      long now = nanoClock.getAsLong();
      if (lastCheck == Long.MIN_VALUE || now - lastCheck >= CHECK_INTERVAL_MILLIS * 1_000_000L) {
        lastCheck = now;
        requested.set(Files.exists(file));
      }
    }
    return requested.get();
  }

  /**
   * Tells whether the stop was requested at some check, without looking for the file.
   *
   * @return true if the run stopped (or will stop) at the user's request
   */
  public boolean wasRequested() {
    return requested.get();
  }
}
