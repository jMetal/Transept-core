package org.transept.run;

import java.util.Iterator;
import java.util.Locale;
import java.util.Map;

/**
 * Minimal JSON serialiser for run metadata: maps with string keys, iterables, strings, numbers,
 * booleans and {@code null}. Map iteration order is preserved, so use ordered maps for stable
 * output.
 */
public final class Json {

  private Json() {}

  /**
   * Serialises a value on a single line.
   *
   * @param value map, iterable, string, number, boolean, enum or {@code null}
   * @return the JSON text
   * @throws IllegalArgumentException for unsupported types or non-finite numbers
   */
  public static String write(Object value) {
    StringBuilder json = new StringBuilder();
    append(json, value);
    return json.toString();
  }

  private static void append(StringBuilder json, Object value) {
    switch (value) {
      case null -> json.append("null");
      case String text -> appendString(json, text);
      case Enum<?> constant -> appendString(json, constant.name());
      case Boolean flag -> json.append(flag);
      case Integer number -> json.append(number);
      case Long number -> json.append(number);
      case Number number -> appendNumber(json, number.doubleValue());
      case Map<?, ?> map -> appendMap(json, map);
      case Iterable<?> items -> appendIterable(json, items);
      default ->
          throw new IllegalArgumentException(
              "Unsupported JSON value type: " + value.getClass().getName());
    }
  }

  private static void appendNumber(StringBuilder json, double number) {
    if (!Double.isFinite(number)) {
      throw new IllegalArgumentException("JSON numbers must be finite: " + number);
    }
    json.append(
        number == Math.rint(number) && Math.abs(number) < 1e15
            ? String.format(Locale.ROOT, "%.1f", number)
            : Double.toString(number));
  }

  private static void appendMap(StringBuilder json, Map<?, ?> map) {
    json.append('{');
    Iterator<? extends Map.Entry<?, ?>> entries = map.entrySet().iterator();
    while (entries.hasNext()) {
      Map.Entry<?, ?> entry = entries.next();
      appendString(json, String.valueOf(entry.getKey()));
      json.append(':');
      append(json, entry.getValue());
      if (entries.hasNext()) {
        json.append(',');
      }
    }
    json.append('}');
  }

  private static void appendIterable(StringBuilder json, Iterable<?> items) {
    json.append('[');
    Iterator<?> iterator = items.iterator();
    while (iterator.hasNext()) {
      append(json, iterator.next());
      if (iterator.hasNext()) {
        json.append(',');
      }
    }
    json.append(']');
  }

  private static void appendString(StringBuilder json, String text) {
    json.append('"');
    for (int index = 0; index < text.length(); index++) {
      char character = text.charAt(index);
      switch (character) {
        case '"' -> json.append("\\\"");
        case '\\' -> json.append("\\\\");
        case '\n' -> json.append("\\n");
        case '\r' -> json.append("\\r");
        case '\t' -> json.append("\\t");
        default -> {
          if (character < 0x20) {
            json.append(String.format(Locale.ROOT, "\\u%04x", (int) character));
          } else {
            json.append(character);
          }
        }
      }
    }
    json.append('"');
  }
}
