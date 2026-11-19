---
name: check
description: Run every check of Transept core (tests, Checkstyle, coverage) before a commit or a pull request, and report what fails. Use before committing, before opening a pull request, or when asked whether the project builds.
---

# Check Transept core

Run the checks that a commit and a pull request must pass. Do not commit, push or change code
here: only run and report.

1. Make sure you are in the root of Transept core (`pom.xml` present) and note the branch
   (`git branch --show-current`). If it is `main`, say so: changes go on a branch from `develop`.
2. Run the full build:

   ```bash
   mvn -q clean verify
   ```

   It runs the unit tests (`*Test`), the integration tests (`*IT`), Checkstyle (Google style) and
   JaCoCo (at least 80 % line coverage).
3. If the benchmark data is next to the repository, also run the integration tests on real data:

   ```bash
   test -d ../Dataset_Structural_GPCRs/GPCRdb && \
     GPCR_DATASET=../Dataset_Structural_GPCRs/GPCRdb mvn -q verify
   ```

   If it is not there, say that these tests were skipped.
4. Report:
   - "All checks pass", with the number of tests run, if everything passed;
   - otherwise, for every failure: the test or rule, the file and line, and the message. Group
     Checkstyle violations by file. Quote the relevant output; do not paraphrase it.
5. Never call the change ready while a check fails. If `BaselineOracleIT` fails after a deliberate
   change to an objective, the snapshot must be regenerated in the same change (see `AGENTS.md`,
   pitfall 4); otherwise the failure is a regression.
