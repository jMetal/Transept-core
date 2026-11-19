# Git and Commit Guidelines

This project follows the [Conventional Commits](https://www.conventionalcommits.org/) specification.

## Branches

| Branch | Role |
|------------|------------------------------------------------------|
| `develop`  | Integration branch. Working branches start here and return here through a pull request |
| `main`     | Published state. Only the maintainers bring `develop` into it |
| `feat/...`, `fix/...`, `docs/...` | Working branches, one per change, started from `develop` |

- Start every change on a working branch from an up-to-date `develop`, and open a pull request
  into `develop`. Never commit directly to `main`.
- Keep a pull request to one topic; describe what changes, why, and how it was checked.
- A change of the command line or the output files of Transept core needs a matching pull request
  in the other repository (Transept core or Transept); link them to each other.
- AI assistants push, open pull requests or merge only when the user asks.

## Message format

```
<type>: <short imperative description>
```

## Allowed types

| Type       | When to use                                          |
|------------|------------------------------------------------------|
| `feat`     | A new feature or public method                       |
| `fix`      | A bug fix                                            |
| `test`     | Adding or correcting tests (no production code)      |
| `refactor` | Code change that is neither a fix nor a new feature  |
| `docs`     | Documentation only (`README.md`, `AGENTS.md`, etc.)  |
| `chore`    | Build config, dependencies, `.gitignore`, etc.       |

## Atomic commits

Each commit must represent **one single logical change**. Guidelines:

- If the commit message needs "and" to describe what it does, split it into two commits.
- Ensure the project builds and all checks pass before committing: `mvn clean verify` in Transept
  core, `make lint` and `make test` in Transept.
- Never mix production code changes with test changes in the same commit.
- Never mix code changes with documentation changes in the same commit.

## Examples

```bash
# Good
git commit -m "feat: add isEven method to MathUtils"
git commit -m "test: add tests for MathUtils.isEven"
git commit -m "fix: throw ArithmeticException on divide by zero"
git commit -m "docs: add README.md"

# Bad — too broad, mixes concerns
git commit -m "add stuff and fix tests and update readme"
```
