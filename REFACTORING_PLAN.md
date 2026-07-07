# THG Protocol Refactoring and Packaging Plan

This document captures the current maintainability issues in the repository and
lays out an implementation plan for turning it into a maintainable Python
package with tests, documentation, and large-file handling through Git LFS.

## Goals

- Make the repository installable as a Python package.
- Make reusable workflow logic available through a documented Python API.
- Keep command-line workflows usable during the migration, but treat CLIs as
  thin wrappers around stable Python APIs rather than the primary interface.
- Move reusable logic out of top-level scripts into importable modules.
- Add tests before and during refactoring so behavior changes are visible.
- Move large model/data artifacts to Git LFS or remove generated artifacts from
  git entirely.
- Improve documentation for installation, usage, development, and data/model
  management.

## Current Issues

### Repository Layout

- Source code, scripts, generated models, logs, reports, notebooks, caches, and
  backups are mixed together.
- There is no package metadata: no `pyproject.toml`, `setup.py`, `setup.cfg`, or
  central pytest configuration.
- Several directories are not valid or ideal Python package names:
  - `generate_data-base` contains a hyphen.
  - `merge_metabolic_netowrks_and_network_consistency` is misspelled.
  - `functions` is too generic for a public package namespace.
- There is a stray `functions/__init.py__` next to `functions/__init__.py`.
- Large generated outputs and platform/editor files are present, including
  `.DS_Store`, notebook checkpoints, log files, and Excel lock files.

### Import and Execution Model

- Many scripts modify `sys.path` to import local modules.
- Several scripts assume repository-relative `models/`, `files/`, and `logs/`
  paths.
- Some modules perform heavy work at import time instead of exposing functions
  and using a guarded `main()` entry point.
- There are many wildcard imports from `functions.*`, making dependencies and
  public APIs hard to understand.
- Many scripts write directly into fixed repository folders, which makes tests
  and repeated runs harder to isolate.
- The repository has some script-level CLI support, but not a cohesive installed
  command-line interface. For example, `gapfill/gapfill.py` and
  `implement_pathway/pathway_implementation.py` expose useful `--help` output,
  while other workflows still run work at import/top level or fail before help
  because heavy dependencies are imported before argument parsing.

### Code Maintainability

- Core modules are very large and combine unrelated responsibilities.
- Several workflows mix configuration, I/O, external API calls, model mutation,
  logging, and report generation in the same script.
- There are many print/debug statements and `pdb` imports in production code.
- Dependency versions are split across subdirectories and sometimes conflict
  with the current Python version.
- External services such as BioCyc, KEGG, PubChem, and Ensembl are called
  directly from workflow code, which makes testing unreliable.

### Tests

- `pytest` was not installed in the active environment during review, so current
  test collection could not be verified.
- Existing tests often import copied helper files from test directories instead
  of testing production modules.
- Some tests depend on live external services.
- Some fixtures are very large SBML/XML files.
- There is no central test marker policy for slow, online, solver-dependent, or
  memote tests.

### Large Files

- The repository currently tracks roughly 963 MB of files.
- About 50 tracked files larger than 1 MB account for roughly 944 MB.
- No files are currently tracked by Git LFS.
- `.gitignore` contains entries for `models`, `files`, and `logs`, but many
  files in those directories are already tracked, so the ignore rules do not
  affect them.

## Target Package Structure

Use a `src/` layout so imports are tested against the installed package rather
than accidentally resolving local files.

```text
.
├── pyproject.toml
├── README.md
├── REFACTORING_PLAN.md
├── docs/
│   ├── index.md
│   ├── installation.md
│   ├── usage.md
│   ├── development.md
│   └── data-and-model-files.md
├── src/
│   └── thg_protocol/
│       ├── __init__.py
│       ├── config.py
│       ├── io/
│       ├── annotation/
│       │   ├── metabolites.py
│       │   └── reactions.py
│       ├── gpr/
│       ├── database/
│       ├── model_build/
│       ├── merge/
│       ├── gapfill/
│       ├── cell_specific/
│       ├── pathway/
│       ├── analysis/
│       ├── figures/
│       └── cli/
├── tests/
│   ├── unit/
│   ├── integration/
│   ├── data/
│   └── conftest.py
└── examples/
```

Suggested mapping from existing folders:

| Current path | Target package area |
| --- | --- |
| `functions/function_metabolite_identification.py` | `thg_protocol.annotation.metabolites` |
| `functions/function_reac_identification.py` | `thg_protocol.annotation.reactions` |
| `functions/gpr/` | `thg_protocol.gpr` |
| `generate_data-base/` | `thg_protocol.database` |
| `build_model/` | `thg_protocol.model_build` |
| `merge_metabolic_netowrks_and_network_consistency/` | `thg_protocol.merge` |
| `gapfill/` | `thg_protocol.gapfill` |
| `cell_type_specific_model/` | `thg_protocol.cell_specific` |
| `implement_pathway/` | `thg_protocol.pathway` |
| `compare_models/` | `thg_protocol.analysis.compare` |
| `generate_figures/` | `thg_protocol.figures` |
| `network_analysis/` | `thg_protocol.analysis.network` |
| `utils/` | `thg_protocol.io` or `thg_protocol.utils` |

## Packaging Plan

### Minimal `pyproject.toml`

Start with a conservative package definition and expand dependency groups as the
code is migrated. Do not add console scripts until the corresponding command has
an import-safe module, a stable Python API, argument parsing, and at least a help
smoke test.

```toml
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "thg-protocol"
version = "0.1.0"
description = "Protocol tools for constructing and analyzing human genome-scale metabolic models."
readme = "README.md"
requires-python = ">=3.10,<3.13"
license = { file = "LICENSE" }
authors = [
  { name = "THG contributors" }
]
dependencies = [
  "cobra>=0.26",
  "numpy",
  "pandas",
  "scipy",
  "requests",
  "python-libsbml",
  "tqdm",
]

[project.optional-dependencies]
dev = [
  "pytest",
  "pytest-cov",
  "ruff",
  "black",
]
memote = [
  "memote",
]
solver = [
  "swiglpk",
]
cell-specific = [
  "troppo",
  "pathos",
]
docs = [
  "mkdocs",
  "mkdocs-material",
]

# Add [project.scripts] entries only as commands are migrated.
# Likely first stable commands:
# thg-gapfill = "thg_protocol.gapfill.cli:main"
# thg-compare = "thg_protocol.analysis.compare_cli:main"
# thg-pathway = "thg_protocol.pathway.cli:main"
```

### Python API Principles

- The Python API is the primary interface for reusable logic.
- Core functions should accept explicit inputs such as models, paths, client
  objects, configuration objects, and output directories.
- Core functions should return structured results or write only to explicitly
  provided paths.
- Importing a module should not load large models, write files, open network
  sessions, or start long-running work.
- External services should be isolated behind client modules and injected into
  workflow functions where practical.
- Existing import paths may remain temporarily as thin compatibility wrappers
  during migration, but new code should import from `thg_protocol`.

### CLI Principles

- Add CLI commands for stable, repeatable workflows, not for every helper or
  exploratory script.
- Prefer Python APIs for lower-level annotation, matching, mass-balance, model
  manipulation, and analysis functions.
- Each workflow should expose a `main(argv: list[str] | None = None) -> int`.
- CLI wrappers should only parse arguments, load configuration, call package
  functions, and return an exit code.
- All input/output paths should be arguments or config values.
- No workflow should write to `models/`, `files/`, or `logs/` unless explicitly
  requested.
- `--help` should work without importing optional heavy dependencies or loading
  models.
- Every installed CLI command should have a smoke test for `--help`.
- Good initial CLI candidates are `thg-gapfill`, `thg-compare`, and
  `thg-pathway`. Defer `thg-build` and `thg-generate-db` until their inputs,
  outputs, credentials, and dependency requirements are cleanly parameterized.

## Git LFS and Artifact Plan

### Decide File Ownership

Use three categories:

1. **Source-controlled text**: Python, Markdown, config templates, small JSON
   examples, small TSV/CSV fixtures.
2. **Git LFS artifacts**: canonical model files and large curated datasets that
   are required for reproducible workflows.
3. **Ignored generated artifacts**: logs, temporary outputs, checkpoints, local
   caches, backup folders, generated reports, and exploratory figures.

### Recommended `.gitattributes`

Do not add broad LFS rules until file ownership has been decided. Prefer a
targeted `.gitattributes` generated from the artifact inventory, with canonical
large model/data files tracked by LFS and small fixtures kept as normal text.

A starting point might look like this after ownership is decided. Use explicit
approved paths from the artifact inventory, not broad directory rules:

```gitattributes
models/<approved-canonical-model>.xml filter=lfs diff=lfs merge=lfs -text
memote_and_task_analysis/models/<approved-memote-model>.xml filter=lfs diff=lfs merge=lfs -text
files/<approved-large-dataset>.gct filter=lfs diff=lfs merge=lfs -text
```

Be careful with broad `*.xml`, `*.json`, `*.csv`, `*.xlsx`, or `*.pdf` LFS rules
because the repository also contains small configs, useful text fixtures,
published reports, and documentation assets. Prefer targeted rules for known
large canonical artifacts. Do not use broad rules such as `models/**` unless the
whole directory has been reviewed and approved as canonical LFS content.

### Recommended `.gitignore` Additions

```gitignore
# Python
__pycache__/
*.py[cod]
.pytest_cache/
.coverage
htmlcov/
*.egg-info/
build/
dist/

# Local environment
.env
.venv/
venv/

# Generated outputs
logs/
backups/
*.log
run_full.log
run_full_db.log
run_subset.log

# Local/editor files
.DS_Store
**/.DS_Store
**/.ipynb_checkpoints/
~$*

# Caches and checkpoints
files/caches/
files/*cache*
files/*checkpoint*
*.tmp
```

### LFS Migration Steps

Run this on a dedicated branch. The history rewrite should be coordinated with
all collaborators.

```bash
git checkout -b chore/package-refactor-plan
git lfs install
# After reviewing the artifact inventory, track only approved canonical artifacts.
# Replace the placeholders below with exact paths from the approved inventory.
git lfs track "models/<approved-canonical-model>.xml"
git lfs track "memote_and_task_analysis/models/<approved-memote-model>.xml"
git lfs track "files/<approved-large-dataset>.gct"
git add .gitattributes
git commit -m "Track large model artifacts with Git LFS"
```

If existing history must be cleaned:

```bash
# Replace this example with the reviewed comma-separated path list.
git lfs migrate import --include="models/<approved-canonical-model>.xml,memote_and_task_analysis/models/<approved-memote-model>.xml,files/<approved-large-dataset>.gct"
```

After migration:

```bash
git lfs ls-files
git status
```

## Testing Strategy

### Test Categories

- **Unit tests**: pure functions, small inputs, no network, no large models.
- **Fixture integration tests**: small COBRA models stored in `tests/data`.
- **Workflow smoke tests**: one tiny end-to-end example per CLI.
- **Slow tests**: full or near-full model workflows.
- **Online tests**: BioCyc, KEGG, PubChem, Ensembl.
- **Solver tests**: tests that require GLPK, Gurobi, or another solver.
- **Memote tests**: custom memote/task analysis.

### Pytest Markers

Add central markers in `pyproject.toml`:

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-ra --strict-markers"
markers = [
  "slow: long-running model or workflow test",
  "online: requires external network service",
  "solver: requires an optimization solver",
  "gurobi: requires Gurobi",
  "memote: requires memote",
]
```

Default CI should run:

```bash
pytest -m "not slow and not online and not gurobi and not memote"
```

### Priority Test Coverage

Start with characterization tests before moving code. These tests lock current
behavior and reduce refactoring risk. Each module-move PR should add or preserve
tests for the behavior it moves.

1. GPR parsing and sanitization.
2. Metabolite identification and annotation helpers.
3. Reaction parsing, Jaccard matching, and annotation helpers.
4. Mass-balance equation helpers.
5. Compartment mapping and ID normalization.
6. Pathway builder functions.
7. JSON/SBML load/save helpers.
8. Gapfill phase functions using tiny toy models.
9. Network analysis functions using toy graphs/models.
10. CLI smoke tests with temporary directories.

### Test Data Rules

- Keep unit test fixtures small and text-readable where possible.
- Create tiny COBRA models in Python fixtures when practical.
- Use large full models only in marked integration or slow tests.
- Do not require online services in default tests.
- Mock external API clients and provide recorded/static response fixtures.

## Documentation Plan

Create documentation in `docs/` and keep the README focused on orientation.

Recommended docs:

- `docs/installation.md`: supported Python versions, solver setup, optional
  extras, editable install.
- `docs/usage.md`: Python API examples first, followed by CLI examples for
  stable installed commands.
- `docs/development.md`: test commands, linting, package layout, release steps.
- `docs/data-and-model-files.md`: what lives in Git LFS, what is ignored, where
  to place local/generated files.
- `docs/workflows/`: one page per major pipeline:
  - metabolite/reaction identification
  - database generation
  - model build
  - merge/network consistency
  - gapfill
  - pathway implementation
  - cell-specific models
  - memote/task analysis
  - comparison and figures

## Implementation Plan

### Phase 0: Stabilize the Current State

- Create a migration branch.
- Record the current large-file inventory.
- Record current dependency files and active Python version.
- Map existing requirement files into planned base, dev, solver, memote,
  cell-specific, and workflow-specific extras.
- Audit current script entry points: which scripts have usable `--help`, which
  scripts execute at import time, and which scripts require hard-coded repo
  paths.
- Decide initial CLI candidates and deferred CLI candidates. Start with stable,
  repeatable workflows such as gapfill, compare, and pathway implementation.
- Decide the supported Python range. A practical starting point is
  `>=3.10,<3.13` unless dependency constraints require otherwise.
- Decide initial artifact ownership: canonical LFS artifacts, ignored generated
  files, and small fixtures that stay in normal git.
- Do not rename or move large directories yet.

Deliverables:

- `REFACTORING_PLAN.md`
- Initial issue list in GitHub or a local checklist.
- Dependency inventory and proposed extras.
- Large-file inventory and proposed LFS/ignore policy.
- Agreement on Python version, solver support, and initial CLI scope.

### Phase 1: Add Packaging and Developer Tooling

- Add `pyproject.toml`.
- Add `src/thg_protocol/__init__.py`.
- Add central pytest config.
- Add `ruff` and `black` configuration.
- Add a `tests/` directory with one trivial import test.
- Add the first characterization tests for the safest pure functions before
  moving their code.
- Add `docs/` skeleton.
- Add updated `.gitignore`.
- Do not add `[project.scripts]` entries yet unless the matching CLI module is
  already import-safe and has a help smoke test.

Validation:

```bash
python -m pip install -e ".[dev]"
pytest -m "not slow and not online and not gurobi and not memote"
python -c "import thg_protocol; print(thg_protocol.__version__)"
```

### Phase 2: Move Pure Utilities First

- For each utility, add characterization coverage before or in the same PR as
  the move.
- Move config loading utilities into `thg_protocol.config`.
- Move GPR parsing/sanitization into `thg_protocol.gpr`.
- Move metabolite annotation helpers into `thg_protocol.annotation.metabolites`.
- Move reaction annotation helpers into `thg_protocol.annotation.reactions`.
- Replace wildcard imports with explicit imports.
- Keep old import paths as thin compatibility wrappers during transition when
  existing scripts still depend on them. Add tests for both old and new imports
  until wrappers are removed.

Validation:

```bash
pytest tests/unit
ruff check src tests
```

### Phase 3: Define Workflow APIs and Refactor Script Entry Points

For each major script:

- Move reusable logic into package modules.
- Expose a Python workflow function with explicit inputs and outputs.
- Add `main(argv=None)`.
- Put execution behind `if __name__ == "__main__": raise SystemExit(main())`.
- Replace hard-coded paths with arguments or config values.
- Use `pathlib.Path`.
- Write outputs to a user-provided output directory.
- Ensure `--help` works without optional heavy dependencies or model files.
- Add a `[project.scripts]` entry only after the command is import-safe,
  parameterized, and covered by a help smoke test.

Initial targets:

1. `gapfill/gapfill.py` -> `thg_protocol.gapfill` API plus `thg-gapfill`.
2. `implement_pathway/pathway_implementation.py` -> `thg_protocol.pathway` API
   plus `thg-pathway`.
3. `compare_models/compare_models.py` -> `thg_protocol.analysis.compare` API
   plus `thg-compare`.

Deferred targets until inputs, outputs, credentials, and dependencies are
cleanly parameterized:

1. `build_model/build_model.py`
2. `build_model/build_model_batch.py`
3. `generate_data-base/generate_db.py`
4. `metabolite_reac_identification/metabolite_reac_identification.py`
5. `generate_figures/create_figure.py`

Validation:

```bash
pytest tests/unit tests/integration
thg-compare --help
thg-gapfill --help
thg-pathway --help
```

### Phase 4: Expand Integration and CLI Tests

- Convert existing `test_algorithms/` tests to import production modules.
- Move reusable fixtures into `tests/conftest.py`.
- Replace large fixture models with small toy models where possible.
- Mark tests that still require large models as `slow`.
- Mark tests requiring network as `online`.
- Add CLI tests using `tmp_path` for installed commands.
- Add workflow API tests that call Python functions directly without shelling
  out.

Validation:

```bash
pytest -m "not slow and not online and not gurobi and not memote"
pytest -m slow
```

### Phase 5: External Service Isolation

- Create client modules for BioCyc, KEGG, PubChem, and Ensembl.
- Inject clients into workflow functions instead of calling services directly.
- Add retry, timeout, and cache behavior in one place.
- Add mock/static-response tests for all client-dependent logic.
- Ensure credentials are read only from environment variables or config files
  outside git.

Validation:

```bash
pytest -m "not online"
pytest -m online
```

### Phase 6: Large File Cleanup and Git LFS

- Decide which model/data artifacts are canonical.
- Move canonical large artifacts to Git LFS.
- Remove generated logs, backups, reports, temporary files, and checkpoints from
  git tracking unless they are intentionally published outputs.
- Replace duplicate model copies with one canonical location or documented
  download/LFS paths.
- Keep small test fixtures in normal git.

Validation:

```bash
git lfs ls-files
git status
git ls-files | wc -l
```

### Phase 7: Documentation and CI

- Update `README.md` with a short overview and links to docs.
- Add install instructions for base, dev, solver, memote, and cell-specific
  extras.
- Add examples for the main Python workflow APIs.
- Add CLI examples only for installed, tested commands.
- Add GitHub Actions or equivalent CI:
  - install package
  - run lint
  - run default tests
  - run CLI help smoke tests for installed commands
  - optionally run slow/online tests on manual schedule

Validation:

```bash
python -m pip install -e ".[dev]"
ruff check src tests
pytest -m "not slow and not online and not gurobi and not memote"
```

## Suggested First Pull Requests

1. Add `REFACTORING_PLAN.md`, dependency inventory, large-file inventory, and an
   improved `.gitignore`.
2. Add `pyproject.toml`, `src/thg_protocol`, pytest configuration, and one
   import test without console scripts.
3. Add first characterization tests for pure utilities.
4. Move GPR parsing utilities with compatibility wrappers and tests.
5. Move metabolite/reaction annotation utilities with compatibility wrappers
   and tests.
6. Convert gapfill into a package API plus `thg-gapfill` once help and smoke
   tests pass.
7. Convert pathway implementation or compare models into a package API plus a
   CLI wrapper.
8. Move selected canonical large files to Git LFS after artifact ownership is
   approved.
9. Convert existing tests to use production modules.

## Definition of Done

The refactor can be considered successful when:

- `pip install -e ".[dev]"` works from a clean checkout.
- Reusable logic is available through documented Python APIs.
- Core modules can be imported without loading models, writing files, or calling
  external services.
- Default tests run without network access and without full-size model files.
- Stable, repeatable major workflows are available as thin CLI commands where a
  CLI adds value.
- Every installed CLI command supports `--help` and has a smoke test.
- Large canonical artifacts are tracked with Git LFS or documented as external
  downloads.
- Generated artifacts are ignored or written outside the source tree.
- Documentation explains installation, development, model/data handling, and the
  main Python and CLI workflows.
