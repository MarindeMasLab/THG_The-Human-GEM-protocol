# Tests

The first tests in this repository are characterization tests. They describe
how selected current production modules behave before the package layout is
changed.

Run them from the repository root without installing the project:

```bash
python -m unittest discover -s tests
```

These tests intentionally import current modules such as `functions.config`.
During the refactor, keep the expected behavior stable while moving the
implementation into `thg_protocol` modules or compatibility wrappers.
