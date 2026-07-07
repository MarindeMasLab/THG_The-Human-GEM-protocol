# Development

Run the default test set with:

```bash
pytest -m "not slow and not online and not gurobi and not memote"
```

Run linting with:

```bash
ruff check src tests
```
