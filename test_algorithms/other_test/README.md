# Unit Tests

This directory contains unit tests for the THG model generation code.

## Setup

Install test dependencies:
```bash
pip install pytest pytest-cov dill
```

## Running Tests

Run all tests:
```bash
pytest test_algorithms/other_test/ -v
```

Run specific test file:
```bash
pytest test_algorithms/other_test/test_checkpoint.py -v
```

Run specific test class:
```bash
pytest test_algorithms/other_test/test_checkpoint.py::TestReactionPickle -v
```

Run specific test:
```bash
pytest test_algorithms/other_test/test_checkpoint.py::TestReactionPickle::test_reaction_pickle_no_lambda_error -v
```

Run with coverage:
```bash
pytest test_algorithms/other_test/ --cov=generate_data-base --cov=functions --cov-report=html
```

## Test Files

### test_checkpoint.py

Tests for checkpoint save/load functionality, including:
- **TestReactionPickle**: Verifies reactions with MethodType methods pickle correctly
- **TestCheckpointDataStructure**: Tests complete checkpoint data structure
- **TestSanitizeLoadedReactions**: Tests backward compatibility sanitization
- **TestCompartmentalizedReactions**: Tests compartmentalized reactions (_c, _n, etc.)
- **TestFailedLocationReactions**: Tests filtering of failed location reactions
- **TestOperatorAddReplacement**: Tests operator.add for list concatenation
- **TestCheckpointBackwardCompatibility**: Tests loading old checkpoint formats

## Test Coverage Goals

- Checkpoint save/load: 100%
- Reaction pickle/unpickle: 100%
- Data integrity after reload: 100%
- Sanitization functions: 100%

## Continuous Integration

These tests should be run automatically on:
- Every commit
- Every pull request
- Before releases

## Adding New Tests

When adding new tests:
1. Follow the existing test structure
2. Use descriptive test names starting with `test_`
3. Include docstrings explaining what is being tested
4. Use `assert` statements with clear error messages
5. Clean up temporary files in `finally` blocks
