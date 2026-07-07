# Usage

Reusable workflow code will be exposed through the `thg_protocol` Python API.

Existing top-level scripts remain available during the migration. New command
line entry points will be added only after their workflow modules are
import-safe, parameterized, and covered by help smoke tests.
