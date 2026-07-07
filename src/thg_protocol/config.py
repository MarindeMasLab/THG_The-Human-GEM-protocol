"""Configuration helpers for THG protocol workflows."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def get_project_root() -> Path:
    """Return the repository root for the current source checkout."""
    return Path(__file__).resolve().parents[2]


def load_config(config_path: str | Path | None = None) -> dict[str, Any]:
    """Load a JSON configuration file.

    Relative paths are resolved first from the current working directory and
    then from the project root. This preserves the behavior of the original
    ``functions.config`` helper while making it importable from the package.
    """
    if config_path is None:
        resolved_path = get_project_root() / "config.json"
    else:
        resolved_path = Path(config_path)
        if not resolved_path.is_absolute():
            cwd_path = Path.cwd() / resolved_path
            root_path = get_project_root() / resolved_path

            if cwd_path.exists():
                resolved_path = cwd_path
            elif root_path.exists():
                resolved_path = root_path

    if not resolved_path.exists():
        raise FileNotFoundError(
            f"Configuration file not found: {resolved_path}\n"
            f"  Searched in:\n"
            f"    - {resolved_path.absolute()}\n"
            f"  Please ensure the config file exists."
        )

    try:
        with resolved_path.open() as file_handle:
            config = json.load(file_handle)
    except json.JSONDecodeError as error:
        raise json.JSONDecodeError(
            f"Invalid JSON in config file: {resolved_path}",
            error.doc,
            error.pos,
        ) from error

    print(f"  Configuration file: {resolved_path.absolute()}")
    return config


def get_model_paths(config: dict[str, Any] | None = None) -> tuple[str, str]:
    """Return base and output model paths from a configuration mapping."""
    if config is None:
        config = load_config()

    project_root = get_project_root()
    base_model = config.get("model", {}).get("base")
    output_model = config.get("model", {}).get("output")

    if not base_model:
        raise ValueError("'model.base' not found in config.json")
    if not output_model:
        raise ValueError("'model.output' not found in config.json")

    base_path = (
        project_root / base_model
        if not Path(base_model).is_absolute()
        else Path(base_model)
    )
    output_path = (
        project_root / output_model
        if not Path(output_model).is_absolute()
        else Path(output_model)
    )

    return str(base_path), str(output_path)


def get_compartments(config: dict[str, Any] | None = None) -> list[dict[str, str]]:
    """Return compartment specifications from a configuration mapping."""
    if config is None:
        config = load_config()

    compartments = config.get("compartments", [])
    if not compartments:
        raise ValueError("'compartments' not found or empty in config.json")

    for compartment in compartments:
        if "name" not in compartment or "abbreviation" not in compartment:
            raise ValueError(
                f"Invalid compartment definition: {compartment}. "
                "Each compartment must have 'name' and 'abbreviation' fields."
            )

    return compartments


def resolve_compartment_abbreviation(
    model: Any,
    compartment_config: dict[str, str],
) -> dict[str, str | bool | None]:
    """Resolve the compartment abbreviation that should be used for a model."""
    desired_name = compartment_config["name"].lower()
    config_abbrev = compartment_config["abbreviation"]

    if hasattr(model, "compartments"):
        compartments = model.compartments
    elif isinstance(model, dict) and "compartments" in model:
        compartments = model["compartments"]
    else:
        raise ValueError("Invalid model format. Expected COBRApy model or JSON dict.")

    for abbrev, name in compartments.items():
        if name.lower() == desired_name or abbrev.lower() == config_abbrev.lower():
            return {
                "abbreviation": abbrev,
                "exists": True,
                "existing_abbrev": abbrev,
                "name": name,
            }

    return {
        "abbreviation": config_abbrev,
        "exists": False,
        "existing_abbrev": None,
        "name": compartment_config["name"],
    }


def resolve_all_compartments(
    model: Any,
    config: dict[str, Any] | None = None,
) -> list[dict[str, str | bool | None]]:
    """Resolve configured compartment abbreviations against a model."""
    if config is None:
        config = load_config()

    return [
        resolve_compartment_abbreviation(model, compartment)
        for compartment in get_compartments(config)
    ]
