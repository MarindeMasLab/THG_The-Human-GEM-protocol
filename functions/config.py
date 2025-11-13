#!/usr/bin/env python3
"""
Configuration loader for cytoskeleton implementation.

Loads settings from config.json in the project root.
"""

import json
import os
from pathlib import Path


def get_project_root():
    """Get the project root directory (where config.json is located)."""
    # This file is in src/, so parent is project root
    return Path(__file__).parent.parent


def load_config(config_path=None):
    """
    Load configuration from a JSON file.
    
    Args:
        config_path: Path to config file (str or Path). 
                    If None, uses 'config.json' in project root.
                    Can be absolute or relative to current directory.
        
    Returns:
        dict: Configuration dictionary
        
    Raises:
        FileNotFoundError: If config file doesn't exist
        json.JSONDecodeError: If config file is invalid JSON
    """
    if config_path is None:
        config_path = get_project_root() / "config.json"
    else:
        config_path = Path(config_path)
        # If relative path, try relative to current directory first
        if not config_path.is_absolute():
            # Try relative to current working directory
            cwd_path = Path.cwd() / config_path
            # Try relative to project root
            root_path = get_project_root() / config_path
            
            if cwd_path.exists():
                config_path = cwd_path
            elif root_path.exists():
                config_path = root_path
            # else keep as-is and let it fail with clear error message
    
    if not config_path.exists():
        raise FileNotFoundError(
            f"Configuration file not found: {config_path}\n"
            f"  Searched in:\n"
            f"    - {config_path.absolute()}\n"
            f"  Please ensure the config file exists."
        )
    
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
    except json.JSONDecodeError as e:
        raise json.JSONDecodeError(
            f"Invalid JSON in config file: {config_path}",
            e.doc, e.pos
        )
    
    print(f"  Configuration file: {config_path.absolute()}")
    return config


def get_model_paths(config=None):
    """
    Get model paths from configuration.
    
    Args:
        config: Optional configuration dict. If None, loads from config.json
        
    Returns:
        tuple: (base_model_path, output_model_path)
    """
    if config is None:
        config = load_config()
    
    project_root = get_project_root()
    
    # Get model paths and resolve relative to project root
    base_model = config.get('model', {}).get('base')
    output_model = config.get('model', {}).get('output')
    
    if not base_model:
        raise ValueError("'model.base' not found in config.json")
    if not output_model:
        raise ValueError("'model.output' not found in config.json")
    
    # Convert to absolute paths if they're relative
    base_path = project_root / base_model if not Path(base_model).is_absolute() else Path(base_model)
    output_path = project_root / output_model if not Path(output_model).is_absolute() else Path(output_model)
    
    return str(base_path), str(output_path)


def get_compartments(config=None):
    """
    Get compartment specifications from configuration.
    
    Args:
        config: Optional configuration dict. If None, loads from config.json
        
    Returns:
        list: List of compartment dicts with 'name' and 'abbreviation' keys
    """
    if config is None:
        config = load_config()
    
    compartments = config.get('compartments', [])
    
    if not compartments:
        raise ValueError("'compartments' not found or empty in config.json")
    
    # Validate compartment structure
    for comp in compartments:
        if 'name' not in comp or 'abbreviation' not in comp:
            raise ValueError(
                f"Invalid compartment definition: {comp}. "
                "Each compartment must have 'name' and 'abbreviation' fields."
            )
    
    return compartments


def resolve_compartment_abbreviation(model, compartment_config):
    """
    Resolve compartment abbreviation for a model.
    
    Checks if a compartment with the desired name already exists in the model.
    If it exists, returns the existing abbreviation.
    If not, returns the abbreviation from config.
    
    Args:
        model: Model dict (JSON format) or COBRApy model
        compartment_config: Compartment dict with 'name' and 'abbreviation'
        
    Returns:
        dict: {
            'abbreviation': str,  # The abbreviation to use
            'exists': bool,       # Whether compartment already exists
            'existing_abbrev': str or None  # Existing abbreviation if found
        }
    """
    desired_name = compartment_config['name'].lower()
    config_abbrev = compartment_config['abbreviation']
    
    # Handle both JSON dict and COBRApy model
    if hasattr(model, 'compartments'):
        # COBRApy model
        compartments = model.compartments
    elif isinstance(model, dict) and 'compartments' in model:
        # JSON dict format
        compartments = model['compartments']
    else:
        raise ValueError("Invalid model format. Expected COBRApy model or JSON dict.")
    
    # Search for existing compartment by name
    for abbrev, name in compartments.items():
        if name.lower() == desired_name or abbrev.lower() == config_abbrev.lower():
            return {
                'abbreviation': abbrev,
                'exists': True,
                'existing_abbrev': abbrev,
                'name': name
            }
    
    # Compartment doesn't exist
    return {
        'abbreviation': config_abbrev,
        'exists': False,
        'existing_abbrev': None,
        'name': compartment_config['name']
    }


def resolve_all_compartments(model, config=None):
    """
    Resolve abbreviations for all compartments in config.
    
    Args:
        model: Model dict (JSON format) or COBRApy model
        config: Optional configuration dict. If None, loads from config.json
        
    Returns:
        list: List of resolved compartment dicts with keys:
            - name: compartment name
            - abbreviation: abbreviation to use
            - exists: whether it already exists in model
            - existing_abbrev: existing abbreviation if found
    """
    if config is None:
        config = load_config()
    
    compartments = get_compartments(config)
    resolved = []
    
    for comp_config in compartments:
        resolved_comp = resolve_compartment_abbreviation(model, comp_config)
        resolved.append(resolved_comp)
    
    return resolved


if __name__ == "__main__":
    # Test config loading
    try:
        config = load_config()
        print("✓ Configuration loaded successfully")
        print(f"\nConfig contents:")
        print(json.dumps(config, indent=2))
        
        base, output = get_model_paths(config)
        print(f"\nModel paths:")
        print(f"  Base model:   {base}")
        print(f"  Output model: {output}")
        
        # Check if files exist
        if Path(base).exists():
            print(f"  ✓ Base model exists")
        else:
            print(f"  ✗ Base model not found")
        
        # Test compartment resolution
        print(f"\nCompartments from config:")
        compartments = get_compartments(config)
        for comp in compartments:
            print(f"  - {comp['name']} [{comp['abbreviation']}]")
        
        # Test compartment resolution with actual model
        if Path(base).exists():
            print(f"\nResolving compartments against base model...")
            with open(base, 'r') as f:
                model = json.load(f)
            
            resolved = resolve_all_compartments(model, config)
            for res in resolved:
                if res['exists']:
                    print(f"  ✓ '{res['name']}' already exists as [{res['abbreviation']}]")
                else:
                    print(f"  + '{res['name']}' will be created as [{res['abbreviation']}]")
            
    except Exception as e:
        print(f"✗ Error: {e}")
