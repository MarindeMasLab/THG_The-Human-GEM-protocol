import unittest
from pathlib import Path

from functions import config


class ConfigCharacterizationTests(unittest.TestCase):
    def test_get_model_paths_resolves_relative_paths_from_project_root(self):
        root = config.get_project_root()
        base, output = config.get_model_paths(
            {
                "model": {
                    "base": "models/input.json",
                    "output": "models/output.json",
                }
            }
        )

        self.assertEqual(base, str(root / "models/input.json"))
        self.assertEqual(output, str(root / "models/output.json"))

    def test_get_model_paths_preserves_absolute_paths(self):
        base_path = Path("/tmp/base.json")
        output_path = Path("/tmp/output.json")

        base, output = config.get_model_paths(
            {
                "model": {
                    "base": str(base_path),
                    "output": str(output_path),
                }
            }
        )

        self.assertEqual(base, str(base_path))
        self.assertEqual(output, str(output_path))

    def test_get_model_paths_requires_base_and_output(self):
        with self.assertRaisesRegex(ValueError, "model.base"):
            config.get_model_paths({"model": {"output": "models/output.json"}})

        with self.assertRaisesRegex(ValueError, "model.output"):
            config.get_model_paths({"model": {"base": "models/input.json"}})

    def test_get_compartments_validates_required_fields(self):
        compartments = config.get_compartments(
            {
                "compartments": [
                    {"name": "cytosol", "abbreviation": "c"},
                    {"name": "mitochondrion", "abbreviation": "m"},
                ]
            }
        )

        self.assertEqual(
            compartments,
            [
                {"name": "cytosol", "abbreviation": "c"},
                {"name": "mitochondrion", "abbreviation": "m"},
            ],
        )

        with self.assertRaisesRegex(ValueError, "Invalid compartment definition"):
            config.get_compartments({"compartments": [{"name": "cytosol"}]})

    def test_resolve_compartment_abbreviation_reuses_existing_compartment(self):
        model = {"compartments": {"c": "cytosol", "m": "mitochondrion"}}

        resolved = config.resolve_compartment_abbreviation(
            model, {"name": "Cytosol", "abbreviation": "cyto"}
        )

        self.assertEqual(
            resolved,
            {
                "abbreviation": "c",
                "exists": True,
                "existing_abbrev": "c",
                "name": "cytosol",
            },
        )

    def test_resolve_compartment_abbreviation_uses_config_for_new_compartment(self):
        model = {"compartments": {"c": "cytosol"}}

        resolved = config.resolve_compartment_abbreviation(
            model, {"name": "extracellular", "abbreviation": "e"}
        )

        self.assertEqual(
            resolved,
            {
                "abbreviation": "e",
                "exists": False,
                "existing_abbrev": None,
                "name": "extracellular",
            },
        )


if __name__ == "__main__":
    unittest.main()
