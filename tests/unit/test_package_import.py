import thg_protocol
from thg_protocol import config


def test_package_exposes_version():
    assert thg_protocol.__version__ == "0.1.0"


def test_config_module_resolves_project_root():
    assert (config.get_project_root() / "REFACTORING_PLAN.md").exists()
