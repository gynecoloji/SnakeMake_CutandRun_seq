import yaml
from pathlib import Path
from jsonschema import validate


def _load(p):
    return yaml.safe_load(Path(p).read_text())


def test_default_config_matches_schema():
    schema = _load("workflow/schemas/config.schema.yaml")
    config = _load("config/config.yaml")
    validate(instance=config, schema=schema)


def test_test_config_matches_schema():
    schema = _load("workflow/schemas/config.schema.yaml")
    config = _load(".test/config/config.yaml")
    validate(instance=config, schema=schema)
