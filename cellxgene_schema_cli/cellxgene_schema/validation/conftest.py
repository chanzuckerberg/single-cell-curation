import pytest

from .. import schema


def pytest_addoption(parser):
    parser.addoption("--dataset", action="store", required=True)
    parser.addoption("--ignore-labels", action="store_true", default=False, required=False)


@pytest.hookimpl(tryfirst=True)
def pytest_configure(config: pytest.Config) -> None:
    pytest.cxg_schema_version = schema.get_current_schema_version()
    pytest.cxg_schema_def = schema.get_schema_definition()
    pytest.cxg_ignore_labels = config.getoption("--ignore-labels")
