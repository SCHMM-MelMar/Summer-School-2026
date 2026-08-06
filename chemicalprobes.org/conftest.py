"""Registers the one custom marker, so the integration tests can be deselected.

    pytest chemicalprobes.org -m "not slow"     # skip the ones that read the export
"""


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: reads the 1247-probe export and builds a database"
    )
