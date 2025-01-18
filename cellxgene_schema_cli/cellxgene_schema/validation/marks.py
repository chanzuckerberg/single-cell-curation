import pytest

ignoring_labels = pytest.mark.skipif(pytest.cxg_ignore_labels, reason="Skipping because --ignore-labels is set")
