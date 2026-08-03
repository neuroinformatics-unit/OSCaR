import logging

import pytest

from oscar_colony.logging_setup import init_logging

logger = logging.getLogger(__name__)


@pytest.fixture
def log_file(tmp_path):
    init_logging(output_dir=tmp_path, filename="TEST_LOG")
    return next(tmp_path.glob("*.log"))


def test_file_name_is_correct(log_file):
    assert log_file.name.startswith("TEST_LOG")
    assert log_file.suffix == ".log"


def test_log_writes_to_file(log_file):
    logger.info("test")
    with open(log_file, "r") as f:
        contents = f.readlines()
    assert "test" in contents[-1]
