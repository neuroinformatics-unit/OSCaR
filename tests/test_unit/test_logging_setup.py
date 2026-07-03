import logging
from pathlib import Path

import pytest

from oscar_colony.logging_setup import init_logging


@pytest.fixture
def log_file(tmp_path):
    init_logging(output_dir=tmp_path, filename=Path(__file__).stem)
    logging.getLogger(__name__)
    return next(tmp_path.glob("*.log"))


def test_log_folder_created(log_file):
    assert log_file.exists()


def test_file_name_is_correct(log_file):
    assert log_file.name.startswith("test_logging_setup")
    assert log_file.suffix == ".log"


def test_file_writes_file(log_file):
    logging.info("test")
    with open(log_file, "r") as f:
        contents = f.readlines()
    assert "test" in contents[-1]
