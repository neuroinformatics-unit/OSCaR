import logging

import pytest

from oscar_colony.logging_setup import init_logging


def init_test_log(tmp_path):
    init_logging(output_dir=tmp_path)
    logging.getLogger(__name__)
    log_file = next(tmp_path.glob("*.log"))
    return log_file


@pytest.fixture
def test_log_folder_created(log_file):
    assert log_file.exists()


def test_file_name_is_correct(): ...
