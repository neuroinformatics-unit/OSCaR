from pathlib import Path

import pooch
import pytest

GIN_REPO = pooch.create(
    path=Path(__file__).parents[1] / "test_data",
    base_url="https://gin.g-node.org/neuroinformatics/oscar-test-data/raw/master/",
    registry=None,
    retry_if_failed=5,
)
GIN_REPO.load_registry(Path(__file__).parent / "pooch_registry.txt")


@pytest.fixture
def single_mutation_csv_path():
    return GIN_REPO.fetch("test-data-single-mutation.csv")
