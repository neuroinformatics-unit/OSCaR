from pathlib import Path

import pooch

GIN_REPO = pooch.create(
    path=Path(__file__).parents[1] / "test_data",
    base_url="https://gin.swc.ucl.ac.uk/neuroinformatics/oscar-test-data/src/hc/specify-offspring-sex-standardisation",
    registry=None,
    retry_if_failed=5,
)
GIN_REPO.load_registry(Path(__file__).parent / "pooch_registry.txt")


def pooch_data_path(filename: str) -> Path:
    """Get path of test data file managed by pooch"""

    return GIN_REPO.fetch(filename)
