import pandas as pd
import pytest

from oscar.pyrat import standardise_pyrat_csv


@pytest.mark.parametrize(
    "pyrat_csv_path, expected_csv_path",
    [
        pytest.param(
            "pyrat_single_mutation_csv_path",
            "standardised_single_mutation_csv_path",
            id="1 mutation",
        ),
        pytest.param(
            "pyrat_2_mutations_csv_path",
            "standardised_2_mutations_csv_path",
            id="2 mutation",
        ),
        pytest.param(
            "pyrat_3_mutations_csv_path",
            "standardised_3_mutations_csv_path",
            id="3 mutation",
        ),
    ],
)
def test_standardise_pyrat_csv(pyrat_csv_path, expected_csv_path, request):
    pyrat_csv = pd.read_csv(request.getfixturevalue(pyrat_csv_path))
    expected_csv = pd.read_csv(request.getfixturevalue(expected_csv_path))

    standard_csv = standardise_pyrat_csv(pyrat_csv)
    pd.testing.assert_frame_equal(
        standard_csv.reset_index(drop=True),
        expected_csv.reset_index(drop=True),
    )
