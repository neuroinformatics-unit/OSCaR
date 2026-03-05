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
@pytest.mark.parametrize("input_type", ["path", "dataframe"])
def test_standardise_pyrat_csv(
    pyrat_csv_path, expected_csv_path, input_type, request
):
    """
    Test standardisation of dataframes containing lines with 1, 2 or
    3 mutations.
    """

    pyrat_csv_path = request.getfixturevalue(pyrat_csv_path)
    expected_csv_path = request.getfixturevalue(expected_csv_path)

    pyrat_csv = pd.read_csv(pyrat_csv_path)
    expected_csv = pd.read_csv(expected_csv_path)

    if input_type == "path":
        standard_csv = standardise_pyrat_csv(pyrat_csv_path)
    else:
        standard_csv = standardise_pyrat_csv(pyrat_csv)

    pd.testing.assert_frame_equal(
        standard_csv.reset_index(drop=True),
        expected_csv.reset_index(drop=True),
    )


def test_standardise_genotypes(
    pyrat_forbidden_genotypes_csv_path,
    standardised_forbidden_genotypes_csv_path,
):
    """
    Test standardisation of a dataframe containing forbidden genotypes e.g.
    +, -, T, Tg, ko/ko.
    """

    pyrat_csv = pd.read_csv(pyrat_forbidden_genotypes_csv_path)
    expected_csv = pd.read_csv(standardised_forbidden_genotypes_csv_path)

    standard_csv = standardise_pyrat_csv(pyrat_csv)
    pd.testing.assert_frame_equal(
        standard_csv.reset_index(drop=True),
        expected_csv.reset_index(drop=True),
    )
