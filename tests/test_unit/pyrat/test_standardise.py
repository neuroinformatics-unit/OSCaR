import pandas as pd
import pytest

from oscar.pyrat.standardise import standardise_pyrat_csv
from tests.pooch_test_data import pooch_data_path


@pytest.mark.parametrize(
    "pyrat_csv_name, expected_csv_name",
    [
        pytest.param(
            "pyrat-data-single-mutation.csv",
            "standardised-data-single-mutation.csv",
            id="1 mutation",
        ),
        pytest.param(
            "pyrat-data-2-mutations.csv",
            "standardised-data-2-mutations.csv",
            id="2 mutation",
        ),
        pytest.param(
            "pyrat-data-3-mutations.csv",
            "standardised-data-3-mutations.csv",
            id="3 mutation",
        ),
    ],
)
@pytest.mark.parametrize("input_type", ["path", "dataframe"])
def test_standardise_pyrat_csv(pyrat_csv_name, expected_csv_name, input_type):
    """
    Test standardisation of dataframes containing lines with 1, 2 or
    3 mutations.
    """

    pyrat_csv_path = pooch_data_path(pyrat_csv_name)
    expected_csv_path = pooch_data_path(expected_csv_name)

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


def test_standardise_genotypes():
    """
    Test standardisation of a dataframe containing forbidden genotypes (e.g.
    +, -, T, Tg, ko/ko), as well as un-genotyped individuals.
    """

    pyrat_csv = pd.read_csv(
        pooch_data_path("pyrat-data-forbidden-genotypes.csv")
    )
    expected_csv = pd.read_csv(
        pooch_data_path("standardised-data-forbidden-genotypes.csv")
    )

    standard_csv = standardise_pyrat_csv(pyrat_csv)
    pd.testing.assert_frame_equal(
        standard_csv.reset_index(drop=True),
        expected_csv.reset_index(drop=True),
    )
