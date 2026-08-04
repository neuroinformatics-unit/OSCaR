import logging

import pandas as pd
import pytest

from oscar_colony.colony_management.pyrat.standardise import (
    standardise_pyrat_csv,
)
from tests.pooch_test_data import pooch_data_path


@pytest.mark.parametrize(
    "pyrat_csv_name, expected_csv_name",
    [
        pytest.param(
            "pyrat-data-1-mutation.csv",
            "standardised-data-1-mutation.csv",
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


@pytest.mark.parametrize(
    "pyrat_csv_name, expected_csv_name, expected_log",
    [
        pytest.param(
            "pyrat-data-forbidden-genotypes.csv",
            "standardised-data-forbidden-genotypes.csv",
            "Filtered out 5 invalid genotype row(s)",
            id="forbidden genotypes",
        ),
        pytest.param(
            "pyrat-data-forbidden-schemes.csv",
            "standardised-data-forbidden-schemes.csv",
            "invalid breeding data",
            id="forbidden schemes",
        ),
        pytest.param(
            "pyrat-data-forbidden-parents.csv",
            "standardised-data-forbidden-parents.csv",
            "invalid breeding data",
            id="forbidden parents",
        ),
    ],
)
@pytest.mark.parametrize("input_type", ["path", "dataframe"])
def test_forbidden_data(
    pyrat_csv_name, expected_csv_name, expected_log, input_type, caplog
):
    """
    Test forbidden data combinations, to make sure they are correctly filtered.

    1. Test standardisation of a dataframe containing forbidden genotypes (e.g.
    +, -, T, Tg, ko/ko), as well as un-genotyped individuals.

    2. Test that impossible breeding schemes are removed from raw data.
    (e.g. hom x hom parents cannot make wt offspring)

    3. Test that impossible parent schemes are removed. Cases where there are
    one parent or no parents.
    """

    pyrat_csv_path = pooch_data_path(pyrat_csv_name)
    expected_csv_path = pooch_data_path(expected_csv_name)

    pyrat_csv = pd.read_csv(pyrat_csv_path)
    expected_csv = pd.read_csv(expected_csv_path)

    with caplog.at_level(logging.INFO):
        if input_type == "path":
            standard_csv = standardise_pyrat_csv(pyrat_csv_path)
        else:
            standard_csv = standardise_pyrat_csv(pyrat_csv)

    pd.testing.assert_frame_equal(
        standard_csv.reset_index(drop=True),
        expected_csv.reset_index(drop=True),
    )

    log_messages = [record.getMessage() for record in caplog.records]
    assert any(expected_log in message for message in log_messages)


@pytest.mark.parametrize(
    "pyrat_csv_name, expected_csv_name",
    [
        pytest.param(
            "pyrat-data-multiple-parents-1-mutation.csv",
            "standardised-data-multiple-parents-1-mutation.csv",
            id="1 mutation",
        ),
        pytest.param(
            "pyrat-data-multiple-parents-2-mutations.csv",
            "standardised-data-multiple-parents-2-mutations.csv",
            id="2 mutation",
        ),
        pytest.param(
            "pyrat-data-multiple-parents-3-mutations.csv",
            "standardised-data-multiple-parents-3-mutations.csv",
            id="3 mutation",
        ),
    ],
)
@pytest.mark.parametrize("input_type", ["path", "dataframe"])
def test_standardise_multiple_parents_pyrat_csv(
    pyrat_csv_name, expected_csv_name, input_type
):
    """
    Test standardisation of dataframes containing multiple parents
    with 1, 2 or 3 mutations.
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


def test_standardise_pyrat_csv_logs(caplog):

    pyrat_csv_path = pooch_data_path(
        "pyrat-data-multiple-parents-3-mutations.csv"
    )

    with caplog.at_level(logging.INFO):
        standardise_pyrat_csv(pyrat_csv_path)

    log_messages = [record.getMessage() for record in caplog.records]
    assert any(
        "Starting standardisation of pyRAT data" in message
        for message in log_messages
    )
    assert any(
        "Standardisation complete:" in message for message in log_messages
    )
