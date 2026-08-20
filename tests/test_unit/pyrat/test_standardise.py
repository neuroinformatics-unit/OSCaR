import logging

import pandas as pd
import pytest

from oscar_colony.colony_management.pyrat.standardise import (
    WildtypeSymbol,
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
    "pyrat_csv_name, expected_csv_name, expected_logs, wt_plus_or_minus",
    [
        pytest.param(
            "pyrat-data-forbidden-genotypes.csv",
            "standardised-data-forbidden-genotypes.csv",
            [
                "Filtered out 6 invalid genotype row(s) for these offspring "
                "IDs : ['ID-002', 'ID-003', 'ID-004', 'ID-005', 'ID-010', "
                "'ID-016'] - "
                "10 remaining",
                "3 offspring have no genotype recorded: "
                "['ID-012', 'ID-013', 'ID-014']",
            ],
            None,
            id="forbidden genotypes",
        ),
        pytest.param(
            "pyrat-data-forbidden-genotypes.csv",
            "standardised-data-forbidden-genotypes-plus-wt.csv",
            [
                "Filtered out 2 invalid genotype row(s) for these offspring "
                "IDs : ['ID-003', 'ID-016'] - "
                "14 remaining",
                "3 offspring have no genotype recorded: "
                "['ID-012', 'ID-013', 'ID-014']",
            ],
            WildtypeSymbol.PLUS,
            id="forbidden genotypes - plus is wt",
        ),
        pytest.param(
            "pyrat-data-forbidden-genotypes.csv",
            "standardised-data-forbidden-genotypes-minus-wt.csv",
            [
                "Filtered out 2 invalid genotype row(s) for these offspring "
                "IDs : ['ID-003', 'ID-016'] - "
                "14 remaining",
                "3 offspring have no genotype recorded: "
                "['ID-012', 'ID-013', 'ID-014']",
            ],
            WildtypeSymbol.MINUS,
            id="forbidden genotypes - minus is wt",
        ),
        pytest.param(
            "pyrat-data-forbidden-schemes.csv",
            "standardised-data-forbidden-schemes.csv",
            [
                "Filtered out 3 row(s) with invalid breeding data for these "
                "offspring IDs: ['ID-008', 'ID-009', 'ID-010']"
                " - 7 remaining",
            ],
            None,
            id="forbidden schemes",
        ),
        pytest.param(
            "pyrat-data-forbidden-parents.csv",
            "standardised-data-forbidden-parents.csv",
            [
                "Filtered out 5 row(s) with invalid breeding data for these "
                "offspring IDs: ['ID-001', 'ID-002', 'ID-003', 'ID-004', "
                "'ID-005'] - 1 remaining",
            ],
            None,
            id="forbidden parents",
        ),
    ],
)
@pytest.mark.parametrize("input_type", ["path", "dataframe"])
def test_forbidden_data(
    pyrat_csv_name,
    expected_csv_name,
    expected_logs,
    wt_plus_or_minus,
    input_type,
    caplog,
):
    """
    Test forbidden data combinations, to make sure they are correctly filtered.

    1. Test standardisation of dataframe containing forbidden genotypes (e.g.
    +, -, T, Tg, ko/ko), as well as un-genotyped individuals.

    2. Test standardisation of dataframe containing forbidden genotypes where
    the user has specified that + is WT (and - is HOM)

    3. Test standardisation of dataframe containing forbidden genotypes where
    the user has specified that - is WT (and + is HOM)

    4. Test that impossible breeding schemes are removed from raw data.
    (e.g. hom x hom parents cannot make wt offspring)

    5. Test that impossible parent schemes are removed. Cases where there are
    one parent or no parents.
    """

    pyrat_csv_path = pooch_data_path(pyrat_csv_name)
    expected_csv_path = pooch_data_path(expected_csv_name)

    pyrat_csv = pd.read_csv(pyrat_csv_path)
    expected_csv = pd.read_csv(expected_csv_path)

    with caplog.at_level(logging.INFO):
        if input_type == "path":
            standard_csv = standardise_pyrat_csv(
                pyrat_csv_path, wt_plus_or_minus=wt_plus_or_minus
            )
        else:
            standard_csv = standardise_pyrat_csv(
                pyrat_csv, wt_plus_or_minus=wt_plus_or_minus
            )

    pd.testing.assert_frame_equal(
        standard_csv.reset_index(drop=True),
        expected_csv.reset_index(drop=True),
    )

    for expected_log in expected_logs:
        assert expected_log in caplog.text


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
    """
    Test log statements are recorded correctly, for a dataframe with no
    invalid data.
    """

    pyrat_csv_path = pooch_data_path("pyrat-data-3-mutations.csv")

    with caplog.at_level(logging.INFO):
        standardise_pyrat_csv(pyrat_csv_path)

    expected_messages = [
        "Starting standardisation of pyRAT data: 29 rows",
        "Standardisation complete: 29 rows",
    ]
    for message in expected_messages:
        assert message in caplog.text
