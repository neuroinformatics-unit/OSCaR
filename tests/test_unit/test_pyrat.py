import pandas as pd

from oscar.pyrat import standardise_pyrat_csv


def test_standardise_single_mutation_csv(
    pyrat_single_mutation_csv_path, standardised_single_mutation_csv_path
):
    pyrat_csv = pd.read_csv(pyrat_single_mutation_csv_path)
    expected_csv = pd.read_csv(standardised_single_mutation_csv_path)

    standard_csv = standardise_pyrat_csv(pyrat_csv)
    pd.testing.assert_frame_equal(
        standard_csv.reset_index(drop=True),
        expected_csv.reset_index(drop=True),
    )


def test_standardise_2_mutations_csv(pyrat_2_mutations_csv_path):
    pyrat_csv = pd.read_csv(pyrat_2_mutations_csv_path)
    # expected_csv = pd.read_csv(standardised_single_mutation_csv_path)

    standardise_pyrat_csv(pyrat_csv)
    # pd.testing.assert_frame_equal(
    #     standard_csv.reset_index(drop=True),
    #     expected_csv.reset_index(drop=True),
    # )
