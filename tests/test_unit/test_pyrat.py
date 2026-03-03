import pandas as pd

from oscar.pyrat import standardise_pyrat_csv


def test_standardise_pyrat_csv(single_mutation_csv_path):
    pyrat_csv = pd.read_csv(single_mutation_csv_path)

    standardise_pyrat_csv(pyrat_csv)
