import pandas as pd

from oscar.pyrat import standardise_pyrat_csv


def test_standardise_pyrat_csv():
    pyrat_csv = pd.read_csv("data/test-data-single-mutation.csv")

    standardise_pyrat_csv(pyrat_csv)
