import pandas as pd
import pytest

from oscar_colony.optimise.surplus_summary import create_surplus_summary
from tests.helpers import assert_dataclass_equal
from tests.pooch_test_data import pooch_data_path


@pytest.mark.parametrize(
    (
        "required_n_per_genotype, n_matings_per_scheme, offspring_per_scheme,"
        "expected_surplus"
    ),
    [
        pytest.param(
            "required_n_per_genotype_1_mutation",
            "n_matings_1_mutation",
            "offspring_per_scheme_1_mutation",
            "surplus_1_mutation",
            id="1 mutation",
        ),
        pytest.param(
            "required_n_per_genotype_2_mutations",
            "n_matings_2_mutations",
            "offspring_per_scheme_2_mutations",
            "surplus_2_mutations",
            id="2 mutations",
        ),
    ],
)
def test_create_surplus_summary(
    required_n_per_genotype,
    n_matings_per_scheme,
    offspring_per_scheme,
    expected_surplus,
    request,
):
    surplus_summary = create_surplus_summary(
        request.getfixturevalue(required_n_per_genotype),
        request.getfixturevalue(n_matings_per_scheme),
        request.getfixturevalue(offspring_per_scheme),
    )

    assert_dataclass_equal(
        surplus_summary, request.getfixturevalue(expected_surplus), abs=1e-4
    )


def test_create_genotype_df(surplus_2_mutations):
    """Test creation of a dataframe from SurplusSummary.surplus_per_genotype"""

    genotype_df = surplus_2_mutations.create_genotype_df(decimal_places=2)
    expected_df = pd.read_csv(
        pooch_data_path("converted-surplus-per-genotype.csv")
    )

    pd.testing.assert_frame_equal(genotype_df, expected_df)
