import pandas as pd
import pytest

from oscar_colony.breeding_scheme import Genotype
from oscar_colony.optimise.surplus_summary import (
    SexSplit,
    create_surplus_summary,
)
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
        pytest.param(
            "required_n_per_genotype_1_mutation_sex_split",
            "n_matings_1_mutation_sex_split",
            "offspring_per_scheme_1_mutation_sex_split",
            "surplus_1_mutation_sex_split",
            id="1 mutation - sex split",
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
        surplus_summary,
        request.getfixturevalue(expected_surplus),
        abs_tolerance=1e-4,
    )


def test_invalid_surplus_summary(
    n_matings_1_mutation_sex_split,
    offspring_per_scheme_1_mutation_sex_split,
):
    """Test that an error is raised for each genotype where the plan doesn't
    reach the required number - in total, and of males / females. Including
    genotypes that can't be produced by any of the breeding schemes (wt here).
    """

    sex_split = SexSplit(n_males=20, n_females=20)

    with pytest.raises(ValueError) as error:
        create_surplus_summary(
            required_n_per_genotype={
                (Genotype.WT,): 100,
                (Genotype.HET,): sex_split,
            },
            n_matings_per_scheme=n_matings_1_mutation_sex_split,
            offspring_per_scheme=offspring_per_scheme_1_mutation_sex_split,
        )

    expected_messages = [
        "wt: plan is 100.0 individuals short of the required number",
        "het: plan is 10.0 individuals short of the required number",
        "het: plan is 2.0 male individuals short of the required number",
        "het: plan is 8.0 female individuals short of the required number",
    ]
    for message in expected_messages:
        assert message in str(error.value)


@pytest.mark.parametrize(
    "surplus_summary, expected_csv",
    [
        pytest.param(
            "surplus_2_mutations",
            "converted-surplus-per-genotype.csv",
            id="no sex split",
        ),
        pytest.param(
            "surplus_2_mutations_sex_split",
            "converted-surplus-per-genotype-sex-split.csv",
            id="sex split",
        ),
    ],
)
def test_create_genotype_df(surplus_summary, expected_csv, request):
    """Test creation of a dataframe from SurplusSummary.surplus_per_genotype.

    The male / female columns should only be present when a sex split was
    required.
    """

    genotype_df = request.getfixturevalue(surplus_summary).create_genotype_df(
        decimal_places=2
    )
    expected_df = pd.read_csv(pooch_data_path(expected_csv))

    pd.testing.assert_frame_equal(genotype_df, expected_df)
