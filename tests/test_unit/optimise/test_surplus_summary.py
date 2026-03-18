import pytest

from oscar.breeding_scheme import Genotype
from oscar.optimise.surplus_summary import (
    GenotypeSurplus,
    SurplusSummary,
    create_surplus_summary,
)


@pytest.mark.parametrize(
    (
        "required_n_per_genotype, n_matings_per_scheme, offspring_per_scheme,"
        "expected_surplus"
    ),
    [
        pytest.param(
            "single_mutation_required_n_per_genotype",
            "single_mutation_n_matings",
            "single_mutation_offspring_per_scheme",
            SurplusSummary(
                total_n=3213,
                total_n_surplus=8,
                surplus_per_genotype={
                    (Genotype.HET,): GenotypeSurplus(
                        total_n=1561,
                        total_n_surplus=3,
                        percent_surplus=pytest.approx((3 / 1561) * 100),
                    ),
                    (Genotype.HOM,): GenotypeSurplus(
                        total_n=473,
                        total_n_surplus=2,
                        percent_surplus=pytest.approx((2 / 473) * 100),
                    ),
                    (Genotype.WT,): GenotypeSurplus(
                        total_n=1179,
                        total_n_surplus=3,
                        percent_surplus=pytest.approx((3 / 1179) * 100),
                    ),
                },
            ),
            id="1 mutation",
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

    assert surplus_summary == expected_surplus
