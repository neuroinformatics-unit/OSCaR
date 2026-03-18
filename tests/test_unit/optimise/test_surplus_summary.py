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
            "required_n_per_genotype_1_mutation",
            "n_matings_1_mutation",
            "offspring_per_scheme_1_mutation",
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
        pytest.param(
            "required_n_per_genotype_2_mutations",
            "n_matings_2_mutations",
            "offspring_per_scheme_2_mutations",
            SurplusSummary(
                total_n=570,
                total_n_surplus=20,
                surplus_per_genotype={
                    (Genotype.WT, Genotype.HET): GenotypeSurplus(
                        total_n=52,
                        total_n_surplus=4,
                        percent_surplus=pytest.approx((4 / 52) * 100),
                    ),
                    (Genotype.WT, Genotype.HOM): GenotypeSurplus(
                        total_n=8,
                        total_n_surplus=1,
                        percent_surplus=pytest.approx((1 / 8) * 100),
                    ),
                    (Genotype.HET, Genotype.HET): GenotypeSurplus(
                        total_n=384,
                        total_n_surplus=4,
                        percent_surplus=pytest.approx((4 / 384) * 100),
                    ),
                    (Genotype.HET, Genotype.WT): GenotypeSurplus(
                        total_n=25,
                        total_n_surplus=2,
                        percent_surplus=pytest.approx((2 / 25) * 100),
                    ),
                    (Genotype.WT, Genotype.WT): GenotypeSurplus(
                        total_n=5,
                        total_n_surplus=1,
                        percent_surplus=pytest.approx((1 / 5) * 100),
                    ),
                    (Genotype.HET, Genotype.HOM): GenotypeSurplus(
                        total_n=20,
                        total_n_surplus=2,
                        percent_surplus=pytest.approx((2 / 20) * 100),
                    ),
                    (Genotype.HOM, Genotype.HET): GenotypeSurplus(
                        total_n=53,
                        total_n_surplus=4,
                        percent_surplus=pytest.approx((4 / 53) * 100),
                    ),
                    (Genotype.HOM, Genotype.HOM): GenotypeSurplus(
                        total_n=2,
                        total_n_surplus=1,
                        percent_surplus=pytest.approx((1 / 2) * 100),
                    ),
                    (Genotype.HOM, Genotype.WT): GenotypeSurplus(
                        total_n=21,
                        total_n_surplus=1,
                        percent_surplus=pytest.approx((1 / 21) * 100),
                    ),
                },
            ),
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

    assert surplus_summary == expected_surplus
