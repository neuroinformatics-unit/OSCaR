import pytest

from oscar.optimise.surplus_summary import create_surplus_summary


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

    assert surplus_summary == request.getfixturevalue(expected_surplus)
