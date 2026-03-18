import pytest

from oscar.breeding_scheme import BreedingScheme, Genotype
from oscar.historical_stats import LineStatistics
from oscar.optimise.estimate_offspring import ExpectedOffspring
from oscar.optimise.optimal_scheme_calculator import (
    _optimise_n_matings,
    calculate_optimal_scheme,
)


@pytest.mark.parametrize(
    (
        "required_n_per_genotype, line_stats, default_litter_size, "
        "expected_n_matings, expected_surplus"
    ),
    [
        pytest.param(
            "required_n_per_genotype_1_mutation",
            LineStatistics(n_mutations=1),
            5.58,
            "n_matings_1_mutation",
            "surplus_1_mutation",
            id="1 mutation",
        ),
        pytest.param(
            "required_n_per_genotype_2_mutations",
            LineStatistics(n_mutations=2),
            5.79,
            "n_matings_2_mutations",
            "surplus_2_mutations",
            id="2 mutations",
        ),
    ],
)
def test_calculate_optimal_scheme(
    required_n_per_genotype,
    line_stats,
    default_litter_size,
    expected_n_matings,
    expected_surplus,
    request,
):
    """Simple test of overall calculate_optimal_scheme function with
    default line_stats.

    This will use the theoretical mendelian ratio for all, with the provided
    default_litter_size.
    """

    n_matings_per_scheme, surplus_summary = calculate_optimal_scheme(
        required_n_per_genotype=request.getfixturevalue(
            required_n_per_genotype
        ),
        line_stats=line_stats,
        default_litter_size=default_litter_size,
    )

    assert n_matings_per_scheme == request.getfixturevalue(expected_n_matings)
    assert surplus_summary == request.getfixturevalue(expected_surplus)


@pytest.mark.parametrize(
    "required_n_per_genotype, offspring_per_scheme, expected_n_matings",
    [
        pytest.param(
            "required_n_per_genotype_1_mutation",
            "offspring_per_scheme_1_mutation",
            "n_matings_1_mutation",
            id="1 mutation",
        ),
        pytest.param(
            "required_n_per_genotype_2_mutations",
            "offspring_per_scheme_2_mutations",
            "n_matings_2_mutations",
            id="2 mutations",
        ),
    ],
)
def test_optimise_n_matings(
    required_n_per_genotype, offspring_per_scheme, expected_n_matings, request
):
    n_matings_per_scheme = _optimise_n_matings(
        request.getfixturevalue(required_n_per_genotype),
        request.getfixturevalue(offspring_per_scheme),
    )

    assert n_matings_per_scheme == request.getfixturevalue(expected_n_matings)


def test_optimise_impossible_scheme():
    """Test situation where the given breeding schemes cannot produce the
    required offspring genotypes."""

    required_n_per_genotype = {(Genotype.WT,): 100, (Genotype.HET,): 50}

    offspring_per_scheme = {
        BreedingScheme("hom", "hom"): ExpectedOffspring(
            total_n=8, n_per_genotype={(Genotype.HOM,): 8}
        )
    }

    error_msg = (
        "Number of matings couldn't be optimised: The problem is infeasible"
    )
    with pytest.raises(ValueError, match=error_msg):
        _optimise_n_matings(
            required_n_per_genotype,
            offspring_per_scheme,
        )
