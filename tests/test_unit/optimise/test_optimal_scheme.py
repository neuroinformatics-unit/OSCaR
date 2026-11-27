from contextlib import nullcontext as does_not_raise

import pytest

from oscar_colony.breeding_scheme import BreedingScheme, Genotype
from oscar_colony.historical_stats import LineStatistics
from oscar_colony.optimise.estimate_offspring import ExpectedOffspring
from oscar_colony.optimise.optimal_scheme_calculator import (
    _optimise_n_matings,
    calculate_optimal_scheme,
)
from oscar_colony.optimise.surplus_summary import SexSplit
from tests.helpers import assert_dataclass_equal


@pytest.mark.parametrize(
    (
        "required_n_per_genotype, line_stats, default_litter_size, "
        "expected_n_matings, expected_surplus"
    ),
    [
        pytest.param(
            "required_n_per_genotype_1_mutation",
            LineStatistics(line_name="Line-A", n_mutations=1),
            5.58,
            "n_matings_1_mutation",
            "surplus_1_mutation",
            id="1 mutation",
        ),
        pytest.param(
            "required_n_per_genotype_2_mutations",
            LineStatistics(line_name="Line-A", n_mutations=2),
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
    assert_dataclass_equal(
        surplus_summary,
        request.getfixturevalue(expected_surplus),
        abs_tolerance=1e-4,
    )


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
        pytest.param(
            "required_n_per_genotype_1_mutation_sex_split",
            "offspring_per_scheme_1_mutation_sex_split",
            "n_matings_1_mutation_sex_split",
            id="1 mutation - sex split",
        ),
    ],
)
def test_optimise_n_matings(
    required_n_per_genotype, offspring_per_scheme, expected_n_matings, request
):
    """Test individual optimisation function with 1 and 2 mutations."""

    n_matings_per_scheme = _optimise_n_matings(
        request.getfixturevalue(required_n_per_genotype),
        request.getfixturevalue(offspring_per_scheme),
    )

    assert n_matings_per_scheme == request.getfixturevalue(expected_n_matings)


@pytest.mark.parametrize(
    (
        "required_n_per_genotype, available_breeding_stock, "
        "expected_n_matings, expectation"
    ),
    [
        pytest.param(
            {(Genotype.HET,): 20, (Genotype.HOM,): 20},
            {
                (Genotype.HET,): SexSplit(0, 10),
                (Genotype.HOM,): SexSplit(10, 0),
            },
            {BreedingScheme("hom", "het"): 8},
            does_not_raise(),
            id="hom males and het females only",
        ),
        pytest.param(
            {(Genotype.HET,): 20, (Genotype.HOM,): 20},
            {(Genotype.HOM,): SexSplit(10, 0)},
            None,
            pytest.raises(
                ValueError,
                match="No breeding scheme can be set up from the available "
                "breeding stock",
            ),
            id="no female breeders",
        ),
        pytest.param(
            {(Genotype.HET,): 20, (Genotype.HOM,): 20},
            {(Genotype.HOM,): SexSplit(10, 10)},
            None,
            pytest.raises(
                ValueError,
                match="Can't breed het from the available breeding stock",
            ),
            id="only hom breeders",
        ),
        pytest.param(
            {(Genotype.HET,): 0, (Genotype.HOM,): 20},
            {(Genotype.HOM,): SexSplit(10, 10)},
            {BreedingScheme("hom", "hom"): 4},
            does_not_raise(),
            id="only hom breeders, 0 het required",
        ),
        pytest.param(
            {(Genotype.HET,): SexSplit(0, 0), (Genotype.HOM,): 20},
            {(Genotype.HOM,): SexSplit(10, 10)},
            {BreedingScheme("hom", "hom"): 4},
            does_not_raise(),
            id="only hom breeders, 0 het of each sex required",
        ),
        pytest.param(
            {(Genotype.HET,): 20, (Genotype.HOM,): 20},
            {
                (Genotype.HET,): SexSplit(10, 10),
                (Genotype.HOM,): SexSplit(6, 0),
            },
            {BreedingScheme("hom", "het"): 6, BreedingScheme("het", "het"): 3},
            does_not_raise(),
            id="six hom males",
        ),
        pytest.param(
            {(Genotype.HET,): 20, (Genotype.HOM,): 20},
            {
                (Genotype.HET,): SexSplit(4, 4),
                (Genotype.HOM,): SexSplit(4, 4),
            },
            {BreedingScheme("hom", "het"): 8},
            does_not_raise(),
            id="four of each, both ways",
        ),
        pytest.param(
            {(Genotype.HET,): 20, (Genotype.HOM,): 20},
            {
                (Genotype.HET,): SexSplit(0, 10),
                (Genotype.HOM,): SexSplit(1, 0),
            },
            None,
            pytest.raises(ValueError, match="Not enough breeders in stock"),
            id="one hom male",
        ),
    ],
)
def test_optimise_n_matings_breeding_stock(
    required_n_per_genotype,
    available_breeding_stock,
    expected_n_matings,
    expectation,
    offspring_per_scheme_1_mutation,
):
    """Test optimisation only uses the available breeding stock when values are
    provided, and raises an error when the stock can't be used."""

    with expectation:
        n_matings_per_scheme = _optimise_n_matings(
            required_n_per_genotype,
            offspring_per_scheme_1_mutation,
            available_breeding_stock,
        )
        assert n_matings_per_scheme == expected_n_matings


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
