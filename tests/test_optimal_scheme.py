import pytest

from oscar.breeding_scheme import BreedingScheme, Genotype
from oscar.optimise.estimate_offspring import ExpectedOffspring
from oscar.optimise.optimal_scheme_calculator import _optimise_n_matings


@pytest.mark.parametrize(
    "required_n_per_genotype, offspring_per_scheme, expected_n_matings",
    [
        pytest.param(
            {
                (Genotype.WT,): 1176,
                (Genotype.HET,): 1558,
                (Genotype.HOM,): 471,
            },
            {
                BreedingScheme("het", "hom"): ExpectedOffspring(
                    total_n=5.58,
                    n_per_genotype={
                        (Genotype.HET,): 2.79,
                        (Genotype.HOM,): 2.79,
                    },
                ),
                BreedingScheme("het", "het"): ExpectedOffspring(
                    total_n=5.58,
                    n_per_genotype={
                        (Genotype.HET,): 2.79,
                        (Genotype.HOM,): 1.395,
                        (Genotype.WT,): 1.395,
                    },
                ),
                BreedingScheme("hom", "hom"): ExpectedOffspring(
                    total_n=5.58,
                    n_per_genotype={
                        (Genotype.HOM,): 5.58,
                    },
                ),
                BreedingScheme("hom", "wt"): ExpectedOffspring(
                    total_n=5.58,
                    n_per_genotype={
                        (Genotype.HET,): 5.58,
                    },
                ),
                BreedingScheme("wt", "het"): ExpectedOffspring(
                    total_n=5.58,
                    n_per_genotype={
                        (Genotype.HET,): 2.79,
                        (Genotype.WT,): 2.79,
                    },
                ),
            },
            {
                BreedingScheme("het", "het"): 274,
                BreedingScheme("hom", "hom"): 16,
                BreedingScheme("wt", "het"): 285,
            },
            id="1 mutation",
        ),
    ],
)
def test_optimise_n_matings(
    required_n_per_genotype, offspring_per_scheme, expected_n_matings
):
    n_matings_per_scheme = _optimise_n_matings(
        required_n_per_genotype,
        offspring_per_scheme,
    )

    assert n_matings_per_scheme == expected_n_matings
