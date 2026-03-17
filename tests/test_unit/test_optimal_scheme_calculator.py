import pytest

from oscar.breeding_scheme import BreedingScheme, Genotype
from oscar.historical_stats import (
    BreedingSchemeStatistics,
    LineStatistics,
)
from oscar.optimal_scheme_calculator import (
    ExpectedOffspring,
    _estimate_n_offspring_per_mating,
)


@pytest.fixture
def example_line_stats():
    return LineStatistics(
        total_n_offspring=100,
        total_n_offspring_per_genotype={
            (Genotype.WT,): 70,
            (Genotype.HET,): 24,
            (Genotype.HOM,): 6,
        },
        stats_per_breeding_scheme={
            BreedingScheme("wt", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=3,
                n_successful_matings=9,
                average_litter_size=6,
                average_n_litters_per_pair=3,
                total_n_offspring=54,
                n_offspring_per_genotype={
                    (Genotype.WT,): 40,
                    (Genotype.HET,): 14,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): 40 / 54,
                    (Genotype.HET,): 14 / 54,
                },
            ),
            BreedingScheme("het", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=5,
                average_litter_size=46 / 5,
                average_n_litters_per_pair=2.5,
                total_n_offspring=46,
                n_offspring_per_genotype={
                    (Genotype.WT,): 30,
                    (Genotype.HET,): 10,
                    (Genotype.HOM,): 6,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): 30 / 46,
                    (Genotype.HET,): 10 / 46,
                    (Genotype.HOM,): 6 / 46,
                },
            ),
        },
    )


def test_estimate_n_offspring_per_mating(example_line_stats):
    offspring_per_scheme = _estimate_n_offspring_per_mating(
        example_line_stats,
        min_n_offspring=0,
        min_n_matings=0,
        default_litter_size=4,
    )

    expected_offspring_per_scheme = {
        BreedingScheme("wt", "hom"): ExpectedOffspring(
            total_n=0, n_per_genotype={}
        ),
        BreedingScheme("wt", "het"): ExpectedOffspring(
            total_n=6,
            n_per_genotype={
                (Genotype.WT,): pytest.approx((40 / 54) * 6),
                (Genotype.HET,): pytest.approx((14 / 54) * 6),
            },
        ),
        BreedingScheme("hom", "hom"): ExpectedOffspring(
            total_n=0, n_per_genotype={}
        ),
        BreedingScheme("hom", "het"): ExpectedOffspring(
            total_n=0, n_per_genotype={}
        ),
        BreedingScheme("het", "het"): ExpectedOffspring(
            total_n=9.2,
            n_per_genotype={
                (Genotype.WT,): pytest.approx((30 / 46) * 9.2),
                (Genotype.HET,): pytest.approx((10 / 46) * 9.2),
                (Genotype.HOM,): pytest.approx((6 / 46) * 9.2),
            },
        ),
    }

    assert offspring_per_scheme == expected_offspring_per_scheme
