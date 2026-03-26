import pytest

from oscar.breeding_scheme import BreedingScheme, Genotype
from oscar.historical_stats import (
    BreedingSchemeStatistics,
    LineStatistics,
)
from oscar.optimise.estimate_offspring import (
    ExpectedOffspring,
    estimate_n_offspring_per_mating,
)


@pytest.fixture
def example_line_stats():
    return LineStatistics(
        n_mutations=1,
        total_n_offspring=140,
        total_n_offspring_per_genotype={
            (Genotype.WT,): 110,
            (Genotype.HET,): 24,
            (Genotype.HOM,): 6,
        },
        total_n_successful_matings=14,
        average_litter_size=10,
        stats_per_breeding_scheme={
            BreedingScheme("wt", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=3,
                n_successful_matings=9,
                average_litter_size=8,
                average_n_litters_per_pair=3,
                total_n_offspring=72,
                n_offspring_per_genotype={
                    (Genotype.WT,): 58,
                    (Genotype.HET,): 14,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): 58 / 72,
                    (Genotype.HET,): 14 / 72,
                },
            ),
            BreedingScheme("het", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=5,
                average_litter_size=13.6,
                average_n_litters_per_pair=2.5,
                total_n_offspring=68,
                n_offspring_per_genotype={
                    (Genotype.WT,): 52,
                    (Genotype.HET,): 10,
                    (Genotype.HOM,): 6,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): 52 / 68,
                    (Genotype.HET,): 10 / 68,
                    (Genotype.HOM,): 6 / 68,
                },
            ),
        },
    )


@pytest.mark.parametrize(
    "min_n_offspring, min_n_matings,expected_offspring_per_scheme",
    [
        pytest.param(
            1,
            1,
            {
                BreedingScheme("wt", "hom"): ExpectedOffspring(
                    total_n=10, n_per_genotype={(Genotype.HET,): 10.0}
                ),
                BreedingScheme("wt", "het"): ExpectedOffspring(
                    total_n=8,
                    n_per_genotype={
                        (Genotype.WT,): pytest.approx((58 / 72) * 8),
                        (Genotype.HET,): pytest.approx((14 / 72) * 8),
                    },
                ),
                BreedingScheme("hom", "hom"): ExpectedOffspring(
                    total_n=10, n_per_genotype={(Genotype.HOM,): 10.0}
                ),
                BreedingScheme("hom", "het"): ExpectedOffspring(
                    total_n=10,
                    n_per_genotype={(Genotype.HOM,): 5, (Genotype.HET,): 5},
                ),
                BreedingScheme("het", "het"): ExpectedOffspring(
                    total_n=13.6,
                    n_per_genotype={
                        (Genotype.WT,): pytest.approx((52 / 68) * 13.6),
                        (Genotype.HET,): pytest.approx((10 / 68) * 13.6),
                        (Genotype.HOM,): pytest.approx((6 / 68) * 13.6),
                    },
                ),
            },
            id="Min values == 1",
        ),
        pytest.param(
            70,
            0,
            {
                BreedingScheme("wt", "hom"): ExpectedOffspring(
                    total_n=10, n_per_genotype={(Genotype.HET,): 10.0}
                ),
                BreedingScheme("wt", "het"): ExpectedOffspring(
                    total_n=8,
                    n_per_genotype={
                        (Genotype.WT,): pytest.approx((58 / 72) * 8),
                        (Genotype.HET,): pytest.approx((14 / 72) * 8),
                    },
                ),
                BreedingScheme("hom", "hom"): ExpectedOffspring(
                    total_n=10, n_per_genotype={(Genotype.HOM,): 10.0}
                ),
                BreedingScheme("hom", "het"): ExpectedOffspring(
                    total_n=10,
                    n_per_genotype={(Genotype.HOM,): 5, (Genotype.HET,): 5},
                ),
                BreedingScheme("het", "het"): ExpectedOffspring(
                    total_n=13.6,
                    n_per_genotype={
                        (Genotype.WT,): pytest.approx(0.25 * 13.6),
                        (Genotype.HET,): pytest.approx(0.5 * 13.6),
                        (Genotype.HOM,): pytest.approx(0.25 * 13.6),
                    },
                ),
            },
            id="min_n_offspring=70. hetxhet should use mendelian ratios now",
        ),
        pytest.param(
            0,
            6,
            {
                BreedingScheme("wt", "hom"): ExpectedOffspring(
                    total_n=10, n_per_genotype={(Genotype.HET,): 10.0}
                ),
                BreedingScheme("wt", "het"): ExpectedOffspring(
                    total_n=8,
                    n_per_genotype={
                        (Genotype.WT,): pytest.approx((58 / 72) * 8),
                        (Genotype.HET,): pytest.approx((14 / 72) * 8),
                    },
                ),
                BreedingScheme("hom", "hom"): ExpectedOffspring(
                    total_n=10, n_per_genotype={(Genotype.HOM,): 10.0}
                ),
                BreedingScheme("hom", "het"): ExpectedOffspring(
                    total_n=10,
                    n_per_genotype={(Genotype.HOM,): 5, (Genotype.HET,): 5},
                ),
                BreedingScheme("het", "het"): ExpectedOffspring(
                    total_n=10.0,
                    n_per_genotype={
                        (Genotype.WT,): pytest.approx((52 / 68) * 10.0),
                        (Genotype.HET,): pytest.approx((10 / 68) * 10.0),
                        (Genotype.HOM,): pytest.approx((6 / 68) * 10.0),
                    },
                ),
            },
            id="min_n_matings=6. hetxhet should use line average_litter_size",
        ),
        pytest.param(
            0,
            20,
            {
                BreedingScheme("wt", "hom"): ExpectedOffspring(
                    total_n=4, n_per_genotype={(Genotype.HET,): 4.0}
                ),
                BreedingScheme("wt", "het"): ExpectedOffspring(
                    total_n=4,
                    n_per_genotype={
                        (Genotype.WT,): pytest.approx((58 / 72) * 4),
                        (Genotype.HET,): pytest.approx((14 / 72) * 4),
                    },
                ),
                BreedingScheme("hom", "hom"): ExpectedOffspring(
                    total_n=4, n_per_genotype={(Genotype.HOM,): 4.0}
                ),
                BreedingScheme("hom", "het"): ExpectedOffspring(
                    total_n=4,
                    n_per_genotype={(Genotype.HOM,): 2, (Genotype.HET,): 2},
                ),
                BreedingScheme("het", "het"): ExpectedOffspring(
                    total_n=4.0,
                    n_per_genotype={
                        (Genotype.WT,): pytest.approx((52 / 68) * 4.0),
                        (Genotype.HET,): pytest.approx((10 / 68) * 4.0),
                        (Genotype.HOM,): pytest.approx((6 / 68) * 4.0),
                    },
                ),
            },
            id="min_n_matings=20. All should use default_litter_size=4",
        ),
    ],
)
def test_estimate_n_offspring_per_mating(
    example_line_stats,
    min_n_offspring,
    min_n_matings,
    expected_offspring_per_scheme,
):
    """Test estimates of the number of offspring per mating, varying
    min_n_offspring and min_n_matings.
    """
    offspring_per_scheme = estimate_n_offspring_per_mating(
        example_line_stats,
        min_n_offspring=min_n_offspring,
        min_n_matings=min_n_matings,
        default_litter_size=4,
    )

    assert offspring_per_scheme == expected_offspring_per_scheme
