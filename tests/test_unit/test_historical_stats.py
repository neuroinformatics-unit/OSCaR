import pandas as pd
import pytest

from oscar.breeding_scheme import BreedingScheme, Genotype
from oscar.historical_stats import (
    BreedingSchemeStatistics,
    LineStatistics,
    calculate_historical_stats_for_line,
)


@pytest.fixture
def expected_stats_single_mutation():
    return LineStatistics(
        total_n_offspring=18,
        total_n_offspring_per_genotype={
            (Genotype.WT,): 6,
            (Genotype.HET,): 10,
            (Genotype.HOM,): 2,
        },
        total_n_successful_matings=9,
        average_litter_size=pytest.approx(2, abs=1e-3),
        stats_per_breeding_scheme={
            BreedingScheme("wt", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=3,
                average_litter_size=pytest.approx(2.666, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.5, abs=1e-3),
                total_n_offspring=8,
                n_offspring_per_genotype={
                    (Genotype.WT,): 4,
                    (Genotype.HET,): 4,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): pytest.approx(0.5, abs=1e-3),
                    (Genotype.HET,): pytest.approx(0.5, abs=1e-3),
                },
            ),
            BreedingScheme("wt", "hom"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=pytest.approx(2, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.0, abs=1e-3),
                total_n_offspring=4,
                n_offspring_per_genotype={(Genotype.HET,): 4},
                proportion_offspring_per_genotype={
                    (Genotype.HET,): pytest.approx(1.0, abs=1e-3)
                },
            ),
            BreedingScheme("het", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=pytest.approx(2, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.0, abs=1e-3),
                total_n_offspring=4,
                n_offspring_per_genotype={
                    (Genotype.WT,): 2,
                    (Genotype.HET,): 1,
                    (Genotype.HOM,): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): pytest.approx(0.5, abs=1e-3),
                    (Genotype.HET,): pytest.approx(0.25, abs=1e-3),
                    (Genotype.HOM,): pytest.approx(0.25, abs=1e-3),
                },
            ),
            BreedingScheme("het", "hom"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=pytest.approx(1.0, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.0, abs=1e-3),
                total_n_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.HET,): 1,
                    (Genotype.HOM,): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET,): pytest.approx(0.5, abs=1e-3),
                    (Genotype.HOM,): pytest.approx(0.5, abs=1e-3),
                },
            ),
        },
    )


@pytest.fixture
def expected_stats_2_mutations():
    return LineStatistics(
        total_n_offspring=20,
        total_n_offspring_per_genotype={
            (Genotype.HOM, Genotype.HOM): 2,
            (Genotype.HET, Genotype.HOM): 6,
            (Genotype.HOM, Genotype.HET): 4,
            (Genotype.HET, Genotype.WT): 2,
            (Genotype.WT, Genotype.HET): 4,
            (Genotype.HET, Genotype.HET): 2,
        },
        total_n_successful_matings=10,
        average_litter_size=pytest.approx(2, abs=1e-3),
        stats_per_breeding_scheme={
            BreedingScheme("het_hom", "hom_het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=3,
                average_litter_size=pytest.approx(2.666, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.5, abs=1e-3),
                total_n_offspring=8,
                n_offspring_per_genotype={
                    (Genotype.HOM, Genotype.HOM): 2,
                    (Genotype.HET, Genotype.HOM): 6,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HOM, Genotype.HOM): pytest.approx(
                        0.25, abs=1e-3
                    ),
                    (Genotype.HET, Genotype.HOM): pytest.approx(
                        0.75, abs=1e-3
                    ),
                },
            ),
            BreedingScheme("hom_wt", "het_het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=pytest.approx(2.0, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.0, abs=1e-3),
                total_n_offspring=4,
                n_offspring_per_genotype={
                    (
                        Genotype.HOM,
                        Genotype.HET,
                    ): 4
                },
                proportion_offspring_per_genotype={
                    (
                        Genotype.HOM,
                        Genotype.HET,
                    ): pytest.approx(1.0, abs=1e-3)
                },
            ),
            BreedingScheme("het_wt", "wt_het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=pytest.approx(1.5, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.0, abs=1e-3),
                total_n_offspring=3,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT): 2,
                    (Genotype.WT, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT): pytest.approx(
                        0.666, abs=1e-3
                    ),
                    (Genotype.WT, Genotype.HET): pytest.approx(
                        0.333, abs=1e-3
                    ),
                },
            ),
            BreedingScheme("het_wt", "wt_hom"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=pytest.approx(1.0, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.0, abs=1e-3),
                total_n_offspring=2,
                n_offspring_per_genotype={
                    (
                        Genotype.WT,
                        Genotype.HET,
                    ): 2,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET): pytest.approx(1.0, abs=1e-3),
                },
            ),
            BreedingScheme("hom_wt", "wt_hom"): BreedingSchemeStatistics(
                n_breeding_pairs=1,
                n_successful_matings=1,
                average_litter_size=pytest.approx(3.0, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.0, abs=1e-3),
                total_n_offspring=3,
                n_offspring_per_genotype={
                    (
                        Genotype.HET,
                        Genotype.HET,
                    ): 2,
                    (Genotype.WT, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET): pytest.approx(
                        0.666, abs=1e-3
                    ),
                    (Genotype.WT, Genotype.HET): pytest.approx(
                        0.333, abs=1e-3
                    ),
                },
            ),
        },
    )


@pytest.fixture
def expected_stats_3_mutations():
    return LineStatistics(
        total_n_offspring=20,
        total_n_offspring_per_genotype={
            (Genotype.HET, Genotype.WT, Genotype.HOM): 8,
            (Genotype.WT, Genotype.HET, Genotype.HOM): 2,
            (Genotype.HOM, Genotype.WT, Genotype.HET): 2,
            (Genotype.WT, Genotype.WT, Genotype.HET): 1,
            (Genotype.HET, Genotype.WT, Genotype.WT): 1,
            (Genotype.WT, Genotype.HET, Genotype.WT): 3,
            (Genotype.HET, Genotype.HET, Genotype.WT): 2,
            (Genotype.WT, Genotype.HET, Genotype.HET): 1,
        },
        total_n_successful_matings=10,
        average_litter_size=pytest.approx(2, abs=1e-3),
        stats_per_breeding_scheme={
            BreedingScheme(
                "wt_wt_het", "het_het_het"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=3,
                average_litter_size=pytest.approx(2.666, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.5, abs=1e-3),
                total_n_offspring=8,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.HOM): 4,
                    (Genotype.WT, Genotype.HET, Genotype.HOM): 2,
                    (Genotype.HOM, Genotype.WT, Genotype.HET): 2,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.HOM): pytest.approx(
                        0.5, abs=1e-3
                    ),
                    (Genotype.WT, Genotype.HET, Genotype.HOM): pytest.approx(
                        0.25, abs=1e-3
                    ),
                    (Genotype.HOM, Genotype.WT, Genotype.HET): pytest.approx(
                        0.25, abs=1e-3
                    ),
                },
            ),
            BreedingScheme(
                "wt_het_het", "het_het_hom"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=pytest.approx(2.0, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.0, abs=1e-3),
                total_n_offspring=4,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.HOM): 4
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.HOM): pytest.approx(
                        1.0, abs=1e-3
                    )
                },
            ),
            BreedingScheme(
                "het_wt_wt", "wt_het_hom"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=pytest.approx(1.5, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.0, abs=1e-3),
                total_n_offspring=3,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.WT, Genotype.HET): 1,
                    (Genotype.HET, Genotype.WT, Genotype.WT): 1,
                    (Genotype.WT, Genotype.HET, Genotype.WT): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.WT, Genotype.HET): pytest.approx(
                        0.333, abs=1e-3
                    ),
                    (Genotype.HET, Genotype.WT, Genotype.WT): pytest.approx(
                        0.333, abs=1e-3
                    ),
                    (Genotype.WT, Genotype.HET, Genotype.WT): pytest.approx(
                        0.333, abs=1e-3
                    ),
                },
            ),
            BreedingScheme(
                "het_wt_wt", "wt_hom_het"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=pytest.approx(1.0, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.0, abs=1e-3),
                total_n_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET, Genotype.WT): 2,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET, Genotype.WT): pytest.approx(
                        1.0, abs=1e-3
                    ),
                },
            ),
            BreedingScheme(
                "wt_hom_het", "hom_wt_wt"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=1,
                n_successful_matings=1,
                average_litter_size=pytest.approx(3.0, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.0, abs=1e-3),
                total_n_offspring=3,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET, Genotype.WT): 2,
                    (Genotype.WT, Genotype.HET, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET, Genotype.WT): pytest.approx(
                        0.666, abs=1e-3
                    ),
                    (Genotype.WT, Genotype.HET, Genotype.HET): pytest.approx(
                        0.333, abs=1e-3
                    ),
                },
            ),
        },
    )


@pytest.mark.parametrize(
    "standardised_csv_path, line_name, expected_stats",
    [
        pytest.param(
            "standardised_single_mutation_csv_path",
            "Line-A",
            "expected_stats_single_mutation",
            id="1 mutation",
        ),
        pytest.param(
            "standardised_2_mutations_csv_path",
            "Line-AB",
            "expected_stats_2_mutations",
            id="2 mutations",
        ),
        pytest.param(
            "standardised_3_mutations_csv_path",
            "Line-ABC",
            "expected_stats_3_mutations",
            id="3 mutations",
        ),
    ],
)
def test_calculate_historical_stats_for_line(
    standardised_csv_path, line_name, expected_stats, request
):
    """
    Test calculation of summary historical statistics for lines with 1, 2 or
    3 mutations.
    """

    standardised_csv = pd.read_csv(
        request.getfixturevalue(standardised_csv_path)
    )
    expected_stats = request.getfixturevalue(expected_stats)

    line_stats = calculate_historical_stats_for_line(
        standardised_csv, line_name
    )
    assert line_stats == expected_stats
