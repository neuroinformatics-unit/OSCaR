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
        total_n_offspring=13,
        total_n_offspring_per_genotype={
            (Genotype.WT,): 4,
            (Genotype.HET,): 7,
            (Genotype.HOM,): 2,
        },
        stats_per_breeding_scheme={
            BreedingScheme("wt", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=3,
                n_successful_matings=4,
                average_litter_size=pytest.approx(1.0, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(1.333, abs=1e-3),
                total_n_offspring=4,
                n_offspring_per_genotype={
                    (Genotype.WT,): 2,
                    (Genotype.HET,): 2,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): pytest.approx(0.5, abs=1e-3),
                    (Genotype.HET,): pytest.approx(0.5, abs=1e-3),
                },
            ),
            BreedingScheme("wt", "hom"): BreedingSchemeStatistics(
                n_breeding_pairs=1,
                n_successful_matings=2,
                average_litter_size=pytest.approx(1.5, abs=1e-3),
                average_n_litters_per_pair=pytest.approx(2.0, abs=1e-3),
                total_n_offspring=3,
                n_offspring_per_genotype={(Genotype.HET,): 3},
                proportion_offspring_per_genotype={
                    (Genotype.HET,): pytest.approx(1.0, abs=1e-3)
                },
            ),
            BreedingScheme("het", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=3,
                n_successful_matings=3,
                average_litter_size=pytest.approx(1.333, abs=1e-3),
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


@pytest.mark.parametrize(
    "standardised_csv_path, line_name, expected_stats",
    [
        pytest.param(
            "standardised_single_mutation_csv_path",
            "Line-A",
            "expected_stats_single_mutation",
            id="1 mutation",
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
