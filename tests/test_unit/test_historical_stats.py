import pandas as pd
import pytest

from oscar.breeding_scheme import BreedingScheme, Genotype
from oscar.historical_stats import (
    BreedingSchemeStatistics,
    LineStatistics,
    calculate_historical_stats_for_line,
)
from tests.pooch_test_data import pooch_data_path


@pytest.fixture
def expected_stats_single_mutation():
    return LineStatistics(
        n_mutations=1,
        total_n_offspring=18,
        total_n_genotyped_offspring=18,
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
                total_n_genotyped_offspring=8,
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
                total_n_genotyped_offspring=4,
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
                total_n_genotyped_offspring=4,
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
                total_n_genotyped_offspring=2,
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
        n_mutations=2,
        total_n_offspring=20,
        total_n_genotyped_offspring=20,
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
                total_n_genotyped_offspring=8,
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
                total_n_genotyped_offspring=4,
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
                total_n_genotyped_offspring=3,
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
                total_n_genotyped_offspring=2,
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
                total_n_genotyped_offspring=3,
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
        n_mutations=3,
        total_n_offspring=20,
        total_n_genotyped_offspring=20,
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
                total_n_genotyped_offspring=8,
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
                total_n_genotyped_offspring=4,
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
                total_n_genotyped_offspring=3,
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
                total_n_genotyped_offspring=2,
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
                total_n_genotyped_offspring=3,
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
    "standardised_csv_name, line_name, expected_stats",
    [
        pytest.param(
            "standardised-data-single-mutation.csv",
            "Line-A",
            "expected_stats_single_mutation",
            id="1 mutation",
        ),
        pytest.param(
            "standardised-data-2-mutations.csv",
            "Line-AB",
            "expected_stats_2_mutations",
            id="2 mutations",
        ),
        pytest.param(
            "standardised-data-3-mutations.csv",
            "Line-ABC",
            "expected_stats_3_mutations",
            id="3 mutations",
        ),
    ],
)
def test_calculate_historical_stats_for_line(
    standardised_csv_name, line_name, expected_stats, request
):
    """
    Test calculation of summary historical statistics for lines with 1, 2 or
    3 mutations.
    """

    standardised_csv = pd.read_csv(pooch_data_path(standardised_csv_name))
    expected_stats = request.getfixturevalue(expected_stats)

    line_stats = calculate_historical_stats_for_line(
        standardised_csv, line_name
    )
    assert line_stats == expected_stats


@pytest.fixture
def expected_stats_ungenotyped():
    return LineStatistics(
        n_mutations=2,
        total_n_offspring=9,
        total_n_genotyped_offspring=6,
        total_n_offspring_per_genotype={
            (Genotype.HET, Genotype.WT): 1,
            (Genotype.HET, Genotype.HOM): 3,
            (Genotype.HET, Genotype.HET): 2,
        },
        total_n_successful_matings=3,
        average_litter_size=3.0,
        stats_per_breeding_scheme={
            BreedingScheme("wt_hom", "hom_het"): BreedingSchemeStatistics(
                n_breeding_pairs=1,
                n_successful_matings=1,
                average_litter_size=6.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=6,
                total_n_genotyped_offspring=4,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT): 1,
                    (Genotype.HET, Genotype.HOM): 3,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT): 0.25,
                    (Genotype.HET, Genotype.HOM): 0.75,
                },
            ),
            BreedingScheme("het_het", "het_het"): BreedingSchemeStatistics(
                n_breeding_pairs=1,
                n_successful_matings=1,
                average_litter_size=2.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={(Genotype.HET, Genotype.HET): 2},
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET): 1.0
                },
            ),
            BreedingScheme("wt_wt", "wt_wt"): BreedingSchemeStatistics(
                n_breeding_pairs=1,
                n_successful_matings=1,
                average_litter_size=1.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=1,
                total_n_genotyped_offspring=0,
                n_offspring_per_genotype={},
                proportion_offspring_per_genotype={},
            ),
        },
    )


def test_handling_ungenotyped_individuals_in_stats(expected_stats_ungenotyped):
    """
    Test that un-genotyped individuals (empty genotype_offspring column)
    are included in totals for litter size calculations, but _excluded_ from
    totals for genotype proportions.
    """

    # This csv contains 3 un-genotyped individuals
    standard_csv = pd.read_csv(
        pooch_data_path("standardised-data-forbidden-genotypes.csv")
    )

    line_stats = calculate_historical_stats_for_line(standard_csv, "Line-AB")
    assert line_stats == expected_stats_ungenotyped
