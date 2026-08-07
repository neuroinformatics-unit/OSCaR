import logging

import pandas as pd
import pytest

from oscar_colony.breeding_scheme import BreedingScheme, Genotype
from oscar_colony.historical_stats import (
    BreedingSchemeStatistics,
    LineStatistics,
    calculate_historical_stats_for_line,
)
from tests.helpers import assert_dataclass_equal
from tests.pooch_test_data import pooch_data_path


@pytest.fixture
def expected_stats_1_mutation():
    return LineStatistics(
        line_name="Line-A",
        mutations=["Mut-A"],
        n_mutations=1,
        total_n_offspring=18,
        total_n_genotyped_offspring=18,
        total_n_offspring_per_genotype={
            (Genotype.WT,): 6,
            (Genotype.HET,): 10,
            (Genotype.HOM,): 2,
        },
        total_n_successful_matings=9,
        average_litter_size=2,
        stats_per_breeding_scheme={
            BreedingScheme("wt", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=3,
                average_litter_size=2.666,
                average_n_litters_per_pair=1.5,
                total_n_offspring=8,
                total_n_genotyped_offspring=8,
                n_offspring_per_genotype={
                    (Genotype.WT,): 4,
                    (Genotype.HET,): 4,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): 0.5,
                    (Genotype.HET,): 0.5,
                },
            ),
            BreedingScheme("wt", "hom"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=2,
                average_n_litters_per_pair=1.0,
                total_n_offspring=4,
                total_n_genotyped_offspring=4,
                n_offspring_per_genotype={(Genotype.HET,): 4},
                proportion_offspring_per_genotype={(Genotype.HET,): 1.0},
            ),
            BreedingScheme("het", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=2,
                average_n_litters_per_pair=1.0,
                total_n_offspring=4,
                total_n_genotyped_offspring=4,
                n_offspring_per_genotype={
                    (Genotype.WT,): 2,
                    (Genotype.HET,): 1,
                    (Genotype.HOM,): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): 0.5,
                    (Genotype.HET,): 0.25,
                    (Genotype.HOM,): 0.25,
                },
            ),
            BreedingScheme("het", "hom"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=1.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.HET,): 1,
                    (Genotype.HOM,): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET,): 0.5,
                    (Genotype.HOM,): 0.5,
                },
            ),
        },
    )


@pytest.fixture
def expected_stats_2_mutations():
    return LineStatistics(
        line_name="Line-AB",
        mutations=["Mut-A", "Mut-B"],
        n_mutations=2,
        total_n_offspring=20,
        total_n_genotyped_offspring=20,
        total_n_offspring_per_genotype={
            (Genotype.HOM, Genotype.HOM): 2,
            (Genotype.HET, Genotype.HOM): 6,
            (Genotype.HOM, Genotype.HET): 4,
            (Genotype.HET, Genotype.WT): 2,
            (Genotype.WT, Genotype.HET): 2,
            (Genotype.HET, Genotype.HET): 4,
        },
        total_n_successful_matings=10,
        average_litter_size=2,
        stats_per_breeding_scheme={
            BreedingScheme("het_hom", "hom_het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=3,
                average_litter_size=2.666,
                average_n_litters_per_pair=1.5,
                total_n_offspring=8,
                total_n_genotyped_offspring=8,
                n_offspring_per_genotype={
                    (Genotype.HOM, Genotype.HOM): 2,
                    (Genotype.HET, Genotype.HOM): 6,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HOM, Genotype.HOM): 0.25,
                    (Genotype.HET, Genotype.HOM): 0.75,
                },
            ),
            BreedingScheme("hom_wt", "het_het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=2.0,
                average_n_litters_per_pair=1.0,
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
                    ): 1.0
                },
            ),
            BreedingScheme("het_wt", "wt_het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=1.5,
                average_n_litters_per_pair=1.0,
                total_n_offspring=3,
                total_n_genotyped_offspring=3,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT): 2,
                    (Genotype.WT, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT): 0.666,
                    (Genotype.WT, Genotype.HET): 0.333,
                },
            ),
            BreedingScheme("het_wt", "wt_hom"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=1.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={
                    (
                        Genotype.WT,
                        Genotype.HET,
                    ): 1,
                    (
                        Genotype.HET,
                        Genotype.HET,
                    ): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET): 0.5,
                    (Genotype.HET, Genotype.HET): 0.5,
                },
            ),
            BreedingScheme("hom_wt", "wt_hom"): BreedingSchemeStatistics(
                n_breeding_pairs=1,
                n_successful_matings=1,
                average_litter_size=3.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=3,
                total_n_genotyped_offspring=3,
                n_offspring_per_genotype={
                    (
                        Genotype.HET,
                        Genotype.HET,
                    ): 3,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET): 1,
                },
            ),
        },
    )


@pytest.fixture
def expected_stats_3_mutations():
    return LineStatistics(
        line_name="Line-ABC",
        mutations=["Mut-A", "Mut-B", "Mut-C"],
        n_mutations=3,
        total_n_offspring=20,
        total_n_genotyped_offspring=20,
        total_n_offspring_per_genotype={
            (Genotype.HET, Genotype.WT, Genotype.HOM): 10,
            (Genotype.WT, Genotype.HET, Genotype.HOM): 2,
            (Genotype.WT, Genotype.WT, Genotype.HET): 1,
            (Genotype.HET, Genotype.WT, Genotype.HET): 1,
            (Genotype.WT, Genotype.HET, Genotype.WT): 2,
            (Genotype.WT, Genotype.HET, Genotype.HET): 1,
            (Genotype.HET, Genotype.HET, Genotype.WT): 2,
            (Genotype.HET, Genotype.HET, Genotype.HET): 1,
        },
        total_n_successful_matings=10,
        average_litter_size=2,
        stats_per_breeding_scheme={
            BreedingScheme(
                "wt_wt_het", "het_het_het"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=3,
                average_litter_size=2.666,
                average_n_litters_per_pair=1.5,
                total_n_offspring=8,
                total_n_genotyped_offspring=8,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.HOM): 6,
                    (Genotype.WT, Genotype.HET, Genotype.HOM): 2,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.HOM): 0.75,
                    (Genotype.WT, Genotype.HET, Genotype.HOM): 0.25,
                },
            ),
            BreedingScheme(
                "wt_het_het", "het_het_hom"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=2.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=4,
                total_n_genotyped_offspring=4,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.HOM): 4
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.HOM): 1.0
                },
            ),
            BreedingScheme(
                "het_wt_wt", "wt_het_hom"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=1.5,
                average_n_litters_per_pair=1.0,
                total_n_offspring=3,
                total_n_genotyped_offspring=3,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.WT, Genotype.HET): 1,
                    (Genotype.HET, Genotype.WT, Genotype.HET): 1,
                    (Genotype.WT, Genotype.HET, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.WT, Genotype.HET): 0.333,
                    (Genotype.HET, Genotype.WT, Genotype.HET): 0.333,
                    (Genotype.WT, Genotype.HET, Genotype.HET): 0.333,
                },
            ),
            BreedingScheme(
                "het_wt_wt", "wt_hom_het"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=1.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET, Genotype.WT): 2,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET, Genotype.WT): 1.0,
                },
            ),
            BreedingScheme(
                "wt_hom_het", "hom_wt_wt"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=1,
                n_successful_matings=1,
                average_litter_size=3.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=3,
                total_n_genotyped_offspring=3,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET, Genotype.WT): 2,
                    (Genotype.HET, Genotype.HET, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET, Genotype.WT): 0.666,
                    (Genotype.HET, Genotype.HET, Genotype.HET): 0.333,
                },
            ),
        },
    )


@pytest.mark.parametrize(
    "standardised_csv_name, line_name, expected_stats",
    [
        pytest.param(
            "standardised-data-1-mutation.csv",
            "Line-A",
            "expected_stats_1_mutation",
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
    assert_dataclass_equal(line_stats, expected_stats, abs_tolerance=1e-3)


def test_calculate_historical_stats_for_line_logs(caplog):
    standardised_csv = pd.read_csv(
        pooch_data_path("standardised-data-1-mutation.csv")
    )

    with caplog.at_level(logging.DEBUG):
        calculate_historical_stats_for_line(standardised_csv, "Line-A")

    expected_messages = [
        "Calculating historical stats for line 'Line-A' with "
        "18 offspring records",
        "Calculating stats for breeding scheme",
    ]
    for message in expected_messages:
        assert message in caplog.text


@pytest.fixture
def expected_stats_ungenotyped():
    return LineStatistics(
        line_name="Line-AB",
        mutations=["Mut-A", "Mut-B"],
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
                average_litter_size=5.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=5,
                total_n_genotyped_offspring=3,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.HOM): 3,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.HOM): 1,
                },
            ),
            BreedingScheme("het_het", "het_het"): BreedingSchemeStatistics(
                n_breeding_pairs=1,
                n_successful_matings=1,
                average_litter_size=3.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=3,
                total_n_genotyped_offspring=3,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT): 1,
                    (Genotype.HET, Genotype.HET): 2,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT): 0.333,
                    (Genotype.HET, Genotype.HET): 0.666,
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
    assert_dataclass_equal(
        line_stats, expected_stats_ungenotyped, abs_tolerance=1e-3
    )


@pytest.fixture
def expected_stats_multi_parent_1_mutation():
    return LineStatistics(
        line_name="Line-A",
        mutations=["Mut-A"],
        n_mutations=1,
        total_n_offspring=16,
        total_n_genotyped_offspring=16,
        total_n_offspring_per_genotype={
            (Genotype.HET,): 8,
            (Genotype.HOM,): 3,
            (Genotype.WT,): 5,
        },
        total_n_successful_matings=23,
        average_litter_size=0.696,
        stats_per_breeding_scheme={
            BreedingScheme("hom", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=6,
                n_successful_matings=4,
                average_litter_size=0.75,
                average_n_litters_per_pair=0.667,
                total_n_offspring=3,
                total_n_genotyped_offspring=3,
                n_offspring_per_genotype={
                    (Genotype.HET,): 2,
                    (Genotype.HOM,): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET,): 0.667,
                    (Genotype.HOM,): 0.333,
                },
            ),
            BreedingScheme("wt", "wt"): BreedingSchemeStatistics(
                n_breeding_pairs=4,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=0.5,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.WT,): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): 1.0,
                },
            ),
            BreedingScheme("het", "wt"): BreedingSchemeStatistics(
                n_breeding_pairs=12,
                n_successful_matings=10,
                average_litter_size=0.7,
                average_n_litters_per_pair=0.833,
                total_n_offspring=7,
                total_n_genotyped_offspring=7,
                n_offspring_per_genotype={
                    (Genotype.WT,): 3,
                    (Genotype.HET,): 4,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): 0.429,
                    (Genotype.HET,): 0.571,
                },
            ),
            BreedingScheme("hom", "hom"): BreedingSchemeStatistics(
                n_breeding_pairs=4,
                n_successful_matings=4,
                average_litter_size=0.5,
                average_n_litters_per_pair=1.0,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.HOM,): 2,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HOM,): 1.0,
                },
            ),
            BreedingScheme("het", "het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=1,
                average_litter_size=2.0,
                average_n_litters_per_pair=0.5,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.WT,): 1,
                    (Genotype.HET,): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT,): 0.5,
                    (Genotype.HET,): 0.5,
                },
            ),
            BreedingScheme("hom", "wt"): BreedingSchemeStatistics(
                n_breeding_pairs=4,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=0.5,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.HET,): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET,): 1.0,
                },
            ),
        },
    )


@pytest.fixture
def expected_stats_multi_parent_2_mutations():
    return LineStatistics(
        line_name="Line-AB",
        mutations=["Mut-A", "Mut-B"],
        n_mutations=2,
        total_n_offspring=15,
        total_n_genotyped_offspring=15,
        total_n_offspring_per_genotype={
            (Genotype.WT, Genotype.HET): 3,
            (Genotype.HET, Genotype.WT): 5,
            (Genotype.HET, Genotype.HET): 5,
            (Genotype.HOM, Genotype.HET): 2,
        },
        total_n_successful_matings=22,
        average_litter_size=0.682,
        stats_per_breeding_scheme={
            BreedingScheme("wt_het", "het_wt"): BreedingSchemeStatistics(
                n_breeding_pairs=10,
                n_successful_matings=8,
                average_litter_size=0.625,
                average_n_litters_per_pair=0.8,
                total_n_offspring=5,
                total_n_genotyped_offspring=5,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET): 1,
                    (Genotype.HET, Genotype.WT): 2,
                    (Genotype.HET, Genotype.HET): 2,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET): 0.2,
                    (Genotype.HET, Genotype.WT): 0.4,
                    (Genotype.HET, Genotype.HET): 0.4,
                },
            ),
            BreedingScheme("hom_wt", "hom_het"): BreedingSchemeStatistics(
                n_breeding_pairs=6,
                n_successful_matings=4,
                average_litter_size=0.5,
                average_n_litters_per_pair=0.667,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.HOM, Genotype.HET): 2,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HOM, Genotype.HET): 1.0,
                },
            ),
            BreedingScheme("hom_wt", "het_het"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=1,
                average_litter_size=1.0,
                average_n_litters_per_pair=0.5,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET): 1.0,
                },
            ),
            BreedingScheme("het_wt", "wt_wt"): BreedingSchemeStatistics(
                n_breeding_pairs=4,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=0.5,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT): 1.0,
                },
            ),
            BreedingScheme("het_het", "het_wt"): BreedingSchemeStatistics(
                n_breeding_pairs=4,
                n_successful_matings=3,
                average_litter_size=1.0,
                average_n_litters_per_pair=0.75,
                total_n_offspring=3,
                total_n_genotyped_offspring=3,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET): 1,
                    (Genotype.HET, Genotype.WT): 1,
                    (Genotype.HET, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET): 0.333,
                    (Genotype.HET, Genotype.WT): 0.333,
                    (Genotype.HET, Genotype.HET): 0.333,
                },
            ),
            BreedingScheme("hom_het", "wt_wt"): BreedingSchemeStatistics(
                n_breeding_pairs=4,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=0.5,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET): 1.0,
                },
            ),
            BreedingScheme("het_het", "wt_wt"): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=1.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET): 1,
                    (Genotype.HET, Genotype.WT): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET): 0.5,
                    (Genotype.HET, Genotype.WT): 0.5,
                },
            ),
        },
    )


@pytest.fixture
def expected_stats_multi_parent_3_mutations():
    return LineStatistics(
        line_name="Line-ABC",
        mutations=["Mut-A", "Mut-B", "Mut-C"],
        n_mutations=3,
        total_n_offspring=16,
        total_n_genotyped_offspring=16,
        total_n_offspring_per_genotype={
            (Genotype.WT, Genotype.WT, Genotype.HOM): 2,
            (Genotype.HOM, Genotype.WT, Genotype.WT): 1,
            (Genotype.HET, Genotype.HET, Genotype.HET): 3,
            (Genotype.HET, Genotype.HOM, Genotype.HET): 1,
            (Genotype.HET, Genotype.WT, Genotype.WT): 1,
            (Genotype.HET, Genotype.WT, Genotype.HOM): 1,
            (Genotype.WT, Genotype.HET, Genotype.HET): 1,
            (Genotype.HET, Genotype.WT, Genotype.HET): 2,
            (Genotype.WT, Genotype.HET, Genotype.WT): 2,
            (Genotype.HET, Genotype.HET, Genotype.HOM): 1,
            (Genotype.HOM, Genotype.HET, Genotype.WT): 1,
        },
        total_n_successful_matings=25,
        average_litter_size=0.64,
        stats_per_breeding_scheme={
            BreedingScheme(
                "wt_wt_het", "het_het_het"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=4,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=0.5,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.WT, Genotype.HOM): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.WT, Genotype.HOM): 1.0,
                },
            ),
            BreedingScheme(
                "het_het_het", "hom_wt_het"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=4,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=0.5,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.HOM, Genotype.WT, Genotype.WT): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HOM, Genotype.WT, Genotype.WT): 1.0,
                },
            ),
            BreedingScheme(
                "het_het_het", "wt_het_het"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=1.0,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET, Genotype.HET): 1.0,
                },
            ),
            BreedingScheme(
                "het_hom_wt", "wt_het_hom"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=1.0,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.HOM, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.HOM, Genotype.HET): 1.0,
                },
            ),
            BreedingScheme("het_wt_wt", "wt_wt_wt"): BreedingSchemeStatistics(
                n_breeding_pairs=4,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=0.5,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.WT): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.WT): 1.0,
                },
            ),
            BreedingScheme(
                "het_wt_hom", "wt_wt_hom"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=1.0,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.HOM): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.WT, Genotype.HOM): 1.0,
                },
            ),
            BreedingScheme(
                "het_wt_wt", "wt_het_hom"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=4,
                average_litter_size=0.5,
                average_n_litters_per_pair=2.0,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET, Genotype.HET): 1,
                    (Genotype.HET, Genotype.WT, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET, Genotype.HET): 0.5,
                    (Genotype.HET, Genotype.WT, Genotype.HET): 0.5,
                },
            ),
            BreedingScheme(
                "wt_het_hom", "het_wt_het"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=4,
                n_successful_matings=2,
                average_litter_size=1.0,
                average_n_litters_per_pair=0.5,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.WT, Genotype.HOM): 1,
                    (Genotype.HET, Genotype.HET, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.WT, Genotype.HOM): 0.5,
                    (Genotype.HET, Genotype.HET, Genotype.HET): 0.5,
                },
            ),
            BreedingScheme(
                "het_wt_het", "het_het_wt"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=1,
                average_litter_size=2.0,
                average_n_litters_per_pair=0.5,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET, Genotype.WT): 1,
                    (Genotype.HET, Genotype.HET, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET, Genotype.WT): 0.5,
                    (Genotype.HET, Genotype.HET, Genotype.HET): 0.5,
                },
            ),
            BreedingScheme(
                "hom_het_het", "wt_wt_het"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=4,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=0.5,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET, Genotype.HOM): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HET, Genotype.HET, Genotype.HOM): 1.0,
                },
            ),
            BreedingScheme(
                "hom_wt_het", "hom_het_wt"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=0.5,
                average_n_litters_per_pair=1.0,
                total_n_offspring=1,
                total_n_genotyped_offspring=1,
                n_offspring_per_genotype={
                    (Genotype.HOM, Genotype.HET, Genotype.WT): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.HOM, Genotype.HET, Genotype.WT): 1.0,
                },
            ),
            BreedingScheme(
                "het_het_het", "wt_wt_wt"
            ): BreedingSchemeStatistics(
                n_breeding_pairs=2,
                n_successful_matings=2,
                average_litter_size=1.0,
                average_n_litters_per_pair=1.0,
                total_n_offspring=2,
                total_n_genotyped_offspring=2,
                n_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET, Genotype.WT): 1,
                    (Genotype.HET, Genotype.WT, Genotype.HET): 1,
                },
                proportion_offspring_per_genotype={
                    (Genotype.WT, Genotype.HET, Genotype.WT): 0.5,
                    (Genotype.HET, Genotype.WT, Genotype.HET): 0.5,
                },
            ),
        },
    )


@pytest.mark.parametrize(
    "standardised_csv_name, line_name, expected_stats",
    [
        pytest.param(
            "standardised-data-multiple-parents-1-mutation.csv",
            "Line-A",
            "expected_stats_multi_parent_1_mutation",
            id="1 mutation",
        ),
        pytest.param(
            "standardised-data-multiple-parents-2-mutations.csv",
            "Line-AB",
            "expected_stats_multi_parent_2_mutations",
            id="2 mutations",
        ),
        pytest.param(
            "standardised-data-multiple-parents-3-mutations.csv",
            "Line-ABC",
            "expected_stats_multi_parent_3_mutations",
            id="3 mutations",
        ),
    ],
)
def test_calculate_historical_stats_line_multiple_parents(
    standardised_csv_name, line_name, expected_stats, request
):
    """
    Test calculation of summary historical statistics for lines with 1, 2 or
    3 mutations account parents.
    """

    standardised_csv = pd.read_csv(pooch_data_path(standardised_csv_name))
    expected_stats = request.getfixturevalue(expected_stats)

    line_stats = calculate_historical_stats_for_line(
        standardised_csv, line_name
    )
    assert_dataclass_equal(line_stats, expected_stats, abs_tolerance=1e-3)


def test_total_number_genotype_df(expected_stats_2_mutations):
    """
    Test creation of a summary dataframe from
    LineStatistics.total_n_offspring_per_genotype
    """

    genotype_df = expected_stats_2_mutations.create_n_per_genotype_df()
    expected_genotype_df = pd.read_csv(
        pooch_data_path("converted-total-number.csv")
    )

    pd.testing.assert_frame_equal(genotype_df, expected_genotype_df)


def test_scheme_summary_df(expected_stats_2_mutations):
    """
    Test creation of a summary dataframe from
    LineStatistics.stats_per_breeding_scheme
    """

    summary_df = expected_stats_2_mutations.create_scheme_summary_df(
        decimal_places=2
    )
    expected_summary_df = pd.read_csv(
        pooch_data_path("converted-scheme-summary.csv")
    )

    # convert summary df Scheme column to string for comparison
    summary_df = summary_df.astype({"Scheme": "str"})

    pd.testing.assert_frame_equal(summary_df, expected_summary_df)


def test_scheme_number_df(expected_stats_2_mutations):
    """
    Test creation of a dataframe from
    LineStatistics.stats_per_breeding_scheme.n_offspring_per_genotype
    """

    number_df = expected_stats_2_mutations.create_scheme_number_df()
    expected_number_df = pd.read_csv(
        pooch_data_path("converted-scheme-number.csv")
    )

    # convert summary df Scheme column to string for comparison
    number_df = number_df.astype({"Scheme": "str"})

    pd.testing.assert_frame_equal(number_df, expected_number_df)


def test_scheme_proportion_df(expected_stats_2_mutations):
    """
    Test creation of a dataframe from
    LineStatistics.stats_per_breeding_scheme.proportion_offspring_per_genotype
    """

    proportion_df = expected_stats_2_mutations.create_scheme_proportion_df(
        decimal_places=2
    )
    expected_proportion_df = pd.read_csv(
        pooch_data_path("converted-scheme-proportion.csv")
    )

    # convert summary df Scheme column to string for comparison
    proportion_df = proportion_df.astype({"Scheme": "str"})

    pd.testing.assert_frame_equal(proportion_df, expected_proportion_df)
