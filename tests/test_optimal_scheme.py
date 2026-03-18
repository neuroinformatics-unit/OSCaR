import pytest

from oscar.breeding_scheme import BreedingScheme, Genotype
from oscar.optimise.estimate_offspring import ExpectedOffspring
from oscar.optimise.optimal_scheme_calculator import (
    _optimise_n_matings,
)

# def test_calculate_optimal_scheme():
#     calculate_optimal_scheme(
#         required_n_per_genotype,
#         line_stats,
#         default_litter_size,
#         min_n_matings,
#         min_n_offspring
#     )


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
        pytest.param(
            {
                (Genotype.HET, Genotype.HET): 380,
                (Genotype.HET, Genotype.HOM): 18,
                (Genotype.HET, Genotype.WT): 23,
                (Genotype.HOM, Genotype.HET): 49,
                (Genotype.HOM, Genotype.HOM): 1,
                (Genotype.HOM, Genotype.WT): 20,
                (Genotype.WT, Genotype.HET): 48,
                (Genotype.WT, Genotype.HOM): 7,
                (Genotype.WT, Genotype.WT): 4,
            },
            {
                BreedingScheme("wt_wt", "het_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.WT): 1.4475,
                        (Genotype.WT, Genotype.HET): 1.4475,
                        (Genotype.WT, Genotype.WT): 1.4475,
                    },
                ),
                BreedingScheme("wt_wt", "het_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 2.895,
                        (Genotype.WT, Genotype.HET): 2.895,
                    },
                ),
                BreedingScheme("wt_wt", "hom_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 2.895,
                        (Genotype.HET, Genotype.WT): 2.895,
                    },
                ),
                BreedingScheme("wt_wt", "hom_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 5.79,
                    },
                ),
                BreedingScheme("wt_het", "wt_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.WT, Genotype.HET): 2.895,
                        (Genotype.WT, Genotype.HOM): 2.895,
                    },
                ),
                BreedingScheme("wt_het", "het_wt"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.WT): 1.4475,
                        (Genotype.WT, Genotype.HET): 1.4475,
                        (Genotype.WT, Genotype.WT): 1.4475,
                    },
                ),
                BreedingScheme("wt_het", "het_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.HOM): 0.72375,
                        (Genotype.HET, Genotype.WT): 0.72375,
                        (Genotype.WT, Genotype.HET): 1.4475,
                        (Genotype.WT, Genotype.HOM): 0.72375,
                        (Genotype.WT, Genotype.WT): 0.72375,
                    },
                ),
                BreedingScheme("wt_het", "het_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.HOM): 1.4475,
                        (Genotype.WT, Genotype.HET): 1.4475,
                        (Genotype.WT, Genotype.HOM): 1.4475,
                    },
                ),
                BreedingScheme("wt_het", "hom_wt"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 2.895,
                        (Genotype.HET, Genotype.WT): 2.985,
                    },
                ),
                BreedingScheme("wt_het", "hom_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 2.895,
                        (Genotype.HET, Genotype.HOM): 1.4475,
                        (Genotype.HET, Genotype.WT): 1.4475,
                    },
                ),
                BreedingScheme("wt_het", "hom_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 2.895,
                        (Genotype.HET, Genotype.HOM): 2.895,
                    },
                ),
                BreedingScheme("wt_hom", "het_wt"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 2.895,
                        (Genotype.WT, Genotype.HET): 2.895,
                    },
                ),
                BreedingScheme("wt_hom", "het_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.HOM): 1.4475,
                        (Genotype.WT, Genotype.HET): 1.4475,
                        (Genotype.WT, Genotype.HOM): 1.4475,
                    },
                ),
                BreedingScheme("wt_hom", "het_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HOM): 2.895,
                        (Genotype.WT, Genotype.HOM): 2.895,
                    },
                ),
                BreedingScheme("wt_hom", "hom_wt"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 5.79,
                    },
                ),
                BreedingScheme("wt_hom", "hom_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 2.895,
                        (Genotype.HET, Genotype.HOM): 2.895,
                    },
                ),
                BreedingScheme("wt_hom", "hom_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HOM): 5.79,
                    },
                ),
                BreedingScheme("het_wt", "het_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.WT): 1.4475,
                        (Genotype.HOM, Genotype.HET): 0.72375,
                        (Genotype.HOM, Genotype.WT): 0.72375,
                        (Genotype.WT, Genotype.HET): 0.72375,
                        (Genotype.WT, Genotype.WT): 0.72375,
                    },
                ),
                BreedingScheme("het_wt", "het_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 2.895,
                        (Genotype.HOM, Genotype.HET): 1.4475,
                        (Genotype.WT, Genotype.HET): 1.4475,
                    },
                ),
                BreedingScheme("het_wt", "hom_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.WT): 1.4475,
                        (Genotype.HOM, Genotype.HET): 1.4475,
                        (Genotype.HOM, Genotype.WT): 1.4475,
                    },
                ),
                BreedingScheme("het_wt", "hom_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 2.895,
                        (Genotype.HOM, Genotype.HET): 2.895,
                    },
                ),
                BreedingScheme("het_het", "het_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.HOM): 0.72375,
                        (Genotype.HET, Genotype.WT): 0.72375,
                        (Genotype.HOM, Genotype.HET): 0.72375,
                        (Genotype.HOM, Genotype.HOM): 0.361875,
                        (Genotype.HOM, Genotype.WT): 0.361875,
                        (Genotype.WT, Genotype.HET): 0.72375,
                        (Genotype.WT, Genotype.HOM): 0.361875,
                        (Genotype.WT, Genotype.WT): 0.361875,
                    },
                ),
                BreedingScheme("het_het", "het_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.HOM): 1.4475,
                        (Genotype.HOM, Genotype.HET): 0.72375,
                        (Genotype.HOM, Genotype.HOM): 0.72375,
                        (Genotype.WT, Genotype.HET): 0.72375,
                        (Genotype.WT, Genotype.HOM): 0.72375,
                    },
                ),
                BreedingScheme("het_het", "hom_wt"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.WT): 1.4475,
                        (Genotype.HOM, Genotype.HET): 1.4475,
                        (Genotype.HOM, Genotype.WT): 1.4475,
                    },
                ),
                BreedingScheme("het_het", "hom_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.HOM): 0.72375,
                        (Genotype.HET, Genotype.WT): 0.72375,
                        (Genotype.HOM, Genotype.HET): 1.4475,
                        (Genotype.HOM, Genotype.HOM): 0.72375,
                        (Genotype.HOM, Genotype.WT): 0.72375,
                    },
                ),
                BreedingScheme("het_het", "hom_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.HOM): 1.4475,
                        (Genotype.HOM, Genotype.HET): 1.4475,
                        (Genotype.HOM, Genotype.HOM): 1.4475,
                    },
                ),
                BreedingScheme("het_hom", "het_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HOM): 2.895,
                        (Genotype.HOM, Genotype.HOM): 1.4475,
                        (Genotype.WT, Genotype.HOM): 1.4475,
                    },
                ),
                BreedingScheme("het_hom", "hom_wt"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 2.895,
                        (Genotype.HOM, Genotype.HET): 2.895,
                    },
                ),
                BreedingScheme("het_hom", "hom_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HET): 1.4475,
                        (Genotype.HET, Genotype.HOM): 1.4475,
                        (Genotype.HOM, Genotype.HET): 1.4475,
                        (Genotype.HOM, Genotype.HOM): 1.4475,
                    },
                ),
                BreedingScheme("het_hom", "hom_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HET, Genotype.HOM): 2.895,
                        (Genotype.HOM, Genotype.HOM): 2.895,
                    },
                ),
                BreedingScheme("hom_wt", "hom_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HOM, Genotype.HET): 2.895,
                        (Genotype.HOM, Genotype.WT): 2.895,
                    },
                ),
                BreedingScheme("hom_wt", "hom_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HOM, Genotype.HET): 5.79,
                    },
                ),
                BreedingScheme("hom_het", "hom_het"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HOM, Genotype.HET): 2.895,
                        (Genotype.HOM, Genotype.HOM): 1.4475,
                        (Genotype.HOM, Genotype.WT): 1.4475,
                    },
                ),
                BreedingScheme("hom_het", "hom_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HOM, Genotype.HET): 2.895,
                        (Genotype.HOM, Genotype.HOM): 2.895,
                    },
                ),
                BreedingScheme("hom_hom", "hom_hom"): ExpectedOffspring(
                    total_n=5.79,
                    n_per_genotype={
                        (Genotype.HOM, Genotype.HOM): 5.79,
                    },
                ),
            },
            {
                BreedingScheme("het", "het"): 274,
                BreedingScheme("hom", "hom"): 16,
                BreedingScheme("wt", "het"): 285,
            },
            id="2 mutations",
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
