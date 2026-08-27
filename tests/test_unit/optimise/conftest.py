import platform

import pytest

from oscar_colony.breeding_scheme import BreedingScheme, Genotype
from oscar_colony.optimise.estimate_offspring import ExpectedOffspring
from oscar_colony.optimise.surplus_summary import (
    GenotypeSurplus,
    SexSurplus,
    SurplusSummary,
)


@pytest.fixture
def required_n_per_genotype_1_mutation():
    return {
        (Genotype.WT,): 1176,
        (Genotype.HET,): 1558,
        (Genotype.HOM,): 471,
    }


@pytest.fixture
def n_matings_1_mutation():
    return {
        BreedingScheme("hom", "het"): 137,
        BreedingScheme("hom", "hom"): 16,
        BreedingScheme("wt", "het"): 422,
    }


@pytest.fixture
def offspring_per_scheme_1_mutation():
    return {
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
    }


@pytest.fixture
def surplus_1_mutation():
    return SurplusSummary(
        total_n=3208.5,
        total_n_surplus=3.5,
        surplus_per_genotype={
            (Genotype.HET,): GenotypeSurplus(
                total_n=1559.6100,
                total_n_surplus=1.610,
                percent_surplus=0.1032,
            ),
            (Genotype.HOM,): GenotypeSurplus(
                total_n=471.51,
                total_n_surplus=0.5100,
                percent_surplus=0.1082,
            ),
            (Genotype.WT,): GenotypeSurplus(
                total_n=1177.38,
                total_n_surplus=1.3800,
                percent_surplus=0.1172,
            ),
        },
    )


@pytest.fixture
def required_n_per_genotype_1_mutation_sex_split():
    return {
        (Genotype.HET,): (18, 6),
        (Genotype.HOM,): 10,
    }


@pytest.fixture
def n_matings_1_mutation_sex_split():
    return {BreedingScheme("het", "hom"): 10}


@pytest.fixture
def offspring_per_scheme_1_mutation_sex_split():
    """A single breeding scheme, where 60% of each litter is male."""

    return {
        BreedingScheme("het", "hom"): ExpectedOffspring(
            total_n=6,
            proportion_male=0.6,
            n_per_genotype={
                (Genotype.HET,): 3,
                (Genotype.HOM,): 3,
            },
        )
    }


@pytest.fixture
def surplus_1_mutation_sex_split():
    """Surplus from 10 matings giving 30 het and 30 hom offspring, 60% of
    them male. Only het was required as a split of (n_males, n_females), so
    only it has a male / female surplus.
    """

    return SurplusSummary(
        total_n=60,
        total_n_surplus=26,
        sex_surplus=SexSurplus(
            total_males=18.0,
            total_male_surplus=0.0,
            male_percent_surplus=0.0,
            total_females=12.0,
            total_female_surplus=6.0,
            female_percent_surplus=50.0,
        ),
        surplus_per_genotype={
            (Genotype.HET,): GenotypeSurplus(
                total_n=30,
                total_n_surplus=6,
                percent_surplus=20.0,
                sex_surplus=SexSurplus(
                    total_males=18.0,
                    total_male_surplus=0.0,
                    male_percent_surplus=0.0,
                    total_females=12.0,
                    total_female_surplus=6.0,
                    female_percent_surplus=50.0,
                ),
            ),
            (Genotype.HOM,): GenotypeSurplus(
                total_n=30,
                total_n_surplus=20,
                percent_surplus=66.6667,
            ),
        },
    )


@pytest.fixture
def required_n_per_genotype_2_mutations():
    return {
        (Genotype.HET, Genotype.HET): 380,
        (Genotype.HET, Genotype.HOM): 18,
        (Genotype.HET, Genotype.WT): 23,
        (Genotype.HOM, Genotype.HET): 49,
        (Genotype.HOM, Genotype.HOM): 1,
        (Genotype.HOM, Genotype.WT): 20,
        (Genotype.WT, Genotype.HET): 48,
        (Genotype.WT, Genotype.HOM): 7,
        (Genotype.WT, Genotype.WT): 4,
    }


@pytest.fixture
def n_matings_2_mutations():
    """The solution returned by _optimise_n_matings differs on Mac vs other
    operating systems.

    Both solutions below are equally optimal (same total number of matings),
    but it's likely small differences in the compiler used etc. on Mac cause
    scipy.milp to make slightly different decisions and return a slightly
    different solution.
    """

    if platform.system() == "Darwin":  # check if running on Mac
        return {
            BreedingScheme("het_het", "hom_het"): 2,
            BreedingScheme("het_wt", "het_het"): 6,
            BreedingScheme("het_wt", "het_hom"): 20,
            BreedingScheme("het_wt", "hom_het"): 8,
            BreedingScheme("hom_wt", "hom_het"): 1,
            BreedingScheme("wt_het", "hom_het"): 1,
            BreedingScheme("wt_hom", "het_het"): 5,
            BreedingScheme("wt_hom", "hom_het"): 3,
            BreedingScheme("wt_wt", "het_hom"): 3,
            BreedingScheme("wt_wt", "hom_hom"): 47,
        }
    else:
        return {
            BreedingScheme("het_hom", "hom_het"): 1,
            BreedingScheme("het_wt", "het_het"): 6,
            BreedingScheme("het_wt", "het_hom"): 20,
            BreedingScheme("het_wt", "hom_het"): 9,
            BreedingScheme("hom_wt", "hom_het"): 1,
            BreedingScheme("wt_het", "hom_het"): 1,
            BreedingScheme("wt_hom", "het_het"): 5,
            BreedingScheme("wt_hom", "hom_het"): 3,
            BreedingScheme("wt_wt", "het_hom"): 3,
            BreedingScheme("wt_wt", "hom_hom"): 47,
        }


@pytest.fixture
def offspring_per_scheme_2_mutations():
    return {
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
                (Genotype.HET, Genotype.WT): 2.895,
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
    }


@pytest.fixture
def surplus_2_mutations():
    return SurplusSummary(
        total_n=555.8400,
        total_n_surplus=5.8400,
        surplus_per_genotype={
            (Genotype.HET, Genotype.HET): GenotypeSurplus(
                total_n=380.6925,
                total_n_surplus=0.6925,
                percent_surplus=0.1819,
            ),
            (Genotype.HET, Genotype.HOM): GenotypeSurplus(
                total_n=18.8175,
                total_n_surplus=0.8175,
                percent_surplus=4.3444,
            ),
            (Genotype.HOM, Genotype.HET): GenotypeSurplus(
                total_n=50.6625,
                total_n_surplus=1.6625,
                percent_surplus=3.2815,
            ),
            (Genotype.HOM, Genotype.HOM): GenotypeSurplus(
                total_n=1.4475,
                total_n_surplus=0.4475,
                percent_surplus=30.9154,
            ),
            (Genotype.HET, Genotype.WT): GenotypeSurplus(
                total_n=23.16,
                total_n_surplus=0.1600,
                percent_surplus=0.6908,
            ),
            (Genotype.HOM, Genotype.WT): GenotypeSurplus(
                total_n=20.265,
                total_n_surplus=0.2650,
                percent_surplus=1.3077,
            ),
            (Genotype.WT, Genotype.HET): GenotypeSurplus(
                total_n=49.2150,
                total_n_surplus=1.2150,
                percent_surplus=2.4688,
            ),
            (Genotype.WT, Genotype.WT): GenotypeSurplus(
                total_n=4.3425,
                total_n_surplus=0.3425,
                percent_surplus=7.8872,
            ),
            (Genotype.WT, Genotype.HOM): GenotypeSurplus(
                total_n=7.2375,
                total_n_surplus=0.2375,
                percent_surplus=3.2815,
            ),
        },
    )


@pytest.fixture
def surplus_2_mutations_sex_split():
    """Surplus where (n_males, n_females) was required for het_het
    (190, 190) and hom_het (25, 25). All other genotypes were required as a
    total only, so have no male / female surplus.
    """

    return SurplusSummary(
        total_n=555.8400,
        total_n_surplus=4.8400,
        sex_surplus=SexSurplus(
            total_males=215.6775,
            total_male_surplus=0.6775,
            male_percent_surplus=0.3141,
            total_females=215.6775,
            total_female_surplus=0.6775,
            female_percent_surplus=0.3141,
        ),
        surplus_per_genotype={
            (Genotype.HET, Genotype.HET): GenotypeSurplus(
                total_n=380.6925,
                total_n_surplus=0.6925,
                percent_surplus=0.1819,
                sex_surplus=SexSurplus(
                    total_males=190.34625,
                    total_male_surplus=0.34625,
                    male_percent_surplus=0.1819,
                    total_females=190.34625,
                    total_female_surplus=0.34625,
                    female_percent_surplus=0.1819,
                ),
            ),
            (Genotype.HET, Genotype.HOM): GenotypeSurplus(
                total_n=18.8175,
                total_n_surplus=0.8175,
                percent_surplus=4.3444,
            ),
            (Genotype.HOM, Genotype.HET): GenotypeSurplus(
                total_n=50.6625,
                total_n_surplus=0.6625,
                percent_surplus=1.3077,
                sex_surplus=SexSurplus(
                    total_males=25.33125,
                    total_male_surplus=0.33125,
                    male_percent_surplus=1.3077,
                    total_females=25.33125,
                    total_female_surplus=0.33125,
                    female_percent_surplus=1.3077,
                ),
            ),
            (Genotype.HOM, Genotype.HOM): GenotypeSurplus(
                total_n=1.4475,
                total_n_surplus=0.4475,
                percent_surplus=30.9154,
            ),
            (Genotype.HET, Genotype.WT): GenotypeSurplus(
                total_n=23.16,
                total_n_surplus=0.1600,
                percent_surplus=0.6908,
            ),
            (Genotype.HOM, Genotype.WT): GenotypeSurplus(
                total_n=20.265,
                total_n_surplus=0.2650,
                percent_surplus=1.3077,
            ),
            (Genotype.WT, Genotype.HET): GenotypeSurplus(
                total_n=49.2150,
                total_n_surplus=1.2150,
                percent_surplus=2.4688,
            ),
            (Genotype.WT, Genotype.WT): GenotypeSurplus(
                total_n=4.3425,
                total_n_surplus=0.3425,
                percent_surplus=7.8872,
            ),
            (Genotype.WT, Genotype.HOM): GenotypeSurplus(
                total_n=7.2375,
                total_n_surplus=0.2375,
                percent_surplus=3.2815,
            ),
        },
    )
