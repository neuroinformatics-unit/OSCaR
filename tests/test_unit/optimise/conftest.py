import pytest

from oscar.breeding_scheme import BreedingScheme, Genotype
from oscar.optimise.estimate_offspring import ExpectedOffspring
from oscar.optimise.surplus_summary import GenotypeSurplus, SurplusSummary


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
        total_n=3212,
        total_n_surplus=7,
        surplus_per_genotype={
            (Genotype.HET,): GenotypeSurplus(
                total_n=1561,
                total_n_surplus=3,
                percent_surplus=pytest.approx((3 / 1561) * 100),
            ),
            (Genotype.HOM,): GenotypeSurplus(
                total_n=473,
                total_n_surplus=2,
                percent_surplus=pytest.approx((2 / 473) * 100),
            ),
            (Genotype.WT,): GenotypeSurplus(
                total_n=1178,
                total_n_surplus=2,
                percent_surplus=pytest.approx((2 / 1178) * 100),
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
        total_n=pytest.approx(555.8400, abs=1e-4),
        total_n_surplus=pytest.approx(5.8400, abs=1e-4),
        surplus_per_genotype={
            (Genotype.HET, Genotype.HET): GenotypeSurplus(
                total_n=pytest.approx(380.6925, abs=1e-4),
                total_n_surplus=pytest.approx(0.6925, abs=1e-4),
                percent_surplus=pytest.approx(0.1819, abs=1e-4),
            ),
            (Genotype.HET, Genotype.HOM): GenotypeSurplus(
                total_n=pytest.approx(18.8175, abs=1e-4),
                total_n_surplus=pytest.approx(0.8175, abs=1e-4),
                percent_surplus=pytest.approx(4.3444, abs=1e-4),
            ),
            (Genotype.HOM, Genotype.HET): GenotypeSurplus(
                total_n=pytest.approx(50.6625, abs=1e-4),
                total_n_surplus=pytest.approx(1.6625, abs=1e-4),
                percent_surplus=pytest.approx(3.2815, abs=1e-4),
            ),
            (Genotype.HOM, Genotype.HOM): GenotypeSurplus(
                total_n=pytest.approx(1.4475, abs=1e-4),
                total_n_surplus=pytest.approx(0.4475, abs=1e-4),
                percent_surplus=pytest.approx(30.9154, abs=1e-4),
            ),
            (Genotype.HET, Genotype.WT): GenotypeSurplus(
                total_n=pytest.approx(23.16, abs=1e-4),
                total_n_surplus=pytest.approx(0.1600, abs=1e-4),
                percent_surplus=pytest.approx(0.6908, abs=1e-4),
            ),
            (Genotype.HOM, Genotype.WT): GenotypeSurplus(
                total_n=pytest.approx(20.265, abs=1e-4),
                total_n_surplus=pytest.approx(0.2650, abs=1e-4),
                percent_surplus=pytest.approx(1.3077, abs=1e-4),
            ),
            (Genotype.WT, Genotype.HET): GenotypeSurplus(
                total_n=pytest.approx(49.2150, abs=1e-4),
                total_n_surplus=pytest.approx(1.2150, abs=1e-4),
                percent_surplus=pytest.approx(2.4688, abs=1e-4),
            ),
            (Genotype.WT, Genotype.WT): GenotypeSurplus(
                total_n=pytest.approx(4.3425, abs=1e-4),
                total_n_surplus=pytest.approx(0.3425, abs=1e-4),
                percent_surplus=pytest.approx(7.8872, abs=1e-4),
            ),
            (Genotype.WT, Genotype.HOM): GenotypeSurplus(
                total_n=pytest.approx(7.2375, abs=1e-4),
                total_n_surplus=pytest.approx(0.2375, abs=1e-4),
                percent_surplus=pytest.approx(3.2815, abs=1e-4),
            ),
        },
    )
