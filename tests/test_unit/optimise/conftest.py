import pytest

from oscar.breeding_scheme import BreedingScheme, Genotype
from oscar.optimise.estimate_offspring import ExpectedOffspring


@pytest.fixture
def single_mutation_required_n_per_genotype():
    return {
        (Genotype.WT,): 1176,
        (Genotype.HET,): 1558,
        (Genotype.HOM,): 471,
    }


@pytest.fixture
def single_mutation_offspring_per_scheme():
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
def single_mutation_n_matings():
    return {
        BreedingScheme("het", "het"): 274,
        BreedingScheme("hom", "hom"): 16,
        BreedingScheme("wt", "het"): 285,
    }
