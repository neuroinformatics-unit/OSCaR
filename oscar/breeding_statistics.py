import itertools
from enum import Enum

import numpy as np


class Genotype(Enum):
    """Genotype status: homozygous (HOM), heterozygous (HET) or WT (wildtype).

    Each mouse will have two copies a particular gene - each being either
    wildtype or mutated. The value of the enum is the number of wildtype copies
    for that genotype.
    """

    HOM = 0
    HET = 1
    WT = 2


class BreedingScheme:
    def __init__(
        self,
        parent_1_genotype: tuple[Genotype, ...],
        parent_2_genotype: tuple[Genotype, ...],
    ):
        if len(parent_1_genotype) != len(parent_2_genotype):
            raise ValueError(
                "Both parents must have a genotype of the same length"
            )

        self.parent_1_genotype = parent_1_genotype
        self.parent_2_genotype = parent_2_genotype
        self.n_mutations = len(self.parent_1_genotype)

    def __eq__(self, other):
        # The order of parent 1 vs parent 2 doesn't matter. Breeding
        # schemes are equal if they are combining the same two genotypes
        # in any order.
        return sorted(
            [self.parent_1_genotype, self.parent_2_genotype]
        ) == sorted([other.parent_1_genotype, other.parent_2_genotype])

    def mendelian_ratio(self) -> dict[Genotype, float]:
        """Calculate the theoretical mendelian ratio for this breeding scheme.

        Returns
        -------
        dict[tuple[Genotype, ...], float]
            Returns a dictionary with keys being the possible genotypes of
            offspring, and values the expected proportion of offspring of
            that genotype.
        """

        # first gene
        parent_1 = self.parent_1_genotype[0]
        parent_2 = self.parent_2_genotype[0]

        # Each parent has 2 copies of that gene. For wt = both true. For
        # het = one true, one false. For hom = both false.
        alleles = []
        for parent in [parent_1, parent_2]:
            parent_alleles = [True] * parent.value
            parent_alleles.extend([False] * (2 - parent.value))
            alleles.append(parent_alleles)

        offspring_combos = itertools.product(alleles[0], alleles[1])

        total_sum = 0
        genotype_sums: dict[Genotype, int] = {}

        for combo in offspring_combos:
            total_sum += 1
            n_wt_copies = np.array(combo).sum()
            genotype = Genotype(n_wt_copies)

            if genotype in genotype_sums:
                genotype_sums[genotype] += 1
            else:
                genotype_sums[genotype] = 1

        mendelian_ratios = {}
        for genotype, n_with_genotype in genotype_sums.items():
            mendelian_ratios[genotype] = n_with_genotype / total_sum

        return mendelian_ratios


def generate_breeding_schemes(
    n_mutations: int,
) -> list[list[tuple[Genotype, ...]]]:
    """Generate all possible combinations of parent1 x parent2 genotype, for
    the given number of mutations.

    Whether mutations come from the father or mother is ignored e.g.
    wt x hom (father x mother) is identical to hom x wt (father x mother), and
    only one of these combinations will be returned.

    Parameters
    ----------
    n_mutations : int
        Number of mutations

    Returns
    -------
    list[list[tuple[Genotype, ...]]]
        A list where each item is a possible breeding scheme. Each item is
        a list of form: [(parent_1_genotype), (parent_2_genotype)]
    """

    breeding_schemes = []

    # First, generate all possible genotypes of a single parent.
    parent_genotypes = generate_genotypes(n_mutations)

    # Then combine two parents, bearing in mind order doesn't matter e.g.
    # wt x hom == hom x wt
    breeding_combos = itertools.combinations_with_replacement(
        parent_genotypes, 2
    )

    for combo in breeding_combos:
        parent_1 = combo[0]
        parent_2 = combo[1]

        # exclude combos that have wt for the same allele in both parents.
        # e.g. for a 2 mutations scenario, we could generate wt, het x wt, hom.
        # Here, both parents are wt for allele 1, so this is actually a 1
        # mutation scenario and should be excluded.
        if not _breeding_scheme_contains_wt_pairs(parent_1, parent_2):
            breeding_schemes.append([parent_1, parent_2])

    return breeding_schemes


def generate_genotypes(n_mutations: int) -> itertools.product:
    """Generate all possible genotypes of a single mouse, with the given
    number of mutations.

    Parameters
    ----------
    n_mutations : int
        Number of mutations

    Returns
    -------
    itertools.product
        An iterable of all possible genotypes
    """
    single_genotypes = [Genotype.WT, Genotype.HOM, Genotype.HET]
    return itertools.product(single_genotypes, repeat=n_mutations)


def _breeding_scheme_contains_wt_pairs(
    parent_1_genotype: tuple[Genotype, ...],
    parent_2_genotype: tuple[Genotype, ...],
) -> bool:
    """Check if a given breeding scheme (parent_1 x parent_2) contains pairs of
    wt alleles.

    E.g. wt, het x wt, hom would return true as both parents are wt for the
    first allele.

    Parameters
    ----------
    parent_1_genotype : tuple[Genotype]
        Genotype of parent 1 - tuple of one or more Genotype
        (depending on number of mutations)
    parent_2_genotype : tuple[Genotype]
        Genotype of parent 2 - tuple of one or more Genotype
        (depending on number of mutations)

    Returns
    -------
    bool
        Whether wt pairs are present.
    """

    for parent_1_allele, parent_2_allele in zip(
        parent_1_genotype, parent_2_genotype
    ):
        if parent_1_allele == Genotype.WT and parent_2_allele == Genotype.WT:
            return True

    return False
