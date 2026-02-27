import itertools
from enum import Enum

import numpy as np


class Genotype(Enum):
    """Genotype status: homozygous (HOM), heterozygous (HET) or WT (wildtype).

    Each animal will have two copies (alleles) of a particular gene - each
    being either wildtype or mutated. The value of the enum is the number of
    wildtype copies for that genotype.
    """

    HOM = 0
    HET = 1
    WT = 2


class BreedingScheme:
    """
    Class representing a particular breeding scheme, i.e. the breeding
    of two specific parent genotypes.
    """

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
        return set([self.parent_1_genotype, self.parent_2_genotype]) == set(
            [other.parent_1_genotype, other.parent_2_genotype]
        )

    def mendelian_ratio(self) -> dict[tuple[Genotype, ...], float]:
        """Calculate the theoretical mendelian ratio for this breeding scheme.

        Returns
        -------
        dict[tuple[Genotype, ...], float]
            Returns a dictionary with keys being the possible genotypes of
            offspring, and values the expected proportion of offspring of
            that genotype.
        """

        # For each parent, determine all the combinations of alleles they
        # could pass on to their offspring. This is equivalent to making the
        # header for that parent in a punnet square.
        parent_1_alleles = self._parent_allele_combos(self.parent_1_genotype)
        parent_2_alleles = self._parent_allele_combos(self.parent_2_genotype)

        # All combinations of parent 1 x parent 2. This is equivalent to the
        # contents of a punnet square.
        offspring_combos = itertools.product(
            parent_1_alleles, parent_2_alleles
        )

        # Calculate total numbers for all offspring genotypes
        total_sum = 0
        genotype_sums: dict[tuple[Genotype, ...], int] = {}

        for combo in offspring_combos:
            # each combo is a tuple of two items: the first being the alleles
            # from parent 1, the second being the alleles from parent 2
            offspring_genotype = self._determine_offspring_genotype(
                combo[0], combo[1]
            )

            total_sum += 1

            if offspring_genotype in genotype_sums:
                genotype_sums[offspring_genotype] += 1
            else:
                genotype_sums[offspring_genotype] = 1

        # calculate overall proportion of offspring of each type - termed
        # the 'mendelian ratio'
        mendelian_ratios = {}
        for genotype, n_with_genotype in genotype_sums.items():
            mendelian_ratios[genotype] = n_with_genotype / total_sum

        return mendelian_ratios

    def _parent_allele_combos(
        self, parent_genotype: tuple[Genotype, ...]
    ) -> itertools.product:
        """For a parent, determine all the combinations of alleles they
        could pass on to their offspring.

        Bear in mind that each parent has 2 alleles for each gene, and will
        pass on 1 to their offspring.

        You can think of this like n separate containers (n == number of
        mutations), each with 2 items inside (the alleles). We must generate
        all combos of picking one item from each container.

        Parameters
        ----------
        parent_genotype : tuple[Genotype, ...]
            The genotype of the parent

        Returns
        -------
        itertools.product
            An iterable of all allele combinations they could pass to their
            offpsring. Each item is a list of bools with length == the number
            of mutations. e.g. an item of [true, false, true] would mean they
            are wildtype for gene 1, mutant for gene 2 and wildtype for gene 3
        """

        alleles = []
        for gene in parent_genotype:
            alleles.append(self._alleles_for_genotype(gene))

        return itertools.product(*alleles)

    def _determine_offspring_genotype(
        self,
        parent_1_alleles: tuple[bool, ...],
        parent_2_alleles: tuple[bool, ...],
    ) -> tuple[Genotype, ...]:
        """Determine the genotype of the offspring, based on the alleles
        it inherited from each parent.

        Parameters
        ----------
        parent_1_alleles : tuple[bool, ...]
            Alleles inherited from parent 1. True = wildtype, False = mutated.
        parent_2_alleles : tuple[bool, ...]
            Alleles inherited from parent 2. True = wildtype, False = mutated.

        Returns
        -------
        tuple[Genotype, ...]
            The genotype of the offspring. The tuple will have length ==
            the number of mutations.
        """

        offspring_genotype = []

        # Loop through each gene, and combine the allele from parent 1 with
        # that from parent 2. The result may be wt, het or hom for each gene.
        for parent_1_allele, parent_2_allele in zip(
            parent_1_alleles, parent_2_alleles
        ):
            n_wt_copies = np.array([parent_1_allele, parent_2_allele]).sum()
            offspring_genotype.append(Genotype(n_wt_copies))

        return tuple(offspring_genotype)

    def _alleles_for_genotype(self, genotype: Genotype) -> list[bool]:
        """Return the 2 alleles for the given genotype.

        Parameters
        ----------
        genotype : Genotype
            The genotype

        Returns
        -------
        list[bool]
            The 2 alleles: True = wildtype, False = mutated.
        """
        alleles = [True] * genotype.value
        alleles.extend([False] * (2 - genotype.value))
        return alleles


def generate_breeding_schemes(
    n_mutations: int,
) -> list[BreedingScheme]:
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
    list[BreedingScheme]
        A list where each item is a possible breeding scheme (a cross between
        two specific parent genotypes)
    """

    breeding_schemes = []

    # First, generate all possible genotypes of a single parent.
    single_genotypes = [Genotype.WT, Genotype.HOM, Genotype.HET]
    parent_genotypes = itertools.product(single_genotypes, repeat=n_mutations)

    # Then combine two parents, bearing in mind order doesn't matter e.g.
    # wt x hom == hom x wt
    breeding_combos = itertools.combinations_with_replacement(
        parent_genotypes, 2
    )

    for combo in breeding_combos:
        parent_1 = combo[0]
        parent_2 = combo[1]

        # exclude combos that have wt for the same gene in both parents.
        # e.g. for a 2 mutations scenario, we could generate wt, het x wt, hom.
        # Here, both parents are wt for gene 1, so this is actually a 1
        # mutation scenario and should be excluded.
        if not _breeding_scheme_contains_wt_pairs(parent_1, parent_2):
            breeding_schemes.append(BreedingScheme(parent_1, parent_2))

    return breeding_schemes


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
