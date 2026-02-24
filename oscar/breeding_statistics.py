import itertools
from enum import StrEnum, auto


class Genotype(StrEnum):
    WT = "wildtype"
    HOM = "homozygous"
    HET = "heterozygous"
    UNGENOTYPED = auto()


def generate_breeding_schemes(
    number_of_mutations: int,
) -> list[list[tuple[Genotype, ...]]]:
    """Generate all possible combinations of parent1 x parent2 genotype, for
    the given number of mutations.

    Whether mutations come from the father or mother is ignored e.g.
    wt x hom (father x mother) is identical to hom x wt (father x mother), and
    only one of these combinations will be returned.

    Parameters
    ----------
    number_of_mutations : int
        Number of mutations
    """

    genotypes = [Genotype.WT, Genotype.HOM, Genotype.HET]
    breeding_schemes = []

    # First, generate all possible genotypes of a single parent.
    parent_genotype = itertools.product(genotypes, repeat=number_of_mutations)

    # Then combine two parents, bearing in mind order doesn't matter e.g.
    # wt x hom == hom x wt
    breeding_combos = itertools.combinations_with_replacement(
        parent_genotype, 2
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
