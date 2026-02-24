from enum import StrEnum, auto
from itertools import combinations_with_replacement


class Genotype(StrEnum):
    WT = "wildtype"
    HOM = "homozygous"
    HET = "heterozygous"
    UNGENOTYPED = auto()


def generate_breeding_schemes(number_of_mutations: int):
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

    # as order doesn't matter e.g. wt x hom == hom x wt, these are
    # mathematically 'combinations'
    breeding_combos = combinations_with_replacement(
        genotypes, 2 * number_of_mutations
    )

    for combo in breeding_combos:
        parent_1 = combo[:number_of_mutations]
        parent_2 = combo[number_of_mutations:]

        if not _breeding_scheme_contains_wt_pairs(parent_1, parent_2):
            breeding_schemes.append([parent_1, parent_2])

    return breeding_schemes


def _breeding_scheme_contains_wt_pairs(
    parent_1_genotype: tuple[Genotype, ...],
    parent_2_genotype: tuple[Genotype, ...],
):
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
