import math
from dataclasses import dataclass, field

from oscar.breeding_scheme import (
    BreedingScheme,
    Genotype,
)
from oscar.optimise.estimate_offspring import ExpectedOffspring


@dataclass
class GenotypeSurplus:
    """Summary of surplus for a single genotype"""

    total_n: int = 0
    total_n_surplus: int = 0
    percent_surplus: float = 0


@dataclass
class SurplusSummary:
    """Summary of surplus across all genotypes"""

    total_n: int = 0
    total_n_surplus: int = 0

    surplus_per_genotype: dict[tuple[Genotype, ...], GenotypeSurplus] = field(
        default_factory=dict
    )


def create_surplus_summary(
    required_n_per_genotype: dict[tuple[Genotype, ...], int],
    n_matings_per_scheme: dict[BreedingScheme, int],
    offspring_per_scheme: dict[BreedingScheme, ExpectedOffspring],
) -> SurplusSummary:
    """Create a summary of the total and surplus numbers for the
    given combination of breeding schemes.

    Parameters
    ----------
    required_n_per_genotype : dict[tuple[Genotype, ...], int]
        Required number of individuals of each genotype
    n_matings_per_scheme : dict[BreedingScheme, int]
        Optimal number of matings per breeding scheme
    offspring_per_scheme : dict[BreedingScheme, ExpectedOffspring]
        The estimated number of offspring produced per mating of each
        breeding scheme

    Returns
    -------
    SurplusSummary
        Summary of total and surplus numbers
    """
    surplus_summary = SurplusSummary()
    surplus_per_genotype = surplus_summary.surplus_per_genotype

    for breeding_scheme, n_matings in n_matings_per_scheme.items():
        n_per_genotype = offspring_per_scheme[breeding_scheme].n_per_genotype

        total_for_scheme = 0
        for genotype, n_per_mating in n_per_genotype.items():
            # round value up to the nearest int - as we can't have
            # partial animals
            total_n = math.ceil(n_per_mating * n_matings)
            total_for_scheme += total_n

            if genotype not in surplus_per_genotype:
                surplus_per_genotype[genotype] = GenotypeSurplus()
            surplus_per_genotype[genotype].total_n += total_n

        surplus_summary.total_n += total_for_scheme

    total_required = sum(
        [required_n for required_n in required_n_per_genotype.values()]
    )
    surplus_summary.total_n_surplus = surplus_summary.total_n - total_required

    for genotype, surplus in surplus_per_genotype.items():
        total_surplus = surplus.total_n - required_n_per_genotype[genotype]
        surplus.total_n_surplus = total_surplus
        surplus.percent_surplus = (total_surplus / surplus.total_n) * 100

    return surplus_summary
