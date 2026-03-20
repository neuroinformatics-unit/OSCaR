from dataclasses import dataclass, field

from oscar.breeding_scheme import (
    BreedingScheme,
    Genotype,
)
from oscar.optimise.estimate_offspring import ExpectedOffspring


@dataclass
class GenotypeSurplus:
    """Summary of surplus for a single genotype"""

    total_n: float = 0
    total_n_surplus: float = 0
    percent_surplus: float = 0


@dataclass
class SurplusSummary:
    """Summary of surplus across all genotypes"""

    total_n: float = 0
    total_n_surplus: float = 0

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

    # Get number of expected offspring overall / per genotype
    for breeding_scheme, n_matings in n_matings_per_scheme.items():
        expected_offspring = offspring_per_scheme[breeding_scheme]
        surplus_summary.total_n += expected_offspring.total_n * n_matings

        n_per_genotype = expected_offspring.n_per_genotype
        for genotype, n_per_mating in n_per_genotype.items():
            if genotype not in surplus_per_genotype:
                surplus_per_genotype[genotype] = GenotypeSurplus()
            surplus_per_genotype[genotype].total_n += n_per_mating * n_matings

    # Calculate total surplus
    total_required = sum(
        [required_n for required_n in required_n_per_genotype.values()]
    )
    surplus_summary.total_n_surplus = surplus_summary.total_n - total_required

    # Calculate surplus per genotype
    for genotype, surplus in surplus_per_genotype.items():
        if genotype in required_n_per_genotype:
            required_n = required_n_per_genotype[genotype]
        else:
            required_n = 0

        surplus.total_n_surplus = surplus.total_n - required_n
        surplus.percent_surplus = (
            surplus.total_n_surplus / surplus.total_n
        ) * 100

    return surplus_summary
