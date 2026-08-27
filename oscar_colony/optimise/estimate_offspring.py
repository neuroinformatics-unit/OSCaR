from dataclasses import dataclass, field

from oscar_colony.breeding_scheme import (
    BreedingScheme,
    Genotype,
    generate_breeding_schemes,
)
from oscar_colony.historical_stats import LineStatistics


@dataclass
class ExpectedOffspring:
    """Summary of expected average number of offspring from a single mating"""

    total_n: float = 0  # litter size
    proportion_male: float = 0.5  # proportion of the litter that is male

    n_per_genotype: dict[tuple[Genotype, ...], float] = field(
        default_factory=dict
    )


def estimate_n_offspring_per_mating(
    line_stats: LineStatistics,
    default_litter_size: int,
    min_n_matings: int = 3,
    min_n_offspring: int = 10,
) -> dict[BreedingScheme, ExpectedOffspring]:
    """For all possible breeding schemes for the given line, estimate the
    number of offspring produced per mating.

    Calculates the total number across all genotypes == the litter size. As
    well as the expected number per offspring genotype, and the proportion
    of the litter expected to be male.

    Parameters
    ----------
    line_stats : LineStatistics
        Statistics from historical data for the line
    default_litter_size: int
        The default value used for average litter size if there isn't enough
        historical data for the line. This should usually be set to the average
        litter size across all available data for all lines.
    min_n_matings : int, optional
        Minimum number of successful matings required to use the measured
        litter size from line_stats. If there aren't enough matings for a
        specific breeding scheme, the average of the whole line will be used
        instead. If the whole line also doesn't have enough matings, then
        default_litter_size is used.
    min_n_offspring: int, optional
        Minimum number of offspring required from a breeding scheme to use
        its measured proportion of each genotype / of males from line_stats.

    Returns
    -------
    dict[BreedingScheme, ExpectedOffspring]
        Returns a dict mapping each breeding scheme to the expected number
        of offspring it will produce per mating.
    """

    # The breeding schemes listed in line_stats.stats_per_breeding_scheme are
    # sparse i.e. only breeding schemes that appeared in the historical data
    # are included.
    # Here we need ALL possible schemes:
    n_mutations = line_stats.n_mutations
    breeding_schemes = generate_breeding_schemes(n_mutations)

    expected_offspring_per_scheme = {}
    for breeding_scheme in breeding_schemes:
        expected_offspring = ExpectedOffspring()
        expected_offspring.total_n = _expected_litter_size(
            breeding_scheme, line_stats, min_n_matings, default_litter_size
        )

        _expected_proportion_of_males(
            breeding_scheme, line_stats, min_n_offspring, expected_offspring
        )

        proportion_per_genotype = _expected_proportion_per_genotype(
            breeding_scheme, line_stats, min_n_offspring
        )

        for genotype, proportion in proportion_per_genotype.items():
            expected_n = proportion * expected_offspring.total_n
            if expected_n > 0:
                expected_offspring.n_per_genotype[genotype] = expected_n

        expected_offspring_per_scheme[breeding_scheme] = expected_offspring

    return expected_offspring_per_scheme


def _expected_litter_size(
    breeding_scheme: BreedingScheme,
    line_stats: LineStatistics,
    min_n_matings: int,
    default_litter_size: float,
) -> float:
    """Create an estimated average litter size (total number of individuals
    produced from one mating).

    Parameters
    ----------
    breeding_scheme : BreedingScheme
        The breeding scheme to return litter size for
    line_stats : LineStatistics
        Summary statistics from historical data for the line
    min_n_matings : int
        Minimum number of successful matings required to use the measured
        litter size from line_stats. If there aren't enough matings for a
        specific breeding scheme, the average of the whole line will be used
        instead. If the whole line also doesn't have enough matings, then
        default_litter_size is used.
    default_litter_size : float
        Default litter size to fallback to if there's not enough historical
        data for the line.

    Returns
    -------
    float
        Expected litter size
    """

    scheme_stats = line_stats.stats_per_breeding_scheme.get(
        breeding_scheme, None
    )

    if (scheme_stats is not None) and (
        scheme_stats.n_successful_matings >= min_n_matings
    ):
        return scheme_stats.average_litter_size
    elif line_stats.total_n_successful_matings >= min_n_matings:
        return line_stats.average_litter_size
    else:
        return default_litter_size


def _expected_proportion_per_genotype(
    breeding_scheme: BreedingScheme,
    line_stats: LineStatistics,
    minimum_n_offspring: int,
) -> dict[tuple[Genotype, ...], float]:
    """Calculate the expected proportion of offspring of each genotype.

    If enough historical data is available in line_stats, the measured
    proportion will be used. Otherwise, the theoretical mendelian ratio
    will be returned.

    Parameters
    ----------
    breeding_scheme: BreedingScheme
        Breeding scheme to return proportions for
    line_stats : LineStatistics
        Summary line statistics from historical data
    minimum_n_offspring: int
        The minimum number of offspring required for this breeding scheme to
        use the genotyping ratio (measured from historical data). Otherwise,
        defaults to theoretical mendelian ratio.

    Returns
    -------
    dict[tuple[Genotype, ...], float]
        A dictionary mapping offspring genotypes to the expected proportion
        of that type
    """

    scheme_stats = line_stats.stats_per_breeding_scheme.get(
        breeding_scheme, None
    )

    # If there's enough recorded offspring, use the observed proportion
    # from historical data
    if (scheme_stats is not None) and (
        scheme_stats.total_n_genotyped_offspring >= minimum_n_offspring
    ):
        return scheme_stats.proportion_offspring_per_genotype

    else:
        return breeding_scheme.mendelian_ratio()


def _expected_proportion_of_males(
    breeding_scheme: BreedingScheme,
    line_stats: LineStatistics,
    min_n_offspring: int,
    expected_offspring: ExpectedOffspring,
) -> None:
    """Update ExpectedOffspring with the expected proportion of male offspring.

    If enough historical data is available in line_stats, the measured
    proportion will be used - first for the specific breeding scheme, then
    for the whole line. Otherwise, defaults to 0.5.

    Parameters
    ----------
    breeding_scheme: BreedingScheme
        Breeding scheme to retrieve summary stats from
    line_stats : LineStatistics
        Summary line statistics from historical data
    min_n_offspring: int
        The minimum number of sexed offspring required to use the measured
        proportion of males (from historical data). Otherwise, defaults
        to 0.5.
    expected_offspring : ExpectedOffspring
        Class to set the expected proportion of males on
    """

    scheme_stats = line_stats.stats_per_breeding_scheme.get(
        breeding_scheme, None
    )
    if (scheme_stats is not None) and (
        scheme_stats.total_n_sexed_offspring >= min_n_offspring
    ):
        expected_offspring.proportion_male = scheme_stats.proportion_male
    elif line_stats.total_n_sexed_offspring >= min_n_offspring:
        expected_offspring.proportion_male = line_stats.proportion_male
    else:
        expected_offspring.proportion_male = 0.5
