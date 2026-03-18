from dataclasses import dataclass, field

from oscar.breeding_scheme import (
    BreedingScheme,
    Genotype,
    generate_breeding_schemes,
)
from oscar.historical_stats import LineStatistics


@dataclass
class ExpectedOffspring:
    """Summary of expected average number of offspring from a single mating"""

    total_n: float = 0  # litter size

    n_per_genotype: dict[tuple[Genotype, ...], float] = field(
        default_factory=dict
    )


def estimate_n_offspring_per_mating(
    line_stats: LineStatistics,
    min_n_offspring: int,
    min_n_matings: int,
    default_litter_size: int,
) -> dict[BreedingScheme, ExpectedOffspring]:
    """For all possible breeding schemes for the given line, estimate the
    number of offspring produced per mating.

    Calculates the total number across all genotypes == the litter size. As
    well as the expected number per offspring genotype.

    Parameters
    ----------
    line_stats : LineStatistics
        Statistics from historical data for the line
    min_n_offspring: int, optional
        Minimum number of offspring required from a breeding scheme to use
        the measured proportion of each genotype from line_stats. If not met,
        the theoretical mendelian ratio will be used instead.
    min_n_matings : int, optional
        Minimum number of successful matings required to use the measured
        litter size from line_stats. If there aren't enough matings for a
        specific breeding scheme, the average of the whole line will be used
        instead. If the whole line also doesn't have enough matings, then
        default_litter_size is used.
    default_litter_size: float
        The default value used for average litter size if there isn't enough
        historical data for the line. This should usually be set to the average
        litter size across all available data for all lines.
    """

    # The breeding schemes listed in line_stats.stats_per_breeding_scheme are
    # sparse. i.e. only breeding schemes that appeared in the historical data
    # are included (rather than all that are possible)
    # Here we need ALL possible schemes:
    n_mutations = len(
        list(line_stats.total_n_offspring_per_genotype.keys())[0]
    )
    breeding_schemes = generate_breeding_schemes(n_mutations)

    expected_offspring_per_scheme = {}
    for breeding_scheme in breeding_schemes:
        expected_offspring = ExpectedOffspring()
        litter_size = _expected_litter_size(
            breeding_scheme, line_stats, min_n_matings, default_litter_size
        )
        expected_offspring.total_n = litter_size

        proportion_per_genotype = _expected_proportion_per_genotype(
            breeding_scheme, line_stats, min_n_offspring
        )

        for genotype, proportion in proportion_per_genotype.items():
            expected_n = proportion * litter_size
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
        Minimum number of successful matings required to use the
        measured litter size from line_stats.
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
        scheme_stats.total_n_offspring >= minimum_n_offspring
    ):
        return scheme_stats.proportion_offspring_per_genotype

    else:
        return breeding_scheme.mendelian_ratio()
