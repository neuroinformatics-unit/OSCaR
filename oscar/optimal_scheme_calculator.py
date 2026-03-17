import math
from dataclasses import dataclass, field

from scipy.optimize import linprog

from oscar.breeding_scheme import (
    BreedingScheme,
    Genotype,
    generate_breeding_schemes,
)
from oscar.historical_stats import LineStatistics


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


@dataclass
class ExpectedOffspring:
    """Summary of expected average number of offspring from a single mating"""

    total_n: float = 0  # litter size

    n_per_genotype: dict[tuple[Genotype, ...], float] = field(
        default_factory=dict
    )


def calculate_optimal_scheme(
    required_n_per_genotype: dict[tuple[Genotype, ...], int],
    line_stats: LineStatistics,
    default_litter_size: int,
    min_n_matings: int = 3,
    min_n_offspring: int = 10,
) -> tuple[dict[BreedingScheme, int], SurplusSummary]:
    """Calculate the optimal combination of breeding schemes to produce
    the required_n_per_genotype.

    Parameters
    ----------
    required_n_per_genotype : dict[tuple[Genotype, ...], int]
        Required number of individuals per genotype
    line_stats : LineStatistics
        Statistics from historical data for the line
    default_litter_size: float
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
        the measured proportion of each genotype from line_stats. If not met,
        the theoretical mendelian ratio will be used instead.

    Returns
    -------
    tuple[dict[BreedingScheme, int], SurplusSummary]
        Returns 2 items:
            - optimal number of matings per breeding scheme
            - a summary of the surplus numbers for this optimal scheme
              combination
    """

    offspring_per_scheme = _estimate_n_offspring_per_mating(
        line_stats, min_n_offspring, min_n_matings, default_litter_size
    )

    n_matings_per_scheme = _optimise_n_matings(
        required_n_per_genotype, offspring_per_scheme
    )

    surplus_summary = _create_surplus_summary(
        required_n_per_genotype, n_matings_per_scheme, offspring_per_scheme
    )

    return n_matings_per_scheme, surplus_summary


def _optimise_n_matings(
    required_n_per_genotype: dict[tuple[Genotype, ...], int],
    offspring_per_scheme: dict[BreedingScheme, ExpectedOffspring],
) -> dict[BreedingScheme, int]:
    """Calculate the optimal number of matings of each breeding scheme to
    meet the required_genotypes.

    This is framed as a linear programming (optimisation) problem with three
    parts: the decision variables (the number of matings for each scheme),
    the objective function (the surplus to minimise) and the constraints
    (that the final result must meet the given required_genotypes).

    The objective function is:
    total_surplus =
      sum(litter_size_per_scheme * n_matings_per_scheme) - total_n_required

    There is one constraint per required offspring genotype of form:
    sum(n_of_genotype_offspring_per_mating * n_matings_per_scheme)
      >= required_n_for_genotype

    Parameters
    ----------
    required_n_per_genotype : dict[tuple[Genotype, ...], int]
        The required number of individuals per genotype
    offspring_per_scheme : dict[BreedingScheme, ExpectedOffspring]
        The estimated number of offspring produced per mating of each
        breeding scheme
    """

    # Extract names of breeding schemes / required genotypes as a list, so we
    # can make sure we always iterate through them in the same order
    breeding_schemes = list(offspring_per_scheme.keys())
    required_genotypes = list(required_n_per_genotype.keys())

    # Coefficients of the objective function we want to minimise. This will
    # be a list with length == number of breeding schemes, containing the
    # expected litter size for each. (see docstring for description of
    # objective function. We can ignore the constant, total _n_required, as
    # it won't affect the position of the solutions)
    objective_coefficients = []

    # Coefficients and constants of our constraints (see docstring for
    # description - coefficients come from the left side of the equation,
    # constants from the right).
    # Both will be a list of length == len(required_n_per_genotype). With
    # each item being:
    # - For constraint_coefficients, a list of the expected number of offspring
    #   of that genotype per mating for all breeding schemes
    # - For constraint_constants, the required number of individuals of that
    #   genotype
    constraint_coefficients: list[list[float]] = []
    constraint_constants = []

    for genotype in required_genotypes:
        constraint_coefficients.append([])
        # Multiply by -1 as linprog only accepts <= constraints, and ours
        # is >=
        constraint_constants.append(required_n_per_genotype[genotype] * -1)

    for breeding_scheme in breeding_schemes:
        expected_offspring = offspring_per_scheme[breeding_scheme]
        objective_coefficients.append(
            expected_offspring.total_n  # litter size
        )

        n_per_genotype = expected_offspring.n_per_genotype

        for i, genotype in enumerate(required_genotypes):
            if genotype in n_per_genotype:
                # again, multiply by -1 as linprog only accepts <= constraints
                constraint_coefficients[i].append(
                    n_per_genotype[genotype] * -1
                )
            else:
                constraint_coefficients[i].append(0)

    # By default linprog has bounds of 0-infinity for variables (i.e. the
    # n_matings) so no need to specify further here.

    optimised_result = linprog(
        objective_coefficients,
        A_ub=constraint_coefficients,
        b_ub=constraint_constants,
        integrality=1,  # returned values for n of matings must be integers
        method="highs",
    )

    if (optimised_result.status != 0) or (not optimised_result.success):
        raise ValueError(
            f"Number of matings couldn't be optimised: "
            f"{optimised_result.message}"
        )

    n_matings_per_scheme = dict(zip(breeding_schemes, optimised_result.x))
    return n_matings_per_scheme


def _create_surplus_summary(
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


def _estimate_n_offspring_per_mating(
    line_stats: LineStatistics,
    min_n_offspring: int,
    min_n_matings: int,
    default_litter_size: int,
) -> dict[BreedingScheme, ExpectedOffspring]:
    """For each breeding scheme, estimate the number of offspring produced
    per mating.

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
            expected_offspring.n_per_genotype[genotype] = (
                proportion * litter_size
            )

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
