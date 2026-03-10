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
class SurplusSummary:
    total_n: int = 0
    total_n_surplus: int = 0

    total_n_per_genotype: dict[tuple[Genotype, ...], int] = field(
        default_factory=dict
    )
    total_surplus_per_genotype: dict[tuple[Genotype, ...], int] = field(
        default_factory=dict
    )
    percent_surplus_per_genotype: dict[tuple[Genotype, ...], float] = field(
        default_factory=dict
    )


def calculate_optimal_scheme(
    required_n_per_genotype: dict[tuple[Genotype, ...], int],
    line_stats: LineStatistics,
    min_n_matings: int = 3,
) -> tuple[dict[BreedingScheme, int], SurplusSummary]:
    """Calculate optimal combination of breeding schemes to produce
    the required_n_per_genotype.

    Parameters
    ----------
    required_n_per_genotype : dict[tuple[Genotype, ...], int]
        Required number of individuals per genotype
    line_stats : LineStatistics
        Statistics from historical data for the line
    min_n_matings : int, optional
        Minimum number of recorded matings required in historical data to use
        the measured litter size. If there aren't enough matings for a specific
        breeding scheme, the average of the whole line will be used instead.
        If the whole line also doesn't have enough matings, then the average
        across all lines will be used.

    Returns
    -------
    tuple[dict[BreedingScheme, int], SurplusSummary]
        Returns 2 items:
            - optimal number of matings per breeding scheme
            - a summary of the surplus numbers for this optimal scheme
              combination
    """

    # TODO - do we want to enable the scenario where there is no historical
    # data for a line?
    combined_ratios = _get_combined_ratios(line_stats)
    litter_size_per_scheme = _get_estimated_litter_sizes(
        combined_ratios.keys(), line_stats, min_n_matings
    )

    # Convert expected ratios into an expected number of offspring of
    # each genotype per litter
    n_per_scheme_genotype = combined_ratios.copy()
    for (
        breeding_scheme,
        genotype_to_ratio_dict,
    ) in n_per_scheme_genotype.items():
        for genotype in genotype_to_ratio_dict:
            genotype_to_ratio_dict[genotype] = (
                genotype_to_ratio_dict[genotype]
                * litter_size_per_scheme[breeding_scheme]
            )

    n_matings_per_scheme = _optimise_n_matings(
        required_n_per_genotype, litter_size_per_scheme, n_per_scheme_genotype
    )

    surplus_summary = _create_surplus_summary(
        required_n_per_genotype, n_matings_per_scheme, n_per_scheme_genotype
    )

    return n_matings_per_scheme, surplus_summary


def _optimise_n_matings(
    required_n_per_genotype: dict[tuple[Genotype, ...], int],
    litter_size_per_scheme: dict[BreedingScheme, float],
    n_per_scheme_genotype: dict[
        BreedingScheme, dict[tuple[Genotype, ...], float]
    ],
) -> dict[BreedingScheme, int]:
    """Calculate the optimal number of matings of each breeding scheme to
    meet the required_genotypes.

    This is framed as a linear programming (optimisation) problem with three
    parts: the decision variables (the number of matings for each scheme),
    the objective function (the surplus to minimise) and the constraints
    (that the final result must meet the given required_genotypes).

    Parameters
    ----------
    required_n_per_genotype : dict[tuple[Genotype, ...], int]
        The required number of individuals per genotype
    litter_size_per_scheme : dict[BreedingScheme, float]
        The estimated litter size per breeding scheme
    n_per_scheme_genotype :
        dict[BreedingScheme, dict[tuple[Genotype, ...], float]]
        The estimated number of individuals per litter of the given
        genotype + breeding scheme
    """

    # Extract names of breeding schemes / required genotypes as a list, so we
    # can make sure we always iterate through them in the same order
    breeding_schemes = list(litter_size_per_scheme.keys())
    required_genotypes = list(required_n_per_genotype.keys())

    # Coefficients of the function defining the amount of surplus (we want to
    # minimise).
    # This is the sum of
    # (estimated litter size for scheme * n matings for scheme)
    # for all schemes minus the total number of individuals required.
    # Therefore, coefficients are the estimated litter sizes (we ignore the
    # constant = total number of individuals for now, as it won't affect the
    # position of the solutions)
    objective_coefficients = []

    # linprog constraints are given in form Ax + By + ... <= C
    # The coefficients on the left side are added to 'constraint_coefficients',
    # and the constant on the right to 'constraint_constants'.

    # Our constraints are that the final number of
    # individuals for each genotype must be >= the required number.
    # e.g. for genotype hom:
    # sum of (n matings for scheme * expected number of hom per litter)
    # across all schemes >= required number of hom
    # So coefficients will be expected_n_per_genotype for that genotype, and
    # the constant the required number for that genotype.
    constraint_coefficients: list[list[float]] = []
    constraint_constants = []

    for genotype in required_genotypes:
        constraint_coefficients.append([])
        # Multiply by -1 as linprog only accepts <= constraints, and ours
        # is >=
        constraint_constants.append(required_n_per_genotype[genotype] * -1)

    for breeding_scheme in breeding_schemes:
        objective_coefficients.append(litter_size_per_scheme[breeding_scheme])
        n_per_genotype = n_per_scheme_genotype[breeding_scheme]

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

    if (not optimised_result.status == 0) or (not optimised_result.success):
        raise ValueError(
            f"Number of matings couldn't be optimised: "
            f"{optimised_result.message}"
        )

    n_matings_per_scheme = dict(zip(breeding_schemes, optimised_result.x))

    return n_matings_per_scheme


def _create_surplus_summary(
    required_n_per_genotype: dict[tuple[Genotype, ...], int],
    n_matings_per_scheme: dict[BreedingScheme, int],
    n_per_scheme_genotype,
):
    surplus_summary = SurplusSummary()

    for breeding_scheme, n_matings in n_matings_per_scheme.items():
        n_per_genotype = n_per_scheme_genotype[breeding_scheme]

        total_for_scheme = 0
        for genotype, n in n_per_genotype.items():
            # round value up to the nearest int - as we can't have
            # partial animals
            total_n = math.ceil(n * n_matings)
            total_for_scheme += total_n

            if genotype in surplus_summary.total_n_per_genotype:
                surplus_summary.total_n_per_genotype[genotype] += total_n
            else:
                surplus_summary.total_n_per_genotype[genotype] = total_n

        surplus_summary.total_n += total_for_scheme

    total_required = sum(
        [required_n for required_n in required_n_per_genotype.values()]
    )
    surplus_summary.total_n_surplus = surplus_summary.total_n - total_required

    for genotype, total_n in surplus_summary.total_n_per_genotype.items():
        total_surplus = total_n - required_n_per_genotype[genotype]
        surplus_summary.total_surplus_per_genotype[genotype] = total_surplus
        surplus_summary.percent_surplus_per_genotype[genotype] = (
            total_surplus / surplus_summary.total_n_per_genotype[genotype]
        ) * 100

    return surplus_summary


def _get_estimated_litter_sizes(breeding_schemes, line_stats, min_n_matings):
    litter_sizes = {}

    for breeding_scheme in breeding_schemes:
        scheme_stats = line_stats.stats_per_breeding_scheme[breeding_scheme]

        if scheme_stats.n_successful_matings >= min_n_matings:
            litter_sizes[breeding_scheme] = scheme_stats.average_litter_size
        elif line_stats.total_n_successful_matings >= min_n_matings:
            litter_sizes[breeding_scheme] = line_stats.average_litter_size
        else:
            # TODO - Fallback to whole institution stats
            litter_sizes[breeding_scheme] = 0

    return litter_sizes


def _get_combined_ratios(
    line_stats: LineStatistics,
) -> dict[BreedingScheme, dict[tuple[Genotype, ...], float]]:
    """Create dictionary combining the observed genotyping ratio (proportion of
    offspring of each genotype from historical data) with the theoretical
    mendelian ratio (when no historical data is present).

    Parameters
    ----------
    line_stats : LineStatistics
        Summary line statistics from historical data

    Returns
    -------
    dict[BreedingScheme, float]
        Expected proportion of offspring ...
    """

    # The breeding schemes listed in line_stats.stats_per_breeding_scheme are
    # sparse. i.e. only breeding schemes that appeared in the historical data
    # are included (rather than all that are possible)
    # Here we need ALL possible schemes:
    n_mutations = len(
        list(line_stats.total_n_offspring_per_genotype.keys())[0]
    )
    breeding_schemes = generate_breeding_schemes(n_mutations)

    combined_ratios = {}

    for breeding_scheme in breeding_schemes:
        # If present, use observed proportion from historical data
        if breeding_scheme in line_stats.stats_per_breeding_scheme:
            combined_ratios[breeding_scheme] = (
                line_stats.stats_per_breeding_scheme[
                    breeding_scheme
                ].proportion_offspring_per_genotype
            )
        # Otherwise, use the theoretical ratio
        else:
            combined_ratios[breeding_scheme] = (
                breeding_scheme.mendelian_ratio()
            )

    return combined_ratios
