from scipy.optimize import linprog

from oscar.breeding_scheme import (
    BreedingScheme,
    Genotype,
)
from oscar.historical_stats import LineStatistics
from oscar.optimise.estimate_offspring import (
    ExpectedOffspring,
    estimate_n_offspring_per_mating,
)
from oscar.optimise.surplus_summary import (
    SurplusSummary,
    create_surplus_summary,
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

    offspring_per_scheme = estimate_n_offspring_per_mating(
        line_stats, min_n_offspring, min_n_matings, default_litter_size
    )

    n_matings_per_scheme = _optimise_n_matings(
        required_n_per_genotype, offspring_per_scheme
    )

    surplus_summary = create_surplus_summary(
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
