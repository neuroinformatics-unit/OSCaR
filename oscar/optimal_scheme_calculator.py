from oscar.breeding_scheme import (
    BreedingScheme,
    Genotype,
    generate_breeding_schemes,
)
from oscar.historical_stats import LineStatistics


def calculate_optimal_scheme(
    required_genotypes: dict[tuple[Genotype, ...], int],
    line_stats: LineStatistics,
    min_n_matings: int = 3,
):
    # TODO - do we want to enable the scenario where there is no historical
    # data for a line?
    combined_ratios = _get_combined_ratios(line_stats)

    # Convert expected ratios into an expected number of offspring of
    # each genotype per litter
    estimated_n_per_litter = combined_ratios.copy()
    for (
        breeding_scheme,
        genotype_to_ratio_dict,
    ) in estimated_n_per_litter.items():
        scheme_stats = line_stats.stats_per_breeding_scheme[breeding_scheme]
        if scheme_stats.n_successful_matings >= min_n_matings:
            expected_litter_size = scheme_stats.average_litter_size
        elif line_stats.total_n_successful_matings >= min_n_matings:
            expected_litter_size = line_stats.average_litter_size
        else:
            # TODO - Fallback to whole institution stats
            expected_litter_size = 0

        for genotype in genotype_to_ratio_dict:
            genotype_to_ratio_dict[genotype] = (
                genotype_to_ratio_dict[genotype] * expected_litter_size
            )


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
