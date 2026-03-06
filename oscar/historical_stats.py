from dataclasses import dataclass, field

import pandas as pd

from oscar.breeding_scheme import (
    BreedingScheme,
    Genotype,
)


@dataclass
class BreedingSchemeStatistics:
    n_breeding_pairs: int = 0
    n_successful_matings: int = 0
    average_litter_size: float = 0
    average_n_litters_per_pair: float = 0
    total_n_offspring: int = 0
    n_offspring_per_genotype: dict[tuple[Genotype, ...], int] = field(
        default_factory=dict
    )
    proportion_offspring_per_genotype: dict[tuple[Genotype, ...], float] = (
        field(default_factory=dict)
    )


@dataclass
class LineStatistics:
    total_n_offspring: int = 0
    total_n_offspring_per_genotype: dict[tuple[Genotype, ...], int] = field(
        default_factory=dict
    )

    stats_per_breeding_scheme: dict[
        BreedingScheme, BreedingSchemeStatistics
    ] = field(default_factory=dict)


def calculate_historical_stats_for_line(
    line_data: pd.DataFrame,
) -> LineStatistics:
    breeding_schemes = line_data.apply(_create_breeding_scheme, axis=1)
    data_with_schemes = line_data.copy()
    data_with_schemes["breeding_scheme"] = breeding_schemes

    line_stats = LineStatistics(total_n_offspring=len(line_data))

    for breeding_scheme in data_with_schemes["breeding_scheme"].unique():
        breeding_scheme_data = data_with_schemes.loc[
            data_with_schemes.breeding_scheme == breeding_scheme, :
        ]
        scheme_stats = _historical_stats_for_breeding_scheme(
            breeding_scheme_data
        )
        line_stats.stats_per_breeding_scheme[breeding_scheme] = scheme_stats

        # Update summary of number of offspring per genotype across entire line
        for (
            genotype,
            n_offspring,
        ) in scheme_stats.n_offspring_per_genotype.items():
            if genotype in line_stats.total_n_offspring_per_genotype:
                line_stats.total_n_offspring_per_genotype[genotype] += (
                    n_offspring
                )
            else:
                line_stats.total_n_offspring_per_genotype[genotype] = (
                    n_offspring
                )

    return line_stats


def _create_breeding_scheme(row: pd.Series) -> BreedingScheme:
    return BreedingScheme(row.genotype_father, row.genotype_mother)


def _historical_stats_for_breeding_scheme(
    scheme_data: pd.DataFrame,
) -> BreedingSchemeStatistics:
    stats = BreedingSchemeStatistics()

    # breeding pairs is unique combos of father ID x mother ID
    stats.n_breeding_pairs = scheme_data.groupby(
        ["ID_father", "ID_mother"]
    ).ngroups

    # Successful matings is unique combos of father ID x mother ID x date
    # (assuming only one per day)
    stats.n_successful_matings = scheme_data.groupby(
        ["ID_father", "ID_mother", "date_of_birth"]
    ).ngroups

    stats.total_n_offspring = len(scheme_data)
    stats.average_litter_size = (
        stats.total_n_offspring / stats.n_successful_matings
    )
    stats.average_n_litters_per_pair = (
        stats.n_successful_matings / stats.n_breeding_pairs
    )

    # convert string representation e.g. wt_hom_het to tuple representation
    # of genotype: (Genotype.WT, Genotype.HOM, Genotype.HET)
    scheme_data["genotype_offspring"] = scheme_data[
        "genotype_offspring"
    ].apply(Genotype.from_string)

    # Number and proportion of offspring per genotype
    stats.n_offspring_per_genotype = (
        scheme_data.groupby("genotype_offspring").size().to_dict()
    )

    for genotype, n_offspring in stats.n_offspring_per_genotype.items():
        proportion = n_offspring / stats.total_n_offspring
        stats.proportion_offspring_per_genotype[genotype] = proportion

    return stats
