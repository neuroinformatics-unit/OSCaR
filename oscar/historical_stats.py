from dataclasses import dataclass, field

import numpy as np
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
    average_litter_per_pair: float = 0
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
    parent_genotypes = line_data[["genotype_father", "genotype_mother"]]

    # Get unique rows - regardless of order e.g. wt_het x het_het; should be
    # treated as a duplicate of het_het x wt_het
    sorted_genotypes = pd.DataFrame(
        np.sort(parent_genotypes, 1), index=parent_genotypes.index
    )
    unique_genotypes = parent_genotypes[~sorted_genotypes.duplicated()]

    # Get list of all unique breeding schemes present in this data
    breeding_schemes = list(
        unique_genotypes.apply(_create_breeding_scheme, axis=1)
    )

    line_stats = LineStatistics()
    for breeding_scheme in breeding_schemes:
        breeding_scheme_stats = _historical_stats_for_breeding_scheme(
            breeding_scheme, line_data
        )
        line_stats.stats_per_breeding_scheme[breeding_scheme] = (
            breeding_scheme_stats
        )

    return line_stats


def _create_breeding_scheme(row: pd.Series) -> BreedingScheme:
    return BreedingScheme(row.genotype_father, row.genotype_mother)


def _historical_stats_for_breeding_scheme(
    scheme: BreedingScheme, line_data: pd.DataFrame
) -> BreedingSchemeStatistics:
    stats = BreedingSchemeStatistics()

    # filter all rows for that breeding scheme - allow it to match in
    # any order
    scheme_f_m = (
        line_data["genotype_father"] + "x" + line_data["genotype_mother"]
    )
    scheme_m_f = (
        line_data["genotype_mother"] + "x" + line_data["genotype_father"]
    )

    scheme_data = line_data[
        (scheme_f_m == str(scheme)) | (scheme_m_f == str(scheme))
    ]
    if len(scheme_data) == 0:
        return stats

    # breeding pairs is unique combos of father ID x mother ID
    stats.n_breeding_pairs = scheme_data.groupby(
        ["ID_father", "ID_mother"]
    ).ngroups

    return stats
