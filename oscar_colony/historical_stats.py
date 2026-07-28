from dataclasses import dataclass, field

import pandas as pd

from oscar_colony.breeding_scheme import (
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
    total_n_genotyped_offspring: int = 0
    n_offspring_per_genotype: dict[tuple[Genotype, ...], int] = field(
        default_factory=dict
    )
    proportion_offspring_per_genotype: dict[tuple[Genotype, ...], float] = (
        field(default_factory=dict)
    )


@dataclass
class LineStatistics:
    line_name: str
    mutations: list[str] = field(default_factory=list)
    n_mutations: int = 0
    total_n_offspring: int = 0
    total_n_genotyped_offspring: int = 0
    total_n_offspring_per_genotype: dict[tuple[Genotype, ...], int] = field(
        default_factory=dict
    )
    total_n_successful_matings: int = 0
    average_litter_size: float = 0

    stats_per_breeding_scheme: dict[
        BreedingScheme, BreedingSchemeStatistics
    ] = field(default_factory=dict)


def calculate_historical_stats_for_line(
    standardised_data: pd.DataFrame, line_name: str
) -> LineStatistics:
    """Calculate summary statistics for a specific line from standardised
    historical data.

    Parameters
    ----------
    standardised_data : pd.DataFrame
        Standardised historical data e.g. from standardise_pyrat_csv
    line_name : str
        Name of line

    Returns
    -------
    LineStatistics
        Summary statistics for the given line
    """

    line_data = standardised_data.loc[
        standardised_data.line_name == line_name, :
    ]
    if len(line_data) == 0:
        raise ValueError(f"No data for {line_name} found")

    breeding_schemes = line_data.apply(_create_breeding_scheme, axis=1)
    data_with_schemes = line_data.copy()
    data_with_schemes["breeding_scheme"] = breeding_schemes

    # Get mutation column names in numeric order e.g. mutation_1, mutation_2...
    mutation_cols = sorted(line_data.filter(regex=r"^mutation_\d+$").columns)
    mutations = line_data[mutation_cols].iloc[0].dropna().tolist()

    line_stats = LineStatistics(
        line_name=line_name,
        mutations=mutations,
        n_mutations=line_data.n_mutations.iloc[0],
        total_n_offspring=len(line_data),
        total_n_genotyped_offspring=sum(~line_data.genotype_offspring.isna()),
    )

    for breeding_scheme in data_with_schemes["breeding_scheme"].unique():
        breeding_scheme_data = data_with_schemes.loc[
            data_with_schemes.breeding_scheme == breeding_scheme, :
        ]
        scheme_stats = _historical_stats_for_breeding_scheme(
            breeding_scheme_data
        )

        line_stats.stats_per_breeding_scheme[breeding_scheme] = scheme_stats
        line_stats.total_n_successful_matings += (
            scheme_stats.n_successful_matings
        )

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

    # Use total_n_offspring for litter size (including un-genotyped)
    line_stats.average_litter_size = (
        line_stats.total_n_offspring / line_stats.total_n_successful_matings
    )

    return line_stats


def _create_breeding_scheme(row: pd.Series) -> BreedingScheme:
    return BreedingScheme(row.genotype_father, row.genotype_mother)


def _historical_stats_for_breeding_scheme(
    scheme_data: pd.DataFrame,
) -> BreedingSchemeStatistics:
    """Calculate summary statistics for an individual breeding scheme
    (within a specific line).

    Accounting for multiple parents by finding every combination of parent
    then removing the duplicates for total breeding pairs. For Successful
    matings unique mothers with unique dates are quantified.

    Parameters
    ----------
    scheme_data : pd.DataFrame
        Dataframe of data for a single breeding scheme and line

    Returns
    -------
    BreedingSchemeStatistics
        Summary statistics for the breeding scheme
    """
    stats = BreedingSchemeStatistics()

    father_cols = [
        c for c in scheme_data.columns if c.startswith("ID_father_")
    ]
    mother_cols = [
        c for c in scheme_data.columns if c.startswith("ID_mother_")
    ]

    # breeding pairs: all unique (father_id, mother_id) combinations, taking
    # the cartesian product of fathers x mothers per row - de-duplicated
    father_combos = (
        scheme_data[father_cols]
        .melt(ignore_index=False, value_name="father_id")
        .dropna()
    )
    mother_combos = (
        scheme_data[mother_cols]
        .melt(ignore_index=False, value_name="mother_id")
        .dropna()
    )
    breeding_pairs = father_combos[["father_id"]].merge(
        mother_combos[["mother_id"]], left_index=True, right_index=True
    )
    stats.n_breeding_pairs = breeding_pairs.drop_duplicates().shape[0]

    # unique (mother_id, date_of_birth) pairs — mothers
    mothers_with_date = (
        scheme_data[mother_cols + ["date_of_birth"]]
        .melt(
            id_vars=["date_of_birth"],
            value_vars=mother_cols,
            value_name="mother_id",
        )
        .dropna(subset=["mother_id"])
    )
    stats.n_successful_matings = (
        mothers_with_date[["mother_id", "date_of_birth"]]
        .drop_duplicates()
        .shape[0]
    )

    stats.total_n_offspring = len(scheme_data)

    genotyped_rows = scheme_data.loc[~scheme_data.genotype_offspring.isna()]
    stats.total_n_genotyped_offspring = len(genotyped_rows)

    # litter size calculations use total_no_offspring (including
    # un-genotyped individuals)
    stats.average_litter_size = (
        stats.total_n_offspring / stats.n_successful_matings
    )
    stats.average_n_litters_per_pair = (
        stats.n_successful_matings / stats.n_breeding_pairs
    )

    # convert string representation e.g. wt_hom_het to tuple representation
    # of genotype: (Genotype.WT, Genotype.HOM, Genotype.HET)
    genotyped_rows["genotype_offspring"] = genotyped_rows[
        "genotype_offspring"
    ].apply(Genotype.from_string)

    # Number and proportion of offspring per genotype
    stats.n_offspring_per_genotype = (
        genotyped_rows.groupby("genotype_offspring").size().to_dict()
    )

    for genotype, n_offspring in stats.n_offspring_per_genotype.items():
        # Proportions use total_n_genotyped_offspring (excluding un-genotyped)
        proportion = n_offspring / stats.total_n_genotyped_offspring
        stats.proportion_offspring_per_genotype[genotype] = proportion

    return stats
