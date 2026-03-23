from dataclasses import dataclass, field
from enum import Enum

import pandas as pd

from oscar.breeding_scheme import (
    BreedingScheme,
    Genotype,
)


# ---------------------------------------------------------------------------
# Surplus type classification
# ---------------------------------------------------------------------------

class SurplusType(Enum):
    """Categories of surplus animals, as defined in issue #3."""

    EXPERIMENTAL = "Experimental"
    UNAVOIDABLE = "Unavoidable"
    AVOIDABLE = "Avoidable"
    OTHER = "Other"

    @classmethod
    def from_sacrifice_reason(cls, reason: str | None) -> "SurplusType":
        """Map a raw sacrifice_reason string to a SurplusType.

        Parameters
        ----------
        reason : str | None
            Raw value from the 'sacrifice_reason' column.

        Returns
        -------
        SurplusType
            Matched surplus category, or SurplusType.OTHER if unrecognised.
        """
        if reason is None or (isinstance(reason, float)):
            # NaN / missing → Other
            return cls.OTHER

        reason_lower = reason.strip().lower()

        # Map known pyRAT sacrifice reasons to surplus categories.
        # Adjust these mappings to match real pyRAT vocabulary as needed.
        _EXPERIMENTAL = {
            "experimental", "experiment", "used for experiment",
            "killed for experiment", "terminal", "scientific procedure",
        }
        _UNAVOIDABLE = {
            "found dead", "died", "natural death", "spontaneous death",
            "illness", "sick", "humane endpoint", "welfare",
            "unexpected death",
        }
        _AVOIDABLE = {
            "surplus", "cull", "culled", "excess", "overproduction",
            "wrong genotype", "not required",
        }

        if reason_lower in _EXPERIMENTAL:
            return cls.EXPERIMENTAL
        if reason_lower in _UNAVOIDABLE:
            return cls.UNAVOIDABLE
        if reason_lower in _AVOIDABLE:
            return cls.AVOIDABLE
        return cls.OTHER


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class SurplusStatistics:
    """Counts and proportions for each surplus type."""

    n_per_surplus_type: dict[SurplusType, int] = field(default_factory=dict)
    proportion_per_surplus_type: dict[SurplusType, float] = field(
        default_factory=dict
    )
    total_n: int = 0

    def _compute_proportions(self) -> None:
        """Derive proportions from counts. Called after all counts are set."""
        for surplus_type, n in self.n_per_surplus_type.items():
            self.proportion_per_surplus_type[surplus_type] = (
                n / self.total_n if self.total_n > 0 else 0.0
            )


@dataclass
class LineUsageStatistics:
    """Animal usage statistics for a single line (issue #3).

    Attributes
    ----------
    surplus_per_genotype : dict
        SurplusStatistics broken down by offspring genotype.
    overall_surplus : SurplusStatistics
        Aggregate SurplusStatistics across all offspring genotypes.
    """

    surplus_per_genotype: dict[tuple[Genotype, ...], SurplusStatistics] = field(
        default_factory=dict
    )
    overall_surplus: SurplusStatistics = field(
        default_factory=SurplusStatistics
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
    usage_statistics: LineUsageStatistics = field(
        default_factory=LineUsageStatistics
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

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
        Summary statistics for the given line, including usage statistics
        (surplus breakdown) as required by issue #3.
    """

    line_data = standardised_data.loc[
        standardised_data.line_name == line_name, :
    ]
    if len(line_data) == 0:
        raise ValueError(f"No data for {line_name} found")

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

    # Issue #3: calculate animal usage (surplus) statistics
    line_stats.usage_statistics = _calculate_usage_statistics(data_with_schemes)

    return line_stats


def _calculate_usage_statistics(line_data: pd.DataFrame) -> LineUsageStatistics:
    """Calculate surplus/usage statistics for a line (issue #3).

    Splits sacrifice_reason into SurplusType categories, then computes
    counts and percentages both per genotype and overall.

    Parameters
    ----------
    line_data : pd.DataFrame
        Standardised data for a single line, must contain
        'genotype_offspring' and 'sacrifice_reason' columns.

    Returns
    -------
    LineUsageStatistics
    """
    usage_stats = LineUsageStatistics()

    # Map each row's sacrifice_reason to a SurplusType
    surplus_series: pd.Series = line_data["sacrifice_reason"].apply(
        SurplusType.from_sacrifice_reason
    )

    # Parse genotype strings once → reuse for grouping
    genotype_series: pd.Series = line_data["genotype_offspring"].apply(
        Genotype.from_string
    )

    # ---- per genotype -------------------------------------------------------
    for genotype in genotype_series.unique():
        mask = genotype_series == genotype
        surplus_for_genotype = surplus_series[mask]
        usage_stats.surplus_per_genotype[genotype] = _surplus_stats_from_series(
            surplus_for_genotype
        )

    # ---- overall (all offspring) --------------------------------------------
    usage_stats.overall_surplus = _surplus_stats_from_series(surplus_series)

    return usage_stats


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _surplus_stats_from_series(surplus_series: pd.Series) -> SurplusStatistics:
    """Build a SurplusStatistics from a Series of SurplusType values.

    Parameters
    ----------
    surplus_series : pd.Series
        Each element is a SurplusType.

    Returns
    -------
    SurplusStatistics
        Counts and proportions for every SurplusType.
    """
    stats = SurplusStatistics()
    stats.total_n = len(surplus_series)

    # Initialise all four types to 0 so every key is always present
    stats.n_per_surplus_type = {st: 0 for st in SurplusType}

    counts = surplus_series.value_counts()
    for surplus_type, count in counts.items():
        stats.n_per_surplus_type[surplus_type] = int(count)

    stats._compute_proportions()
    return stats


def _create_breeding_scheme(row: pd.Series) -> BreedingScheme:
    return BreedingScheme(row.genotype_father, row.genotype_mother)


def _historical_stats_for_breeding_scheme(
    scheme_data: pd.DataFrame,
) -> BreedingSchemeStatistics:
    """Calculate summary statistics for an individual breeding scheme
    (within a specific line).

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
