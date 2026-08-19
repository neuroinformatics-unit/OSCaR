import logging
from dataclasses import dataclass, field

import pandas as pd

from oscar_colony.breeding_scheme import (
    BreedingScheme,
    Genotype,
)
from oscar_colony.optimise.estimate_offspring import ExpectedOffspring

logger = logging.getLogger(__name__)


@dataclass
class GenotypeSurplus:
    """Summary of surplus for a single genotype.

    Male / female surplus is None where a split of (n_males, n_females)
    wasn't requested for this genotype.
    """

    total_n: float = 0
    total_n_surplus: float = 0
    total_m: float = 0
    total_m_surplus: float | None = None
    m_percent_surplus: float | None = None
    total_f: float = 0
    total_f_surplus: float | None = None
    f_percent_surplus: float | None = None
    percent_surplus: float = 0


@dataclass
class SurplusSummary:
    """Summary of surplus across all genotypes.

    Male / female numbers are None where no genotype was requested as a
    split of (n_males, n_females). Where some were, they only cover those
    genotypes.
    """

    total_n: float = 0
    total_n_surplus: float = 0
    total_m: float | None = None
    total_m_surplus: float | None = None
    total_f: float | None = None
    total_f_surplus: float | None = None

    surplus_per_genotype: dict[tuple[Genotype, ...], GenotypeSurplus] = field(
        default_factory=dict
    )

    def create_genotype_df(
        self, decimal_places: int | None = None
    ) -> pd.DataFrame:
        """Create a pandas dataframe from surplus_per_genotype.

        Columns are: 'Genotype', 'Required N', 'Total N', 'Total N Surplus',
        'Total M', 'Total M Surplus', 'Percent M Surplus', 'Total F',
        'Total F Surplus', 'Percent F Surplus' and 'Percent Surplus'

        Parameters
        ----------
        decimal_places : int | None, optional
            Number of decimal places to round float values to

        Returns
        -------
        pd.DataFrame
            DataFrame summarising totals and surplus per genotype
        """
        rows = []
        for genotype, surplus in self.surplus_per_genotype.items():
            required_n = surplus.total_n - surplus.total_n_surplus

            # Male / female numbers are only shown where a split of
            # (n_males, n_females) was requested for this genotype
            sex_requested = surplus.total_m_surplus is not None

            rows.append(
                (
                    Genotype.to_string(genotype),
                    required_n,
                    surplus.total_n,
                    surplus.total_n_surplus,
                    surplus.total_m if sex_requested else None,
                    surplus.total_m_surplus,
                    surplus.m_percent_surplus,
                    surplus.total_f if sex_requested else None,
                    surplus.total_f_surplus,
                    surplus.f_percent_surplus,
                    surplus.percent_surplus,
                )
            )

        genotype_df = pd.DataFrame(
            rows,
            columns=[
                "Genotype",
                "Required N",
                "Total N",
                "Total N Surplus",
                "Total M",
                "Total M Surplus",
                "Percent M Surplus",
                "Total F",
                "Total F Surplus",
                "Percent F Surplus",
                "Percent Surplus",
            ],
        )

        if decimal_places is not None:
            genotype_df = genotype_df.round(decimals=decimal_places)

        return genotype_df


def create_surplus_summary(
    required_n_per_genotype: dict[tuple[Genotype, ...], int | tuple[int, int]],
    n_matings_per_scheme: dict[BreedingScheme, int],
    offspring_per_scheme: dict[BreedingScheme, ExpectedOffspring],
) -> SurplusSummary:
    """Create a summary of the total and surplus numbers for the
    given combination of breeding schemes.

    Parameters
    ----------
    required_n_per_genotype : dict
        Required number of individuals of each genotype. Each value is
        either an int (where male / female doesn't matter), or a tuple of
        (n_males, n_females).
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

    # Total required numbers - males / females are only counted where they
    # were specifically asked for
    total_required = 0
    n_males_required = 0
    n_females_required = 0
    sex_requested = False
    for required_n in required_n_per_genotype.values():
        if isinstance(required_n, tuple):
            sex_requested = True
            n_males_required += required_n[0]
            n_females_required += required_n[1]
            total_required += sum(required_n)
        else:
            total_required += required_n

    # Get number of expected offspring overall / per genotype
    for breeding_scheme, n_matings in n_matings_per_scheme.items():
        expected_offspring = offspring_per_scheme[breeding_scheme]
        total_n = expected_offspring.total_n * n_matings
        proportion_male = expected_offspring.proportion_male

        surplus_summary.total_n += total_n

        n_per_genotype = expected_offspring.n_per_genotype
        for genotype, n_per_mating in n_per_genotype.items():
            if genotype not in surplus_per_genotype:
                surplus_per_genotype[genotype] = GenotypeSurplus()

            n_of_genotype = n_per_mating * n_matings
            genotype_surplus = surplus_per_genotype[genotype]
            genotype_surplus.total_n += n_of_genotype
            genotype_surplus.total_m += n_of_genotype * proportion_male
            genotype_surplus.total_f += n_of_genotype * (1 - proportion_male)

    # Calculate total surplus
    surplus_summary.total_n_surplus = surplus_summary.total_n - total_required

    # Calculate surplus per genotype. Males / females are only reported
    # where a split was specifically asked for.
    for genotype, surplus in surplus_per_genotype.items():
        required_n = required_n_per_genotype.get(genotype, 0)
        if isinstance(required_n, tuple):
            required_m, required_f = required_n
            required_n = required_m + required_f
            surplus.total_m_surplus = surplus.total_m - required_m
            surplus.total_f_surplus = surplus.total_f - required_f
            if surplus.total_m > 0:
                surplus.m_percent_surplus = (
                    surplus.total_m_surplus / surplus.total_m
                ) * 100
            if surplus.total_f > 0:
                surplus.f_percent_surplus = (
                    surplus.total_f_surplus / surplus.total_f
                ) * 100

        surplus.total_n_surplus = surplus.total_n - required_n
        surplus.percent_surplus = (
            surplus.total_n_surplus / surplus.total_n
        ) * 100

    # Summary male / female numbers only cover the genotypes that were
    # requested as a split, so they can be compared to the required numbers
    if sex_requested:
        total_m = 0.0
        total_f = 0.0
        for surplus in surplus_per_genotype.values():
            if surplus.total_m_surplus is not None:
                total_m += surplus.total_m
                total_f += surplus.total_f

        surplus_summary.total_m = total_m
        surplus_summary.total_f = total_f
        surplus_summary.total_m_surplus = total_m - n_males_required
        surplus_summary.total_f_surplus = total_f - n_females_required

    # Warn where the plan doesn't reach the required numbers
    for genotype, surplus in surplus_per_genotype.items():
        for label, n_surplus in (
            ("", surplus.total_n_surplus),
            ("male ", surplus.total_m_surplus),
            ("female ", surplus.total_f_surplus),
        ):
            if (n_surplus is not None) and (n_surplus < 0):
                logger.warning(
                    f"{Genotype.to_string(genotype)}: plan is "
                    f"{abs(n_surplus):.1f} {label}individuals short of the "
                    "required number"
                )

    return surplus_summary
