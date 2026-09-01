from dataclasses import dataclass, field
from typing import NamedTuple

import pandas as pd

from oscar_colony.breeding_scheme import (
    BreedingScheme,
    Genotype,
)
from oscar_colony.optimise.estimate_offspring import ExpectedOffspring


class SexSplit(NamedTuple):
    n_males: int
    n_females: int


@dataclass
class SexSurplus:
    """Summary of surplus split by male / female."""

    total_males: float = 0
    total_male_surplus: float = 0
    male_percent_surplus: float = 0
    total_females: float = 0
    total_female_surplus: float = 0
    female_percent_surplus: float = 0

    def set_surplus(self, required_males: int, required_females: int) -> None:
        """Set the surplus / percent surplus, by comparing the required
        numbers to the total males / females.
        """
        self.total_male_surplus = self.total_males - required_males
        self.total_female_surplus = self.total_females - required_females

        if self.total_males > 0:
            self.male_percent_surplus = (
                self.total_male_surplus / self.total_males
            ) * 100
        if self.total_females > 0:
            self.female_percent_surplus = (
                self.total_female_surplus / self.total_females
            ) * 100


@dataclass
class GenotypeSurplus:
    """Summary of surplus for a single genotype.

    sex_surplus is None where a SexSplit of (n_males, n_females) wasn't
    requested for this genotype.
    """

    total_n: float = 0
    total_n_surplus: float = 0
    percent_surplus: float = 0
    sex_surplus: SexSurplus | None = None

    def set_surplus(self, required_n: int) -> None:
        """Set the surplus / percent surplus, by comparing the required
        number to the total.
        """
        self.total_n_surplus = self.total_n - required_n
        if self.total_n > 0:
            self.percent_surplus = (self.total_n_surplus / self.total_n) * 100


@dataclass
class SurplusSummary:
    """Summary of surplus across all genotypes.

    sex_surplus is None when no genotypes requested a SexSplit of
    (n_males, n_females). Covers cases in which there is a mix of sex-split
    and non-sex-split genotypes.
    """

    total_n: float = 0
    total_n_surplus: float = 0
    sex_surplus: SexSurplus | None = None

    surplus_per_genotype: dict[tuple[Genotype, ...], GenotypeSurplus] = field(
        default_factory=dict
    )

    def create_genotype_df(
        self, decimal_places: int | None = None
    ) -> pd.DataFrame:
        """Create a pandas dataframe from surplus_per_genotype.

        Columns are: 'Genotype', 'Required N', 'Total N',
        'Total N Surplus', 'Percent Surplus', 'Required Males',
        'Total Males', 'Total Male Surplus', 'Percent Male Surplus',
        'Required Females', 'Total Females', 'Total Female Surplus' and
        'Percent Female Surplus'

        The male / female columns are omitted where no genotype was
        requested as a SexSplit of (n_males, n_females).

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
            row = {
                "Genotype": Genotype.to_string(genotype),
                "Required N": surplus.total_n - surplus.total_n_surplus,
                "Total N": surplus.total_n,
                "Total N Surplus": surplus.total_n_surplus,
                "Percent Surplus": surplus.percent_surplus,
            }

            sex_surplus = surplus.sex_surplus
            if sex_surplus is not None:
                row.update(
                    {
                        "Required Males": (
                            sex_surplus.total_males
                            - sex_surplus.total_male_surplus
                        ),
                        "Total Males": sex_surplus.total_males,
                        "Total Male Surplus": sex_surplus.total_male_surplus,
                        "Percent Male Surplus": (
                            sex_surplus.male_percent_surplus
                        ),
                        "Required Females": (
                            sex_surplus.total_females
                            - sex_surplus.total_female_surplus
                        ),
                        "Total Females": sex_surplus.total_females,
                        "Total Female Surplus": (
                            sex_surplus.total_female_surplus
                        ),
                        "Percent Female Surplus": (
                            sex_surplus.female_percent_surplus
                        ),
                    }
                )

            rows.append(row)

        genotype_df = pd.DataFrame(rows)

        if decimal_places is not None:
            genotype_df = genotype_df.round(decimals=decimal_places)

        return genotype_df


def create_surplus_summary(
    required_n_per_genotype: dict[tuple[Genotype, ...], int | SexSplit],
    n_matings_per_scheme: dict[BreedingScheme, int],
    offspring_per_scheme: dict[BreedingScheme, ExpectedOffspring],
) -> SurplusSummary:
    """Create a summary of the total and surplus numbers for the
    given combination of breeding schemes.

    Parameters
    ----------
    required_n_per_genotype : dict[tuple[Genotype, ...], int | SexSplit]
        Required number of individuals per genotype. Each value is either an
        int (where male / female isn't requested), or a SexSplit of
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

    Raises
    ------
    ValueError
        If the plan doesn't reach the required number of individuals for any
        genotype - in total, or of males / females
    """
    surplus_summary = SurplusSummary()
    surplus_per_genotype = surplus_summary.surplus_per_genotype

    # Total required numbers - males / females are only counted when requested
    total_required = 0
    n_males_required = 0
    n_females_required = 0
    sex_requested = False
    for required_n in required_n_per_genotype.values():
        if isinstance(required_n, SexSplit):
            sex_requested = True
            n_males_required += required_n.n_males
            n_females_required += required_n.n_females
            total_required += required_n.n_males + required_n.n_females
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
                sex_split = isinstance(
                    required_n_per_genotype.get(genotype), SexSplit
                )
                surplus_per_genotype[genotype] = GenotypeSurplus(
                    sex_surplus=SexSurplus() if sex_split else None
                )

            n_of_genotype = n_per_mating * n_matings
            genotype_surplus = surplus_per_genotype[genotype]
            genotype_surplus.total_n += n_of_genotype

            sex_surplus = genotype_surplus.sex_surplus
            if sex_surplus is not None:
                sex_surplus.total_males += n_of_genotype * proportion_male
                sex_surplus.total_females += n_of_genotype * (
                    1 - proportion_male
                )

    # Calculate total surplus
    surplus_summary.total_n_surplus = surplus_summary.total_n - total_required

    # Required genotypes that can't be produced by any of the chosen schemes
    # still need a summary, so their shortfall is visible
    for genotype, required_n in required_n_per_genotype.items():
        if genotype not in surplus_per_genotype:
            surplus_per_genotype[genotype] = GenotypeSurplus(
                sex_surplus=SexSurplus()
                if isinstance(required_n, SexSplit)
                else None
            )

    for genotype, surplus in surplus_per_genotype.items():
        required_n = required_n_per_genotype.get(genotype, 0)
        sex_surplus = surplus.sex_surplus
        if isinstance(required_n, SexSplit):
            if sex_surplus is not None:
                sex_surplus.set_surplus(
                    required_n.n_males, required_n.n_females
                )
            required_n = required_n.n_males + required_n.n_females

        surplus.set_surplus(required_n)

    # Summary male / female numbers only cover the genotypes that were
    # requested as a split, so they can be compared to the required numbers
    if sex_requested:
        summary_sex_surplus = SexSurplus()
        for surplus in surplus_per_genotype.values():
            if surplus.sex_surplus is not None:
                summary_sex_surplus.total_males += (
                    surplus.sex_surplus.total_males
                )
                summary_sex_surplus.total_females += (
                    surplus.sex_surplus.total_females
                )
        summary_sex_surplus.set_surplus(n_males_required, n_females_required)
        surplus_summary.sex_surplus = summary_sex_surplus

    # Error where the plan doesn't reach the required numbers
    error_messages = []
    for genotype in required_n_per_genotype:
        surplus = surplus_per_genotype[genotype]
        surplus_per_label = [("", surplus.total_n_surplus)]
        if surplus.sex_surplus is not None:
            surplus_per_label += [
                ("male ", surplus.sex_surplus.total_male_surplus),
                ("female ", surplus.sex_surplus.total_female_surplus),
            ]

        for label, n_surplus in surplus_per_label:
            if n_surplus < 0:
                error_messages.append(
                    f"{Genotype.to_string(genotype)}: plan is "
                    f"{abs(n_surplus):.1f} {label}individuals short of the "
                    "required number"
                )

    if error_messages:
        raise ValueError("\n".join(error_messages))

    return surplus_summary
