import pandas as pd

from oscar_colony.breeding_scheme import Genotype
from oscar_colony.optimise.surplus_summary import SexSplit


def count_breeding_stock(
    standardised_data: pd.DataFrame, line_name: str
) -> dict[tuple[Genotype, ...], SexSplit]:
    """Count the males and females of each genotype in a line.

    Parameters
    ----------
    standardised_data : pd.DataFrame
        Standardised data of live animals e.g. from standardise_pyrat_csv
    line_name : str
        Name of line

    Returns
    -------
    dict[tuple[Genotype, ...], SexSplit]
        Number of (n_males, n_females) per genotype, for use as
        available_breeding_stock. Un-genotyped animals are left out.
    """

    line_data = standardised_data.loc[
        standardised_data.line_name == line_name, :
    ]

    available_breeding_stock = {}
    for genotype_string, genotype_data in line_data.groupby(
        "genotype_offspring"
    ):
        genotype = Genotype.from_string(genotype_string)
        sex_counts = genotype_data["offspring_sex"].value_counts()
        available_breeding_stock[genotype] = SexSplit(
            n_males=int(sex_counts.get("m", 0)),
            n_females=int(sex_counts.get("f", 0)),
        )

    return available_breeding_stock
