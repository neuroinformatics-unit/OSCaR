from enum import Enum
from pathlib import Path

import pandas as pd

from oscar.breeding_scheme import Genotype


class Identifier(Enum):
    """
    Identifier for who a mutation / genotype refers to.
    """

    OFFSPRING = 0
    FATHER = 1
    MOTHER = 2


def standardise_pyrat_csv(input_csv: pd.DataFrame | Path) -> pd.DataFrame:
    """Standardise a csv file exported from pyRAT.

    Parameters
    ----------
    input_csv : pd.DataFrame | Path
        Csv file exported from pyRAT.

    Returns
    -------
    pd.DataFrame
        Standardised dataframe
    """
    if isinstance(input_csv, Path):
        input_csv = pd.read_csv(input_csv)

    mutation_cols, genotype_cols = _create_column_name_dicts(input_csv)
    all_mutation_cols_list = sum(mutation_cols.values(), [])
    all_genotype_cols_list = sum(genotype_cols.values(), [])

    required_cols = (
        [
            "ID",
            "Line / Strain (Name)",
            "DOB",
            "Father",
            "Mother",
            "Sacrifice reason",
        ]
        + all_mutation_cols_list
        + all_genotype_cols_list
    )

    # Get rid of any additional columns + rename to standard names
    standard_csv = input_csv[required_cols]
    standard_csv = standard_csv.rename(
        columns={
            "ID": "ID_offspring",
            "Line / Strain (Name)": "line_name",
            "DOB": "date_of_birth",
            "Father": "ID_father",
            "Mother": "ID_mother",
            "Sacrifice reason": "sacrifice_reason",
        }
    )

    standard_csv = _filter_allowed_genotypes(
        standard_csv, all_genotype_cols_list
    )
    standard_csv = _filter_ungenotyped(
        standard_csv, genotype_cols[Identifier.OFFSPRING]
    )

    standard_csv["n_mutations"] = (
        standard_csv.loc[:, genotype_cols[Identifier.OFFSPRING]]
        .notna()
        .sum(axis=1)
    )

    # make sure number of mutations is the same throughout each line -
    # use the max.
    # Sometimes particular individuals are ungenotyped (n_mutations = 0) or a
    # genotype value is omitted to mean wt.
    standard_csv["n_mutations"] = standard_csv.groupby("line_name")[
        "n_mutations"
    ].transform("max")

    standard_csv = standard_csv.groupby("line_name").apply(
        _make_combined_genotype_columns_for_line, mutation_cols, genotype_cols
    )

    standard_csv = standard_csv.reset_index().drop(
        ["level_1"] + all_genotype_cols_list + all_mutation_cols_list,
        axis=1,
    )

    # for readability, make sure ID_offspring is first
    id_offspring_col = standard_csv.pop("ID_offspring")
    standard_csv.insert(0, "ID_offspring", id_offspring_col)

    return standard_csv


def _create_column_name_dicts(
    input_csv: pd.DataFrame,
) -> tuple[dict[Identifier, list[str]], dict[Identifier, list[str]]]:
    """Create a dict of mutation / genotype column names for all identifiers
    (offspring, father, mother).

    Parameters
    ----------
    input_csv : pd.DataFrame
        Dataframe to extract column names from

    Returns
    -------
    tuple[dict[Identifier, list[str]], dict[Identifier, list[str]]]
        Returns (mutation column dict, genotype column dict). Both dictionaries
        have Identifier as the keys, and a list of column names as values.
    """

    prefixes = {
        Identifier.OFFSPRING: "",
        Identifier.FATHER: "Father: ",
        Identifier.MOTHER: "Mother: ",
    }

    mutation_dict = {}
    genotype_dict = {}

    for identifier, prefix in prefixes.items():
        # columns of form 'PREFIXMutation NUMBER'
        mutation_cols = list(
            input_csv.columns[
                input_csv.columns.str.contains(rf"^{prefix}Mutation \d$")
            ]
        )

        # columns of form 'PREFIXGrade NUMBER'
        genotype_cols = list(
            input_csv.columns[
                input_csv.columns.str.contains(rf"^{prefix}Grade \d$")
            ]
        )

        # Each mutation must have a corresponding genotype
        if len(mutation_cols) != len(genotype_cols):
            raise ValueError(
                f"Not all {identifier} mutation columns have a corresponding "
                f"genotype column."
            )

        # Make sure lists are in numeric order e.g. Grade 1, Grade 2, Grade 3
        mutation_dict[identifier] = sorted(mutation_cols)
        genotype_dict[identifier] = sorted(genotype_cols)

    return mutation_dict, genotype_dict


def _filter_allowed_genotypes(
    standard_csv: pd.DataFrame, genotype_cols: list[str]
) -> pd.DataFrame:
    """Only keep allowed genotypes of wt, het or hom.

    Remove others such as +/-, Tg, ko/ko

    Parameters
    ----------
    standard_csv : pd.DataFrame
        Dataframe to filter
    genotype_cols : list[str]
        Names of all genotype columns

    Returns
    -------
    pd.DataFrame
        Dataframe with forbidden genotypes removed
    """

    genotype_data = standard_csv.loc[:, genotype_cols]
    allowed_genotypes = (
        genotype_data.isin([genotype.name.lower() for genotype in Genotype])
        | genotype_data.isna()
    ).all(axis=1)
    filtered_data = standard_csv.loc[allowed_genotypes, :]

    return filtered_data


def _filter_ungenotyped(
    standard_csv: pd.DataFrame, offspring_genotype_cols: list[str]
) -> pd.DataFrame:
    """Remove rows for ungenotyped individuals.

    These offspring have na for all offspring genotype columns.

    Parameters
    ----------
    standard_csv : pd.DataFrame
        Dataframe to filter
    offspring_genotype_cols : list[str]
        Names of offspring genotype columns

    Returns
    -------
    pd.DataFrame
        Dataframe with ungenotyped individuals removed
    """
    ungenotyped = (
        standard_csv.loc[:, offspring_genotype_cols].isna().all(axis=1)
    )
    filtered_data = standard_csv.loc[~ungenotyped, :]

    return filtered_data


def _make_combined_genotype_columns_for_line(
    line_data: pd.DataFrame,
    mutation_cols: dict[Identifier, list[str]],
    genotype_cols: dict[Identifier, list[str]],
) -> pd.DataFrame:
    """For data on a single line, add columns for 'mutations',
    'genotype_offspring', 'genotype_father' and 'genotype_mother'.

    All genotype columns list genotypes in the same order as given in
    'mutations'. Any missing genotypes are assumed to be wildtype.

    Parameters
    ----------
    line_data : pd.DataFrame
        Data for a single line
    mutation_cols : dict[str, list[str]]
        Mutations columns grouped by identifier
    genotype_cols : dict[str, list[str]]
        Genotype columns grouped by identifier

    Returns
    -------
    pd.DataFrame
        Line data with added columns summarising mutations and genotypes
    """

    # get unique offspring mutations for this line
    unique_mutations = pd.unique(
        line_data[mutation_cols[Identifier.OFFSPRING]].values.ravel("K")
    )
    unique_mutations = list(pd.Series(unique_mutations).dropna())

    # Copy so we don't edit the original dataframe (this can cause issues
    # with apply)
    line_data_with_combined_cols = line_data.copy()
    line_data_with_combined_cols["mutations"] = "_".join(unique_mutations)

    for identifier in Identifier:
        _make_combined_genotype_column_for_identifier(
            line_data_with_combined_cols,
            identifier,
            unique_mutations,
            mutation_cols[identifier],
            genotype_cols[identifier],
        )

    return line_data_with_combined_cols


def _make_combined_genotype_column_for_identifier(
    line_data: pd.DataFrame,
    identifier: Identifier,
    unique_mutations: list[str],
    mutation_cols: list[str],
    genotype_cols: list[str],
) -> None:
    """Add a genotype_IDENTIFIER column summarising all genotype columns.

    E.g. combining Grade 1 / 2 / 3 into a single genotype_offspring column.
    Any missing genotypes are assumed to be wildtype.

    Parameters
    ----------
    line_data : pd.DataFrame
        Data for a single line.
    identifier : Identifier
        The identifier to summarise.
    unique_mutations : list[str]
        The unique mutations for this line. Genotypes in genotype_IDENTIFIER
        will have length equal to this, and be returned in this order.
    mutation_cols : list[str]
        Mutation columns for the given identifier.
    genotype_cols : list[str]
        Genotype columns for the given identifier.
    """

    pivoted_mutations = pd.DataFrame()

    # pivot each pair of mutation / genotype columns. E.g. if Mutation 1 /
    # Grade 1 had rows with a mix of Mutation-A and Mutation-B: this would
    # produce two columns named 'Mutation-A' and 'Mutation-B', with the
    # genotypes as the column values.
    for mutation_col, genotype_col in zip(mutation_cols, genotype_cols):
        pivoted_cols = line_data.pivot(
            columns=mutation_col, values=genotype_col
        )
        # drop columns named NaN
        pivoted_cols = pivoted_cols.loc[:, pivoted_cols.columns.notna()]

        # If all values were NaN for this Mutation/Grade combo
        if pivoted_cols.empty:
            continue

        if pivoted_mutations.empty:
            pivoted_mutations = pivoted_cols
        else:
            # If there are matching column names, use the new pivoted_col to
            # fill na values
            pivoted_mutations = pivoted_mutations.fillna(pivoted_cols)

            # Merge any new column names
            common_cols = list(
                set(pivoted_mutations.columns).intersection(
                    pivoted_cols.columns
                )
            )
            pivoted_cols = pivoted_cols.drop(common_cols, axis=1)
            pivoted_mutations = pivoted_mutations.join(pivoted_cols)

    # Any remaining NaN values are assumed to be wildtype
    pivoted_mutations = pivoted_mutations.fillna("wt")

    # If a mutation name doesn't have a corresponding column -> assume all wt
    for mutation in unique_mutations:
        if mutation not in pivoted_mutations:
            pivoted_mutations[mutation] = "wt"

    # Combine pivoted mutations into a single summary column
    line_data[f"genotype_{identifier.name.lower()}"] = pivoted_mutations[
        unique_mutations
    ].agg("_".join, axis=1)
