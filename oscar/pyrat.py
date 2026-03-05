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
    if isinstance(input_csv, Path):
        input_csv = pd.read_csv(input_csv)

    mutation_cols, genotype_cols = _create_column_name_dicts(input_csv)

    required_cols = (
        [
            "ID",
            "Line / Strain (Name)",
            "DOB",
            "Father",
            "Mother",
            "Sacrifice reason",
        ]
        + sum(mutation_cols.values(), [])
        + sum(genotype_cols.values(), [])
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
        standard_csv, sum(genotype_cols.values(), [])
    )
    standard_csv = _filter_ungenotyped(
        standard_csv, genotype_cols["offspring"]
    )

    standard_csv["n_mutations"] = (
        standard_csv.loc[:, genotype_cols["offspring"]].notna().sum(axis=1)
    )

    # make sure number of mutations is the same throughout each line -
    # use the max.
    # Sometimes particular individuals are ungenotyped (n_mutations = 0) or a
    # genotype value is omitted to mean wt.
    standard_csv["n_mutations"] = standard_csv.groupby("line_name")[
        "n_mutations"
    ].transform("max")

    standard_csv.groupby("line_name").apply(
        _make_combined_genotype_columns_for_line, mutation_cols, genotype_cols
    )

    _make_combined_genotype_columns(standard_csv)

    # # Combine grade 1/2/3 columns into an overall genotype
    # _make_combined_genotype_column(
    #     standard_csv, offspring_genotype_cols, "genotype_offspring"
    # )
    # _make_combined_genotype_column(
    #     standard_csv, father_genotype_cols, "genotype_father"
    # )
    # _make_combined_genotype_column(
    #     standard_csv, mother_genotype_cols, "genotype_mother"
    # )

    # standard_csv = standard_csv.drop(
    #     all_genotype_cols,
    #     axis=1,
    # )
    return standard_csv


def _create_column_name_dicts(
    input_csv: pd.DataFrame,
) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    keys = ["offspring", "father", "mother"]
    prefixes = ["", "Father: ", "Mother : "]

    mutation_dict = {}
    genotype_dict = {}

    for key, prefix in zip(keys, prefixes):
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
                f"Not all {key} mutation columns have a corresponding "
                f"genotype column."
            )

        # Make sure lists are in numeric order e.g. Grade 1, Grade 2, Grade 3
        mutation_dict[key] = sorted(mutation_cols)
        genotype_dict[key] = sorted(genotype_cols)

    return mutation_dict, genotype_dict


def _filter_allowed_genotypes(
    standard_csv: pd.DataFrame, genotype_cols: list[str]
) -> pd.DataFrame:
    """Only keep allowed genotypes of wt, het or hom.

    Remove others such as +/-, Tg, ko, as well as ungenotyped
    individuals (na for grade 1/2/3)

    Parameters
    ----------
    line_data : pd.DataFrame
        Data for a single line
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

    These offspring have na for all genotype columns.

    Parameters
    ----------
    standard_csv : pd.DataFrame
        _description_
    offspring_genotype_cols : list[str]
        _description_

    Returns
    -------
    pd.DataFrame
        _description_
    """
    ungenotyped = (
        standard_csv.loc[:, offspring_genotype_cols].isna().all(axis=1)
    )
    filtered_data = standard_csv.loc[~ungenotyped, :]

    return filtered_data


def _make_combined_genotype_columns(standard_csv: pd.DataFrame):
    # Combine columns for mutation/grade positions into one column.
    # i.e. combine Mutation 1/2/3 -> Mutation; Father: Grade 1/2/3 into
    # Father: Grade and so on...
    merged_columns = [
        "Mutation",
        "Father: Mutation",
        "Mother: Mutation",
        "Grade",
        "Father: Grade",
        "Mother: Grade",
    ]
    long_format = (
        pd.wide_to_long(
            standard_csv,
            stubnames=merged_columns,
            i="ID_offspring",
            j="mutation_number",
            sep=" ",
        )
        .reset_index()
        .sort_values("ID_offspring")
    )

    # Each should have one row per mutation name

    # Rename columns so that the individual the mutation / grade belongs to
    # (offspring, father or mother) is the suffix
    long_format = long_format.rename(
        columns={
            "Mutation": "Mutation_offspring",
            "Grade": "Grade_offspring",
            "Father: Mutation": "Mutation_father",
            "Father: Grade": "Grade_father",
            "Mother: Mutation": "Mutation_mother",
            "Mother: Grade": "Grade_mother",
        }
    )

    # Combine offspring / father / mother columns into one column.
    # i.e. Mutation_offspring/father/mother -> Mutation;
    # Grade_offspring/father/mother -> Grade
    long_format = (
        pd.wide_to_long(
            long_format,
            stubnames=["Mutation", "Grade"],
            i=["ID_offspring", "mutation_number"],
            j="mutation_belongs_to",
            sep="_",
            suffix=r"\w+",
        )
        .reset_index()
        .sort_values("ID_offspring")
    )

    long_format = long_format.dropna(
        axis=0, how="any", subset=["Mutation", "Grade"]
    )

    long_format["mutations"] = long_format.groupby(["line_name"])[
        "Mutation"
    ].transform(lambda x: [x.unique()] * len(x))
    # long_format.groupby(["line_name"])["Mutation"].transform(lambda x:
    # '_'.join(x.unique()))

    long_format.groupby("ID_offspring")

    # required_rows = long_format.groupby(["ID_offspring", "mutations"])

    genotype_offspring_col = pd.Series(dtype="str")
    for i, mut in enumerate(["Mut-A", "Mut-B"]):
        genotypes = long_format.loc[
            (long_format.Mutation == mut)
            & (long_format.mutation_belongs_to == "offspring"),
            "Grade",
        ]

        if i == 0:
            genotype_offspring_col = genotypes
        else:
            genotype_offspring_col = genotype_offspring_col + "_" + genotypes

    # mutations = long_format.groupby("line_name")["Mutation"].unique()
    long_format.groupby("line_name", "offspring")

    pass


def _make_combined_genotype_column(
    standard_csv: pd.DataFrame, cols_to_combine: list[str], new_col_name: str
) -> None:
    """Combine genotype columns for individual genes into one column.

    Parameters
    ----------
    standard_csv : pd.DataFrame
        Csv to read from / add column to
    cols_to_combine : list[str]
        List of columns to combine in gene order e.g.
        [grade 1, grade 2, grade 3]
    new_col_name : str
        Name of new combined column to add
    """

    standard_csv[new_col_name] = pd.Series(dtype="str")

    for i, col_name in enumerate(cols_to_combine):
        filled_col = standard_csv[col_name].fillna("wt")

        if i == 0:
            standard_csv.loc[standard_csv.n_mutations > i, new_col_name] = (
                filled_col
            )
        else:
            standard_csv.loc[standard_csv.n_mutations > i, new_col_name] = (
                standard_csv[new_col_name]
                + "_"
                + standard_csv[col_name].fillna("wt")
            )


def _make_combined_genotype_columns_for_line(
    line_data: pd.DataFrame,
    mutation_cols: dict[str, list[str]],
    genotype_cols: dict[str, list[str]],
):
    # get unique offspring mutations for this line
    unique_mutations = pd.unique(
        line_data[["Mutation 1", "Mutation 2", "Mutation 3"]].values.ravel("K")
    )
    unique_mutations = pd.Series(unique_mutations).dropna()

    # offspring first
    pivoted_mutations = pd.DataFrame()
    # pivot each pair of mutation / genotype columns. E.g. if Mutation 1 /
    # Grade 1 had rows with a mix of Mutation-A and Mutation-B: this would
    # produce two columns named 'Mutation-A' and 'Mutation-B', with the
    # genotype as the column values
    for mutation_col, genotype_col in zip(
        mutation_cols["father"], genotype_cols["father"]
    ):
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

    # If a mutation name doesn't have a corresponding column > assume all wt
    for mutation in unique_mutations:
        if mutation not in pivoted_mutations:
            pivoted_mutations[mutation] = "wt"

    line_data["combined_genotype"] = pivoted_mutations[unique_mutations].agg(
        "_".join, axis=1
    )
