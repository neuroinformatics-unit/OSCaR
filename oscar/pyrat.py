from pathlib import Path

import pandas as pd

from oscar.breeding_scheme import Genotype


def standardise_pyrat_csv(input_csv: pd.DataFrame | Path) -> pd.DataFrame:
    if isinstance(input_csv, Path):
        input_csv = pd.read_csv(input_csv)

    offspring_mutation_cols, offspring_genotype_cols = (
        _get_mutation_and_genotype_columns(input_csv)
    )
    father_mutation_cols, father_genotype_cols = (
        _get_mutation_and_genotype_columns(input_csv, prefix="Father: ")
    )
    mother_mutation_cols, mother_genotype_cols = (
        _get_mutation_and_genotype_columns(input_csv, prefix="Mother: ")
    )

    all_genotype_cols = (
        offspring_genotype_cols + father_genotype_cols + mother_genotype_cols
    )
    all_mutation_cols = (
        offspring_mutation_cols + father_mutation_cols + mother_mutation_cols
    )

    required_cols = (
        [
            "ID",
            "Line / Strain (Name)",
            "DOB",
            "Father",
            "Mother",
            "Sacrifice reason",
        ]
        + all_mutation_cols
        + all_genotype_cols
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

    standard_csv = _filter_allowed_genotypes(standard_csv, all_genotype_cols)
    standard_csv = _filter_ungenotyped(standard_csv, offspring_genotype_cols)

    standard_csv["n_mutations"] = (
        standard_csv.loc[:, offspring_genotype_cols].notna().sum(axis=1)
    )

    # make sure number of mutations is the same throughout each line -
    # use the max.
    # Sometimes particular individuals are ungenotyped (n_mutations = 0) or a
    # genotype value is omitted to mean wt.
    standard_csv["n_mutations"] = standard_csv.groupby("line_name")[
        "n_mutations"
    ].transform("max")

    _make_combined_genotype_columns(standard_csv)

    # Combine grade 1/2/3 columns into an overall genotype
    _make_combined_genotype_column(
        standard_csv, offspring_genotype_cols, "genotype_offspring"
    )
    _make_combined_genotype_column(
        standard_csv, father_genotype_cols, "genotype_father"
    )
    _make_combined_genotype_column(
        standard_csv, mother_genotype_cols, "genotype_mother"
    )

    standard_csv = standard_csv.drop(
        all_genotype_cols,
        axis=1,
    )
    return standard_csv


def _get_mutation_and_genotype_columns(
    input_csv: pd.DataFrame, prefix: str = ""
) -> tuple[list[str], list[str]]:
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

    return mutation_cols, genotype_cols


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
    long_format = pd.wide_to_long(
        standard_csv,
        stubnames=[
            "Mutation",
            "Father: Mutation",
            "Mother: Mutation",
            "Grade",
            "Father: Grade",
            "Mother: Grade",
        ],
        i="ID_offspring",
        j="mutation_number",
        sep=" ",
    ).reset_index()

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


def _make_combined_genotype_column_for_line(
    line_data: pd.DataFrame, mutation_cols: list[str], genotype_cols: list[str]
):
    # convert dataframe to 'long' format i.e. combine all mutation cols into
    # one, and all genotype cols into one. Each ID_offspring will have multiple
    # rows - one per mutation.
    print(line_data)
    print(line_data.columns)
