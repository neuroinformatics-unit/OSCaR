from pathlib import Path

import pandas as pd


def standardise_pyrat_csv(input_csv: pd.DataFrame | Path) -> pd.DataFrame:
    if isinstance(input_csv, Path):
        input_csv = pd.read_csv(input_csv)

    offspring_genotype_cols = ["Grade 1", "Grade 2", "Grade 3"]
    father_genotype_cols = [
        "Father: Grade 1",
        "Father: Grade 2",
        "Father: Grade 3",
    ]
    mother_genotype_cols = [
        "Mother: Grade 1",
        "Mother: Grade 2",
        "Mother: Grade 3",
    ]

    required_cols = (
        ["ID", "Line / Strain (Name)", "DOB", "Sacrifice reason"]
        + offspring_genotype_cols
        + father_genotype_cols
        + mother_genotype_cols
    )

    # Get rid of any additional columns + rename to standard names
    standard_csv = input_csv[required_cols]
    standard_csv = standard_csv.rename(
        columns={
            "Line / Strain (Name)": "line_name",
            "DOB": "date_of_birth",
            "Sacrifice reason": "sacrifice_reason",
        }
    )

    standard_csv["n_mutations"] = (
        standard_csv["Grade 1"].notna().astype(int)
        + standard_csv["Grade 2"].notna().astype(int)
        + standard_csv["Grade 3"].notna().astype(int)
    )
    # make sure number of mutations is the same throughout each line -
    # use the max.
    # Sometimes particular individuals are ungenotyped (n_mutations = 0) or a
    # grade 1/2/3 value is omitted to mean wt.
    standard_csv["n_mutations"] = standard_csv.groupby("line_name")[
        "n_mutations"
    ].transform("max")

    # Combine grade 1/2/3 columns into an overall genotype
    _make_combined_genotype_column(
        standard_csv, offspring_genotype_cols, "offspring_genotype"
    )
    _make_combined_genotype_column(
        standard_csv, father_genotype_cols, "father_genotype"
    )
    _make_combined_genotype_column(
        standard_csv, mother_genotype_cols, "mother_genotype"
    )

    standard_csv = standard_csv.drop(
        offspring_genotype_cols + father_genotype_cols + mother_genotype_cols,
        axis=1,
    )
    return standard_csv


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
