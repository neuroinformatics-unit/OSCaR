from enum import Enum
from pathlib import Path

import pandas as pd

from oscar_colony.breeding_scheme import BreedingScheme, Genotype


class Identifier(Enum):
    """
    Identifier for who a mutation / genotype refers to.
    """

    OFFSPRING = 0
    FATHER = 1
    MOTHER = 2


def standardise_pyrat_csv(
    input_csv: pd.DataFrame | Path | str,
) -> pd.DataFrame:
    """Standardise a csv file exported from pyRAT.

    Processing steps include:
    - Correcting or removing forbidden genotypes like +/-, Tg, ko/ko
    - adding columns for the number of mutations per line (n_mutations) and
    a summary of the mutation names (mutations)
    - adding summary columns for 'genotype_offspring', 'genotype_father_n' and
    'genotype_mother_n' that match the order of 'mutations'.
    - marking ungenotyped-offspring as NaN in the 'genotype_offspring' column
    - filling any missing genotypes with wildtype
    - removing columns that aren't needed for further processing steps

    Parameters
    ----------
    input_csv : pd.DataFrame | Path | str
        Csv file exported from pyRAT.

    Returns
    -------
    pd.DataFrame
        Standardised dataframe, ready for further processing
    """
    if isinstance(input_csv, (Path, str)):
        input_csv = pd.read_csv(input_csv)

    father_cols, mother_cols, mutation_cols, genotype_cols = (
        _create_column_name_dicts(input_csv)
    )

    all_mutation_cols_list = sum(mutation_cols.values(), [])
    all_genotype_cols_list = sum(genotype_cols.values(), [])

    required_cols = (
        ["ID", "Line / Strain (Name)", "DOB"]
        + list(father_cols.keys())
        + list(mother_cols.keys())
        + ["Sacrifice reason"]
        + all_mutation_cols_list
        + all_genotype_cols_list
    )

    rename_dict = {
        "ID": "ID_offspring",
        "Line / Strain (Name)": "line_name",
        "DOB": "date_of_birth",
        "Sacrifice reason": "sacrifice_reason",
        **father_cols,
        **mother_cols,
    }

    standard_csv = input_csv[required_cols].rename(columns=rename_dict)

    standard_csv = _filter_or_correct_genotypes(
        standard_csv, all_genotype_cols_list
    )

    standard_csv = _add_n_mutations_column(
        standard_csv, genotype_cols[(Identifier.OFFSPRING, 0)]
    )
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

    impossible_breeding_schemes = standard_csv.apply(
        _check_data_input_reliability, axis=1
    )
    standard_csv = standard_csv[~impossible_breeding_schemes]

    return standard_csv


def _add_n_mutations_column(
    standard_csv: pd.DataFrame, offspring_genotype_cols: list[str]
) -> pd.DataFrame:
    """Add column with number of mutations per line.

    Parameters
    ----------
    standard_csv : pd.DataFrame
        Dataframe to add column to
    offspring_genotype_cols : list[str]
        Offspring genotype columns e.g. Grade 1, Grade 2, Grade 3

    Returns
    -------
    pd.DataFrame
        Dataframe with n_mutations column added
    """
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

    return standard_csv


def _create_column_name_dicts(
    input_csv: pd.DataFrame,
) -> tuple[
    dict[str, str],
    dict[str, str],
    dict[tuple[Identifier, int], list[str]],
    dict[tuple[Identifier, int], list[str]],
]:
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

    n_mothers = len(input_csv.filter(regex=r"^Mother \d+$").columns)
    n_fathers = len(input_csv.filter(regex=r"^Father \d+$").columns)

    father_cols = {
        f"Father {i}": f"ID_father_{i}" for i in range(1, n_fathers + 1)
    }
    mother_cols = {
        f"Mother {i}": f"ID_mother_{i}" for i in range(1, n_mothers + 1)
    }

    mutation_dict: dict = {}
    genotype_dict: dict = {}

    _prefix_unpacking(
        mutation_dict, genotype_dict, input_csv, (Identifier.OFFSPRING, 0), ""
    )

    father_cols = {}
    mother_cols = {}
    for i in range(1, n_fathers + 1):
        father_cols[f"Father {i}"] = f"ID_father_{i}"
        _prefix_unpacking(
            mutation_dict,
            genotype_dict,
            input_csv,
            (Identifier.FATHER, i),
            f"Father {i}: ",
        )

    for i in range(1, n_mothers + 1):
        mother_cols[f"Mother {i}"] = f"ID_mother_{i}"
        _prefix_unpacking(
            mutation_dict,
            genotype_dict,
            input_csv,
            (Identifier.MOTHER, i),
            f"Mother {i}: ",
        )

    return father_cols, mother_cols, mutation_dict, genotype_dict


def _filter_or_correct_genotypes(
    standard_csv: pd.DataFrame, genotype_cols: list[str]
) -> pd.DataFrame:
    """Filter or correct rows so that only genotypes of wt, het or hom remain.

    Where possible, this will convert alternative forms to wt/het/hom e.g.
    ko/ko -> hom. If an un-ambiguous conversion isn't possible
    (like T, Tg, N, +, -), rows that contain these will be removed.

    Parameters
    ----------
    standard_csv : pd.DataFrame
        Dataframe to filter
    genotype_cols : list[str]
        Names of all genotype columns including offspring, father and mother

    Returns
    -------
    pd.DataFrame
        Dataframe with only wt, het or hom in genotype columns
    """

    genotype_conversions = {
        "ko/ko": Genotype.HOM,
        "ko/+": Genotype.HET,
        "ko/-": Genotype.HET,
        "+/ko": Genotype.HET,
        "-/ko": Genotype.HET,
        "ki/ki": Genotype.HOM,
        "ki/+": Genotype.HET,
        "ki/-": Genotype.HET,
        "+/ki": Genotype.HET,
        "-/ki": Genotype.HET,
    }

    # convert genotypes where possible
    genotype_data = standard_csv.loc[:, genotype_cols]
    for old_genotype, new_genotype in genotype_conversions.items():
        genotype_data = genotype_data.replace(
            to_replace=old_genotype, value=new_genotype.name.lower()
        )

    filtered_data = standard_csv.copy()
    filtered_data.loc[:, genotype_cols] = genotype_data

    # remove rows where any of the genotype values aren't in the allowed set:
    # wt, het, hom or empty
    allowed_genotypes = (
        genotype_data.isin([genotype.name.lower() for genotype in Genotype])
        | genotype_data.isna()
    ).all(axis=1)
    filtered_data = filtered_data.loc[allowed_genotypes, :]

    return filtered_data


def _make_combined_genotype_columns_for_line(
    line_data: pd.DataFrame,
    mutation_cols: dict[tuple[Identifier, int], list[str]],
    genotype_cols: dict[tuple[Identifier, int], list[str]],
) -> pd.DataFrame:
    """For data from a single line, add columns for 'mutations',
    'genotype_offspring', 'genotype_father' and 'genotype_mother'.

    All genotype columns list genotypes in the same order as given in
    'mutations'. If all the offspring genotype columns are empty, they
    are assumed to be un-genotyped (i.e. their genotype was never checked,
    and is unknown) - in these cases, the 'genotype_offspring' value will
    be left empty. In all other cases, individual missing genotypes are
    assumed to be wildtype.

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
        line_data[mutation_cols[(Identifier.OFFSPRING, 0)]].values.ravel("K")
    )
    unique_mutations = list(pd.Series(unique_mutations).dropna().astype(str))

    # Copy so we don't edit the original dataframe (this can cause issues
    # with apply)
    line_data_with_combined_cols = line_data.copy()
    line_data_with_combined_cols["mutations"] = "_".join(unique_mutations)

    for identifier_key in mutation_cols:
        _make_combined_genotype_column_for_identifier(
            line_data_with_combined_cols,
            identifier_key,
            unique_mutations,
            mutation_cols[identifier_key],
            genotype_cols[identifier_key],
        )

    return line_data_with_combined_cols


def _make_combined_genotype_column_for_identifier(
    line_data: pd.DataFrame,
    identifier_key: tuple[Identifier, int],
    unique_mutations: list[str],
    mutation_cols: list[str],
    genotype_cols: list[str],
) -> None:
    """Add a genotype_IDENTIFIER column summarising all genotype columns.

    E.g. combining Grade 1 / 2 / 3 into a single genotype_offspring column.

    All individual missing genotypes are assumed to be wildtype, except in
    the case of un-genotyped offspring (identifier == OFFSPRING and all
    genotype columns empty) - these are left empty.

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

    identifier, parent_num = identifier_key

    pivoted_mutations = pd.DataFrame(index=line_data.index)
    wildtype_str = Genotype.WT.name.lower()

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

    # Add columns for any missing mutation names
    for mutation in unique_mutations:
        if mutation not in pivoted_mutations:
            pivoted_mutations[mutation] = pd.Series(dtype=str)

    if identifier == Identifier.OFFSPRING:
        # If all offspring mutations in a row are NaN, leave as-is -> these are
        # un-genotyped individuals.
        # If only some are NaN, then fill with wt
        genotyped_rows = ~pivoted_mutations.isna().all(axis=1)
        pivoted_mutations.loc[genotyped_rows, :] = pivoted_mutations.loc[
            genotyped_rows, :
        ].fillna(wildtype_str)
    else:
        # Fill wildtype for rows where a parent is actually recorded.
        if identifier == Identifier.FATHER:
            parent_id_col = f"ID_father_{parent_num}"
        else:
            parent_id_col = f"ID_mother_{parent_num}"

        parent_recorded = line_data[parent_id_col].notna()
        pivoted_mutations.loc[parent_recorded, :] = pivoted_mutations.loc[
            parent_recorded, :
        ].fillna(wildtype_str)

    # Combine pivoted mutations into a single summary column
    if identifier == Identifier.OFFSPRING:
        new_col_name = "genotype_offspring"
    else:
        new_col_name = f"genotype_{identifier.name.lower()}_{parent_num}"

    line_data[new_col_name] = pd.Series(dtype=str)
    genotyped_rows = ~pivoted_mutations.isna().all(axis=1)
    line_data.loc[genotyped_rows, new_col_name] = pivoted_mutations.loc[
        genotyped_rows, unique_mutations
    ].agg("_".join, axis=1)


def _check_data_input_reliability(
    standardised_df_row: pd.Series,
) -> bool:

    ambigious_parentage = _is_ambigious_parentage(standardised_df_row)
    impossible_breeding_schemes = _is_impossible_breeding_scheme(
        standardised_df_row
    )

    if impossible_breeding_schemes or ambigious_parentage:
        return True
    return False


def _is_impossible_breeding_scheme(
    standardised_df_row: pd.Series,
) -> bool:
    """Checks whether the given row contains an impossible breeding scheme.

    Retrieves parent genotypes and pulls the mendelian ratios from
    BreedingScheme. Compares offspring to these ratios, returning True for
    those which are not possible.
    e.g. hom x hom parents cannot produce wt offspring.

    Parameters
    ----------
    standardised_df_row : pd.Series
        row from standardised_dataframe (pd.DataFrame): standardised PyRAT df

    Returns
    -------
    bool
        bool of whether or not that row contains an impossible breeding scheme
    """

    genotype_father = standardised_df_row["genotype_father_1"]
    genotype_mother = standardised_df_row["genotype_mother_1"]
    genotype_offspring = standardised_df_row["genotype_offspring"]

    # Only processes when offspring is assigned a genotype
    if not pd.isna(genotype_offspring):
        typed_offspring = Genotype.from_string(genotype_offspring)
        scheme = BreedingScheme(genotype_father, genotype_mother)
        ratio = scheme.mendelian_ratio()

        if typed_offspring not in ratio:
            return True
        elif ratio[typed_offspring] == 0:
            return True

    return False


def _is_ambigious_parentage(standardised_df_row: pd.Series) -> bool:
    """_summary_

    Parameters
    ----------
    standardised_df_row : pd.Series
        _description_

    Returns
    -------
    bool
        _description_
    """

    mother = r"genotype_mother_\d+$"
    father = r"genotype_father_\d+$"

    for parent in [mother, father]:
        n_parent = len(standardised_df_row.filter(regex=parent).index)

        if n_parent == 1:
            continue

        parent_genotypes: list = (
            standardised_df_row.filter(regex=parent).dropna().values.tolist()
        )

        standard_genotype = sorted(parent_genotypes[0].split("_"))

        for i, parent_genotype in enumerate(parent_genotypes):
            if i == 0:
                standard_genotype = sorted(parent_genotypes[0].split("_"))
                continue

            parent_genotype = sorted(parent_genotype.split("_"))
            if parent_genotype != standard_genotype:
                return True

    return False


def _prefix_unpacking(
    mutation_dict: dict,
    genotype_dict: dict,
    input_csv: pd.DataFrame,
    identifier: tuple,
    prefix: str,
):
    # columns of form 'PREFIXMutation NUMBER'
    mutation_cols = sorted(
        input_csv.filter(regex=rf"^{prefix}Mutation \d$").columns.tolist()
    )

    # columns of form 'PREFIXGrade NUMBER'
    genotype_cols = sorted(
        input_csv.filter(regex=rf"^{prefix}Grade \d$").columns.tolist()
    )

    # Each mutation must have a corresponding genotype
    if len(mutation_cols) != len(genotype_cols):
        raise ValueError(
            f"Not all {identifier} mutation columns have a corresponding "
            f"genotype column."
        )

    # Make sure lists are in numeric order e.g. Grade 1, Grade 2, Grade 3
    mutation_dict[identifier] = mutation_cols
    genotype_dict[identifier] = genotype_cols
