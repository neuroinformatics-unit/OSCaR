import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from oscar_colony.breeding_scheme import BreedingScheme, Genotype

logger = logging.getLogger(__name__)


def standardise_pyrat_csv(
    input_df: pd.DataFrame | Path | str,
) -> pd.DataFrame:
    """Standardise a csv file exported from pyRAT.

    Processing steps include:
    - standardising column names with a dynamic dict
    - adding columns for the number of mutations per line (n_mutations) and
    a summary of the mutation names (mutations)
    - Correcting or removing forbidden genotypes like +/-, Tg, ko/ko
    - adding summary columns for 'genotype_offspring', 'genotype_father' and
    'genotype_mother' that match the order of 'mutations'.
    - marking ungenotyped-offspring as NaN in the 'genotype_offspring' column
    - filling any missing genotypes with wildtype
    - removing columns that aren't needed for further processing steps
    - checking data input validity and removing impossible input data

    Parameters
    ----------
    input_csv : pd.DataFrame | Path | str
        Csv file exported from pyRAT.

    Returns
    -------
    pd.DataFrame
        Standardised dataframe, ready for further processing
    """

    if isinstance(input_df, (Path, str)):
        input_df = pd.read_csv(input_df)

    logger.info(
        f"Starting standardisation of pyRAT data: {len(input_df)} rows"
    )

    rename_col_dict = _create_rename_dict(input_df)
    mutation_cols, genotype_cols = _create_mutation_genotype_dicts(input_df)

    all_mutation_cols_list = sum(mutation_cols.values(), [])
    all_genotype_cols_list = sum(genotype_cols.values(), [])

    # Uses keys as previous column name
    required_cols = (
        list(rename_col_dict.keys())
        + all_mutation_cols_list
        + all_genotype_cols_list
    )

    standard_df = input_df[required_cols].rename(columns=rename_col_dict)

    standard_df = _filter_or_correct_genotypes(
        standard_df, all_genotype_cols_list
    )

    standard_df = _add_n_mutations_column(
        standard_df, genotype_cols["offspring"]
    )

    standard_df = standard_df.groupby("line_name").apply(
        _make_combined_genotype_columns_for_line, mutation_cols, genotype_cols
    )

    standard_df = _filter_data_input_validity(standard_df)

    standard_df = _collapse_parent_genotype(standard_df)
    standard_df = standard_df.reset_index().drop(
        "level_1", axis=1, errors="ignore"
    )

    # for readability, make sure ID_offspring is first
    id_offspring_col = standard_df.pop("ID_offspring")
    standard_df.insert(0, "ID_offspring", id_offspring_col)

    _log_ungenotyped_animals(standard_df)

    logger.info(f"Standardisation complete: {len(standard_df)} rows")

    return standard_df


def _log_ungenotyped_animals(standard_df: pd.DataFrame) -> None:
    """Log the count and IDs of offspring with no genotype recorded.

    Parameters
    ----------
    standard_df : pd.DataFrame
        Standardised dataframe to check for ungenotyped offspring.
    """
    ungenotyped = standard_df["genotype_offspring"].isna()
    n_ungenotyped = ungenotyped.sum()
    if n_ungenotyped > 0:
        ungenotyped_ids = standard_df.loc[ungenotyped, "ID_offspring"].tolist()
        logger.info(f"{n_ungenotyped} offspring have no genotype recorded")
        logger.debug(f"Offspring with no genotypes are : {ungenotyped_ids}")


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


def _create_rename_dict(input_csv: pd.DataFrame) -> dict[str, str]:
    """Create renaming dict for generic columns and parent_x ID columns.

    parents_x being any number of mother and father columns, labelled as
    Mother x or Father x and assigning them to the standardised form of
    ID_mother_x or ID_father_x

    Parameters
    ----------
    input_csv : pd.DataFrame
        Dataframe to extract column names from.

    Returns
    -------
    dict[str, str]
        Returns a dictionary where the key is the old column names to be
        replaced with the value in future processing.
    """

    # re-name standard columns
    rename_dict = {
        "ID": "ID_offspring",
        "Line / Strain (Name)": "line_name",
        "DOB": "date_of_birth",
    }

    # re-name any number of parent columns
    parent_cols = []
    for col_name in input_csv.columns:
        m = re.match(r"^(Mother|Father) (\d+)$", col_name)
        if m:
            new_name = f"ID_{m.group(1).lower()}_{m.group(2)}"
            parent_cols.append((col_name, new_name))

    # sorts columns before concatenating so they are in the correct order
    parent_cols.sort()
    rename_dict = rename_dict | dict(parent_cols)
    rename_dict["Sacrifice reason"] = "sacrifice_reason"

    return rename_dict


def _create_mutation_genotype_dicts(
    input_df: pd.DataFrame,
) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Create dicts of mutation / genotype column names for all identifiers.

    Uses _sort_and_name_columns_by_prefix to assign each identifier
    (offspring, father_n, mother_n) to its mutation and genotype columns.

    Parameters
    ----------
    input_csv : pd.DataFrame
        Dataframe to extract column names from.

    Returns
    -------
    tuple[dict[str, list[str]], dict[str, list[str]]]
        mutation/genotype dictionaries, keyed by identifier string
        (offspring, father_n, mother_n). Value is the list of
        genotype/mutation column names
    """

    n_fathers = len(input_df.filter(regex=r"^Father \d+$").columns)
    n_mothers = len(input_df.filter(regex=r"^Mother \d+$").columns)

    mutation_dict: dict = {}
    genotype_dict: dict = {}

    _sort_and_name_columns_by_prefix(
        mutation_dict, genotype_dict, input_df, "offspring", ""
    )

    for i in range(1, n_fathers + 1):
        _sort_and_name_columns_by_prefix(
            mutation_dict,
            genotype_dict,
            input_df,
            f"father_{i}",
            f"Father {i}: ",
        )

    for i in range(1, n_mothers + 1):
        _sort_and_name_columns_by_prefix(
            mutation_dict,
            genotype_dict,
            input_df,
            f"mother_{i}",
            f"Mother {i}: ",
        )

    return mutation_dict, genotype_dict


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
    filtered_count = len(filtered_data)
    removed_count = len(standard_csv) - filtered_count

    if removed_count > 0:
        dropped_rows = standard_csv.loc[~allowed_genotypes]

        # Creates a dictionary of {linename: number of dropped rows}
        lines_affected = dropped_rows["line_name"].value_counts().to_dict()

        # Finds and extracts which genotype caused the dropped values.
        # Deduplicate the raw values first, since a filtered slice can have
        # far more cells than distinct invalid genotype values.
        allowed_genotype_names = {
            genotype.name.lower() for genotype in Genotype
        }
        dropped_values = genotype_data.loc[~allowed_genotypes].to_numpy()
        invalid_values = sorted(
            {
                str(value)
                for value in set(dropped_values.ravel())
                if pd.notna(value) and str(value) not in allowed_genotype_names
            }
        )

        logger.info(
            f"Filtered out {removed_count} invalid genotype row(s) - "
            f"{filtered_count} row(s) remaining. "
            f"Invalid genotype value(s): {invalid_values}. "
            f"Lines affected: {lines_affected}."
        )
        logger.debug(
            f"Filtered offspring IDs: {dropped_rows['ID_offspring'].tolist()}"
        )

    return filtered_data


def _make_combined_genotype_columns_for_line(
    line_data: pd.DataFrame,
    mutation_cols: dict[str, list[str]],
    genotype_cols: dict[str, list[str]],
) -> pd.DataFrame:
    """For data from a single line, standardise the mutation order and
    create summary columns for 'genotype_offspring', 'genotype_father' and
    'genotype_mother'.

    All existing mutation / genotype columns will be removed. New mutation
    columns (with mutations in a consistent order across the line) will be
    added, numbered like mutation_1, mutation_2... All summary columns
    list genotypes in the same order e.g. wt_het is wt for mutation_1 and
    het for mutation_2.

    If all the offspring genotype columns are empty, they
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
        Line data with standard mutation and summary genotype columns
    """

    # get unique offspring mutations for this line
    unique_mutations = pd.unique(
        line_data[mutation_cols["offspring"]].values.ravel("K")
    )
    unique_mutations = sorted(pd.Series(unique_mutations).dropna().astype(str))

    # Copy so we don't edit the original dataframe (this can cause issues
    # with apply)
    line_data_with_combined_cols = line_data.copy()

    for identifier_key in mutation_cols:
        _make_combined_genotype_column_for_identifier(
            line_data_with_combined_cols,
            identifier_key,
            unique_mutations,
            mutation_cols[identifier_key],
            genotype_cols[identifier_key],
        )

    # Add column for each mutation IN ORDER (for the sake of readability,
    # add next to the n_mutations column)
    n_mutations_index = line_data_with_combined_cols.columns.get_loc(
        "n_mutations"
    )
    for i, mutation in enumerate(unique_mutations):
        line_data_with_combined_cols.insert(
            n_mutations_index + (i + 1), f"mutation_{i + 1}", mutation
        )

    return line_data_with_combined_cols


def _make_combined_genotype_column_for_identifier(
    line_data: pd.DataFrame,
    identifier_key: str,
    unique_mutations: list[str],
    mutation_cols: list[str],
    genotype_cols: list[str],
) -> None:
    """Combine all mutation / genotype columns for an identifier, into a single
    summary genotype_IDENTIFIER column.

    E.g. removing Mutation 1 / 2 / 3 and Grade 1 / 2 / 3 columns, and
    adding a single combined genotype_offspring column.

    All individual missing genotypes are assumed to be wildtype, except in
    the case of un-genotyped offspring, these are left empty.

    Parameters
    ----------
    line_data : pd.DataFrame
        Data for a single line.
    identifier_key : str
        The identifier to summarise: "offspring", "father_n" or "mother_n".
    unique_mutations : list[str]
        The unique mutations for this line. Genotypes will have length equal
        to this, and be returned in this order.
    mutation_cols : list[str]
        Mutation columns for the given identifier_key.
    genotype_cols : list[str]
        Genotype columns for the given identifier_key.
    """

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

    if identifier_key == "offspring":
        # If all offspring mutations in a row are NaN, leave as-is -> these are
        # un-genotyped individuals.
        # If only some are NaN, then fill with wt
        genotyped_rows = ~pivoted_mutations.isna().all(axis=1)
        pivoted_mutations.loc[genotyped_rows, :] = pivoted_mutations.loc[
            genotyped_rows, :
        ].fillna(wildtype_str)
    else:
        # Fill wildtype for rows where a parent is actually recorded.
        parent_id_col = f"ID_{identifier_key}"
        parent_recorded = line_data[parent_id_col].notna()
        pivoted_mutations.loc[parent_recorded, :] = pivoted_mutations.loc[
            parent_recorded, :
        ].fillna(wildtype_str)

    # Combine pivoted mutations into a single summary column
    new_col_name = f"genotype_{identifier_key}"

    line_data[new_col_name] = pd.Series(dtype=str)
    genotyped_rows = ~pivoted_mutations.isna().all(axis=1)
    line_data.loc[genotyped_rows, new_col_name] = pivoted_mutations.loc[
        genotyped_rows, unique_mutations
    ].agg("_".join, axis=1)

    # Drop old mutation / grade columns
    line_data.drop(columns=mutation_cols + genotype_cols, inplace=True)


def _filter_data_input_validity(standard_df: pd.DataFrame) -> pd.DataFrame:
    """Removes rows containing invalid data.

    Runs _check_data_input_validity on each row, if that returns an issue for
    a particular row, then that row is removed from the final DataFrame.

    Parameters
    ----------
    standard_df : pd.DataFrame
        DataFrame to filter

    Returns
    -------
    pd.DataFrame
        filtered DataFrame
    """
    mother_col_names = standard_df.filter(
        regex=r"genotype_mother_\d+$"
    ).columns.tolist()
    father_col_names = standard_df.filter(
        regex=r"genotype_father_\d+$"
    ).columns.tolist()

    filter_reasons = standard_df.apply(
        _check_data_input_validity,
        axis=1,
        mother_col_names=mother_col_names,
        father_col_names=father_col_names,
    )

    # Returns True if filter reasons has a value
    impossible_input_data = filter_reasons.notna()
    removed_count = impossible_input_data.sum()
    filtered_df = standard_df[~impossible_input_data]

    if removed_count > 0:
        dropped_rows = standard_df.loc[impossible_input_data]

        # Creates a dictionary of {linename: number of dropped rows}
        lines_affected = (
            dropped_rows.index.get_level_values("line_name")
            .value_counts()
            .to_dict()
        )
        reason_counts = (
            filter_reasons.loc[impossible_input_data].value_counts().to_dict()
        )
        logger.info(
            f"Filtered out {removed_count} row(s) with invalid breeding "
            f"data - {len(filtered_df)} row(s) remaining. "
            f"Reasons: {reason_counts}. "
            f"Lines affected: {lines_affected}."
        )
        logger.debug(
            f"Filtered offspring IDs: {dropped_rows['ID_offspring'].tolist()}"
        )

    return filtered_df


def _check_data_input_validity(
    standardised_df_row: pd.Series,
    mother_col_names: list[str],
    father_col_names: list[str],
) -> str | None:
    """Checks a Dataframe row for common data input errors.

    Takes a row from the standardised df, and runs two functions that test the
    validity of recorded data. Whether each sex of parent have the same
    genotype, or whether the breeding scheme is possible. If either of these
    detect an issue, then this function will flag for removal.

    Parameters
    ----------
    standardised_df_row : pd.Series
        row from standardised_dataframe (pd.DataFrame): standardised PyRAT df

    mother_col_names: list[str]
        a list of genotype column names for any number of mothers in the
        standardised_df

    father_col_names: list[str]
        a list of genotype column names for any number of fathers in the
        standardised_df

    Returns
    -------
    str | None
        The reason to filter this row out ("ambiguous parentage" or
        "impossible breeding scheme"), or None if the row is valid.
    """

    offspring_genotype = standardised_df_row["genotype_offspring"]
    mother_genotypes = (
        standardised_df_row[mother_col_names].dropna().to_numpy(dtype=object)
    )
    father_genotypes = (
        standardised_df_row[father_col_names].dropna().to_numpy(dtype=object)
    )

    if _is_ambiguous_parentage(mother_genotypes, father_genotypes):
        return "ambiguous parentage"

    mother_genotype = mother_genotypes[0]
    father_genotype = father_genotypes[0]

    if _is_impossible_breeding_scheme(
        offspring_genotype, mother_genotype, father_genotype
    ):
        return "impossible breeding scheme"

    return None


def _is_impossible_breeding_scheme(
    offspring_genotype: str,
    mother_genotype: str,
    father_genotype: str,
) -> bool:
    """Checks whether the given row contains an impossible breeding scheme.

    Retrieves parent genotypes and pulls the mendelian ratios from
    BreedingScheme. Compares offspring to these ratios, returning True for
    those which are not possible.
    e.g. hom x hom parents cannot produce wt offspring.

    Parameters
    ----------
    offspring_genotype: str
        string of the offspring genotype from standardised dataframe row

    mother_genotype: str
        the first mother genotype in the standardised dataframe row

    father_genotype: str
        the first father genotype in the standardised dataframe row

    Returns
    -------
    bool
        bool of whether or not that row contains an impossible breeding scheme
    """

    # Only processes when offspring is assigned a genotype
    if not pd.isna(offspring_genotype):
        typed_offspring = Genotype.from_string(offspring_genotype)
        scheme = BreedingScheme(father_genotype, mother_genotype)
        ratio = scheme.mendelian_ratio()

        if typed_offspring not in ratio:
            logger.debug(
                f"Genotype {typed_offspring} is not possible "
                f"for parents {father_genotype} x {mother_genotype} "
                f"using the Mendelian ratio"
            )
            return True
        elif ratio[typed_offspring] == 0:
            logger.debug(
                f"Possibility of genotype {typed_offspring} is 0% "
                f"for parents {father_genotype} x {mother_genotype}"
            )
            return True

    return False


def _is_ambiguous_parentage(
    mother_genotypes: np.ndarray,
    father_genotypes: np.ndarray,
) -> bool:
    """checks if parent exists and that all same sex genotypes are equal

    Parameters
    ----------
    mother_genotypes: np.ndarray
        an array of genotypes for any number of mothers in the standardised_df

    father_genotypes: np.ndarray
        an array of genotypes for any number of fathers in the standardised_df

    Returns
    -------
    bool
        True if parent ambiguity detected, False if not.
    """

    for parent_genotypes in [mother_genotypes, father_genotypes]:
        if len(set(parent_genotypes)) != 1:
            return True

    return False


def _sort_and_name_columns_by_prefix(
    mutation_dict: dict,
    genotype_dict: dict,
    input_csv: pd.DataFrame,
    identifier: str,
    prefix: str,
):
    """Assigns a given identifier to the corresponding mutation and grade.

    Called for offspring, mothers and fathers. Uses the chosen prefix to sort
    through DataFrame columns and retrieve a sorted list. Then populated two
    dictionaries for both mutation and genotype, with identifier as the key.

    Parameters
    ----------
    mutation_dict : dict
        dictionary to append [identifier] = mutation column names
    genotype_dict : dict
        dictionary to append [identifier] = genotype column names
    input_csv : pd.DataFrame
        Dataframe to filter through the columns of
    identifier : str
        animal identifier: "offspring", "father_n" or "mother_n".
    prefix : str
        prefix of column names to select from input_csv
    """
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

    mutation_dict[identifier] = mutation_cols
    genotype_dict[identifier] = genotype_cols


def _collapse_parent_genotype(standardised_df: pd.DataFrame) -> pd.DataFrame:
    """Collapses multiple same sex parent genotypes into just one"""

    for parent in ["mother", "father"]:
        genotype_col_name = f"genotype_{parent}_1"

        genotype_columns = standardised_df.filter(
            regex=rf"^genotype_{parent}_\d+$"
        ).columns.tolist()

        for column in genotype_columns:
            if column != genotype_col_name:
                standardised_df = standardised_df.drop(columns=column)

        standardised_df = standardised_df.rename(
            columns={genotype_col_name: f"genotype_{parent}"}
        )

    return standardised_df
