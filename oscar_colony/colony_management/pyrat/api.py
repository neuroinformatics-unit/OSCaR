import datetime
import os
from typing import Any, Iterator

import pandas as pd
import requests


def get_pyrat_data(
    line_name: str | None = None,
    species_name: str | None = None,
    birth_date_from: datetime.date | None = None,
    birth_date_to: datetime.date | None = None,
    max_n_rows: int = 10000,
) -> Iterator[pd.DataFrame]:
    """Fetch animal data directly from the pyRAT api.

    To handle the potentially large number of animals returned from pyRAT,
    this function returns a generator of pandas dataframes (each with no
    more than max_n_rows).

    This expects PYRAT_URL, PYRAT_CLIENT_TOKEN and PYRAT_USER_TOKEN to
    be set as environment variables.

    Parameters
    ----------
    line_name : str | None, optional
        Name of line to fetch
    species_name : str | None, optional
        Name of species to fetch
    birth_date_from : datetime.date | None, optional
        Earliest birth date to include
    birth_date_to : datetime.date | None, optional
        Latest birth date to include
    max_no_rows : int, optional
        Maximum number of items in each returned dataframe (and therefore
        returned per request to the pyRAT api)

    Returns
    -------
    Iterator[pd.DataFrame]
        Generator of dataframes of returned animal data, in format matching
        that exported via the pyRAT UI. If no data is available for the query,
        the dataframe will be empty.
    """
    if (birth_date_to is not None and birth_date_from is not None) and (
        birth_date_to < birth_date_from
    ):
        raise ValueError("birth_date_to must be after birth_date_from")

    params = {
        "k": [
            "animalid",
            "eartag_or_id",
            "species_name",
            "strain_name",
            "dateborn",
            "mutations",
            "parents",
            "sacrifice_reason_name",
        ],
        "s": ["eartag_or_id:asc"],
        "state": ["live", "sacrificed", "exported"],
        "l": max_n_rows,
    }

    if line_name is not None:
        params["strain_name_with_id_like"] = line_name

    if species_name is not None:
        params["species"] = _get_species_id(species_name)

    if birth_date_from is not None:
        params["birth_date_from"] = birth_date_from.isoformat()

    if birth_date_to is not None:
        params["birth_date_to"] = birth_date_to.isoformat()

    # Make one request to determine how many results there are
    animals_response = _make_pyrat_request("animals", params)
    yield _convert_animals_to_df(animals_response.json())
    headers = animals_response.headers
    total_n = int(headers["x-total-count"])

    # If more results than max_n_rows, keep making requests and yielding result
    for start_n in range(max_n_rows, total_n, max_n_rows):
        params["o"] = start_n
        animals_response = _make_pyrat_request("animals", params)
        yield _convert_animals_to_df(animals_response.json())


def _make_pyrat_request(
    endpoint_name: str, params: dict[str, Any]
) -> requests.Response:
    """Make request to the pyRAT api.

    This expects PYRAT_URL, PYRAT_CLIENT_TOKEN and PYRAT_USER_TOKEN to
    be set as environment variables.

    Parameters
    ----------
    endpoint_name : str
        Name of endpoint e.g. 'species'
    params : dict[str, Any]
        Extra parameters to pass to the endpoint

    Returns
    -------
    requests.Response
        The requests response object, containing data from pyRAT
    """
    response = requests.get(
        url=f"{os.environ['PYRAT_URL']}/api/v3/{endpoint_name}",
        auth=(
            os.environ["PYRAT_CLIENT_TOKEN"],
            os.environ["PYRAT_USER_TOKEN"],
        ),
        params=params,
        timeout=5,  # number of seconds before timeout
    )

    # If the request didn't succeed, raise an error containing the status
    # code
    response.raise_for_status()

    return response


def _get_species_id(species_name: str) -> int:
    """Get pyRAT database ID for named species"""

    params = {
        "k": ["id", "name"],
        "s": ["name:asc"],
    }
    species_ids = _make_pyrat_request("species", params).json()

    available_names = []
    for species in species_ids:
        if species["name"] == species_name:
            return species["id"]
        else:
            available_names.append(species["name"])

    raise ValueError(
        f"No ID found for species {species_name}: available values "
        f"are {available_names}"
    )


def _get_parent_mutations_with_eartags(eartags: list[str]) -> pd.DataFrame:
    """Get parent mutation information for the given animal eartags"""

    params = {
        "k": ["animalid", "eartag_or_id", "mutations"],
        "s": ["eartag_or_id:asc"],
        "state": ["live", "sacrificed", "exported"],
        "eartag": eartags,
        "l": len(eartags),
    }
    mutation_data = _make_pyrat_request("animals", params).json()
    if len(mutation_data) != len(eartags):
        raise ValueError(
            f"{len(mutation_data)} animals returned for "
            f"{len(eartags)} eartags: {eartags}"
        )

    return pd.DataFrame(mutation_data)


def _convert_animals_to_df(animals_data: list[dict[str, Any]]) -> pd.DataFrame:
    """Convert animal data fetched from the pyRAT api to a pandas DataFrame.

    The structure / column names are matched to that exported from the pyRAT
    UI.
    """

    animals_df = pd.DataFrame(animals_data)
    if animals_df.empty:
        return animals_df

    # Convert dateborn to Year-Month-Day format (removing time info)
    new_dates = pd.to_datetime(animals_df.dateborn).dt.strftime("%Y-%m-%d")
    animals_df.dateborn = new_dates

    # Expand column with information for multiple mutations into their own
    # columns
    animals_df = _expand_mutations_data(animals_df)

    # Expand column with information about multiple parents into a new
    # dataframe, including their mutation info
    parents_df = _expand_parents_data(animals_df)
    animals_df = animals_df.drop(["parents"], axis=1)
    animals_df = animals_df.merge(parents_df, on="animalid", how="left")

    animals_df = animals_df.drop(["animalid"], axis=1)

    # re-name to match data exported via the pyRAT UI, to make downstream
    # analysis easier
    animals_df = animals_df.rename(
        columns={
            "eartag_or_id": "ID",
            "sacrifice_reason_name": "Sacrifice reason",
            "dateborn": "DOB",
            "strain_name": "Line / Strain (Name)",
            "species_name": "Species",
        }
    )

    return animals_df


def _expand_mutations_data(selected_df: pd.DataFrame) -> pd.DataFrame:
    """Expand a mutations column into a full dataframe.

    Each row of a mutations column contains a list of dictionaries
    (one per mutation for the animal). This function expands these into
    their own columns labelled Mutation 1, 2... and Grade 1, 2...

    Parameters
    ----------
    selected_df : pd.DataFrame
        DataFrame of animals data or expanded df, with raw mutations column

    column_prefix: str
        Prefix to add to expanded column names e.g. a prefix of 'Father: '
        would result in columns Father: Mutation 1, Father: Grade 1 etc.

    Returns
    -------
    pd.DataFrame
        Dataframe with separate Mutation and Grade columns
    """

    exploded_mutations_col = selected_df.mutations.explode()
    mutations_df = pd.DataFrame(
        exploded_mutations_col[~exploded_mutations_col.isna()].tolist()
    )

    mutation_col_name = "Mutation"
    grade_col_name = "Grade"

    # If no mutations are listed for any animals, return an empty Mutation 1 /
    # Grade 1 column
    if mutations_df.empty:
        selected_df = selected_df.drop(["mutations"], axis=1)
        selected_df[f"{mutation_col_name} 1"] = pd.Series(dtype=str)
        selected_df[f"{grade_col_name} 1"] = pd.Series(dtype=str)
        return selected_df

    mutations_df = mutations_df[["animalid", "mutationname", "mutationgrade"]]
    mutations_df = mutations_df.rename(
        columns={
            "mutationname": mutation_col_name,
            "mutationgrade": grade_col_name,
        }
    )

    # Adds a counter for the number of mutation rows per animal id
    mutations_df["count"] = (
        mutations_df.groupby("animalid").cumcount() + 1
    ).astype("string")

    # Create one row per animalid, with separate columns for
    # Mutation 1 / Grade 1, Mutation 2 / Grade 2 ...
    pivoted_mutations = mutations_df.pivot(
        columns="count",
        index="animalid",
        values=[mutation_col_name, grade_col_name],
    ).reset_index()
    pivoted_mutations.columns = [
        " ".join(column_names).strip()
        for column_names in pivoted_mutations.columns.to_flat_index()
    ]

    # merge into the original animals_df, so animalids are in the same order,
    # and any animals with no mutations appear with NaN in the correct slots
    merged_df = selected_df.drop(["mutations"], axis=1)
    merged_df = merged_df.merge(pivoted_mutations, on="animalid", how="left")

    return merged_df


def _add_empty_parent_cols(df: pd.DataFrame, parent: str) -> None:
    """Add empty columns for parent mutation and grade"""

    df[parent] = pd.Series(dtype=str)
    df[f"{parent}: Mutation 1"] = pd.Series(dtype=str)
    df[f"{parent}: Grade 1"] = pd.Series(dtype=str)


def _expand_parents_data(animals_df: pd.DataFrame) -> pd.DataFrame:
    """Expand column containing multiple parents' information into separate
    columns.

    This adds columns for Mother / Father IDs, as well as their respective
    mutations and grades.
    """

    exploded_parents_col = animals_df.parents.explode()
    parents_df = pd.DataFrame(
        exploded_parents_col[~exploded_parents_col.isna()].tolist()
    )

    # If no parents are listed for any animals, return empty mother / father
    # columns, with empty mutation / grade
    if parents_df.empty:
        animals_df = animals_df.loc[:, ["animalid"]]
        _add_empty_parent_cols(animals_df, "Mother 1")
        _add_empty_parent_cols(animals_df, "Father 1")
        return animals_df

    # Rename column names as father or mother
    expanded_df = parents_df[["animalid", "parent_eartag", "parent_sex"]]
    expanded_df = expanded_df.rename(columns={"parent_sex": "parent"})
    expanded_df.loc[expanded_df["parent"] == "m", "parent"] = "Father"
    expanded_df.loc[expanded_df["parent"] == "f", "parent"] = "Mother"

    # Number each consecutive parent appearance and append num to column name
    expanded_df["parent_id"] = (
        expanded_df.groupby(["animalid", "parent"]).cumcount() + 1
    )
    expanded_df["parent_id"] = (
        expanded_df["parent"] + " " + expanded_df["parent_id"].astype(str)
    )

    expanded_df, missing_parent = _merge_parent_mutations(expanded_df)

    merged_df = _dynamic_column_renaming(expanded_df, missing_parent)

    return merged_df


def _merge_parent_mutations(
    expanded_df: pd.DataFrame,
) -> tuple[pd.DataFrame, str]:
    """Separate parents by sex - essenital for the testing format
    retrieve mutations using parent eartags. Then update the df
    with their mutations"""

    # creates new dataframe for each sex
    mother_df = expanded_df[expanded_df["parent"] == "Mother"]
    father_df = expanded_df[expanded_df["parent"] == "Father"]

    parent_df_to_concat = []
    missing_parent = "None"
    # Obtain and format mutation data for mother - if present
    if not mother_df.empty:
        mother_mutations_df = _get_parent_mutations_with_eartags(
            mother_df["parent_eartag"].dropna().unique().tolist()
        )
        parent_df_to_concat.append(_expand_mutations_data(mother_mutations_df))
    else:
        missing_parent = "Mother 1"

    # Obtain and format mutation data for father - if present
    if not father_df.empty:
        father_mutations_df = _get_parent_mutations_with_eartags(
            father_df["parent_eartag"].dropna().unique().tolist()
        )
        parent_df_to_concat.append(_expand_mutations_data(father_mutations_df))
    else:
        missing_parent = "Father 1"

    if parent_df_to_concat:
        mutations_df = pd.concat(parent_df_to_concat, ignore_index=True)

        expanded_df = expanded_df.merge(
            mutations_df,
            left_on="parent_eartag",
            right_on="eartag_or_id",
            how="left",
        )

        expanded_df = expanded_df.drop(
            columns=["parent", "eartag_or_id", "animalid_y"]
        )

        expanded_df = expanded_df.rename(columns={"animalid_x": "animalid"})

    return expanded_df, missing_parent


def _dynamic_column_renaming(expanded_df: pd.DataFrame, missing_parent: str):
    """Pivots on multiple values creating a tuple (value, parent_id)
    to account for all parents.

    Iterates through, renaming with parent names and formatting mutation and
    grades accordingly.
    """

    # pull all columns except those used for the index and columns
    # accounts for multiple mutations and grades
    value_columns = [
        col
        for col in expanded_df.columns
        if col not in ["animalid", "parent_id"]
    ]

    # passing multiple values creates a tuple of (values[x], parent_id)
    tuple_columns_df = expanded_df.pivot(
        index="animalid", columns="parent_id", values=value_columns
    )

    new_columns = []
    for col in tuple_columns_df.columns:
        if col[0] == "parent_eartag":
            new_columns.append(col[1])
        else:
            new_columns.append(f"{col[1]}: {col[0]}")

    tuple_columns_df.columns = new_columns
    merged_df = tuple_columns_df.reset_index()

    if missing_parent != "None":
        _add_empty_parent_cols(merged_df, missing_parent)

    return merged_df
