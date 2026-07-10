import datetime
import logging
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
    log_params = {}
    if line_name is not None:
        params["strain_name_with_id_like"] = line_name
        log_params["strain name with id like"] = line_name

    if species_name is not None:
        params["species"] = _get_species_id(species_name)
        log_params["species"] = species_name

    if birth_date_from is not None:
        params["birth_date_from"] = birth_date_from.isoformat()
        log_params["birth date from"] = birth_date_from.isoformat()

    if birth_date_to is not None:
        params["birth_date_to"] = birth_date_to.isoformat()
        log_params["birth date to"] = birth_date_to.isoformat()

    # Make one request to determine how many results there are
    logging.info(f"searching PyRAT using custom parameters: {log_params}")
    animals_response = _make_pyrat_request("animals", params)
    yield _convert_animals_to_df(animals_response.json())
    headers = animals_response.headers
    total_n = int(headers["x-total-count"])
    logging.info(f"{total_n} animals found in PyRAT database")

    # If more results than max_n_rows, keep making requests and yielding result
    for start_n in range(max_n_rows, total_n, max_n_rows):
        params["o"] = start_n
        animals_response = _make_pyrat_request("animals", params)
        yield _convert_animals_to_df(animals_response.json())


def get_pyrat_lines(max_n_rows: int = 10000) -> Iterator[pd.DataFrame]:
    """Fetch available lines directly from the pyRAT api.

    To handle the potentially large number of lines returned from pyRAT,
    this function returns a generator of pandas dataframes (each with no
    more than max_n_rows).

    This expects PYRAT_URL, PYRAT_CLIENT_TOKEN and PYRAT_USER_TOKEN to
    be set as environment variables.

    Parameters
    ----------
    max_n_rows : int, optional
        Maximum number of lines in each returned dataframe (and therefore
        returned per request to the pyRAT api)

    Returns
    -------
    Iterator[pd.DataFrame]
        Generator of dataframes of returned line data, each with columns:
        name and id. Names are in alphabetical order. Id is useful for
        fetching line mutations via get_pyrat_line_mutations.
        If no data is available for the query, the dataframe will be empty.
    """

    params = {
        "k": ["name", "id"],
        "s": ["name:asc"],
        "status": ["available"],
        "l": max_n_rows,
        "o": 0,
    }

    # Make one request to determine how many results there are
    lines_response = _make_pyrat_request("strains", params)
    yield pd.DataFrame(lines_response.json())
    headers = lines_response.headers
    total_n = int(headers["x-total-count"])

    # If more results than max_n_rows, keep making requests and yielding result
    for start_n in range(max_n_rows, total_n, max_n_rows):
        params["o"] = start_n
        lines_response = _make_pyrat_request("strains", params)
        yield pd.DataFrame(lines_response.json())


def get_pyrat_line_mutations(line_id: int) -> list[str]:
    """Get mutation names for the given line id.

    This expects PYRAT_URL, PYRAT_CLIENT_TOKEN and PYRAT_USER_TOKEN to
    be set as environment variables.

    Parameters
    ----------
    line_id : int
        Id of the line (e.g. as returned from get_pyrat_lines)

    Returns
    -------
    list[str]
        List of mutation names in alphabetical order.
    """

    # We use the line id here (rather than the line name) as / characters in
    # line names were causing 404 responses - even when escaped.
    mutations_response = _make_pyrat_request(f"strains/{line_id}/mutations")
    mutations_list = [
        mutations["name"] for mutations in mutations_response.json()
    ]

    return sorted(mutations_list)


def _make_pyrat_request(
    endpoint_name: str, params: dict[str, Any] | None = None
) -> requests.Response:
    """Make request to the pyRAT api.

    This expects PYRAT_URL, PYRAT_CLIENT_TOKEN and PYRAT_USER_TOKEN to
    be set as environment variables.

    Parameters
    ----------
    endpoint_name : str
        Name of endpoint e.g. 'species'
    params : dict[str, Any] | None, optional
        Extra parameters to pass to the endpoint

    Returns
    -------
    requests.Response
        The requests response object, containing data from pyRAT
    """
    logging.info(f"endpoint_name : {endpoint_name}")
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
    response.raise_for_status
    logging.debug(response)

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


def _get_parent_mutations_with_eartags(
    eartags: list[str], batch_size: int = 400
) -> pd.DataFrame:
    """Get parent mutation information for the given animal eartags.

    Since eartags are appended to the URL, it can exceed the maximum request
    size. It is processed in batches, to remain below the threshold.

    Parameters
    ----------
    eartags : list[str]
        all the unique parent eartags
    batch_size : int, optional
        The number of eartags to process per request to the pyRAT api, by
        default 400. To prevent exceeding maximum characters.

    Returns
    -------
    pd.DataFrame
        df containing the parent eartag along with their assigned mutations.
    """

    all_mutation_data = []

    for start in range(0, len(eartags), batch_size):
        eartag_batch = eartags[start : start + batch_size]
        params = {
            "k": ["animalid", "eartag_or_id", "mutations"],
            "s": ["eartag_or_id:asc"],
            "state": ["live", "sacrificed", "exported"],
            "eartag": eartag_batch,
            "l": len(eartag_batch),
        }
        batch_data = _make_pyrat_request("animals", params).json()
        all_mutation_data.extend(batch_data)

    if len(all_mutation_data) != len(eartags):
        raise ValueError(
            f"{len(all_mutation_data)} animals returned for "
            f"{len(eartags)} eartags: {eartags}"
        )

    return pd.DataFrame(all_mutation_data)


def _convert_animals_to_df(animals_data: list[dict[str, Any]]) -> pd.DataFrame:
    """Convert animal data fetched from the pyRAT api to a pandas DataFrame.

    The structure / column names are matched to that exported from the pyRAT
    UI.
    """

    animals_df = pd.DataFrame(animals_data)
    if animals_df.empty:
        logging.info("no animals collected for this search populated")
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
    logging.debug(animals_df)
    return animals_df


def _expand_mutations_data(selected_df: pd.DataFrame) -> pd.DataFrame:
    """Expand a mutations column into a full dataframe.

    Each row of a mutations column contains a list of dictionaries
    (one per mutation for the animal). This function expands these into
    their own columns labelled Mutation 1, 2... and Grade 1, 2...

    Parameters
    ----------
    selected_df : pd.DataFrame
        DataFrame of pyRAT data with raw mutations column

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
        logging.info(
            f"no mutation(s) found for these animal_ids : "
            f"{selected_df['animalid'].explode().unique().tolist()}"
        )
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
    logging.debug(merged_df)

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

    # If no parents are listed for ANY animals, return empty mother / father
    # columns, with empty mutation / grade
    if parents_df.empty:
        logging.info(
            f"no parent(s) found for these animal_ids : "
            f"{animals_df['animalid'].explode().unique().tolist()}"
        )
        animals_df = animals_df.loc[:, ["animalid"]]
        _add_empty_parent_cols(animals_df, "Mother 1")
        _add_empty_parent_cols(animals_df, "Father 1")
        return animals_df

    # Rename column names as father or mother
    parents_df = parents_df[["animalid", "parent_eartag", "parent_sex"]]
    parents_df = parents_df.rename(columns={"parent_sex": "parent"})
    parents_df.loc[parents_df["parent"] == "m", "parent"] = "Father"
    parents_df.loc[parents_df["parent"] == "f", "parent"] = "Mother"

    # Number each consecutive parent appearance and append num to column name
    parents_df["parent_id"] = (
        parents_df.groupby(["animalid", "parent"]).cumcount() + 1
    )
    parents_df["parent_id"] = (
        parents_df["parent"] + " " + parents_df["parent_id"].astype(str)
    )

    parents_df_with_mutations = _merge_parent_mutations(parents_df)
    clean_parents_df = _parent_column_renaming(parents_df_with_mutations)

    return clean_parents_df


def _merge_parent_mutations(parents_df: pd.DataFrame) -> pd.DataFrame:
    """
    Fetch parent mutations from the pyRAT api using their eartag, and creates
    numbered 'Mutation' and 'Grade' columns for each parent.


    Parameters
    ----------
    parents_df : pd.DataFrame
        dataframe containing animalid, parent_eartag, parent and parent_id

    Returns
    -------
    pd.DataFrame
        parent_df with corresponding mutation and grade appended.
    """

    mutations_df = _get_parent_mutations_with_eartags(
        parents_df["parent_eartag"].dropna().unique().tolist()
    )

    mutations_df = _expand_mutations_data(mutations_df)
    mutations_df = mutations_df.drop(columns=["animalid"])

    parents_df = parents_df.merge(
        mutations_df,
        left_on="parent_eartag",
        right_on="eartag_or_id",
        how="left",
    )

    parents_df = parents_df.drop(columns=["eartag_or_id"])

    return parents_df


def _parent_column_renaming(expanded_df: pd.DataFrame):
    """
    Create columns for each unique parent_id.

    This function removes the parent column in favour of parent_id, then it
    collapses all rows with the same animalid into one row. Each unique
    parent_id is given its own column, and Mutation / Grade columns are
    re-named to include the relevant parent_id as a prefix.
    """

    # pivoting multiple values creates column names which are a tuple of
    # (old_column_name, parent_id)
    expanded_df = expanded_df.drop(columns=["parent"])
    tuple_columns_df = expanded_df.pivot(index="animalid", columns="parent_id")

    new_col_names = []
    for col_name, parent_id in tuple_columns_df.columns:
        if col_name == "parent_eartag":
            new_col_names.append(parent_id)
        else:
            new_col_names.append(f"{parent_id}: {col_name}")

    tuple_columns_df.columns = new_col_names
    merged_df = tuple_columns_df.reset_index()

    for parent in ["Mother 1", "Father 1"]:
        if f"{parent}: Mutation 1" not in merged_df.columns:
            logging.info(
                f"{parent} not recorded for animalid(s) : "
                f"{merged_df['animalid'].unique().tolist()}"
            )
            _add_empty_parent_cols(merged_df, parent)

    return merged_df
