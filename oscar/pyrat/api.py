import datetime
import os
from typing import Any

import pandas as pd
import requests


def get_pyrat_data(
    line_name: str | None = None,
    species_name: str | None = None,
    birth_date_from: datetime.date | None = None,
    birth_date_to: datetime.date | None = None,
) -> pd.DataFrame:
    """Fetch animal data directly from the pyRAT api.

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

    Returns
    -------
    pd.DataFrame
        Dataframe of returned animal data, in format matching that
        exported via the pyRAT UI.
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
        "l": 5,
    }

    if line_name is not None:
        params["strain_name_with_id_like"] = line_name

    if species_name is not None:
        params["species"] = _get_species_id(species_name)

    if birth_date_from is not None:
        params["birth_date_from"] = birth_date_from.isoformat()

    if birth_date_to is not None:
        params["birth_date_to"] = birth_date_to.isoformat()

    animals_data = _make_pyrat_request("animals", params)

    # Convert to a pandas dataframe to make processing easier
    animals_df = _convert_animals_to_df(animals_data)

    return animals_df


def _make_pyrat_request(endpoint_name: str, params: dict[str, Any]) -> Any:
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
    Any
        Response decoded from json into a python object
    """
    response = requests.get(
        url=f"{os.environ['PYRAT_URL']}/api/v3/{endpoint_name}",
        auth=(
            os.environ["PYRAT_CLIENT_TOKEN"],
            os.environ["PYRAT_USER_TOKEN"],
        ),
        params=params,
        timeout=2,  # number of seconds before timeout
    )

    # If the request didn't succeed, raise an error containing the status
    # code
    response.raise_for_status()

    return response.json()


def _get_species_id(species_name: str) -> int:
    """Get pyRAT database ID for named species"""

    params = {
        "k": ["id", "name"],
        "s": ["name:asc"],
    }
    species_ids = _make_pyrat_request("species", params)

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


def _get_mutations_for_eartags(eartags: list[str]) -> pd.DataFrame:
    """Get mutation information for the given animal eartags"""

    params = {
        "k": ["animalid", "eartag_or_id", "mutations"],
        "s": ["eartag_or_id:asc"],
        "state": ["live", "sacrificed", "exported"],
        "eartag": eartags,
    }
    mutation_data = _make_pyrat_request("animals", params)

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
        }
    )

    return animals_df


def _expand_mutations_data(
    animals_df: pd.DataFrame, column_prefix: str = ""
) -> pd.DataFrame:
    """Expand a mutations column into a full dataframe.

    Each row of a mutations column contains a list of dictionaries
    (one per mutation for the animal). This function expands these into
    their own columns labelled Mutation 1, 2... and Grade 1, 2...

    Parameters
    ----------
    animals_df : pd.DataFrame
        DataFrame of animals data, with raw mutations column
    column_prefix: str
        Prefix to add to expanded column names e.g. a prefix of 'Father: '
        would result in columns Father: Mutation 1, Father: Grade 1 etc.

    Returns
    -------
    pd.DataFrame
        Dataframe with separate Mutation and Grade columns
    """

    exploded_mutations_col = animals_df.mutations.explode()
    mutations_df = pd.DataFrame(
        exploded_mutations_col[~exploded_mutations_col.isna()].tolist()
    )

    mutation_col_name = f"{column_prefix}Mutation"
    grade_col_name = f"{column_prefix}Grade"

    # If no mutations are listed for any animals, return an empty Mutation 1 /
    # Grade 1 column
    if mutations_df.empty:
        animals_df = animals_df.drop(["mutations"], axis=1)
        animals_df[f"{mutation_col_name} 1"] = pd.Series(dtype=str)
        animals_df[f"{grade_col_name} 1"] = pd.Series(dtype=str)
        return animals_df

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
    merged_df = animals_df.drop(["mutations"], axis=1)
    merged_df = merged_df.merge(pivoted_mutations, on="animalid", how="left")

    return merged_df


def _expand_parents_data(animals_df: pd.DataFrame) -> pd.DataFrame:
    """Expand column containing multiple parents' information into separate
    columns.

    This adds columns for Mother / Father ID, as well as their respective
    mutations.
    """
    parents_df = pd.DataFrame(animals_df.parents.explode().tolist())

    # Create dataframe with one row per animalid, and one column each for
    # ID of mother and father
    expanded_df = parents_df[["animalid", "parent_eartag", "parent_sex"]]
    expanded_df = expanded_df.pivot(
        columns="parent_sex", values="parent_eartag", index="animalid"
    )
    expanded_df = expanded_df.reset_index().rename_axis(None, axis=1)
    expanded_df = expanded_df.rename(columns={"f": "Mother", "m": "Father"})

    # Fetch mutation info for all parents and merge
    for parent in ["Mother", "Father"]:
        parent_df = _get_mutations_for_parent(expanded_df, parent)
        expanded_df = expanded_df.merge(parent_df, on=parent, how="left")

    return expanded_df


def _get_mutations_for_parent(
    parents_df: pd.DataFrame, parent: str
) -> pd.DataFrame:
    """Return a dataframe with mutations for all unique parent IDs.

    Parameters
    ----------
    parents_df : pd.DataFrame
        Dataframe with animalid, Mother and Father columns
    parent : str
        Parent to return mutations for: "Mother" or "Father"

    Returns
    -------
    pd.DataFrame
        Dataframe with parent IDs and mutation / grade columns
    """

    mutations_df = _get_mutations_for_eartags(
        parents_df[parent].unique().tolist()
    )
    mutations_df = _expand_mutations_data(
        mutations_df, column_prefix=f"{parent}: "
    )
    mutations_df = mutations_df.drop(["animalid"], axis=1)
    mutations_df = mutations_df.rename(columns={"eartag_or_id": parent})

    return mutations_df
