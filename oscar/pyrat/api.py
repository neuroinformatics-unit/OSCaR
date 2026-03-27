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
):
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


def _convert_animals_to_df(animals_data: list[dict[str, Any]]) -> pd.DataFrame:
    animals_df = pd.DataFrame(animals_data)

    mutations_df = _expand_mutations_data(animals_df)

    parents_df = pd.DataFrame(animals_df.parents.explode().tolist())

    mother_df = parents_df.loc[
        parents_df.parent_sex == "f", ["animalid", "parent_eartag"]
    ]
    mother_df = mother_df.rename(columns={"parent_eartag": "Mother"})

    # Deal with potential of multiple fathers / mothers

    father_df = parents_df.loc[
        parents_df.parent_sex == "m", ["animalid", "parent_eartag"]
    ]
    father_df = father_df.rename(columns={"parent_eartag": "Father"})

    animals_df = animals_df.drop(["mutations", "parents"], axis=1)
    animals_df = animals_df.merge(mother_df, on="animalid", how="left")
    animals_df = animals_df.merge(father_df, on="animalid", how="left")

    animals_df = animals_df.merge(mutations_df, on="animalid", how="left")
    animals_df = animals_df.drop(["Mutation", "Grade", "count"], axis=1)

    return animals_df


def _expand_mutations_data(animals_df: pd.DataFrame) -> pd.DataFrame:
    """Expand the mutations column into a full dataframe.

    Each row of the mutations column contains a list of dictionaries
    (one per mutation for the animal). This function expands these into
    their own columns labelled Mutation 1, 2... and Grade 1, 2...

    Parameters
    ----------
    animals_df : pd.DataFrame
        DataFrame of animals data, with raw mutations column

    Returns
    -------
    pd.DataFrame
        Dataframe with separate Mutation and Grade columns
    """

    mutations_df = pd.DataFrame(animals_df.mutations.explode().tolist())
    mutations_df = mutations_df[["animalid", "mutationname", "mutationgrade"]]
    mutations_df = mutations_df.rename(
        columns={"mutationname": "Mutation", "mutationgrade": "Grade"}
    )

    # Adds a counter for the number of mutation rows per animal id
    mutations_df["count"] = (
        mutations_df.groupby("animalid").cumcount() + 1
    ).astype("string")

    # Create one row per animalid, with separate columns for
    # Mutation 1 / Grade 1, Mutation 2 / Grade 2 ...
    pivoted_mutations = mutations_df.pivot(
        columns="count", index="animalid", values=["Mutation", "Grade"]
    ).reset_index()
    pivoted_mutations.columns = [
        " ".join(column_names).strip()
        for column_names in pivoted_mutations.columns.to_flat_index()
    ]

    return pivoted_mutations
