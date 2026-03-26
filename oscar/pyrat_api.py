import os
from typing import Any

import requests


def get_pyrat_data(
    line_name: str | None = None,
    species_name: str | None = None,
    start_date: None = None,
    end_date: None = None,
):
    payload = {
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
        payload["strain_name_with_id_like"] = line_name

    if species_name is not None:
        payload["species"] = _get_species_id(species_name)

    animals = _make_pyrat_request("animals", payload)

    return animals


def _make_pyrat_request(endpoint_name: str, payload: dict[str, Any]) -> Any:
    response = requests.get(
        url=f"{os.environ['PYRAT_URL']}/api/v3/{endpoint_name}",
        auth=(
            os.environ["PYRAT_CLIENT_TOKEN"],
            os.environ["PYRAT_USER_TOKEN"],
        ),
        params=payload,
        timeout=2,  # number of seconds before timeout
    )

    # If the request didn't succeed, raise an error containing the status
    # code
    response.raise_for_status()

    return response.json()


def _get_species_id(species_name: str) -> int:
    """Get pyRAT database ID for named species"""

    payload = {
        "k": ["id", "name"],
        "s": ["name:asc"],
    }
    species_ids = _make_pyrat_request("species", payload)

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
