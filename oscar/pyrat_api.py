import os

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
    }

    response = requests.get(
        url=f"{os.environ['PYRAT_URL']}/api/v3/animals",
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
