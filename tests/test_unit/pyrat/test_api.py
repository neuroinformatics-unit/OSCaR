import datetime
import os
from typing import Any

import pytest
import responses
from responses import matchers

from oscar.pyrat.api import get_pyrat_data


def create_animal_response(
    json: list[Any], query_params: dict[str, Any] | None = None
) -> responses.Response:
    pyrat_url = os.environ["PYRAT_URL"]

    response = responses.Response(
        method="GET",
        url=f"{pyrat_url}/api/v3/animals",
        json=json,
        headers={"x-count": str(len(json)), "x-total-count": str(len(json))},
    )

    if query_params is not None:
        response.match = [
            matchers.query_param_matcher(query_params, strict_match=False)
        ]

    return response


@pytest.mark.parametrize(
    "line_name, father_response, mother_response, offspring_response",
    [
        pytest.param(
            "Line-A",
            create_animal_response(
                json=[
                    {
                        "animalid": 100,
                        "eartag_or_id": "ID-100",
                        "mutations": [],
                    }
                ],
                query_params={"eartag": "ID-100"},
            ),
            create_animal_response(
                json=[
                    {
                        "animalid": 101,
                        "eartag_or_id": "ID-101",
                        "mutations": [
                            {
                                "animalid": 101,
                                "mutation_id": 5,
                                "mutationname": "Mut-A",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "het",
                            }
                        ],
                    },
                ],
                query_params={"eartag": "ID-101"},
            ),
            create_animal_response(
                json=[
                    {
                        "animalid": 1,
                        "eartag_or_id": "ID-001",
                        "species_name": "Mouse",
                        "sacrifice_reason_name": "Not usable age",
                        "dateborn": "2026-02-03T00:00:00",
                        "strain_name": "Line-A",
                        "mutations": [
                            {
                                "animalid": 1,
                                "mutation_id": 5,
                                "mutationname": "Mut-A",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "het",
                            }
                        ],
                        "parents": [
                            {
                                "animalid": 1,
                                "parent_id": 100,
                                "parent_eartag": "ID-100",
                                "parent_sex": "m",
                                "parent_labid": "B",
                                "parent_cagenumber": "B2",
                            },
                            {
                                "animalid": 1,
                                "parent_id": 101,
                                "parent_eartag": "ID-101",
                                "parent_sex": "f",
                                "parent_labid": "B",
                                "parent_cagenumber": "B2",
                            },
                        ],
                    }
                ]
            ),
            id="Single item returned",
        ),
    ],
)
@responses.activate
def test_get_pyrat_data(
    line_name, father_response, mother_response, offspring_response
):
    # mock species response
    responses.get(
        f"{os.environ['PYRAT_URL']}/api/v3/species",
        json=[{"id": 1, "name": "Mouse"}, {"id": 2, "name": "Rat"}],
    )

    # add mock animal responses
    for response in [father_response, mother_response, offspring_response]:
        responses.add(response)

    pyrat_dfs = get_pyrat_data(
        line_name=line_name,
        species_name="Mouse",
        birth_date_from=datetime.date(2026, 2, 1),
        birth_date_to=datetime.date(2026, 2, 6),
    )
    pyrat_dfs = list(pyrat_dfs)

    assert len(pyrat_dfs) == 1
