import datetime
import os
from typing import Any

import pandas as pd
import pytest
import responses
from responses import matchers

from oscar.pyrat.api import get_pyrat_data
from tests.pooch_test_data import GIN_REPO, pooch_data_path


@pytest.fixture
def species_response():
    return responses.Response(
        method="GET",
        url=f"{os.environ['PYRAT_URL']}/api/v3/species",
        json=[{"id": 1, "name": "Mouse"}, {"id": 2, "name": "Rat"}],
    )


def create_animal_response(
    json: list[Any], query_params: dict[str, Any] | None = None
) -> responses.Response:
    response = responses.Response(
        method="GET",
        url=f"{os.environ['PYRAT_URL']}/api/v3/animals",
        json=json,
        headers={"x-count": str(len(json)), "x-total-count": str(len(json))},
    )

    if query_params is not None:
        response.match = [
            matchers.query_param_matcher(query_params, strict_match=False)
        ]

    return response


@pytest.mark.parametrize(
    "father_response, mother_response, offspring_response, expected_csv_name",
    [
        pytest.param(
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
            "pyrat-api-single-response.csv",
            id="Single item returned",
        ),
        pytest.param(
            create_animal_response(
                json=[
                    {
                        "animalid": 100,
                        "eartag_or_id": "ID-100",
                        "mutations": [
                            {
                                "animalid": 100,
                                "mutation_id": 5,
                                "mutationname": "Mut-A",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "wt",
                            }
                        ],
                    },
                    {
                        "animalid": 102,
                        "eartag_or_id": "ID-102",
                        "mutations": [
                            {
                                "animalid": 102,
                                "mutation_id": 5,
                                "mutationname": "Mut-A",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "wt",
                            },
                            {
                                "animalid": 102,
                                "mutation_id": 5,
                                "mutationname": "Mut-B",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "hom",
                            },
                        ],
                    },
                    {
                        "animalid": 104,
                        "eartag_or_id": "ID-104",
                        "mutations": [
                            {
                                "animalid": 104,
                                "mutation_id": 5,
                                "mutationname": "Mut-B",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "hom",
                            }
                        ],
                    },
                    {
                        "animalid": 106,
                        "eartag_or_id": "ID-106",
                        "mutations": [],
                    },
                ],
                query_params={
                    "eartag": ["ID-100", "ID-102", "ID-104", "ID-106"]
                },
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
                            },
                        ],
                    },
                    {
                        "animalid": 103,
                        "eartag_or_id": "ID-103",
                        "mutations": [
                            {
                                "animalid": 103,
                                "mutation_id": 5,
                                "mutationname": "Mut-A",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "het",
                            },
                            {
                                "animalid": 103,
                                "mutation_id": 5,
                                "mutationname": "Mut-B",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "wt",
                            },
                        ],
                    },
                    {
                        "animalid": 105,
                        "eartag_or_id": "ID-105",
                        "mutations": [
                            {
                                "animalid": 105,
                                "mutation_id": 5,
                                "mutationname": "Mut-A",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "het",
                            },
                        ],
                    },
                    {
                        "animalid": 107,
                        "eartag_or_id": "ID-107",
                        "mutations": [
                            {
                                "animalid": 107,
                                "mutation_id": 5,
                                "mutationname": "Mut-A",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "het",
                            },
                            {
                                "animalid": 107,
                                "mutation_id": 5,
                                "mutationname": "Mut-B",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "wt",
                            },
                        ],
                    },
                ],
                query_params={
                    "eartag": ["ID-101", "ID-103", "ID-105", "ID-107"]
                },
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
                    },
                    {
                        "animalid": 2,
                        "eartag_or_id": "ID-002",
                        "species_name": "Mouse",
                        "sacrifice_reason_name": "Not usable age",
                        "dateborn": "2026-02-04T00:00:00",
                        "strain_name": "Line-A",
                        "mutations": [
                            {
                                "animalid": 2,
                                "mutation_id": 5,
                                "mutationname": "Mut-A",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "hom",
                            }
                        ],
                        "parents": [
                            {
                                "animalid": 2,
                                "parent_id": 100,
                                "parent_eartag": "ID-100",
                                "parent_sex": "m",
                                "parent_labid": "B",
                                "parent_cagenumber": "B2",
                            },
                            {
                                "animalid": 2,
                                "parent_id": 101,
                                "parent_eartag": "ID-101",
                                "parent_sex": "f",
                                "parent_labid": "B",
                                "parent_cagenumber": "B2",
                            },
                        ],
                    },
                    {
                        "animalid": 3,
                        "eartag_or_id": "ID-003",
                        "species_name": "Mouse",
                        "sacrifice_reason_name": "Not usable age",
                        "dateborn": "2026-02-05T00:00:00",
                        "strain_name": "Line-AB",
                        "mutations": [
                            {
                                "animalid": 3,
                                "mutation_id": 5,
                                "mutationname": "Mut-A",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "het",
                            },
                            {
                                "animalid": 3,
                                "mutation_id": 5,
                                "mutationname": "Mut-B",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "hom",
                            },
                        ],
                        "parents": [
                            {
                                "animalid": 3,
                                "parent_id": 102,
                                "parent_eartag": "ID-102",
                                "parent_sex": "m",
                                "parent_labid": "B",
                                "parent_cagenumber": "B2",
                            },
                            {
                                "animalid": 3,
                                "parent_id": 103,
                                "parent_eartag": "ID-103",
                                "parent_sex": "f",
                                "parent_labid": "B",
                                "parent_cagenumber": "B2",
                            },
                        ],
                    },
                    {
                        "animalid": 4,
                        "eartag_or_id": "ID-004",
                        "species_name": "Mouse",
                        "sacrifice_reason_name": "Not usable age",
                        "dateborn": "2026-02-06T00:00:00",
                        "strain_name": "Line-AB",
                        "mutations": [
                            {
                                "animalid": 4,
                                "mutation_id": 5,
                                "mutationname": "Mut-A",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "het",
                            },
                            {
                                "animalid": 4,
                                "mutation_id": 5,
                                "mutationname": "Mut-B",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "het",
                            },
                        ],
                        "parents": [
                            {
                                "animalid": 4,
                                "parent_id": 104,
                                "parent_eartag": "ID-104",
                                "parent_sex": "m",
                                "parent_labid": "B",
                                "parent_cagenumber": "B2",
                            },
                            {
                                "animalid": 4,
                                "parent_id": 105,
                                "parent_eartag": "ID-105",
                                "parent_sex": "f",
                                "parent_labid": "B",
                                "parent_cagenumber": "B2",
                            },
                        ],
                    },
                    {
                        "animalid": 5,
                        "eartag_or_id": "ID-005",
                        "species_name": "Mouse",
                        "sacrifice_reason_name": "Not usable age",
                        "dateborn": "2026-02-07T00:00:00",
                        "strain_name": "Line-AB",
                        "mutations": [
                            {
                                "animalid": 5,
                                "mutation_id": 5,
                                "mutationname": "Mut-A",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "het",
                            },
                            {
                                "animalid": 5,
                                "mutation_id": 5,
                                "mutationname": "Mut-B",
                                "mutationtype": "BAC transgenic",
                                "grade_id": 5,
                                "mutationgrade": "het",
                            },
                        ],
                        "parents": [
                            {
                                "animalid": 5,
                                "parent_id": 106,
                                "parent_eartag": "ID-106",
                                "parent_sex": "m",
                                "parent_labid": "B",
                                "parent_cagenumber": "B2",
                            },
                            {
                                "animalid": 5,
                                "parent_id": 107,
                                "parent_eartag": "ID-107",
                                "parent_sex": "f",
                                "parent_labid": "B",
                                "parent_cagenumber": "B2",
                            },
                        ],
                    },
                ]
            ),
            "pyrat-api-multiple-responses.csv",
            id="Multiple items returned",
        ),
    ],
)
@responses.activate
def test_get_pyrat_data(
    father_response,
    mother_response,
    offspring_response,
    expected_csv_name,
    species_response,
):
    # stop responses library interfering with pooch requests
    responses.add_passthru(GIN_REPO.base_url)

    # add mock responses
    responses.add(species_response)
    for response in [father_response, mother_response, offspring_response]:
        responses.add(response)

    pyrat_dfs = get_pyrat_data(
        line_name="Line-A",
        species_name="Mouse",
        birth_date_from=datetime.date(2026, 2, 1),
        birth_date_to=datetime.date(2026, 2, 6),
    )
    pyrat_dfs = list(pyrat_dfs)

    expected_csv = pd.read_csv(pooch_data_path(expected_csv_name), dtype=str)
    assert len(pyrat_dfs) == 1
    pd.testing.assert_frame_equal(pyrat_dfs[0], expected_csv)


@responses.activate
def test_no_pyrat_data_exists(species_response):
    """Test fetching pyrat api data when no animals are available"""

    # stop responses library interfering with pooch requests
    responses.add_passthru(GIN_REPO.base_url)

    # add mock responses
    responses.add(species_response)
    responses.add(
        create_animal_response(json=[])  # All animal responses empty
    )

    pyrat_dfs = get_pyrat_data(
        line_name="Line-A",
        species_name="Mouse",
        birth_date_from=datetime.date(2026, 2, 1),
        birth_date_to=datetime.date(2026, 2, 6),
    )
    pyrat_dfs = list(pyrat_dfs)

    assert len(pyrat_dfs) == 1
    assert pyrat_dfs[0].empty
