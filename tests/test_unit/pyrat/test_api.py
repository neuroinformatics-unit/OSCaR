import datetime
import json
import os
import re
from typing import Any

import pandas as pd
import pytest
import responses
from responses import matchers

from oscar_colony.colony_management.pyrat.api import (
    get_pyrat_data,
    get_pyrat_lines,
)
from tests.pooch_test_data import GIN_REPO, pooch_data_path


@pytest.fixture
def species_response():
    """Default response to give for a request to the /species endpoint"""

    return responses.Response(
        method="GET",
        url=f"{os.environ['PYRAT_URL']}/api/v3/species",
        json=[{"id": 1, "name": "Mouse"}, {"id": 2, "name": "Rat"}],
    )


def create_pyrat_response(
    endpoint_name: str = "animals",
    json_filename: str | None = None,
    query_params: dict[str, Any] | None = None,
    total_count: int | None = None,
) -> responses.Response:
    """Create a mock response for a request to the pyRAT api.

    Parameters
    ----------
    endpoint_name: str, optional
        The name of the pyRAT endpoint. Defaults to 'animals'
    json_filename : str | None, optional
        Json filename of response data. If not given, the json data will
        be empty.
    query_params : dict[str, Any] | None, optional
        The required query parameters to return this response. Note: matching
        isn't strict, so while the request must contain the given parameters,
        it may also have additional ones that aren't listed here.
    total_count: int | None, optional
        The total count to provide in the response header. If None,
        will be set equal to the json length.

    Returns
    -------
    responses.Response
        The mock responses object
    """
    if json_filename is not None:
        with open(pooch_data_path(json_filename)) as f:
            json_data = json.load(f)
    else:
        json_data = []

    if total_count is None:
        total_count = len(json_data)

    # x-count is the total number of items returned in this response
    # x-total-count is the total number of available items
    response = responses.Response(
        method="GET",
        url=f"{os.environ['PYRAT_URL']}/api/v3/{endpoint_name}",
        json=json_data,
        headers={
            "x-count": str(len(json_data)),
            "x-total-count": str(total_count),
        },
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
            create_pyrat_response(
                json_filename="pyrat-api-single-response-father.json",
                query_params={"eartag": "ID-100"},
            ),
            create_pyrat_response(
                json_filename="pyrat-api-single-response-mother.json",
                query_params={"eartag": "ID-101"},
            ),
            create_pyrat_response(
                json_filename="pyrat-api-single-response-offspring.json"
            ),
            "pyrat-api-single-response.csv",
            id="Single item returned",
        ),
        pytest.param(
            create_pyrat_response(
                json_filename="pyrat-api-multiple-responses-father.json",
                query_params={
                    "eartag": ["ID-100", "ID-102", "ID-104", "ID-106"]
                },
            ),
            create_pyrat_response(
                json_filename="pyrat-api-multiple-responses-mother.json",
                query_params={
                    "eartag": ["ID-101", "ID-103", "ID-105", "ID-107"]
                },
            ),
            create_pyrat_response(
                json_filename="pyrat-api-multiple-responses-offspring.json",
            ),
            "pyrat-api-multiple-responses.csv",
            id="Multiple items returned",
        ),
        pytest.param(
            None,
            create_pyrat_response(
                json_filename="pyrat-api-single-parent-mother.json",
                query_params={"eartag": "ID-101"},
            ),
            create_pyrat_response(
                json_filename="pyrat-api-single-parent-offspring.json",
            ),
            "pyrat-api-single-parent.csv",
            id="All missing one parent",
        ),
        pytest.param(
            None,
            None,
            create_pyrat_response(
                json_filename="pyrat-api-no-parents-mutations-offspring.json",
            ),
            "pyrat-api-no-parents-mutations.csv",
            id="No listed mutations or parents",
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
        if response is not None:
            responses.add(response)

    pyrat_dfs = get_pyrat_data(
        species_name="Mouse",
        birth_date_from=datetime.date(2026, 2, 1),
        birth_date_to=datetime.date(2026, 3, 1),
    )
    pyrat_dfs = list(pyrat_dfs)

    expected_csv = pd.read_csv(pooch_data_path(expected_csv_name), dtype=str)
    assert len(pyrat_dfs) == 1
    # use check_like=True to ignore column order
    pd.testing.assert_frame_equal(pyrat_dfs[0], expected_csv, check_like=True)


@responses.activate
def test_no_pyrat_data_exists(species_response):
    """Test fetching pyrat api data when no animals are available"""

    # add mock responses
    responses.add(species_response)
    responses.add(
        create_pyrat_response()  # All animal responses empty []
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


def test_invalid_date_range():
    """
    Test fetching pyrat api data when birth_date_from is after
    birth_date_to
    """

    error_msg = "birth_date_to must be after birth_date_from"
    with pytest.raises(ValueError, match=error_msg):
        pyrat_dfs = get_pyrat_data(
            line_name="Line-A",
            species_name="Mouse",
            birth_date_from=datetime.date(2026, 2, 6),
            birth_date_to=datetime.date(2026, 2, 1),
        )
        pyrat_dfs = list(pyrat_dfs)


@responses.activate
def test_invalid_species_name(species_response):
    """
    Test fetching pyrat api data when given species name is invalid
    """

    # add mock species response
    responses.add(species_response)

    error_msg = re.escape(
        "No ID found for species Fish: available values are ['Mouse', 'Rat']"
    )
    with pytest.raises(ValueError, match=error_msg):
        pyrat_dfs = get_pyrat_data(
            line_name="Line-A",
            species_name="Fish",
            birth_date_from=datetime.date(2026, 2, 1),
            birth_date_to=datetime.date(2026, 2, 6),
        )
        pyrat_dfs = list(pyrat_dfs)


@responses.activate
def test_get_pyrat_lines():
    """
    Test fetching of line names and ids from the pyRAT api.

    This test uses 3 available lines with a max_n_rows of 2 - so there should
    be two requests, and two final dataframes.
    """
    # stop responses library interfering with pooch requests
    responses.add_passthru(GIN_REPO.base_url)

    # create mock line responses. We have two here, as there's a total of
    # 3 lines available, but we are limiting max_n_rows to 2 per request.
    max_n_rows = 2
    offset = 0
    for json_filename in [
        "pyrat-api-lines-response-1.json",
        "pyrat-api-lines-response-2.json",
    ]:
        lines_response = create_pyrat_response(
            "strains",
            json_filename=json_filename,
            query_params={
                "k": ["name", "id"],
                "s": "name:asc",
                "status": "available",
                "l": max_n_rows,
                "o": offset,
            },
            total_count=3,
        )
        offset += max_n_rows
        responses.add(lines_response)

    pyrat_dfs = get_pyrat_lines(max_n_rows=max_n_rows)
    pyrat_dfs = list(pyrat_dfs)

    assert len(pyrat_dfs) == 2
    assert len(pyrat_dfs[0]) == 2
    assert len(pyrat_dfs[1]) == 1

    expected_csv = pd.read_csv(pooch_data_path("pyrat-api-lines.csv"))
    pd.testing.assert_frame_equal(
        pd.concat(pyrat_dfs).reset_index(drop=True), expected_csv
    )
