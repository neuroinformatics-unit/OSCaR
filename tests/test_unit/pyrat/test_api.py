import datetime
import json
import logging
import os
import re
from typing import Any

import pandas as pd
import pytest
import responses
from responses import matchers

from oscar_colony.colony_management.pyrat.api import (
    get_pyrat_data,
    get_pyrat_line_mutations,
    get_pyrat_line_name,
    get_pyrat_lines,
)
from tests.pooch_test_data import GIN_REPO, pooch_data_path


@pytest.fixture
def exclude_gin():
    """stop responses library interfering with pooch requests"""
    responses.add_passthru(GIN_REPO.base_url)


@pytest.fixture
def species_response():
    """Default response to give for a request to the /species endpoint"""

    return create_pyrat_response(
        endpoint_name="species",
        json_data=[{"id": 1, "name": "Mouse"}, {"id": 2, "name": "Rat"}],
    )


def create_pyrat_response(
    endpoint_name: str = "animals",
    json_filename: str | None = None,
    json_data: list[Any] | None = None,
    query_params: dict[str, Any] | None = None,
    total_count: int | None = None,
) -> responses.Response:
    """Create a mock response for a request to the pyRAT api.

    If neither json_filename or json_data are provided, then the json data
    will be empty.

    Parameters
    ----------
    endpoint_name: str, optional
        The name of the pyRAT endpoint. Defaults to 'animals'
    json_filename : str | None, optional
        Json filename of response data.
    json_data: list[Any] | None, optional
        Json data for the response. Ignored if json_filename is also provided.
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
            json_response = json.load(f)
    elif json_data is not None:
        json_response = json_data
    else:
        json_response = []

    if total_count is None:
        total_count = len(json_response)

    # x-count is the total number of items returned in this response
    # x-total-count is the total number of available items
    response = responses.Response(
        method="GET",
        url=f"{os.environ['PYRAT_URL']}/api/v3/{endpoint_name}",
        json=json_response,
        headers={
            "x-count": str(len(json_response)),
            "x-total-count": str(total_count),
        },
    )

    if query_params is not None:
        response.match = [
            matchers.query_param_matcher(query_params, strict_match=False)
        ]

    return response


@pytest.mark.parametrize(
    "parent_response, offspring_response, expected_csv_name",
    [
        pytest.param(
            create_pyrat_response(
                json_filename="pyrat-api-single-response-parents.json",
                query_params={"eartag": ["ID-100", "ID-101"]},
            ),
            create_pyrat_response(
                json_filename="pyrat-api-single-response-offspring.json"
            ),
            "pyrat-api-single-response.csv",
            id="Single item returned",
        ),
        pytest.param(
            create_pyrat_response(
                json_filename="pyrat-api-multiple-responses-parents.json",
                query_params={
                    "eartag": [
                        "ID-100",
                        "ID-101",
                        "ID-102",
                        "ID-103",
                        "ID-104",
                        "ID-105",
                        "ID-106",
                        "ID-107",
                    ]
                },
            ),
            create_pyrat_response(
                json_filename="pyrat-api-multiple-responses-offspring.json",
            ),
            "pyrat-api-multiple-responses.csv",
            id="Multiple items returned",
        ),
        pytest.param(
            create_pyrat_response(
                json_filename=(
                    "pyrat-api-multiple-response-and-parents-both-parents.json"
                ),
                query_params={
                    "eartag": [
                        "ID-100",
                        "ID-101",
                        "ID-102",
                        "ID-103",
                        "ID-105",
                    ]
                },
            ),
            create_pyrat_response(
                json_filename=(
                    "pyrat-api-multiple-response-and-parents-offspring.json"
                )
            ),
            "pyrat-api-multiple-response-multiple-parents.csv",
            id="Multiple items with multiple parents returned",
        ),
        pytest.param(
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
            create_pyrat_response(
                json_filename="pyrat-api-no-parents-mutations-offspring.json",
            ),
            "pyrat-api-no-parents-mutations.csv",
            id="No listed mutations or parents",
        ),
    ],
)
@responses.activate
@pytest.mark.usefixtures("exclude_gin")
def test_get_pyrat_data(
    parent_response,
    offspring_response,
    expected_csv_name,
    species_response,
):

    # add mock responses
    responses.add(species_response)
    if parent_response is not None:
        responses.add(parent_response)
    responses.add(offspring_response)

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
@pytest.mark.usefixtures("exclude_gin")
def test_get_pyrat_data_logs(species_response, caplog):
    """Test log statements are recorded correctly when fetching pyRAT data"""

    # add mock responses
    responses.add(species_response)
    responses.add(
        create_pyrat_response(
            json_filename="pyrat-api-single-response-parents.json",
            query_params={"eartag": ["ID-100", "ID-101"]},
        )
    )
    responses.add(
        create_pyrat_response(
            json_filename="pyrat-api-single-response-offspring.json"
        )
    )

    with caplog.at_level(logging.INFO):
        pyrat_dfs = get_pyrat_data(
            species_name="Mouse",
            birth_date_from=datetime.date(2026, 2, 1),
            birth_date_to=datetime.date(2026, 3, 1),
        )
        pyrat_dfs = list(pyrat_dfs)

    expected_messages = [
        "searching PyRAT using custom parameters: {'species': 'Mouse', "
        "'birth date from': '2026-02-01', 'birth date to': '2026-03-01'}",
        "1 animals found in PyRAT database",
        "Converted animal response to dataframe with 1 rows",
    ]
    for message in expected_messages:
        assert message in caplog.text


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
@pytest.mark.usefixtures("exclude_gin")
def test_get_pyrat_lines():
    """
    Test fetching of line names and ids from the pyRAT api.

    This test uses 3 available lines with a max_n_rows of 2 - so there should
    be two requests, and two final dataframes.
    """

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


@responses.activate
@pytest.mark.usefixtures("exclude_gin")
def test_get_pyrat_line_mutations():
    """
    Test fetching of an alphabetical list of line mutations from the
    pyRAT api.
    """

    # create mock line mutations response
    line_id = 111
    line_mutations_response = create_pyrat_response(
        f"strains/{line_id}/mutations",
        json_filename="pyrat-api-line-mutations-response.json",
    )
    responses.add(line_mutations_response)

    line_mutations = get_pyrat_line_mutations(line_id)
    assert line_mutations == ["Mut-A", "Mut-B", "Mut-C"]


@responses.activate
def test_get_pyrat_line_name():
    """Test fetching a line name from a line id."""

    # create mock line name response with one entry
    line_id = 12
    lines_response = create_pyrat_response(
        "strains",
        json_data=[{"id": line_id, "name": "Line-AB"}],
        query_params={
            "k": ["name", "id"],
            "id": 12,
        },
    )
    responses.add(lines_response)

    line_name = get_pyrat_line_name(line_id)
    assert line_name == "Line-AB"


@responses.activate
def test_get_pyrat_line_name_invalid():
    """Test returning multiple line names throws an error."""

    # create mock line name response with two entries
    line_id = 12
    lines_response = create_pyrat_response(
        "strains",
        json_data=[
            {"id": line_id, "name": "Line-AB"},
            {"id": line_id, "name": "Line-BC"},
        ],
        query_params={
            "k": ["name", "id"],
            "id": line_id,
        },
    )
    responses.add(lines_response)

    error_msg = "Multiple lines returned for id: 12"
    with pytest.raises(ValueError, match=error_msg):
        get_pyrat_line_name(line_id)
