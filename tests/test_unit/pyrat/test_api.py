import datetime
import os

import responses
from responses import matchers

from oscar.pyrat.api import get_pyrat_data


@responses.activate
def test_get_pyrat_data():
    pyrat_url = os.environ["PYRAT_URL"]
    line_name = "Line-A"

    # mock species response
    responses.get(
        f"{pyrat_url}/api/v3/species",
        json=[{"id": 1, "name": "Mouse"}, {"id": 2, "name": "Rat"}],
    )

    # mock mother response
    params = {"eartag": "ID-100"}
    responses.get(
        f"{pyrat_url}/api/v3/animals",
        match=[matchers.query_param_matcher(params, strict_match=False)],
        json=[
            {
                "animalid": 100,
                "eartag_or_id": "ID-100",
                "mutations": [],
            }
        ],
    )

    # mock father response
    params = {"eartag": "ID-101"}
    responses.get(
        f"{pyrat_url}/api/v3/animals",
        match=[matchers.query_param_matcher(params, strict_match=False)],
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
    )

    # mock offspring response
    responses.get(
        f"{pyrat_url}/api/v3/animals",
        json=[
            {
                "animalid": 1,
                "eartag_or_id": "ID-001",
                "species_name": "Mouse",
                "sacrificice_reason_name": "Not usable age",
                "dateborn": "2026-02-03T00:00:00",
                "strain_name": line_name,
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
        ],
        headers={"x-count": "1", "x-total-count": "1"},
    )

    pyrat_dfs = get_pyrat_data(
        line_name=line_name,
        species_name="Mouse",
        birth_date_from=datetime.date(2026, 2, 1),
        birth_date_to=datetime.date(2026, 2, 6),
    )
    pyrat_dfs = list(pyrat_dfs)

    assert len(pyrat_dfs) == 1
