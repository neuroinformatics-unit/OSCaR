import datetime
import os

import responses

from oscar.pyrat.api import get_pyrat_data


@responses.activate
def test_get_pyrat_data():
    pyrat_url = os.environ["PYRAT_URL"]
    line_name = "Line-AB"

    # mock species response
    responses.get(
        f"{pyrat_url}/api/v3/species",
        json=[{"id": 1, "name": "Mouse"}, {"id": 2, "name": "Rat"}],
    )

    # mock animals response
    responses.get(
        f"{pyrat_url}/api/v3/animals",
        json=[
            {
                "animalid": 1,
                "eartag_or_id": "A1",
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
                        "animalid": 100,
                        "parent_id": 1000,
                        "parent_eartag": "A100",
                        "parent_sex": "m",
                        "parent_labid": "B",
                        "parent_cagenumber": "B2",
                    },
                    {
                        "animalid": 101,
                        "parent_id": 1001,
                        "parent_eartag": "A101",
                        "parent_sex": "f",
                        "parent_labid": "B",
                        "parent_cagenumber": "B2",
                    },
                ],
            }
        ],
    )

    pyrat_dfs = get_pyrat_data(
        line_name=line_name,
        species_name="Mouse",
        birth_date_from=datetime.date(2026, 2, 1),
        birth_date_to=datetime.date(2026, 2, 6),
    )
    pyrat_dfs = list(pyrat_dfs)
