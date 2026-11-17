from datetime import date, timedelta

import pytest

from oscar_colony.optimise.breeding_timeline import (
    BreedingTimings,
    calculate_breeding_timeline,
)

TODAY = date.today()


@pytest.mark.parametrize(
    "timings, target_date, expected_pairing_date, expected_ready_date",
    [
        pytest.param(
            BreedingTimings(),
            None,
            TODAY,
            TODAY + timedelta(days=50),
            id="asap",
        ),
        pytest.param(
            BreedingTimings(),
            TODAY + timedelta(days=100),
            TODAY + timedelta(days=50),
            TODAY + timedelta(days=100),
            id="target date",
        ),
        pytest.param(
            BreedingTimings(gestation_days=5),
            None,
            TODAY,
            TODAY + timedelta(days=35),
            id="custom gestation",
        ),
    ],
)
def test_calculate_breeding_timeline(
    timings, target_date, expected_pairing_date, expected_ready_date
):
    assert calculate_breeding_timeline(30, target_date, timings) == (
        expected_pairing_date,
        expected_ready_date,
    )


@pytest.mark.parametrize(
    "target_age_days, target_date, timings",
    [
        pytest.param(
            30, TODAY + timedelta(days=49), BreedingTimings(), id="too soon"
        ),
        pytest.param(-10, None, BreedingTimings(), id="negative age"),
        pytest.param(
            30,
            None,
            BreedingTimings(gestation_days=-10),
            id="negative gestation",
        ),
    ],
)
def test_calculate_breeding_timeline_invalid(
    target_age_days, target_date, timings
):
    with pytest.raises(ValueError):
        calculate_breeding_timeline(target_age_days, target_date, timings)
