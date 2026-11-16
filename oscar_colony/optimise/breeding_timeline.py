from dataclasses import dataclass
from datetime import date, timedelta


@dataclass
class BreedingTimings:
    gestation_days: int = 20


def calculate_breeding_timeline(
    target_age_days: int,
    target_date: date | None = None,
    timings: BreedingTimings = BreedingTimings(),
) -> tuple[date, date]:
    """
    Calculate when to pair the parents so their offspring are
    target_age_days old on the target_date.

    If target_date isn't given, pairing starts asap.

    Parameters
    ----------
    target_age_days : int
        Days old the animals need to be.
    target_date : date | None, optional
        Date the animals are needed by. If None - asap.
    timings : BreedingTimings, optional
        Breeding timings to use.

    Returns
    -------
    tuple[date, date]
        The date to pair the parents, and the date the animals reach
        target_age_days.

    Raises
    ------
    ValueError
        If there isn't enough time between today and target_date.
    """

    calculated_days_deadline = timedelta(
        days=timings.gestation_days + target_age_days
    )
    if target_date is None:
        target_date = date.today() + calculated_days_deadline

    breeding_start_date = target_date - calculated_days_deadline
    if breeding_start_date < date.today():
        raise ValueError(f"Can't be ready by {target_date}")

    return breeding_start_date, target_date
