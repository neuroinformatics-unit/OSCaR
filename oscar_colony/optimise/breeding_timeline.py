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
    calculated_days_deadline = timedelta(
        days=timings.gestation_days + target_age_days
    )
    if target_date is None:
        target_date = date.today() + calculated_days_deadline

    breeding_start_date = target_date - calculated_days_deadline
    if breeding_start_date < date.today():
        raise ValueError(f"Can't be ready by {target_date}")

    return breeding_start_date, target_date
