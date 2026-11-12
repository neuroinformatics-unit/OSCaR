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
    days_to_ready = timedelta(days=timings.gestation_days + target_age_days)
    if target_date is None:
        target_date = date.today() + days_to_ready

    pairing_date = target_date - days_to_ready
    if pairing_date < date.today():
        raise ValueError(f"Can't be ready by {target_date}")

    return pairing_date, target_date
