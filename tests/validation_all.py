import logging
from collections import Counter

import pandas as pd
from dotenv import load_dotenv

from oscar_colony.colony_management.pyrat.api import (
    get_pyrat_data,
)
from oscar_colony.colony_management.pyrat.standardise import (
    standardise_pyrat_csv,
)
from oscar_colony.logging_setup import init_logging

logger = logging.getLogger(__name__)

# Leave these values as None to scrape all values. Specify if needed.
SPECIES_NAME = "Mouse"
LINE_NAME = None
BIRTH_DATE_FROM = None
BIRTH_DATE_TO = None


def _tally_offspring(
    df: pd.DataFrame,
    standard_df: pd.DataFrame,
    offspring_counts: Counter,
    filtered_counts: Counter,
) -> None:
    """Update the running per-line counters in place for one page of data."""
    offspring_counts.update(
        df["Line / Strain (Name)"].value_counts().to_dict()
    )

    filtered_ids = set(df["ID"]) - set(standard_df["ID_offspring"])
    raw_row_by_id = df.set_index("ID")

    for offspring_id in filtered_ids:
        line_name = raw_row_by_id.loc[offspring_id, "Line / Strain (Name)"]
        filtered_counts[line_name] += 1


def _log_summary(offspring_counts: Counter, filtered_counts: Counter) -> None:
    """Log one summary line per line/strain, e.g. 'LineX: 120 offspring,
    4 filtered'.
    """

    for line_name in sorted(offspring_counts):
        message = f"'{line_name}': {offspring_counts[line_name]} offspring"

        n_filtered = filtered_counts.get(line_name, 0)
        if n_filtered:
            message += f", {n_filtered} filtered"
            message += (
                f", {offspring_counts[line_name] - n_filtered} remaining"
            )
        logger.info(message)


def main():
    """Query pyRAT, standardise the result, and log why rows were dropped in
    the process."""
    load_dotenv()
    init_logging(output_dir="tmp/logs/", file_log_level="INFO")

    data = get_pyrat_data(
        species_name=SPECIES_NAME,
        line_name=LINE_NAME,
        birth_date_from=BIRTH_DATE_FROM,
        birth_date_to=BIRTH_DATE_TO,
    )

    offspring_counts: Counter = Counter()
    filtered_counts: Counter = Counter()

    for df in data:
        if df.empty:
            continue

        standard_df = standardise_pyrat_csv(df)
        _tally_offspring(df, standard_df, offspring_counts, filtered_counts)

    logger.info(f"{len(offspring_counts)} strains found matching this query")

    _log_summary(offspring_counts, filtered_counts)


if __name__ == "__main__":
    main()
