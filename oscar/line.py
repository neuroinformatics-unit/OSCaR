from pathlib import Path

import pandas as pd

from oscar.breeding_scheme import BreedingScheme, generate_breeding_schemes


class Line:
    def __init__(self, historical_data: Path | pd.DataFrame, line_name: str):
        if isinstance(historical_data, Path):
            historical_data = pd.read_csv(historical_data)

        # TODO - check fits expected structure / column names. Also > 0
        # mutations and same number for whole line

        line_data = historical_data[historical_data["line_name"] == line_name]
        if len(line_data) == 0:
            raise ValueError(f"No data for line {line_name} exists.")

        self._summary_stats_for_line(line_data)

    def _summary_stats_for_line(self, line_data: pd.DataFrame):
        self.breeding_schemes = generate_breeding_schemes(
            line_data.n_mutations.iloc[0]
        )

        for breeding_scheme in self.breeding_schemes:
            self._summary_stats_for_breeding_scheme(breeding_scheme, line_data)

    def _summary_stats_for_breeding_scheme(
        self, scheme: BreedingScheme, line_data: pd.DataFrame
    ):
        pass
