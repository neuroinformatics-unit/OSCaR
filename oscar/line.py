import logging
from pathlib import Path

import pandas as pd

from oscar.breeding_scheme import (
    BreedingScheme,
    generate_breeding_schemes,
)

logger = logging.getLogger(__name__)


class Line:
    def __init__(self, historical_data: Path | pd.DataFrame, line_name: str):
        if isinstance(historical_data, Path):
            historical_data = pd.read_csv(historical_data)

        # TODO - check fits expected structure / column names. Also > 0
        # mutations and same number for whole line

        line_data = historical_data[historical_data["line_name"] == line_name]
        if len(line_data) == 0:
            raise ValueError(f"No data for line {line_name} exists.")

        self._historical_stats: dict[BreedingScheme, float | None] = {}
        self._historical_stats_for_line(line_data)

    def _historical_stats_for_line(self, line_data: pd.DataFrame):
        self.breeding_schemes = generate_breeding_schemes(
            line_data.n_mutations.iloc[0]
        )

        for breeding_scheme in self.breeding_schemes:
            self._historical_stats_for_breeding_scheme(
                breeding_scheme, line_data
            )

    def _historical_stats_for_breeding_scheme(
        self, scheme: BreedingScheme, line_data: pd.DataFrame
    ):
        # filter all rows for that breeding scheme - allow it to match in
        # any order
        scheme_f_m = (
            line_data["father_genotype"] + "x" + line_data["mother_genotype"]
        )
        scheme_m_f = (
            line_data["mother_genotype"] + "x" + line_data["father_genotype"]
        )

        scheme_data = line_data[
            (scheme_f_m == str(scheme)) | (scheme_m_f == str(scheme))
        ]
        if len(scheme_data) == 0:
            msg = (
                f"No data available for scheme {str(scheme)} with line "
                f"{line_data.line_name.iloc[0]}"
            )
            logger.info(msg)
            self._historical_stats[scheme] = None
            return

        # breeding pairs is unique combos of father ID x mother ID
        n_breeding_pairs = scheme_data.groupby(
            ["ID_father", "ID_mother"]
        ).ngroups

        return n_breeding_pairs
