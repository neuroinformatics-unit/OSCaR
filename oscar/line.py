from pathlib import Path

import pandas as pd

from oscar.breeding_scheme import generate_breeding_schemes


class Line:
    def __init__(self, historical_data: Path | pd.DataFrame, line_name: str):
        if isinstance(historical_data, Path):
            historical_data = pd.read_csv(historical_data)

        # TODO - check fits expected structure / column names

        line_data = historical_data[
            historical_data["Line / Strain (Name)"] == line_name
        ]
        if len(line_data) == 0:
            raise ValueError(f"No data for line {line_name} exists.")

        line_data = self._filter_allowed_genotypes(line_data)
        if len(line_data) == 0:
            raise ValueError(
                f"No valid genotypes (grade 1 / 2 / 3) exist for {line_name}."
            )

        self.n_mutations = 0
        for col in ["Grade 1", "Grade 2", "Grade 3"]:
            if not line_data[col].isna().all():
                self.n_mutations += 1

        if self.n_mutations == 0:
            raise ValueError(f"No mutations recorded for line {line_name}")

        self._calculate_summary_statistics()

    def _filter_allowed_genotypes(self, line_data: pd.DataFrame):
        """Only keep allowed genotypes of wt, het or hom.

        Remove others such as +/-, Tg, ko, as well as ungenotyped
        individuals (na for grade 1/2/3)

        Parameters
        ----------
        line_data : pd.DataFrame
            Data for a single line
        """

        genotype_columns = [
            "Grade 1",
            "Grade 2",
            "Grade 3",
            "Father: Grade 1",
            "Father: Grade 2",
            "Father: Grade 3",
            "Mother: Grade 1",
            "Mother: Grade 2",
            "Mother: Grade 3",
        ]

        # remove rows containing forbidden genotypes
        genotype_data = line_data.loc[:, genotype_columns]
        allowed_genotypes = (
            genotype_data.isin(["wt", "het", "hom"]) | genotype_data.isna()
        ).all(axis=1)
        filtered_data = line_data.loc[allowed_genotypes, :]

        # remove rows for ungenotyped individuals (na for grade 1/2/3)
        ungenotyped = (
            line_data.loc[:, ["Grade 1", "Grade 2", "Grade 3"]]
            .isna()
            .all(axis=1)
        )
        filtered_data = filtered_data.loc[~ungenotyped, :]

        return filtered_data

    def _calculate_summary_statistics(self):
        self.breeding_schemes = generate_breeding_schemes(self.n_mutations)
