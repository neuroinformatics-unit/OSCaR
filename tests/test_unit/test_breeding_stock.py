import numpy as np
import pandas as pd

from oscar_colony.breeding_scheme import Genotype
from oscar_colony.breeding_stock import count_breeding_stock
from oscar_colony.optimise.surplus_summary import SexSplit


def test_count_breeding_stock():
    """Test males and females are counted per genotype for the given line
    only, leaving out un-genotyped animals."""

    standardised_data = pd.DataFrame(
        {
            "line_name": ["Line-A"] * 5 + ["Line-B"],
            "genotype_offspring": ["het", "het", "het", "hom", np.nan, "wt"],
            "offspring_sex": ["m", "f", "f", "m", "f", "m"],
        }
    )

    available_breeding_stock = count_breeding_stock(
        standardised_data, "Line-A"
    )

    assert available_breeding_stock == {
        (Genotype.HET,): SexSplit(n_males=1, n_females=2),
        (Genotype.HOM,): SexSplit(n_males=1, n_females=0),
    }
