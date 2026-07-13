from datetime import date

from dotenv import load_dotenv

from oscar_colony.breeding_scheme import BreedingScheme
from oscar_colony.colony_management.pyrat.api import get_pyrat_data
from oscar_colony.colony_management.pyrat.standardise import (
    standardise_pyrat_csv,
)
from oscar_colony.historical_stats import calculate_historical_stats_for_line
from oscar_colony.optimise.estimate_offspring import (
    estimate_n_offspring_per_mating,
)
from oscar_colony.optimise.surplus_summary import (
    create_surplus_summary,
)

load_dotenv()

line_name = "BF/Vgat-cre"
birth_date_from = date(2023, 7, 1)
birth_date_to = date(2023, 7, 28)

data = get_pyrat_data(
    species_name="Mouse",
    line_name=line_name,
    birth_date_from=birth_date_from,
    birth_date_to=birth_date_to,
)

for df in data:
    standard_df = standardise_pyrat_csv(df)
    stats = calculate_historical_stats_for_line(standard_df, line_name)

    proposed_n_matings_per_schemes = {
        BreedingScheme("het", "het"): 11,
        BreedingScheme("wt", "hom"): 25,
        BreedingScheme("het", "hom"): 8,
        BreedingScheme("wt", "het"): 1,
    }

    offspring_per_scheme = estimate_n_offspring_per_mating(
        stats, int(round(stats.average_litter_size))
    )

    required_n_per_genotype: dict = {}
    surplus_summary = create_surplus_summary(
        required_n_per_genotype,
        proposed_n_matings_per_schemes,
        offspring_per_scheme,
    )

    print("breeding scheme:")
    print(f"{birth_date_from} - {birth_date_to}")
    for scheme, n_matings in proposed_n_matings_per_schemes.items():
        print(f"{scheme}: {n_matings}")

    print(
        "\nPredicted total animals (all genotypes):",
        int(round(surplus_summary.total_n)),
    )
    print("Predicted totals per genotype:")
    for genotype, surplus in surplus_summary.surplus_per_genotype.items():
        print(f"{genotype}: {int(round(surplus.total_n))}")
