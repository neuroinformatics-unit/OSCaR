from collections import Counter, defaultdict

import pandas as pd
from dotenv import load_dotenv

from oscar_colony.colony_management.pyrat.api import (
    get_pyrat_data,
    get_pyrat_line_mutations,
    get_pyrat_lines,
)
from oscar_colony.colony_management.pyrat.standardise import (
    standardise_pyrat_csv,
)

load_dotenv()

# Leave these values as None to scrape all values. Specify if needed.
SPECIES_NAME = "Mouse"
LINE_NAME = None
BIRTH_DATE_FROM = None
BIRTH_DATE_TO = None

OUTPUT_CSV_PATH = "tmp/validation_report.csv"

pyrat_lines_df = pd.concat(get_pyrat_lines(), ignore_index=True)
line_id_by_name = dict(zip(pyrat_lines_df["name"], pyrat_lines_df["id"]))

found_count_by_line: Counter = Counter()
filtered_count_by_line: Counter = Counter()
mutations_seen_by_line: dict[str, set[str]] = defaultdict(set)

pyrat_data_pages = get_pyrat_data(
    species_name=SPECIES_NAME,
    line_name=LINE_NAME,
    birth_date_from=BIRTH_DATE_FROM,
    birth_date_to=BIRTH_DATE_TO,
)

for raw_df in pyrat_data_pages:
    if raw_df.empty:
        continue

    standardised_df = standardise_pyrat_csv(raw_df)

    found_count_by_line.update(
        raw_df["Line / Strain (Name)"].value_counts().to_dict()
    )

    filtered_out_ids = set(raw_df["ID"]) - set(standardised_df["ID_offspring"])
    filtered_out_df = raw_df[raw_df["ID"].isin(filtered_out_ids)]
    filtered_count_by_line.update(
        filtered_out_df["Line / Strain (Name)"].value_counts().to_dict()
    )

    offspring_mutation_cols = raw_df.filter(
        regex=r"^Mutation \d+$"
    ).columns.tolist()
    for line_name, line_offspring_df in raw_df.groupby("Line / Strain (Name)"):
        raw_mutation_values = pd.unique(
            line_offspring_df[offspring_mutation_cols].values.ravel("K")
        )
        mutations_seen_by_line[line_name].update(
            pd.Series(raw_mutation_values).dropna()
        )

expected_mutations_by_line = {
    line_name: set(get_pyrat_line_mutations(line_id))
    for line_name, line_id in line_id_by_name.items()
}

report_rows = []
for line_name in sorted(set(line_id_by_name) | set(found_count_by_line)):
    found_count = found_count_by_line.get(line_name, 0)
    filtered_count = filtered_count_by_line.get(line_name, 0)

    expected_mutations_for_line = expected_mutations_by_line.get(line_name)
    seen_mutations_for_line = mutations_seen_by_line.get(line_name, set())

    failures = []
    if found_count == 0:
        failures.append("no offspring found")

    if (
        expected_mutations_for_line is not None
        and not expected_mutations_for_line
    ):
        failures.append("wildtype background stock (no mutations in pyRAT)")

    if expected_mutations_for_line is not None and found_count:
        missing_mutations = (
            expected_mutations_for_line - seen_mutations_for_line
        )
        unexpected_mutations = (
            seen_mutations_for_line - expected_mutations_for_line
        )

        capitalisation_mismatches = []
        for unexpected_mutation in list(unexpected_mutations):
            for expected_mutation in expected_mutations_for_line:
                if unexpected_mutation.lower() == expected_mutation.lower():
                    capitalisation_mismatches.append(
                        f"'{unexpected_mutation}' vs '{expected_mutation}'"
                    )
                    unexpected_mutations.discard(unexpected_mutation)
                    missing_mutations.discard(expected_mutation)
                    break

        if capitalisation_mismatches:
            failures.append(
                "inconsistent capitalisation: "
                f"{sorted(capitalisation_mismatches)}"
            )
        if missing_mutations:
            failures.append(f"missing mutations: {sorted(missing_mutations)}")
        if unexpected_mutations:
            failures.append(
                f"unexpected mutations: {sorted(unexpected_mutations)}"
            )

    if found_count and expected_mutations_for_line is None:
        failures.append("line name not active in pyrat")

    report_rows.append(
        {
            "line_name": line_name,
            "found": found_count,
            "filtered": filtered_count,
            "remaining": found_count - filtered_count,
            "failure": "; ".join(failures),
        }
    )

report_df = pd.DataFrame(report_rows).set_index("line_name")
with open(OUTPUT_CSV_PATH, "w") as csv_file:
    report_df.to_csv(csv_file)
print(f"Validation report written to {OUTPUT_CSV_PATH}")

lines_with_remaining_offspring = report_df["remaining"] > 0
lines_fully_filtered = (report_df["found"] > 0) & (report_df["remaining"] == 0)

print(f"\n{len(report_df)} lines checked")
print("-----------------")
print("")
print(
    f"{lines_with_remaining_offspring.sum()} : lines with offspring "
    "remaining after filtering"
)
print(f"{(report_df['found'] == 0).sum()} : lines with no offspring")

print(
    f"{lines_fully_filtered.sum()} : lines where all offspring were filtered"
)
print(
    f"{report_df['failure'].str.contains('wildtype background').sum()} : "
    "Wildtype stock like BL6"
)
print(
    f"{report_df['failure'].str.contains('inconsistent capitalisation').sum()}"
    " : with inconsistent mutation capitalisation"
)
print(
    f"{report_df['failure'].str.contains('missing mutations').sum()} : "
    "with missing mutation(s)"
)
print(
    f"{report_df['failure'].str.contains('unexpected mutations').sum()} : "
    "with unexpected mutation(s)"
)
print(
    f"{report_df['failure'].str.contains('not active in pyrat').sum()} : "
    "with a line name not active in pyrat"
)
