# Standard table structure

All steps beyond those in: `oscar_colony.colony_management` rely on having data in a standard table format.

Example of format for a line with 2 mutations:

| ID_offspring | line_name | offspring_sex | date_of_birth | ID_father_1 | ID_mother_1 | sacrifice_reason | n_mutations | mutation_1 | mutation_2 | genotype_offspring | genotype_father | genotype_mother |
| ------ | ------  | - | ---------- | ----- | ----- | ----------------- | - | ------|-------| --------| --------| --------|
| ID-001 | Line-AB | m | 15/12/2025 | ID-F1 | ID-M1 | End of experiment | 2 | Mut-A | Mut-B | hom_hom | het_hom | hom_het |
| ID-002 | Line-AB | f | 15/12/2025 | ID-F1 | ID-M1 | End of experiment | 2 | Mut-A | Mut-B | wt_wt   | wt_wt   | wt_wt   |
| ID-003 | Line-AB | f | 02/01/2026 | ID-F1 | ID-M1 | End of experiment | 2 | Mut-A | Mut-B | wt_het  | wt_hom  | wt_wt   |

Each row represents one animal with columns:
- `ID_offspring`: the ID of the animal
- `line_name`: the name of the line
- `offspring_sex`: the sex of the offspring
- `date_of_birth`: date of birth
- `ID_father_1`: the ID of the animal's father
- `ID_mother_1`: the ID of the animal's mother
- `sacrifice_reason`: a description of why the animal was sacrificed
- `n_mutations`: the number of mutations the line has. This should match the number of `mutation_NUMBER` columns (see below)
- `mutation_1` / `mutation_2`...: the names of the mutations for this line (these should be identical, and in the same order, across all animals from the same line)

    For all genotype columns below, the number of values is equal to `n_mutations` and is given in the order of the mutation names `mutation_1`, `mutation_2`...:
- `genotype_offspring`: the genotype of the animal. `wt`, `het` or `hom` _only_, separated by underscores.
- `genotype_father`: the genotype of the animal's father. `wt`, `het` or `hom` _only_, separated by underscores.
- `genotype_mother`: the genotype of the animal's mother. `wt`, `het` or `hom` _only_, separated by underscores.

## Mutations

Any number of mutations is supported, just make sure it is consistent throughout all animals belonging to a particular line.

For example, to update the table above for a line with 3 mutations, you would:
- Add another column for `mutation_3`
- Update `n_mutations` to 3
- Update all genotype columns to contain 3 values e.g. wt_het_hom, het_hom_hom...

## Parents

In the example above, all animals had one father and one mother - but OSCaR does support multi-parent scenarios. For example, you may have put two female animals together with one male for a particular mating. This can be represented by adding further columns for the extra parent ids e.g.:

| ID_offspring | line_name | date_of_birth | ID_father_1 | ID_mother_1 | ID_mother_2 | sacrifice_reason | n_mutations | mutation_1 | mutation_2 | genotype_offspring | genotype_father | genotype_mother |
| ------ | ------  | ---------- | ----- | ----- | ------|------------------ | - | ------|------ | --------| --------| --------|
| ID-001 | Line-AB | 15/12/2025 | ID-F1 | ID-M1 | ID-M2 | End of experiment | 2 | Mut-A | Mut-B | hom_hom | het_hom | hom_het |

Here we have added one extra column for `ID_mother_2`.

Note: there should still only be one `genotype_mother` and one `genotype_father` column, as the genotypes of all parents of the same sex should be identical. If they are not, then there's no unambiguous way of determining the breeding scheme that animal came from.

All animals must have at least one mother and one father listed.

## Un-genotyped animals

Un-genotyped animals can also be included by leaving the `genotype_offspring` value empty for a particular row. These animals will still be included in appropriate parts of the {class}`~oscar_colony.historical_stats.LineStatistics`.
