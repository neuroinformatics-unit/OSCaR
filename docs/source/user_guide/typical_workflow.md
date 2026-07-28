# Typical workflow

Here, we'll walk through a typical workflow for using OSCaR.

## Fetch and standardise data from a colony management system

First, we fetch animal data via the colony management system's api, and convert it to OSCaR's [standard table format](./standard_table.md). At the moment, only [PyRAT](https://www.scionics.com/pyrat.html) is supported - but we hope to expand to more systems in future.

- See the [PyRAT docs]() for details of how to setup and use pyRAT for this step
- Alternatively, you can load your animal data via any other means (e.g. a custom script / from an exported csv file), and convert it to OSCaR's [standard table format](./standard_table.md) manually.

## Calculate historical stats

Next, we calculate summary statistics based on this standardised data.

```python
from oscar_colony.historical_stats import calculate_historical_stats_for_line

line_stats = calculate_historical_stats_for_line(
    standardised_data, # a pandas DataFrame in OSCaR's standard csv format
    line_name="MY-LINE"  # The name of the line we want to process
)
```
The :func:`calculate_historical_stats_for_line` function produces a :class:`LineStatistics` object, with summary statistics for the specified line.

For example:
```python
print(line_stats.total_n_offspring)
print(line_stats.total_n_successful_matings)
```

## Define the required animals

Now we can define the number and genotypes of animals that we need:

```python
from oscar_colony.breeding_scheme import Genotype

# Here we are asking for 20 animals of genotype wt_het, 53 of het_het and
# 27 of hom_hom
required_n_per_genotype = {
    (Genotype.WT, Genotype.HET): 20,
    (Genotype.HET, Genotype.HET): 53,
    (Genotype.HOM, Genotype.HOM): 27
}
```
You can include any number of genotypes here. Just make sure that your Genotype tuples:
- Only contain `Genotype.WT`, `Genotype.HET` or `Genotype.HOM`
- Have length matching the `n_mutations` in your standard table for the line of interest. So e.g. if `n_mutations=1` for this line, your tuples may look like: `(Genotype.WT,)`, or if `n_mutations=3` it may be: `(Genotype.HET, Genotype.HOM, Genotype.WT)`. All the tuples you provide should have the same length.
- The order of values matches the order of `mutations` in your standard table. For example, if `mutations` is `Mut-A_Mut-B` and you provide a genotype of: `(Genotype.HET, Genotype.WT)`, you are asking for animals that are HET for `Mut-A` and WT for `Mut-B`.

## Optimise schemes

Now we have defined the animals we need, OSCaR can calculate an optimal combination of breeding schemes based on the historical stats for that line.

```python
from oscar_colony.optimise.optimal_scheme_calculator import calculate_optimal_scheme

breeding_schemes, surplus = calculate_optimal_scheme(
    required_n_per_genotype,
    line_stats=line_stats,
    default_litter_size=6
)
```
Make sure you set the `default_litter_size` to an appropriate value. For example, you could set it to the average litter size across all your institution's data.

```python
print(breeding_schemes)
```

```python
print(surplus)
```
