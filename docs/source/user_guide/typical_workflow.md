# Typical workflow

Here, we'll walk through a typical workflow for using OSCaR.

## Fetch and standardise data from a colony management system

First, we fetch animal data via a colony management system's API, and convert it to OSCaR's [standard table format](./standard_table). At the moment, only [PyRAT](https://www.scionics.com/pyrat.html) is supported - but we hope to expand to more systems in future.

- See the [PyRAT docs](./pyrat) for details of how to setup and use PyRAT for this step
- Alternatively, you can load your animal data via any other means (e.g. a custom script / from an exported csv file), and convert it to OSCaR's [standard table format](./standard_table) manually.

## Calculate historical stats

Next, we calculate summary statistics based on this standardised data.

```python
from oscar_colony.historical_stats import calculate_historical_stats_for_line

line_stats = calculate_historical_stats_for_line(
    standard_df, # a pandas DataFrame in OSCaR's standard csv format
    line_name="MY-LINE"  # The name of the line we want to process
)
```
The {func}`~oscar_colony.historical_stats.calculate_historical_stats_for_line` function produces a {class}`~oscar_colony.historical_stats.LineStatistics` object, with summary statistics for the specified line.

For example:
```python
# Access the total number of offspring for this line
print(line_stats.total_n_offspring)

# Access the total number of successful matings for this line
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
- The order of values matches the order of `mutations` in your standard table. For example, if your table contains columns for `mutation_1` and `mutation_2` containing `Mut-A` and `Mut-B` respectively, then a tuple of: `(Genotype.HET, Genotype.WT)` means HET for `Mut-A` and WT for `Mut-B`.

### Requesting a specific number of males and females

If you require a certain number of each sex, input an int into the named tuple SexSplits `(n_males, n_females)` instead of a single number:

```python
from oscar_colony.breeding_scheme import Genotype

# Here we are asking for 5 male and 10 female of genotype wt_het, 33 male and 40 female of het_het and
# 52 unspecified animals of hom_hom
required_n_per_genotype = {
    (Genotype.WT, Genotype.HET): (5, 10),
    (Genotype.HET, Genotype.HET): (33, 40),
    (Genotype.HOM, Genotype.HOM): 52
}
```
You can mix the two methods, so you only have to specify a sex split where it matters as shown above.

The expected proportion of males comes from the historical stats for the line or for the specific breeding scheme, where there's enough data. Otherwise it uses a 50:50 split.

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

`breeding_schemes` will contain the recommended schemes, along with the number of matings proposed for each:
```python
print(breeding_schemes)
```

`surplus` will contain a summary of expected totals and surplus (i.e. how many extra animals are expected beyond the `required_n_per_genotype` you asked for). You can see all included values in the {class}`~oscar_colony.optimise.surplus_summary.SurplusSummary` docs. Where you asked for a sex split, the summary will also include the expected male and female totals and surplus.

```python
print(surplus)
```
