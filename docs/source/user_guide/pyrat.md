# PyRAT

[PyRAT](https://www.scionics.com/pyrat.html) is a popular animal facility management software.

To access data in a PyRAT instance via OSCaR, there are some initial setup steps:

## Client tokens

This step should be done by an administrator of your institution's PyRAT instance, and only needs to be done once.

You will need to contact Scionics (the creators of PyRAT) to approve access for OSCaR:

- Ask them to setup a new API client called 'OSCaR' on your PyRAT server; they should send you a 'client token'.
- The new API client will appear under `ADMINISTRATION > API` in PyRAT, and an administrator will need to enable it from there.

## User tokens

Once the client is enabled, any user can request a user token by logging into PyRAT, going to `ADMINISTRATION > API` and clicking the 'Request access' button on the `OSCaR` row.

## Providing tokens to OSCaR

Once all tokens have been created (as above), they can be provided to OSCaR by setting the following environment variables:

- `PYRAT_URL`: the url of the PyRAT instance you want to access
- `PYRAT_CLIENT_TOKEN`: the client token sent to you by Scionics
- `PYRAT_USER_TOKEN`: the user token created via 'Request access' in PyRAT

Set these variables however you prefer - just make sure they are _never_ checked in to version control (tokens should be kept secret). One convenient method is to create a `.env` file and use [`python-dotenv`](https://github.com/theskumar/python-dotenv) to read them. We'll show an example of that below.

## Pulling data from the PyRAT API

See the docs for the {mod}`~oscar_colony.colony_management.pyrat.api` module, for full details of how OSCaR can pull data from the PyRAT API.
Below, we'll walk through an example of retrieving data for a named line.

First, let's set our environment variables. We'll use `python-dotenv` for this, but feel free to use whatever method you like.

We can create a `.env` file containing:
```
PYRAT_URL=https://my-pyrat-instance
PYRAT_CLIENT_TOKEN=my-client-token
PYRAT_USER_TOKEN=my-user-token
```
You would need to update the values to match your PyRAT instance.

Then load it with:
```python
from dotenv import load_dotenv

load_dotenv()
```

Now we can fetch data via OSCaR. Let's retrieve all available animals from a given line:
```
from oscar_colony.colony_management.pyrat.api import get_pyrat_data
animal_data = get_pyrat_data(
    species_name="Mouse",
    line_name="MY-LINE",
)
```

As this could be a very large amount of data, it is returned as a generator of
[pandas](https://pandas.pydata.org/) dataframes that we can loop through to process one-by-one.

## Standardising data from the PyRAT API

Now we have the data, it must be converted into OSCaR's [standard table format](./standard_table).

We can do this with:
```python
from oscar_colony.colony_management.pyrat.standardise import standardise_pyrat_csv
import pandas as pd

standardised_dfs = []

# Loop through the returned animal data, and standardise each
for animal_df in animal_data:
    standard_df = standardise_pyrat_csv(animal_df)
    standardised_dfs.append(standard_df)

# Join all the standardised dataframes into one
standard_df = pd.concat(standardised_dfs)
```

See the docs for the {mod}`~oscar_colony.colony_management.pyrat.standardise` module, for full details of how data is standardised.
