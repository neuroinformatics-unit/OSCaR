# OSCaR

OSCaR (the Open Source Colony-assessment Resource) is a python package for data-driven optimisation of breeding strategies.

## Development setup

First, install the package locally - here we provide instructions for both `uv` and `conda` (use whichever you prefer).

### uv

Install uv following the [instructions on uv's website](https://docs.astral.sh/uv/getting-started/installation/) then run:
```
git clone git@github.com:neuroinformatics-unit/OSCaR.git
cd OSCaR
uv sync --extra dev
```

### Conda

Install conda e.g. [via miniforge](https://conda-forge.org/download/) then run:

```
git clone git@github.com:neuroinformatics-unit/OSCaR.git

conda create -n oscar-dev python=3.12
conda activate oscar-dev
pip install -e ".[dev]"
```

### Pre-commit

We use [pre-commit](https://pre-commit.com/) for automated linting / formatting.
```
# setup pre-commit to run on every commit
pre-commit install
```

### PyRAT access

#### Token creation

To access a [PyRAT](https://www.scionics.com/pyrat.html) instance via the code in `oscar.colony_management.pyrat`, you will first need to contact Scionics (the creators of PyRAT) to approve access.

- Ask them to setup a new API client called 'OSCaR' on your pyRAT server; they should send you a 'client token'.
- The new API client will appear under `ADMINISTRATION > API` in pyRAT, and an administrator will need to enable it from there.

This approval process only needs to be done once per pyRAT instance.

Once the client is enabled, any user can request a client token by logging into PyRAT, going to `ADMINISTRATION > API` and clicking the 'Request access' button on the `OSCaR` row.

#### Providing tokens to OSCaR

Once all tokens have been created (as above), they can be provided to OSCaR by setting the following environment variables:

- `PYRAT_URL`: the url of the pyRAT instance you want to access
- `PYRAT_CLIENT_TOKEN`: the client token sent to you by Scionics
- `PYRAT_USER_TOKEN`: the user token created via 'Request access' in PyRAT

Set these variables however you prefer - just make sure they are _never_ checked in to version control (tokens should be kept secret). One convenient method is to create a `.env` file at the top level of the repository, then use [`python-dotenv`](https://github.com/theskumar/python-dotenv) to read them.
