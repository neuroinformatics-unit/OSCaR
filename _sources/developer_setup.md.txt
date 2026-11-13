# Developer guide

## Installation

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

## Pre-commit

We use [pre-commit](https://pre-commit.com/) for automated linting / formatting.
```
# setup pre-commit to run on every commit
pre-commit install
```

## PyRAT access

See the [PyRAT docs](./user_guide/pyrat), for information on how to set up client and user tokens.

## Test data

All test data is stored in the [oscar-test-data GIN repository](https://gin.swc.ucl.ac.uk/neuroinformatics/oscar-test-data), and fetched using [`pooch`](https://www.fatiando.org/pooch/latest/).

If you add / update a test data file, you will need to update the file names and hashes in the pooch registry at `tests/pooch_registry.txt`. Hashes can be generated using [the instructions in poochs' docs](https://www.fatiando.org/pooch/latest/hashes.html#calculating-hashes).

## Building the docs locally

To build the documentation locally, you will need to install some additional dependencies, then run `sphinx-build` (as below).

### docs install with uv
```
uv sync --all-extras
```

### docs install with conda
```
conda activate oscar-dev
pip install -e ".[dev,docs]"
```

### docs build command
```
sphinx-build docs/source docs/build
```
Then open the generated `docs/build/index.html` file.

To re-build the documentation after making changes, remove the docs/build
folder and re-run the above command:
```
rm -rf docs/build
sphinx-build docs/source docs/build
```

## Deploying the docs

The documentation is deployed automatically from `main` when a new tag is created on GitHub (usually when making a new release).
