# Installation

The `oscar-colony` package can be installed from PyPI using pip:
```
pip install oscar-colony
```

## Using conda

If you want to use conda, you can create an environment first like:
```
conda create -n oscar-env python=3.14
conda activate oscar-env
```

Then install the package into the environment:
```
pip install oscar-colony
```

## Using uv

You can create a new directory with a minimal uv project inside with:
```
uv init oscar-project
```

Then enter the project and add oscar-colony:
```
cd oscar-project
uv add oscar-colony
```

This will create a virtual environment at `.venv` in this directory.
