# Installation

Use conda-forge for the core scientific stack:

```bash
conda env create -f environment.yml
conda activate dominant-structural-variant-ranker
python -m pip install -e ".[dev]"
```

Install optional tools only when the corresponding workflow step is required.

