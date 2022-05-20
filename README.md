# Avatar paper

## Description

Avatar paper is an analysis project that allows researchers to see, re do, or explore figure and data of the "Patient-centric synthetic data generation, no reason to risk re-identification the analysis of biomedical  pseudonymised data".
The aim of the paper is to present the avatarization method.

We really recommend to have a look to the scientific paper (https://doi.org/10.21203/rs.3.rs-1674043/v1) before to explore the repository.

## Prerequisites

Install [Poetry](https://python-poetry.org/) (a Python package and dependency management tool).

This repo uses Git LFS (Large File Storage) to handle large datasets. Make sure to install it for your platform, following the instructions at https://git-lfs.github.com/

## Install

Run:

```bash
make install
```

## R packages

All required packages will be installed using the R package `librarian`.

## How to use it

This code mainly consists in Jupyter notebooks and datasets.

It is written in Python and R.

- `git lfs pull` download large file using `git lfs`, can take around 30min.
- Run the following command in a terminal.  `make notebook`
- Open notebook - such as `0.main.ipynb`.
- Run the cells.

## Structure

```
datasets/
   AIDS/          # original and avatarized AIDS datasets
   WBCD/          # original and avatarized WBCD datasets
   results_df/    # computationally expensive analysis results.

notebooks/        # analysis and graph generation

metrics/
   privacy_metrics/       # multiple function used to compute avatarization metrics

figures/           # figures presented in the article
color.csv         # colors for the figures
```

## Contributing

Feel free to do the analysis again or to explore avatarized datasets.

## License

This code is licensed under the Apache-2.0 license.

## Contributors

- [Morgan Guillaudeux](https://github.com/mguillaudeux) morgan@octopize.io
- [Olivia Rousseau](https://github.com/oliviarousseau) Olivia.Rousseau@univ-nantes.fr 
- [Julien Petot](https://github.com/jpetot) julien@octopize.io
- [Pierre-Antoine Gourraud](https://github.com/gourraud) pierre-antoine.gourraud@univ-nantes.fr, corresponding author.
