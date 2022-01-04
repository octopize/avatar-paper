# Avatar paper

## Description

Avatar paper is an analysis project that allows researchers to see, re do, or explore figure and data of the "PAPER_NAME".
The aim of the paper is to present the avatarization method.

We really recommend to have a look to the scientific paper (link) before to explore the repository.

## Prerequisites

Install [Poetry](https://python-poetry.org/) (a Python package and dependency management tool).

This repo uses Git LFS (Large File Storage) to handle large datasets. Make sure to install it for your platform, following the instructions at https://git-lfs.github.com/

## Install

Run:

```bash
make install
```

## R packages

Install all required packages on the top of each R notebooks with the command `install.packages('PACKAGE_NAME')`.  They are also presented in the `sessioninfo.txt` with their dependencies.

## How to use it

This code mainly consists in Jupyter notebooks and datasets.

It is written in Python and R.

- `git lfs pull` download large file using `git lfs`, can take around 30min.
- Run the following command in a terminal.  `jupyter notebook # open a notebook`
- Open notebook - such as `messageAB_aids.ipynb`.
- For R kernel: install required packages with the command line `install.packages('PACKAGE_NAME')` if it was not done previously.
- Run the cells.

## Structure

```
datasets/
   AIDS/          # original and avatarized AIDS datasets
   WBCD/          # original and avatarized WBCD datasets
   results_df/    # computationally expensive analysis results.
   messageD/      # dataset for figure 5 message

notebooks/
   final_figure/  # analysis and graph generation

lsg/
   dimension/     # multidimensionnal functions such as PCA, MCA or FAMD
   metrics/       # multiple function used to compute avatarization metrics

figure/           # figures presented in the article
color.csv         # colors for the figures
sessioninfo.txt   # R session information.
```

## Contributing

Feel free to do the analysis again or to explore avatarized datasets.

## License

This code is licensed under the MIT license.

## Contributors

- [Julien Petot](https://github.com/jpetot)
- [Morgan Guillaudeux](https://github.com/mguillaudeux)
- [Olivia Rousseau](https://github.com/oliviarousseau)
- [Pierre-Antoine Gourraud](https://github.com/gourraud)
