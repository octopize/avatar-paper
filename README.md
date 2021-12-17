# Description

Avatar paper is an analysis project that allows researchers to see, re do, or explore figure and data of the "PAPER_NAME".
The aim of the paper is to present the avatarization method.

We really recommend to have a look to the scientific paper (link) before to explore the repository.

## Prerequisites

This repo uses Git LFS (Large File Storage) to handle large datasets. Make sure to install it for your platform, following the instructions at https://git-lfs.github.com/

## Python package

```
matplotlib==3.4.2
matplotlib-inline==0.1.2
numpy==1.20.3
pandas==1.2.4
scikit-learn==0.24.2
scipy==1.6.3
seaborn==0.11.1
```

There is also a `requirements.txt`.

## R packages

Install all required packages on the top of each R notebooks with the command `install.packages('PACKAGE_NAME')`.
They are also presented in the `sessioninfo.txt` with their dependencies.

# How to use it

This code is mainly jupyter notebooks and dataset.
It is written in Python and R.

- Run the following command in a terminal.  `jupyter notebook # open a notebook`
- Open notebook - such as `messageAB_aids.ipynb`.
- For R kernel: Install required packages with the command line `install.packages('PACKAGE_NAME')` if it was not done previously.
- Run the cells.

# Structure

```
avatar_paper
|
└─── README.md
|
└─── datasets
        └─── AIDS: original and avatarized AIDS datasets
        └─── WBCD: original and avatarized WBCD datasets
        └─── results_df: computationally expensive analysis results.
        └─── messageD: dataset for figure 5 message
|
└─── notebooks
        └─── final_figure
                | analysis and graphs generation.
|
└─── lsg
        └─── dimension
                | Multidimensionnal functions such as PCA, MCA or FAMD
        └─── metrics
                | Multiple function used to compute avatarization metrics
|
└─── figure
        │ Figures presented in the article
|
└─── color.csv : A csv file that contains color for the figures
|
└─── requirements.txt : requirements for python packages.
|
└─── sessioninfo.txt : R session information.

```

# Contributing / use data

Feel free to redo our analyses or to explore avatarized datasets.

# Contributors

@Jpetot
@mguillaudeux
@oliviarousseau
@pagchun
