# Description 

Avatar paper is an analysis project that allows researchers to see, re do, or explore figure and data of the "PAPER_NAME".
The aim of the paper is to present the avatarization method.


We really recommand to have a look to the scientific paper (link) before to explore the repository.

# How to use it (a developper)

This code is mainly jupyter notebooks and dataset.
It is written in Python and R. 

## Python package : 
```
matplotlib==3.4.2
matplotlib-inline==0.1.2
numpy==1.20.3
pandas==1.2.4
scikit-learn==0.24.2
scipy==1.6.3
seaborn==0.11.1
```
there is also a `requirements.txt`.

## R packages :   
Install all required packages on the top of each R notebooks with the command `install.packages('PACKAGE_NAME')`, nothing special.

# Strucutre  

```
avatar_paper
|
└─── README.md
|
└─── datasets
        └─── AIDS : all AIDS datasets, original and avatarized 
        └─── WBCD : all WBCD datasets, original and avatarized 
        └─── results_df : result dataframes of computationally expensive analysis.
|
└─── notebooks 
        └─── final_figure  
                | analysis and graphs generation.               
|
└─── dimension
        | Multidimensionnal functions such as PCA, MCA or FAMD
|
└─── metrics
        | Multiple function used to compute avatarization metrics 
|
└─── figure
        │ Figures presented in the article
|
└─── color.csv : A csv file that contains color for the figures
|
└─── requirements.txt : requirements for python packages.
```


# Contributing / use data 

Feel free to redo our analyses or to explore avatarized datasets.

# Contributors 

@Jpetot
@mguillaudeux
@oliviarousseau
@pagchun
