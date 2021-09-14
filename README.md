# Description 

Avatar paper is an analysis project that allows researchers to see, re do, or explore figure and data of the "PAPER_NAME".
The aim of the paper is to present the avatarization method.


We really recommand to have a look to the scientific paper (link) before to explore the repository.

# How to use it

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
there is alos a `requirements.txt`.

## R packages :   
Install all required packages on the top of each R notebooks with the command `install.packages('PACKAGE_NAME')`, nothing special.


# Strucutre  

```
avatar_paper
|
└─── README.md
|
└─── datasets
        └─── AIDS : all aids datasets, original and avatarized 
        └─── WBCD : all WBCD datasets, original and avatarized 
        └─── results_df : dataframe of the results of computationally expensive analysis.
|
└─── notebooks 
        └─── final_figure  
                |
                | Folder that contains all the notebook that analyse and generate the article graphs.
                |

        └─── avatarization 
                |
                | Folder that contains all the notebooks that have avatarised the AIDS and WBCD data set. These notebooks
                | require the avatarisation package, which is a private package of octopize's own.
                
|
└─── dimension
        |
        | Multidimensionnal function such as PCA, MCA or FAMD
        |
|
└─── metrics
        |
        | Multiple function used to compute avatarization metrics 
        |
|
└─── figure  : folder that contains all the figure of the article
        │   
        |   multiple figures presented in the article
        |
|
└─── color.csv : A csv file that contains color for the figures
|
└─── requirements.txt : requirements that contains python package of the project.
```


# Contributing / use data 

Feel free to redo our analyses or to explore avatarized datasets.

# Contributors 

@Jpetot
@mguillaudeux
@oliviarousseau
@pagchun
