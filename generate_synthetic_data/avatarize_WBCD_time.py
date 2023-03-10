import pandas as pd 
from lsg import avatarize
import time 
import numpy as np 

NB_ITERATION = 100

if __name__ == "__main__":
    wbcd = pd.read_csv("datasets/WBCD/breast_cancer_wisconsin.csv", sep=",")
    wbcd = wbcd.drop(columns="Sample_code_number")

    wbcd = wbcd.astype("category")
    times = []
    for repeate in range(NB_ITERATION):
        start = time.time()
        avatarize(wbcd, k=20, ncp=2, drop_duplicates=False)
        duration = time.time() - start
        times.append(duration)

    print("duration for ", NB_ITERATION, " avatarizations : ", sum(times))
    print("mean : ", np.mean(times))
    print("std : ", np.std(times))