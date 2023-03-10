import pandas as pd 
from lsg import avatarize
import time 
import numpy as np 

NB_ITERATION = 100

if __name__ == "__main__":
    aids = pd.read_csv("../../paper/avatar-paper/datasets/AIDS/aids_original_data.csv", sep=";")
    categorical_val = []
    continous_val = []
    for column in aids.columns:
        if len(aids[column].unique()) <= 10:
            categorical_val.append(column)
        else:
            continous_val.append(column)

    aids[categorical_val] = aids[categorical_val].astype("category")

    aids = aids.drop(["pidnum"], axis=1)

    times = []
    for repeate in range(NB_ITERATION):
        start = time.time()
        avatarize(
            aids, k=20, ncp=5, drop_duplicates=False, distance_metric="minkowski"
        )
        duration = time.time() - start
        times.append(duration)

    print("duration for ", NB_ITERATION, " avatarizations : ", sum(times))
    print("mean : ", np.mean(times))
    print("std : ", np.std(times))