import pandas as pd 
import numpy as np
from sdv.tabular import CTGAN
import time

NB_ITERATION = 100

if __name__ == "__main__":
    aids = pd.read_csv(
        "datasets/AIDS/aids_original_data.csv", sep=";"
    )
    aids = aids.drop(["pidnum"], axis=1)
    categorical_val = []
    continous_val = []
    for column in aids.columns:
        if len(aids[column].unique()) <= 10:
            categorical_val.append(column)
        else:
            continous_val.append(column)

    aids[categorical_val] = aids[categorical_val].astype("object")

    times = []
    for repeate in range(NB_ITERATION):
        start = time.time()

        model = CTGAN(epochs=1000, batch_size=600, field_transformers={
            'integer': 'age',
            'float': 'wtkg',
            'categorical': 'hemo',
            'categorical': 'homo',
            'categorical': 'drugs',
            'categorical': 'karnof',
            'categorical': 'oprior',
            'categorical': 'z30',
            'categorical': 'zprior',
            'integer': 'preanti',
            'categorical': 'race',
            'categorical': 'gender',
            'categorical': 'str2',
            'categorical': 'strat',
            'categorical': 'symptom',
            'categorical': 'treat',
            'categorical': 'offtrt',
            'integer': 'cd40',
            'integer': 'cd420',
            'float': 'cd496',
            'categorical': 'r',
            'categorical': 'treat',
            'integer': 'cd80',
            'integer': 'cd820',
            'categorical': 'cens',
            'integer': 'days',
            'categorical': 'arms',})
        model.fit(aids)
        model.sample(num_rows=(aids.shape[0]))
        duration = time.time() - start
        times.append(duration)

    print("duration for ", NB_ITERATION, " CTGAN : ", sum(times))
    print("mean : ", np.mean(times))
    print("std : ", np.std(times))