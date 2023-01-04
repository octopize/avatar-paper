
import pandas as pd 
import numpy as np
from sdv.tabular import CTGAN
import time

NB_ITERATION = 100

if __name__ == "__main__":
    aids = pd.read_csv(
    "datasets/WBCD/breast_cancer_wisconsin.csv", sep=","
)
    wbcd = wbcd.drop(["Sample_code_number"], axis=1)
    wbcd = wbcd.astype('object')
    wbcd['Class'] = wbcd['Class'].replace({2: 0, 4: 1})
    wbcd['Class'] = wbcd['Class'].astype('bool')

    times = []
    for repeate in range(NB_ITERATION):
        start = time.time()
        model = CTGAN(epochs=1000, batch_size=600, field_transformers={
                'integer': 'Clump_Thickness',
                'integer': 'Uniformity_of_Cell_Size',
                'integer': 'Uniformity_of_Cell_Shape',
                'integer': 'Marginal_Adhesion',
                'integer': 'Single_Epithelial_Cell_Size',
                'integer': 'Bare_Nuclei',
                'integer': 'Bland_Chromatin',
                'integer': 'Normal_Nucleoli',
                'integer': 'Mitoses',
                'boolean': 'Class',})
        model.fit(wbcd)

        ctgan_data = model.sample(num_rows=(wbcd.shape[0]))

        duration = time.time() - start
        times.append(duration)

    print("duration for ", NB_ITERATION, " CTGAN : ", sum(times))
    print("mean : ", np.mean(times))
    print("std : ", np.std(times))