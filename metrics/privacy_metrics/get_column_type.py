import pandas as pd

def get_categorical_continuous(data: pd.DataFrame):
    categorical_val = []
    continous_val = []
    for column in data.columns:
        if len(data[column].unique()) <= 10:
            categorical_val.append(column)
        else:
            continous_val.append(column)
    
    return categorical_val, continous_val
