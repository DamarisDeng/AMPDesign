from deepforest import CascadeForestRegressor
import numpy as np
import pandas as pd

Ec_test_feature = pd.read_csv("/home/yaolt/pythoncode/AMPActiPred/data/regression/test/Ec.csv").iloc[:,2:].values
Ec_test_value = pd.read_csv("/home/yaolt/pythoncode/AMPActiPred/data/regression/test/Ec.csv")["value"].values
model_Ec = CascadeForestRegressor()
model_Ec.load("/home/yaolt/pythoncode/AMPActiPred/DeepForest-AMP-Regression/model/Ec/")
y_predict_Ec = model_Ec.predict(Ec_test_feature).squeeze()
