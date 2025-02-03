import numpy as np
from pysr import PySRRegressor
import matplotlib.pyplot as plt

data = np.load('depth_data.npy', allow_pickle=True).item()

X_train = np.vstack([data['accret_rate'], data['z']]).T
y_train = data['depth']
y_err_train = np.mean(data['depth_err'], axis=0)

model_sr = PySRRegressor(
    niterations=1000,
    binary_operators=["+", "-", "*", "/", "^"],
    unary_operators=["sin", "cos", "exp", "log"],
    populations=30,  
    elementwise_loss="loss(x, y, y_err) = sum(((x - y) / y_err)^2) / length(x)"
)

model_sr.fit(X_train, y_train, weights=1 / y_err_train**2) 
print(model_sr.equation_)