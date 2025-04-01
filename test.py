#%%
from lablib102 import normalize_to_01
import numpy as np
import matplotlib.pyplot as plt

dataset = np.random.randn(101)
data_normed, _ = normalize_to_01(dataset)
plt.plot(data_normed)