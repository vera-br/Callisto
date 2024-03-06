import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

array = [[0,0,0], [1,1,1], [2,2,2], [3,3,3]]
rolled = np.roll(array, 1, axis=0)
print(rolled)