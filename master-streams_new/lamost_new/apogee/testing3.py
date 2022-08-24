import numpy as np
u = np.array([1,2,3])
b = np.array([4,5,6])
c = np.array([7,8,9])
print(np.concatenate([u,b,c], axis=0))