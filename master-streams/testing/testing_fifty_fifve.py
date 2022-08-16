import numpy as np

uuu = np.array([1,2,3,4,5,6,7,7,7,7])
vvv = np.where(uuu==7)[0]
print(vvv)