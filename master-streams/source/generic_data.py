import numpy as np

# Jorge Gauss Mix Paper
mus = np.array([[426, -4950, 1436],
                [-0.7, 11, 3],
                [144, -278, 122],
                [-3905, -2322, -4664]])
covtrix_1 = [[655**2, 0, 0],
             [0, 1255**2, 0],
             [0, 0, 659**2]]
covtrix_2 = [[11**2, 0, 0],
             [0, 12**2, 0],
             [0, 0, 11**2]]
covtrix_3 = [[1755**2, 0, 0],
             [0, 1926**2, 0],
             [0, 0, 1733**2]]
covtrix_4 = [[13**2, 0, 0],
             [0, 12**2, 0],
             [0, 0, 13**2]]
covtrices = np.array([covtrix_1, covtrix_2, covtrix_3, covtrix_4])