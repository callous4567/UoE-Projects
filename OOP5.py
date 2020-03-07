from OOP1 import utilities
import numpy as np
import matplotlib.pyplot as plt

class mandelbrotter(object):
    def __init__(self):
        self.xMin = float(-2.025)
        self.xMax = float(0.6)
        self.yMin = -1.125
        self.yMax = 1.125
        self.iteration_limit = int(350) # Full iterations excluding zeroth
        self.matrix_parameter = int(600) # Defines matrix length and height
        self.threshold = int(2)
        self.z_0 = 0
    # Default parameters from the website. Choose parameter or whatever you want, asshole. Slap me harder.

    def array_former(self):
        x_range = np.linspace(self.xMin, self.xMax, self.matrix_parameter)
        yi_range = np.linspace(self.yMin, self.yMax, self.matrix_parameter)
        array = utilities.array_formerComplex("null argument", self.matrix_parameter, self.matrix_parameter)
        for n in range(0, self.matrix_parameter):
            for m in range(0, self.matrix_parameter):
                hey = complex(x_range[m], yi_range[int(self.matrix_parameter) - int(1) - n])
                print(hey)
                array[n,m] = hey
        return array
    # This uses the complex array former: np.zeros((rows, cols), dtype=np.complex_) to generate our array. Utilities has it.
    # np.zeros can return an array. Use this in future. np.zeros((dimension rows,cols), dtype = optional for datatype)

    def iterator(self, array):
        newarray = np.zeros((self.matrix_parameter - int(1), self.matrix_parameter-int(1)), dtype=int)
        funct = lambda zn, C: zn**2 + C
        for n in range(0, self.matrix_parameter - int(1)):
            for m in range(0, self.matrix_parameter - int(1)):
                N = 0
                z = self.z_0
                C = array[n, m]
                while True:
                    z = funct(z, C)
                    N += int(1)
                    if np.abs(z) > int(2):
                        newarray[n, m] = N
                        break
                    if N>self.iteration_limit:
                        break
                    else:
                        continue
                continue
        return newarray
    def jul_iterator(self, array):
        newarray = np.zeros((self.matrix_parameter - int(1), self.matrix_parameter-int(1)), dtype=int)
        funct = lambda zn, C: zn**2 + C
        for n in range(0, self.matrix_parameter - int(1)):
            for m in range(0, self.matrix_parameter - int(1)):
                N = 0
                z = array[n,m]
                C = complex(0, -1)
                while True:
                    z = funct(z, C)
                    N += int(1)
                    if np.abs(z) > int(2):
                        newarray[n, m] = N
                        break
                    if N>self.iteration_limit:
                        break
                    else:
                        continue
                continue
        return newarray



hello = mandelbrotter()
hello_array = hello.array_former()
print(hello_array)
newarray = hello.jul_iterator(hello_array)
print(newarray)
plt.matshow(newarray)
plt.show()
print("done")
