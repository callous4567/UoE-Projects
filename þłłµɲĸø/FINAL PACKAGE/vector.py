import numpy as np

# Vector object class, with associated operations. Call as vector(x, y, z) for cartesian coordinates x,y,z desired.
class vector():
    # Set arguments of vector object instance called, i.e. the x,y,z.
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    # Print vector in matrix representation [x,y,z]
    def __str__(self):
        return(str([self.x, self.y, self.z]))
    # Addition operator redefined for vector instances, on a component-by-component basis, i.e. new_x = (x - x_argument) for new_vector = (vector - vector_argument)
    def __add__(self, other):
        xn = self.x + other.x
        yn = self.y + other.y
        zn = self.z + other.z
        return vector(xn, yn, zn)
    # Subtraction operator redfined for vector instances.
    def __sub__(self, other):
        xn = self.x - other.x
        yn = self.y - other.y
        zn = self.z - other.z
        return vector(xn, yn, zn)
    # Scalar division of vector on a component-by-component basis. new_x = (x / argument) for new_vector = vector / argument.
    def __truediv__(self, other):
        xn = self.x / other
        yn = self.y / other
        zn = self.z / other
        return vector(xn, yn, zn)
    # Returns squared modulus... squares components x,y,z and sums them.
    def mag2(self):
        return (self.x**2 + self.y**2 + self.z**2)
    # Returns modulus... the positive root of the squared modulus.
    def mag(self):
        return np.sqrt(self.mag2())
    # Redefine multiplication operator.
    # Example: vector_ex * argument
    # If argument is scalar, isinstance detects scalar and applies scalar multiplication (i.e. [1,2,3] * 4 = [4,8,12]
    # If argument is a vector, does dot product.
    def __mul__(self, other):
        if isinstance(other, vector) == False:
            xn = self.x * other
            yn = self.y * other
            zn = self.z * other
            return vector(xn, yn, zn)
        else:
            xn = self.x * other.x
            yn = self.y * other.y
            zn = self.z * other.z
            return (xn + yn + zn)
    # Useful for iteration in vector comparison and allows subscription of vector objects in terms of simple list, i.e. vector[0] = vector.x... allows listing.
    def __getitem__(self, item):
        if item == int(0):
            return self.x
        if item == int(1):
            return self.y
        if item == int(2):
            return self.z
    # Cross Product in 3D only, using simple matrix formula that can be had via determinant of 3D matrix of [[ex, ey, ez], [self.x, self.y, self.z], [other...]] (you get idea)
    def cross(self, other):
        xn = self.y*other.z - self.z*other.y
        yn = self.z*other.x - self.x*other.z
        zn = self.x*other.y - self.y*other.x
        return vector(xn, yn, zn)
    # Vector comparison.
    # Forms an "index"
    # For each component of the two vectors being compared, is_equal takes this component, turns it into scientific notation for 4 significant figures, and then does a string comparison.
    # If the strings are equal, index increases by 1. If not, index does not increase and the for loop running this process ceases.
    # If final index equals 2, i.e. all three vector components are equal, you get a "True" and otherwise a "False"
    def is_equal(self, other):
        index = int(0)
        for d in range(0, 2):
            str_self = ("{0:.4E}").format(self[index])
            str_othe = ("{0:.4E}").format(other[index])
            if str_self == str_othe:
                index += int(1)
            else:
                break
        if index == int(2):
            return("True")
        else:
            return("False")

