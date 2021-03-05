from OOP1 import utilities
import numpy as np

class polynomial:
    def __init__(self, coefs, powmax):
        self.coefs = []
        self.powmax = []

class polyTaker(object):
    def __init__(self, coefs, powmax):
        self.coefs = []
        self.powmax = []
    def fullstringTaker(self):
        utilities.time_delay("Print polynomial. 4x^2 + 2x + 9.")
        fullstring = str(input())
        return fullstring
    def spaceSplitter(self, fullstring):
        spacestring = fullstring.split(" ")

        newspacestring = []
        newspacestring.append(spacestring[0])
        index = 1
        while index < len(spacestring):
            part1 = spacestring[index]
            part2 = spacestring[index + int(1)]
            part3 = str("{}{}").format(part1, part2)
            newspacestring.append(part3)
            index = index + 2
        return newspacestring
    def powerSplitter(self, spacestring_element):
        powerstring = spacestring_element.split("x^")
        if len(powerstring) == 1:
            x_split = powerstring[0].split("x")
            if len(x_split) == 1:
                powerstring.append(int(0))
            if len(x_split) == 2:
                powerstring = [x_split[0], int(1)]
        return powerstring
    def poly_writer(self):
        fullstring = polyTaker.fullstringTaker(self)
        polystring = polyTaker.spaceSplitter(self, fullstring)

        powerstring = []
        for d in polystring:
            coef_power = polyTaker.powerSplitter(self, d)
            print(coef_power)
            powerstring.append(int(coef_power[1]))
        powermax = max(powerstring)

        coefstring = [0 for d in range(int(powermax) + int(1))]
        for d in polystring:
            coef_power = polyTaker.powerSplitter(self, d)
            coefstring[int(coef_power[1])] = float(coef_power[0])
        return coefstring, powermax
    def array_adder(self, arr1, arr2):
        arr3 = []
        index = 0
        print(arr1)
        print(arr2)
        try:
            while index < len(arr1):
                value = float(arr1[index]) + float(arr2[index])
                arr3.append(value)
                index = index + 1
            return arr3
        except:
            print("Error code #00X4C")
    def poly_adder(self, pol1, pol2):
        if len(pol1) > len(pol2):
            dif = len(pol1) - len(pol2)
            poly2 = pol2
            index = 0
            while index < dif:
                poly2.append(0)
                index = index + 1
            pol3 = polyTaker.array_adder(self, poly2, pol1)
            return pol3
        elif len(pol1) == len(pol2):
            pol3 = polyTaker.array_adder(self, pol1, pol2)
            return pol3
        else:
            dif = len(pol2) - len(pol1)
            poly1 = pol1
            index = 0
            while index < dif:
                poly1.append(0)
                index = index + 1
            pol3 = polyTaker.array_adder(self, poly1, pol2)
            return pol3
    def derivative(self, pol):
        index = 1
        newlist = []
        while index < len(pol):
            valued = pol[index] * index
            newlist.append(valued)
            index = index + 1
        return newlist
    def integral(self, pol, c):
        index = 0
        newpol = [0 for d in range(len(pol) + int(1))]
        while index < len(pol):
            newpol[0] = c
            value = pol[index]/(index + int(1))
            newpol[index + 1] = value
            index = index + 1
        return newpol
    def polyPrinter(self, pol):
        list_0 = str(pol[0]) + " + "
        list_1 = str(pol[1]) + "x"
        printstring = list_0 + list_1
        index = int(2)
        while index < len(pol):
            stringer = str("{}x^{}").format(pol[index], index)
            printstring += " + " + stringer
            index += int(1)
        print(printstring)
    # First one returns fullstring
    # Second one returns split based on spaces and removes the + and the -
    # Third one returns [coefficient, power] including cases for the constant c.
    # Fourth one returns [coefficient list, powermax], with the list for the powers respectively. I.e. [4,2,1] = 4 + 2x^1 + 1x^2... !!!!!!!!!!!!!!!!!!!!!!!!!! [coef_list , AND , powermax]
    # Fifth one adds two same-sized arrays in float
    # Sixth one adds two polynomials of different size
    # Seventh one prints the polynomial based on coef list, i.e. [0 2 3 4] = 0 + 2x^2 + 3x^3 + 4x^4. The [0 0 2] shows up as 0 + 0x + 0x^2.