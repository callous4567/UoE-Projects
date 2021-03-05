from OOP1 import utilities
import numpy as np
import random
import time
import sys

class poly(object):
    def __init__(self, coef):
        self.coef = []

class poly_taker(object):
    def indi_splitter(self, full_string):
        polystringg = full_string.split(" ")
        # Returns [4x^2, + , 2x , - , 4]

        newstring = []
        index = 0
        while index < len(polystringg):
            print("hello")
            if polystringg[index] == "+":
                index = index + 1
                continue
            elif polystringg[index] == "-":
                index = index + 1
                continue
            else:
                newstring.append(polystringg[index])
                index = index + 1
                continue

        return newstring
    # This splits up the full input string of 4x^2 + ...

    def power_splitter(self, polystringer):
        powerstring = polystringer.split("^")
        if powerstring[1] == "0":
            newpowerstring = [powerstring[0], "^0"]
            return newpowerstring
        else:
            return powerstring
    # This splits each element of the list up between the power and etc, ie. 4x^2 ->->-> [4x, 2]

    def coef_splitter(self, powerstringg):
        coefstring = powerstringg.split("x")
        if len(coefstring) == int(1):
            coefstring.append(int(0))
        return coefstring
    # This splits to the 4x into [4, " "]


    def input(self):
        inputtext = "Please print polynomial, i.e. 4x^2 + 2"
        utilities.time_delay(inputtext)
        fullstring = str(input())
        return fullstring

def main():
    null = " "
    fullstring = poly_taker.input("self")
    polystring = poly_taker.indi_splitter("self", fullstring) # Now we have [4x^2, 2x] with example no zeroth coefficient
    powerstring = [poly_taker.coef_splitter("self", d)[1] for d in polystring]
    #powerstringg = [poly_taker.power_splitter("self", d)[1] for d in powerstring]
    print(powerstring)
    #print(powerstringg)


    """
    coef_list = [0 for d in range(powerMax[1])]
    for d in range(0, len(polystring) - int(1)):
        coefstring = poly_taker.coef_splitter("self", polystring[d])
        coef_list[coefstring[1]] = coefstring[0]
    print(coef_list)
    """








main()


