import numpy as np
import time
import random
from OOP1 import utilities
import matplotlib.pyplot as plt


# noinspection PyTypeChecker
class road(object):
    def __init__(self):
        self.length = 50
        self.number = 30
        self.clump_yesno = "y"
        self.max_size = int(10)
    def parameter_getter(self):
        utilities.time_delay("Please provide density of cars and length of road. Density is number of cars per X cells... will be rounded if not appropriate.")
        density = utilities.value_getter("Car density")
        length = utilities.value_getter("Road length")
        number = np.round(density * length)
        self.length = int(np.round(length))
        self.number = int(np.round(number))

        utilities.time_delay("Do you want the road to have clumps or not? y/n")
        t1 = str(input())
        if t1 == "y":
            self.clump_yesno = "y"
        else:
            utilities.time_delay("No clumps it is.")
    # Run this to get the number and length of the road. It sets self.length and self.number. Also gets clumps.
    def road_maker(self):
        road = [0 for d in range(self.length)]# Blank road.
        print(road)
        if self.clump_yesno == "n":
            number = self.number
            N = int(0)
            while True:
                random = np.random.randint(0, self.length)
                if road[random] == int(0):
                    road[random] = int(1)
                    N += int(1)
                else:
                    print("Populated.")
                if N == number:
                    break
        if self.clump_yesno == "y":
            while True:
                base = np.random.randint(0, self.length) # Start position of road
                size = np.random.randint(0, self.max_size)
                numberino = int(0)
                token = int(0)
                maxi = self.number

                ###-------------------------------
                for d in range(0, self.length):
                    if road[d] == int(1):
                        token+= int(1)
                if token == maxi:
                    break
                else:
                    token = int(0)
                if numberino >= size:
                    continue
                ###-------------------------------

                if road[base] == int(1):
                    print("Skipping generation point.")
                    continue
                else:
                    road[base] = int(1)
                    numberino += int(1)

                ###-------------------------------
                for d in range(0, self.length):
                    if road[d] == int(1):
                        token+= int(1)
                if token == maxi:
                    break
                else:
                    token = int(0)
                if numberino >= size:
                    continue
                ###-------------------------------

                for d in range(1, size):
                    if int(base + d) >= self.length:
                        road[base + d - self.length] = int(1)
                        numberino += int(1)
                        ###-------------------------------
                        for d in range(0, self.length):
                            if road[d] == int(1):
                                token += int(1)
                        if token == maxi:
                            break
                        else:
                            token = int(0)
                            pass

                        if numberino >= size:
                            break
                        else:
                            continue
                        ###-------------------------------
                    else:
                        road[base + d] = int(1)
                        numberino += int(1)
                        ###-------------------------------
                        for d in range(0, self.length):
                            if road[d] == int(1):
                                token += int(1)
                        if token == maxi:
                            break
                        else:
                            token = int(0)
                            pass

                        if numberino >= size:
                            break
                        else:
                            continue
                        ###-------------------------------
                ###-------------------------------
                for d in range(0, self.length):
                    if road[d] == int(1):
                        token+= int(1)
                if token == maxi:
                    break
                else:
                    token = int(0)
                    continue
                ###-------------------------------

        return road
    # Returns a road either with no clumps or with randomly placed clumps up to a size you decide (hard coded)
    def cell_iterator(self, road):
        newroad = [0 for d in road]
        motion_number = int(0)
        for d in range(0, self.length):
            if road[d] == int(1):
                if int(d + int(1)) == self.length:
                    if road[0] == int(1):
                        newroad[d] = int(1)
                    else:
                        if road[0] == int(0):
                            newroad[d] = int(0)
                            motion_number += int(1)
                else:
                    if road[d + int(1)] == int(1):
                        newroad[d] = int(1)
                    else:
                        if road[d + int(1)] == int(0):
                            newroad[d] = int(0)
                            motion_number += int(1)
            if road[d] == int(0):
                if d == int(0):
                    if road[self.length - int(1)] == int(1):
                        newroad[d] = int(1)
                        motion_number += int(1)
                    else:
                        if road[self.length - int(1)] == int(0):
                            newroad[d] = int(0)
                else:
                    if road[d - int(1)] == int(1):
                        newroad[d] = int(1)
                        motion_number += int(1)
                    else:
                        if road[d - int(1)] == int(0):
                            newroad[d] = int(0)
        return newroad, motion_number
    # Iterates the road by one step and returns [newroad, motion_number]
    def step_control(self):
        N = int(0)
        road = self.road_maker()

        while True:
            N += int(1)
            roader = self.cell_iterator(road)


            speed = roader[1]/self.number
            print(0.5*speed)
            road = roader[0]
            roadnew = np.expand_dims(road, axis=0)
            plt.imshow(roadnew)
            plt.show()
            time.sleep(int(1))
            if N > int(200):
                break
    # Basic step control to iterate over and print/control shit.



    # Iterates for the road and returns the newroad.

hey = road()
# hey.parameter_getter()
road = hey.road_maker()
print(road)
hey.step_control()