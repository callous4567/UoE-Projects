from OOP1 import utilities
import numpy as np
import random

class twod_arrays(object):
    def __init__(self, rows, cols, decay_1, step):
        self.rows = int(0)
        self.cols = int(0)
        self.decay_1 = float(0)
        self.step = float(0)
    def array_former(self):
        newlist = [3 for d in range(int(self.cols))]
        newarraylist = [newlist for d in range(int(self.rows))]
        array = np.array(newarraylist)
        return array
    # array_former gives an empty matrix with rows ROWS and columns COLS, all initially 1 (for active element example)
    def probaratorOne(self, element):
        probability = self.decay_1 * self.step  # Of decay
        probCount = random.uniform(int(0), int(1))
        if element == int(1) and probCount > probability:
            return int(1)
        elif element == int(1) and probCount <= probability:
            return int(0)
        else:
            return int(0)
    # probaratorOne iterates for one step taking the initial value (1 or 0) and then for one step, returning the result of decay.
    def probaratorThree(self, element, decay_two, decay_three):
        if element == int(3):
            probability = decay_three * step
            probCount = random.uniform(int(0), int(1))
            if probCount > probability:
                return int(3)
            elif probCount <= probability:
                return int(2)
            else:
                print("Error code IX00CX")
        elif element == int(2):
            probability2 = decay_two * step
            probCount2 = random.uniform(int(0), int(1))
            if probCount2 > probability2:
                return int(2)
            elif probCount2 <= probability2:
                return int(1)
            else:
                print("Error code IX00CY")
        else:
            newvalue = self.probaratorOne(element)
            return newvalue
    # probaratorThree has two extra decays which must be specified, like the checkpoint. Uses probaratorOne for the final decay.
    def array_iterator(self, array):
        new_array = array
        for n in range(0, int(self.rows)):
            for m in range(0, int(self.cols)):
                new_element = twod_arrays.probaratorOne(self, array[n, m])
                array[n, m] = new_element
        return new_array
    # Iterates probaratorOne over one single step for an array that you input... or any method for that fact. Just copy the code for others or edit for extra levels.
    def array_iteratorThree(self, array, decay_two, decay_three):
        newarray = array
        for n in range(0, int(self.rows)):
            for m in range(0, int(self.cols)):
                value = newarray[n, m]
                newvalue = self.probaratorThree(int(value), decay_two, decay_three)
                newarray[n, m] = newvalue
        return newarray
    # Iterates except it uses probaratorThree instead of probaratorOne
    def step_control(self, array):
        time = 0
        stepnumber = 0
        truelife = np.log(2)/time
        total_N = self.rows * self.cols
        newarray = array
        while True:
            newarray = twod_arrays.array_iterator(self, newarray)
            N = 0
            time += self.step
            stepnumber += int(1)
            for n in range(0, int(self.rows)):
                for m in range(0, int(self.cols)):
                    if array[n, m] == int(0):
                        N += int(1)
                    else:
                        continue
            if N >= (1/2)*total_N:
                print("Total number has halved.")
                print(total_N)
                print(N)
                break
        print(array)
        print(newarray)
        print(truelife)
        print(time)
    # Takes array and iterates and such, prints the final array.
    def step_controlThree(self, array, decay_two, decay_three):
        initial_array = array
        N = 0
        time = 0
        utilities.time_delay("How long do we iterate for?")
        total_t = float(input())
        desired_N = total_t/self.step
        while True:
            time += step
            initial_array = twod_arrays.array_iteratorThree(self, initial_array, decay_two, decay_three)
            N += int(1)
            print(initial_array)
            if N >= desired_N:
                break
        print(initial_array)
        print(N)
        print(time)


rows = utilities.value_getter("number of rows")
cols = utilities.value_getter("number of columns")
decay = utilities.value_getter("decay constant")
step = utilities.value_getter("timestep")
decay_2 = utilities.value_getter("decay 2")
decay_3 = utilities.value_getter("decay 3")
arrayinstant = twod_arrays(rows, cols, decay, step)
arrayinstant.rows = rows
arrayinstant.cols = cols
arrayinstant.decay_1 = decay
arrayinstant.step = step
array = arrayinstant.array_former()
print(array)
newarray = arrayinstant.step_controlThree(array, decay_2, decay_3)
print(newarray)







