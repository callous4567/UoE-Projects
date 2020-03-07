import numpy as np
import matplotlib.pyplot as plt
import time
import random
import math
import sys
import scipy

def time_delay(text):
    for c in text:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform(0.001, 0.005)
        time.sleep(numbertime)
    print()
# Typical <3

def file_reader(filedesired):
    file = open(filedesired, "r")
    total_list = []
    total_list_2 = []
    time_delay("Which line does the data stream begin?")
    linestart = float(input()) - 1
    for line in file:
        newline = line.split(" ")
        if len(newline) == 3:
            total_list.append(float(newline[0]))
            total_list_2.append(float(newline[2]))
        else:
            print("Hah this is bullshit.")
    file.close()
    return total_list, total_list_2
# 25 KHz with the first list being the voltages, the second list being the currents.


def main():
    time_delay("What is the filename and extension? Try writing hello.txt or fuckyou.pdf<3")
    filetotal = str(input())
    lists = file_reader(filetotal)
    voltage = lists[0]
    current = lists[1]
    frequency = 25000
    interval = 1/frequency
    timelist = np.arange(0, len(lists[0])*interval, interval)
    p_list = []
    for d in range(0, len(lists[0])):
        p = np.log((voltage[d])*(current[d]))
        p_list.append(p)
    # p_list is p(t) for timelist
    plt.plot(timelist, p_list)
    plt.show()




main()
