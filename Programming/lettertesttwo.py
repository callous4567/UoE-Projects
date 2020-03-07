# First let's import math.

import math
import cmath
import random
import sys
import time

def timeDelay(t):
    for c in t:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform( 0.01, 0.1)
        time.sleep(numbertime)
    print(t)

oya = "Oya user! This is the volume to radius converter! ^.^"
timeDelay(oya)




