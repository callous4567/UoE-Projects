"""
Print out a series of integers
"""

nmax = 10
for n in range(nmax + 1):
  print(("The number is {}").format(int(n)))
# This assumes you actually want the value of 0... if you don't, run:

"""
from OOP1 import utilities
import time
for d in range(1, nmax + int(1)):
  utilities.time_delay(str(int(d)))
  time.sleep(0.69)
  # The sleep just makes it more annoying to run >.<
"""