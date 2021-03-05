# Quadratic Solver Checkpoint
import math
import cmath
import time
import sys
import random

# This just lets you print stuff in a more natural manner...
def time_delay(t):
    print()
    for c in t:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform( 0.001, 0.002)
        time.sleep(numbertime)
# This is just to get you the variables you might want...
def value(name):
    value_text = ("Please give functions {0} coefficient...: ").format(name)
    time_delay(value_text)
    try:
        coefficient = float(input())
        return coefficient
    except:
        time_delay("Preferably a number... Rerun if you want to try again")
        print()
        quit()
# Next up is the solve. There are some cases to consider.
def quadratic_solver(a, b, c):
    if a == float(0):
        if b == float(0):
            time_delay("Sorry... no solution. Discriminant = 0")
        else:
            solution = -c/b
            discrim = b**2
            time_delay(("Your singular solution is {0:.2e} and your discriminant is {1:.2e}").format(solution, discrim))
    else:
        # Function is (-b +- root(b^2 - 4ac))/2a...
        discriminant = b**2 - float(4)*a*c
        solution_one = (-b + cmath.sqrt(discriminant))/(float(2)*a)
        solution_two = (-b - cmath.sqrt(discriminant))/(float(2)*a)
        print_text = ("Your solutions are {0:.2e} and {1:.2e}").format(solution_one, solution_two)
        time_delay(print_text)
# Okidoke. Now the actual program...
time_delay("Please provide the values of the coefficients of your quadratic! ax^2 + bx + c.")
a = value("a")
b = value("b")
c = value("c")
quadratic_solver(a, b, c)
