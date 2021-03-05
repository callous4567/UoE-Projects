"""
Okidoke. This will be weirder. Messy too...
Let's let A = coefficient of X^2, B for X, and C for the constant.
Then just set up equations to solve for the roots. Should be able to handle complex  variables too.
"""

import math
import cmath

print("Heyo! I'll need some variables.")
A = float(input("Provide coefficient of x^2 please:"))
B = float(input("Provide coefficient of x please:"))
C = float(input("Also provide the constant, signed (+-):"))

def discriminant(a, b, c):
    return b**2 - float(4)*a*c

# There's the function to get the discriminant... now... quickie issue.
# Need to validate whether or not the variables would work, i.e...
# If a = 0, we get strange division, etc... So... I'll build the solver first.

def quadraticSolver(disc, a, b, c) :
    firstRoot = (-b + cmath.sqrt(disc))/(float(2)*a)
    secondRoot = (-b - cmath.sqrt(disc))/(float(2)*a)
    # These two determine the actual roots!
    realFirst = firstRoot.real
    realSecond = secondRoot.real
    imagFirst = float(firstRoot.imag)
    imagSecond = float(secondRoot.imag)
    # Sets up the variables for the comparitive part/formatting/etc
    # Exponential formatting lets us calculate any kind of scale, big or small,
    # hence why I find it preferential.
    if imagFirst == 0 :
        print(("Your first root for x will be {0:.2e}").format(realFirst))
    else :
        if realFirst == 0 :
            print(("Your first root is {0:.2e}{1}").format(imagFirst, "j"))
        else :
            print(("Your first root is {0:.2e}").format(firstRoot))
    if imagSecond == 0 :
        print(("Your second root for x will be {0:.2e}").format(realSecond))
    else :
        if realSecond == 0 :
            print(("Your second root is {0:.2e}{1}").format(imagSecond, "j"))
        else :
            print(("Your second root is {0:.2e}").format(secondRoot))
    # This rather bulky section of comparison makes it easier on the eyes for
    # the final solution, i.e. giving just the real or imaginary component when
    # the other component is equal to zero.
    
def validityChecker(disc, e, f, g) :
    disc = discriminant(A, B, C)
    a = A
    b = B
    c = C
    if a == 0 :
        if c == 0 :
            print("Sorry, no ax^2 and no constant? x = 0.")
        else :
            solution = float(-1)*c/b
            print(("Your solution for x will be {0:.2e}").format(solution))
    else :
        quadraticSolver(disc, a, b, c)
    # This component/function pushes out the non-quadratic solution, i.e. if
    # there is a solution that would invalidate the function I used above via
    # infinite division (1/0) or there exist no non-trivial solutions
    
validityChecker(discriminant(A, B, C), A, B, C)
# This just runs and starts the program. The checker automatically calls up
# the actual solver once it figures out whether non-trivial solutions are even
# possible.
