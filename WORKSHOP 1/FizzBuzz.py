from OOP1 import utilities
import numpy

def main():
    newvalue = int(utilities.value_getter("integer value")) # Gets the value that you want using utilities
    for d in range(1, newvalue + int(1)): # Runs a chain of if and elif arguments (with an else at the end) to decide Fizz (multiple 3) Buzz (multiple 5) or Fizzbuzz (multiple 3 and 5). Prints int if neither 3 satisfied.
        value = d
        if value%15 == 0:
            print("Fizzbuzz")
        elif value%3 == 0:
            print("Fizz")
        elif value%5 == 0:
            print("Buzz")
        else:
            print(value)

main() # Execute function