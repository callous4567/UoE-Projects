"""
Check whether a given number 
appears as entry in a list of integers.
"""

mynum = 3
numbers = list(range(10))
print(numbers)

for n in numbers:
    if n == mynum:
        print(("Found my number: {}").format(mynum))
