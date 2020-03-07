"""
Return the square of a user-given number
"""

mynum = input("An integer number, please: ")

print("The square of your number is:")
print(mynum**2)

# Error is in casting mynum. Should be int or float (heck I'd even take a complex) if you want it to work as a number at all.
