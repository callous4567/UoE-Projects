"""
Create a list of random integers, then delete the even entries,
then print the remaining entries.
"""

import random

# Create and print a list of 20 integers
# where all entries are 0<=n<=100
numbers = [random.randint(0,100) for n in range(20)]
print(("Original list: {}").format(numbers))


# Not sure what this is trying to do...

new_numbers = []
for i in range(len(numbers)):
    # if remainder of division by 2 is zero, delete list entry
    if numbers[i]%2 != 0:
        new_numbers.append(numbers[i])

print("Odd Entries")
print(new_numbers)

"""
numbers_new = []
for d in range(len(numbers)):
    frac = numbers[d]/int(2)
    string_frac = str(frac)
    split_frac = string_frac.split(".")
    if split_frac[1] != str(0):
        numbers_new.append(numbers[d])
    elif split_frac[1] == str(0):
        continue


# Print the remaining list with all odd numbers
print("Odd entries:")
print(numbers_new)
"""