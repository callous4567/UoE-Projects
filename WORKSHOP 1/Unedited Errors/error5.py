"""
Create a list of random integers, then delete the even entries,
then print the remaining entries.
"""

import random

# Create and print a list of 20 integers
# where all entries are 0<=n<=100
numbers = [random.randint(0,100) for n in range(20)]
print("Original list:")
print(numbers)

for i in range(len(numbers)):
    # if remainder of division by 2 is zero, delete list entry
    if numbers[i]%2 == 0:
        del numbers[i]

# Print the remaining list with all odd numbers
print("Odd entries:")
print(numbers)

# No space when you're doing the if argument... plus that dreadful error whereby if you run a "d in range" loop you end up removing an element,
"""
example:
numbers = [1,3,4,5,6,7]
for i in range(len(numbers)): 
    if numbers[i]%2 == int(0):
        del numbers[i]
        
i = 0
list[0] = 1 and hence you trim this element.
[3,4,5,6,7]
i = 1
list[1] = 4 and hence you continue
i = 2
list[2] = 5 and hence you trim
[3, 4, 6, 7]
i = 3
list[3] = 7 and hence you trim
[3, 4, 6]
i = 4 but list is only 3 long, you're ****ed. The original "i in range" thing was passed an argument len(numbers) using the INITIAL length of the list.
This initial length of the list far exceeded that pathetic little 4 that you're trying to force down its throat,
Thus you get an indexerror or something like that. 
"""