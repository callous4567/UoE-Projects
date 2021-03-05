import math

maximumN = 40
for i in range(0,maximumN):
    theta = 2.0*math.pi*i/maximumN
    value = math.cos(theta)
    print("Theta is: ",theta," and the cos value is",value)

print("LOOPED BABY!")

printloop = 6
for i in range( 0, printloop):
    spacing = "_"*i
    print(spacing)
print("Another loop!")

# Use i for the integer signification. Elem for element. So... example..

newList = []
print("Input List Variables")
for i in range(0,5):
    elementIn = float(input())
    elementOut = elementIn + float(3)
    newList.append(elementOut)
print(newList)

# That's how you form a list based on user input. If we need to now sort the list...
newList.sort()
print(newList)
    
