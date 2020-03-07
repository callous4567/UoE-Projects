# This will test out the listing system and printing all the lists, and multiplying values of a list.
# So, let's go.


def getlistValues(valueNumber):
    print(("Provide a value for list value {0:d}").format(valueNumber + int(1)))
    the_value = float(input())
    return(the_value)
print("Please define the values of your six valued list")
zero = getlistValues(0)
one = getlistValues(1)
two = getlistValues(2)
three = getlistValues(3)
four = getlistValues(4)
five = getlistValues(5)
listValues = [zero,one,two,three,four,five]
print(listValues)
print(("And, the length of our value list is {}".format(str(len(listValues)))))
# Let's go for another quick test.
# Sorting the variables...
# So, variable sorting!

newList = listValues
def sortingFunction(list):
    list.sort()
    x = float(0)
    list.append(x)
    print(list)

sortingFunction(newList)

# So, list.sort() automatically sorts and updates the list variable, list.append() adds the variable
# within the brackets automatically and appends it.
# To select a value from our list...
