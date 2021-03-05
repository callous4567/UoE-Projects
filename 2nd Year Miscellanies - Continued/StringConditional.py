# Okidoke this will be a string conditional test...
# So, let's go ahead and test this out! <3 ^.^
print("Welcome. This will be a good test.")

import time
import sys
import random

# Two valid imports that are useful as heck


Name = str(input("What's your name, Young Padojuan? "))
print("This will be a hard test... you may not survive.")
print("Do you wish to continue?... ")
yesorno = str(input("y/n: "))

def main(yesorno):
    if yesorno == "y" :
        print("Fantastical! Let's continue...")
        songster(Name)
    elif yesorno == "n" :
        print("Okay. Well...")
        quit()
    else :
        print("Well... that wasn't really a valid answer. So screw you.")
        quit()

def songster(Name):
    song = ("Damn it's a rock lobster! You a damn rock lobster! Oh Baby!\
    What if I told you that {} Was a rocking god damn Lobster! Oh baby...")\
    .format(Name)
    for char in song:
        sys.stdout.write(char)
        sys.stdout.flush()
        numbertime = random.uniform( 0.05, 0.1)
        time.sleep(numbertime)
    print(song)

main(yesorno)
# Splendid that works great! Wonderful, SUBARASHIII! Okidoke. Let's try out the # next step that we must pass
# This is just a note though... For conditional statements, there is a space
# associated with the :, i.e. if a == b_:, the _ is a space.
# For a function, the : follows (), i.e. main():, same for for a_in_b:, :
# follows immediately... <3

