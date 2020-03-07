import sys
import time
import random
import numpy as np

def time_delay(text):
    for c in text:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform(0.005, 0.01)
        time.sleep(numbertime)
    print()

document = open("file_name.file_extension","r") #r or w. r is read, w is write
# We can split the files
column_1 = []
column_2 = []
# For two columns.
line_array = document.readlines() # This gives you an array of strings, i.e. [(1,2,3),(1,2,3)...] that can be split up later. Basically all the lines of the document.
# We need to split up each line of the document: each line is inputted as a STRING. A STRING! So, you split up by a "token."
for line in line_array: # Use "Line" to show it's a line operation.
    # This shows that it's an operation that takes each element of the formed matrix as a line/string that is to be split up.
    line_tokens = line.split(",") # Splits the line, cutting out the , and taking each element split by the , as an element of the list line_tokens.
    # I.e. (hello,i,like,corn) goes to [hello, i, like, corn]. A list.
    column_1.append(line_tokens[0])
    column_2.append(line_tokens[1])
    # Appends to the rows.
    # If you specify float(line_tokens[1]), it takes that string, converts to a float, and appends that. If not a float, you get an error.
    # So yeah.
document.close() closes the document, needs to be done before you can carry on.


# document.readlines() reads the lines and forms an array of them which can be called as "line" for use in the token-splitting, using line.split(","), splitting each line string by the tokens specified.
# 
