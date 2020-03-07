f = open("test.txt", 'w+')
for i in range(10):
    f.write("Line number " + str(i) +"..., \n")
f.close()