import pandas as pd
import pickle
import windows_directories

with open(windows_directories.duplimontedir + "\\" + "full_raw" + "\\" + "LE_157.txt", 'rb') as f:
    file = pickle.load(f)

print(file[0])