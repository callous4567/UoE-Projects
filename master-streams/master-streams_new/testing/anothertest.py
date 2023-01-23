from codecs import decode
import struct
import numpy as np
test = "adc_data_roundrobin_test_uint16.txt"
def bits(f):
    bytes = (ord(b) for b in f.read())
    for b in bytes:
        for i in range(8):
            yield (b >> i) & 1

u = ""
for b in bits(open(test, 'r')):
    u += str(b)
print(u[::-1])
print(len(u))
