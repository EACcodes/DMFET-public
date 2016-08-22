#!/usr/bin/python
import sys
from numpy import *

ifile = file( sys.argv[1], 'r' )
ndim = int( sys.argv[2] )

data = array([[0.0 for i in range(ndim)] for j in range(ndim)])

i = 0
j = 0
for line in ifile:
    words = line.split()
    for word in words:
        data[i,j] = float(word)
        j += 1
        if j%ndim == 0:
            j = 0
            i += 1

for i in range(ndim):
    for j in range(ndim):
        print '%5d%5d%20.8e'%(i,j,data[i,j])
    print ''
