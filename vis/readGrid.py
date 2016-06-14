import sys
import numpy as np

def readGrid(filename):
    ext = filename.split('.')[-1]
    if ext == "txt":
        return readGridASCII(filename)
    else:
        print("Unknown extension: {0:s}\nTrying ASCII...".format(ext))
        return readGridASCII(filename)

def readGridASCII(filename):

    f = open(filename, "r")

    #Version String
    version = f.readline()

    #Time
    line = f.readline()
    t = float(line.split()[-1])

    #Faces
    line = f.readline()
    xf = np.array([float(x) for x in line.split()[1:]])
    x = 0.5*(xf[:-1]+xf[1:])
    f.close()

    rho, P, v1, v2 = np.loadtxt(filename, unpack=True, usecols=[1,2,3,4],
                                skiprows=3)

    return version, t, xf, x, rho, P, v1, v2

def readParfile(filename):
# Loads the parameter file into a dictionary for easy use in plotting scripts.

    f = open(filename, "r")

    pars = {}
    for line in f:
        words = line.split()
        if len(words) < 2:
            continue

        key = words[0]
        sval = words[1]

        if key[0] == '/':
            continue
        try:
            val = int(sval)
        except ValueError:
            try:
                val = float(sval)
            except ValueError:
                val = sval
        pars[key] = val

    f.close()

    return pars

