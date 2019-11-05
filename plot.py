import numpy as np
import matplotlib.pyplot as pp

for i in range(0, 6):
    txt_str = "/home/dcrush/CLionProjects/mpi_test/cmake-build-release/data{}.txt".format(i)
    fmt_str = "/home/dcrush/CLionProjects/mpi_test/cmake-build-release/data{}.fmt".format(i)
    shape = np.loadtxt(fmt_str).astype(int)
    half = int(shape[0]/2)
    data = np.loadtxt(txt_str).reshape(shape)
    plane = data[half]
    del data
    if i == 5:
        line = plane[int(shape[1]/2)]
        pp.plot(line)
    else:
        pp.pcolor(plane)
    pp.title("data{}".format(i))
    pp.show()
