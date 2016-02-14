import multiprocessing
from multiprocessing import Pool
import numpy
from matplotlib import pyplot as plt
import ctypes



class gridholder(object):
    def __init__(self):
            grid = multiprocessing.Array(ctypes.c_double, 10*10)
            grid = numpy.ctypeslib.as_array(grid.get_obj())
            grid = grid.reshape(10,10)
            self.grid = grid
            print grid

def func(grid, n):
    #print inp
    #grid = inp[0]
    #n = inp[1]
    print n
    for i in range(n):
        size = grid.shape
        c1 = numpy.random.randint(0, size[0])
        c2 = numpy.random.randint(0, size[1])
        grid[c1,c2] += 1.0
        print "wrote a point"


if __name__=="__main__":
    plt.clf()
    g = gridholder()

    p = Pool(4)
    func( g.grid, 100000)
    p.close()
    p.join()
    plt.imshow(g.grid)
    plt.colorbar()
    plt.show()
    print g.grid.sum()
    
