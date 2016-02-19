import multiprocessing
from multiprocessing import Pool
import numpy
from matplotlib import pyplot as plt
import ctypes
from functools import partial



class gridholder(object):
    def __init__(self):
            grid = multiprocessing.Array(ctypes.c_double, 10*10)
            grid = numpy.ctypeslib.as_array(grid.get_obj())
            grid = grid.reshape(10,10)
            print type(grid)
            self.grid = grid
            #print grid

def func(n):
    global g
    
    print n
    for i in range(n):
        size = g.grid.shape
        c1 = numpy.random.randint(0, size[0])
        c2 = numpy.random.randint(0, size[1])
        print c1, c2
        g.grid[c1,c2] += 1.0
        #print "wrote a point"
    
def callbackfunc(result):
    global g
    
    g.grid += result

if __name__=="__main__":
    plt.clf()
    g = gridholder()
    results = []
    #func1 = partial(func, grid = g.grid)
    p = Pool(4)
    for i in range(1000):
        p.apply_async(func, args=(i,))
    
    p.close()
    p.join()
    
    
    plt.imshow(g.grid)
    plt.colorbar()
    plt.show()
    print g.grid.sum()
    print type(g.grid)
    
