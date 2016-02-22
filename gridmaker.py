import numpy
from multiprocessing import Pool

from otfregrid import convolve_mp
from otfregrid import LMTOTFRegrid_mp

def gridmaker(xmin, xmax, ymin, ymax, filelist):
    ### going to move the 'make_grid' function outside of the class
    ### for multiprocessing purposes? we'll see if it works
    g = LMTOTFRegrid_mp(xmin, xmax, ymin, ymax, filelist)
    print "Starting Grid Making Process..."
    T = g.T
    TSYS = g.TSYS
    WT = g.WT
    

if __name__ == "__main__":
    files = ['35065305.nc']#,'35065304.nc','35065306.nc','35065300.nc','35065301.nc','35065302.nc','35065303.nc','35065307.nc','35065308.nc','35065309.nc']
    gridholder = gridmaker(525.0,535.0,645.0,655.0, files)
