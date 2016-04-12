import numpy
from multiprocessing import Pool

from otfregrid import convolve_mp
from otfregrid import LMTOTFRegrid_mp
from file_compatibility.LMTOTFFile import LMTOTFNetCDFFile
from convolve_cy import convolve_cy

# this version of gridmaker is designed to handle a large-ish (1-30?)
# filelist and make a grid using multiprocessing on one computer.
# As it stands currently, it can process approximately 4 files per 12
# minute increment.

def gridmaker_1c(xmin, xmax, ymin, ymax, filelist, cython=True):
    ### going to move the 'make_grid' function outside of the class
    ### for multiprocessing purposes? we'll see if it works
    #g = LMTOTFRegrid_mp(xmin, xmax, ymin, ymax, filelist)
    initialize_regrid(xmin, xmax, ymin, ymax, filelist)
    print "Starting Grid Making Process..."
    p = Pool(4)
    if (cython==True):
        for i in range(len(g.filelist)):
            p.apply_async(convolve_wrapper_cy, args=(g.filelist[i], g.biased, g.sigmaweight, g.tsysweight, g.RMAX, g.crval2, g.crval3, g.weights, g.naxes0, g.naxes1, g.naxes2, g.theta_n,), callback = callback_update)
    if (cython==False):
            for i in range(len(g.filelist)):
                p.apply_async(convolve_wrapper, args=(g.filelist[i], g.biased, g.sigmaweight, g.tsysweight, g.RMAX, g.crval2, g.crval3, g.weights, g.naxes0, g.naxes1, g.naxes2, g.theta_n,), callback = callback_update)
    p.close()
    p.join()
    g.normalize_grid()
    g.T = g.T.reshape((g.naxes2, g.naxes1, g.naxes0))
    print "Grid normalized and reshaped!"
    print "The median value of the grid is:"
    print g.T.var()**.5
    g.create_netcdf()
    
def initialize_regrid(xmin, xmax, ymin, ymax, filelist):
    global g
    g = LMTOTFRegrid_mp(xmin, xmax, ymin, ymax, filelist)

def callback_update(results):
    global g
    g.T += results[0]
    g.WT += results[1]
    g.TSYS += results[2]

def convolve_wrapper(filename, biased, sigmaweight,tsysweight,RMAX, crval2, crval3, weights, naxes0, naxes1, naxes2, theta_n):
    #global g
    t, wt, tsys = convolve_mp(filename, biased, sigmaweight,tsysweight,RMAX, crval2, crval3, weights, naxes0, naxes1, naxes2, theta_n)
    #g.T += t
    #g.WT += wt
    #g.TSYS += tsys
    return (t, wt, tsys)

def convolve_wrapper_cy(filename, biased, sigmaweight,tsysweight,RMAX, crval2, crval3, weights, naxes0, naxes1, naxes2, theta_n):
    #global g
    t, wt, tsys = convolve_cy(filename, biased, sigmaweight,tsysweight,RMAX, crval2, crval3, weights, naxes0, naxes1, naxes2, theta_n)
    #g.T += t
    #g.WT += wt
    #g.TSYS += tsys
    return (t, wt, tsys)




if __name__ == "__main__":
    files = ['35065305.nc','35065304.nc','35065306.nc','35065300.nc']#,'35065301.nc','35065302.nc','35065303.nc','35065307.nc','35065308.nc','35065309.nc']
    grid = gridmaker_1c(525.0,535.0,645.0,655.0, files)
    
