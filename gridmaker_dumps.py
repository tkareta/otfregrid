import numpy
from multiprocessing import Pool, cpu_count

from otfregrid import convolve_mp
from otfregrid import LMTOTFRegrid_mp
from file_compatibility.LMTOTFFile import LMTOTFNetCDFFile
print "If importing 'convolve_dump' fails, compile the Cython module on your system as the wiki describes"
from convolve_dump import convolve_dump

# this version of gridmaker is designed to handle a large-ish (1-30?)
# filelist and make a grid using multiprocessing on one computer.
# As it stands currently, it can process approximately 4 files per 12
# minute increment.

def gridmaker_dumps(xmin, xmax, ymin, ymax, filelist, dataloc="null", writeloc="null",cython=True, normalize=True):
    ### going to move the 'make_grid' function outside of the class
    ### for multiprocessing purposes? we'll see if it works
    #g = LMTOTFRegrid_mp(xmin, xmax, ymin, ymax, filelist)
    backupfilelist = filelist
    if(writeloc=="null"):
        print "Assuming that all needed files are in local directory!"
    else:
        for i in range(len(filelist)):
            filelist[i] = dataloc+filelist[i] #this pre-appends the location of the data
    g = initialize_regrid(xmin, xmax, ymin, ymax, filelist)
    print "Starting Grid Making Process..."
    
    p = Pool(cpu_count())
    if (cython==True):
        print "number of files = ", len(g.filelist)
        for i in range(len(g.filelist)):
            filename = LMTOTFNetCDFFile(g.filelist[i])
            filename._reduce_data(g.biased)
            if (g.sigmaweight):
                filename.baseline()
            if (g.tsysweight):
                wt1 = numpy.zeros(g.nhorns)
                filename.hdu.data.sigma = 0.0
                for ih in range(g.nhorns):
                    if (filename/hdu.header.tsys[ih] != 0.0):
                        wt1[ih] = 1.0 / (filename.hdu.header.tsys[ih] * filename.hdu.header.tsys[ih])
                    else:
                        wt1[ih] = 1.0
            else:
                wt1 = 1.0
            
            print "number of dumps: ", filename.hdu.header.nsample - 2, " in file ", g.filelist[i]
            for idmp in range(filename.hdu.header.nsample - 2):
                #print idmp
                p.apply_async(convolve_wrapper_dump, args=(g.biased, g.sigmaweight, g.tsysweight, g.RMAX, g.crval2, g.crval3, g.weights, g.naxes0, g.naxes1, g.naxes2, g.theta_n, filename.hdu.header.nhorns, g.nchan, filename.hdu.header.fsky, idmp, filename.hdu.data.XPOS, filename.hdu.data.YPOS, filename.hdu.data.reduced, filename.hdu.data.sigma, wt1, filename.hdu.header.tsys), callback = callback_update)
                #print "process added to Pool"
    if (cython==False):
        for i in range(len(g.filelist)):
            p.apply_async(convolve_wrapper, args=(g.filelist[i], g.biased, g.sigmaweight, g.tsysweight, g.RMAX, g.crval2, g.crval3, g.weights, g.naxes0, g.naxes1, g.naxes2, g.theta_n,), callback = callback_update)
    p.close()
    p.join()
    if (normalize==True):
        g.normalize_grid()
        g.T = g.T.reshape((g.naxes2, g.naxes1, g.naxes0))
        print "Grid normalized and reshaped!"
        print "The stdev value of the grid is:"
        print g.T.var()**.5
        g.create_netcdf()
    if (normalize==False):
        print "Assuming the filelist was delegated to multiple computers:"
        if ((writeloc == "null")or(dataloc == "null")):
            print "Need to know where the file is to be written!"
        else:
            file0 = filelist[0]
            print file0
            file0 = file0[len(writeloc):len(file0)]
            filenew = file0.strip(".nc")
            print "The T, WT, TSYS, MAX_WT and naxes arrays are stored in that order in the file:"
            print writeloc+filenew
            naxes = numpy.array([g.naxes2, g.naxes1, g.naxes0])
            numpy.savez(writeloc+filenew, g.T, g.WT, g.TSYS, g.MAX_WT, naxes)
        
    
####
def initialize_regrid(xmin, xmax, ymin, ymax, filelist):
    global g
    g = LMTOTFRegrid_mp(xmin, xmax, ymin, ymax, filelist)
    return g

def callback_update(results):
    global g
    g.T += results[0]
    #print results[0]
    g.WT += results[1]
    g.TSYS += results[2]
    g.MAX_WT = numpy.maximum(g.MAX_WT, results[3])
    #print results[0]

def convolve_wrapper(filename, biased, sigmaweight,tsysweight,RMAX, crval2, crval3, weights, naxes0, naxes1, naxes2, theta_n):
    #global g
    t, wt, tsys = convolve_mp(filename, biased, sigmaweight,tsysweight,RMAX, crval2, crval3, weights, naxes0, naxes1, naxes2, theta_n)
    #g.T += t
    #g.WT += wt
    #g.TSYS += tsys
    return (t, wt, tsys)

def convolve_wrapper_dump(biased, sigmaweight,tsysweight,RMAX, crval2, crval3, weights,  naxes0,  naxes1, naxes2, theta_n, nhorns, nchan, fsky, idmp, XPOS, YPOS, reduced, sigma, wt1, tsys):
    #global g
    #print "convolve wrapper dump"
    t, wt, tsys, max_wt = convolve_dump(biased, sigmaweight,tsysweight,RMAX, crval2, crval3, weights,  naxes0,  naxes1, naxes2, theta_n, nhorns, nchan, fsky, idmp, XPOS, YPOS, reduced, sigma, wt1, tsys)
    #g.T += t
    #g.WT += wt
    #g.TSYS += tsys
    return (t, wt, tsys, max_wt)




if __name__ == "__main__":
    files = ['35065304.nc','35065306.nc','35065300.nc','35065301.nc','35065302.nc','35065303.nc','35065307.nc','35065308.nc','35065305.nc', '35065309.nc']
    #import sys
    import os
    cwd = os.getcwd() #gets the current working directory
    #sys.path.append(cwd + "/test_data")
    os.chdir(cwd + "/test_data")
    grid = gridmaker_dumps(525.0,535.0,645.0,655.0, files)
    os.chdir(cwd)
    
