cimport numpy
import numpy
#from file_compatibility.LMTOTFFile import LMTOTFNetCDFFile
cimport cython
import time

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef convolve_dump(biased, sigmaweight,tsysweight,RMAX, crval2, crval3, weights,  naxes0,  naxes1, naxes2, theta_n, nhorns, nchan, fsky, idmp, XPOS, YPOS, reduced, sigma, wt1, tsys):
    print "starting a loop through a data dump"
    #filename = LMTOTFNetCDFFile(filename)

    cpdef numpy.ndarray[double, ndim=1] T
    cpdef numpy.ndarray[double, ndim=1] WT
    cpdef numpy.ndarray[double, ndim=1] TSYS

    # these arrays need to be shared in memory between
    # all of the processes writing to the output grid
    T = numpy.zeros(naxes0*naxes1*naxes2)    
    WT = numpy.zeros(naxes1*naxes2)
    TSYS = numpy.zeros(naxes1*naxes2)
    #global T
    #global WT
    #global TSYS
    
    cpdef numpy.ndarray[double, ndim=1] MAX_WT
    cpdef numpy.ndarray[double, ndim=1] INT_TIME
    #cpdef int nchan
    MAX_WT = numpy.ones(naxes1*naxes2) * (-1.0 * 10**30)
    INT_TIME = numpy.zeros(naxes1*naxes2)
    #nchan = filename.hdu.header.nchan

    #filename._reduce_data(biased)
    
    #if (sigmaweight):
        #filename.baseline()
    #if (tsysweight):
        #wt1 = numpy.zeros(filename.hdu.header.nhorns)
        #for ih in range(filename.hdu.header.nhorns):
            #if (filename.hdu.header.tsys[ih] != 0.0):
                #wt1[ih] = 1.0 / (filename.hdu.header.tsys[ih] * filename.hdu.header.tsys[ih])
            #else:
                #wt1[ih] = 1.0
    #else:
	
        #wt1 = 1.0


    cpdef double clight
    clight = 2.99792458e10 #cm/s

    #in the future, we can just check what telescope it was from
    #and then assign a diameter value based on that
    cpdef double D
    D = 1370.0 #FCRAO diameter in cm

    cpdef int ic
    ic = int(nchan / 2.0)
    if (nchan > 128):
        ic = 505
	
    cpdef double dx
    dx = RMAX / 256.0

    # for all uses of XPOS and YPOS, I'm going to bump up the idmp index
    # by 1 - this is to compensate for them being also having the reference
    # positions as their 0th and last indices as well as the actual data in between
    cpdef double XMAX
    cpdef double YMIN
    XMAX = 60.0 * crval2 # in arcseconds? (crvals must be in arcminutes?)
    YMIN = 60.0 * crval3 # in arcseconds?

    cpdef double lambda_D
    lambda_D = 206264.81 * (clight / fsky / ((1.0e9))/ D)
    step = theta_n / lambda_D

    #cpdef double weight1, weight0, wt
    cpdef double YBEG
    cpdef int iybeg

    cpdef double deltax, dx2
    cpdef double deltay, rad2
    cpdef int IS, ii
    print ""
    time1 = time.time()
    for fake in range(1):
        #print idmp, "of", filename.hdu.header.nsample - 2
        for ih in range(nhorns):
            if (tsysweight):
                weight1 = wt1[ih]
            if (sigmaweight):
                weight1 = 1.0 / (sigma[idmp,ih] * sigma[idmp,ih])
            else:
                weight1 = wt1
            YBEG = 60.0 * (YPOS[ih, idmp + 1] - YMIN) / lambda_D
            iybeg = int((YBEG - RMAX) / step  + 1)
            if (iybeg < 0):
                iybeg = 0
            iyend = int((YBEG + RMAX) / step + 1)
            if (iyend > naxes2):
                iyend = naxes2

            XBEG = 60.0 * (XMAX - XPOS[ih, idmp + 1]) / lambda_D
            ixbeg = int((XBEG - RMAX) / step + 1)
            if (ixbeg < 0):
                ixbeg = 0
            ixend = int((XBEG + RMAX) / step + 1)
            if (ixend > naxes1):
                ixend = naxes1

            #iy = iybeg
            #ix = ixbeg
            #T = self.T
            for ix_l in range(ixend - ixbeg):
                ix = ix_l + ixbeg
                deltax = XBEG - (ix * step)
                dx2 = deltax * deltax
                #print deltay
                #dely2 = deltay * deltay
                for iy_l in range(iyend - iybeg):
                    iy = iy_l + iybeg
                    deltay = (iy * step) - YBEG
                    #print deltax
                    rad2 = dx2 + deltay * deltay
                    #print rad2
                    if (rad2 < (RMAX * RMAX)):
                        IS = (numpy.sqrt(rad2) / dx)
                        #print IS
                        weight0 = weights[IS]
                        wt = weight0 * weight1
                        #wt = 1.0 # temporary
                        ii = (ix + iy * naxes1)
                        #print ii
                        WT[ii] += wt
                        TSYS[ii] += tsys[ih] * wt
                        #INT_TIME[ii] += wt * filename.hdu.header.tdump
                        if (weight0 > MAX_WT[ii]):
                            MAX_WT[ii] = weight0
                        jj = int(ii * naxes0)
                        #print jj
                        #print (ix, iy)
                        for k in range(naxes0):
                            if (reduced[idmp, ih, k] < (10.0)**30):
                                T[jj + k] += (wt*reduced[idmp, ih, k])
    time2 = time.time()
    print "took", (time2-time1), "seconds for this data dump"
    return T, WT, TSYS

    
