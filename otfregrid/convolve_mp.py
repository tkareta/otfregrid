import numpy
from file_compatibility.LMTOTFFile import LMTOTFNetCDFFile

def convolve_mp(filename, biased, sigmaweight,tsysweight,RMAX, crval2, crval3, weights, naxes0, naxes1, naxes2, theta_n):
    filename = LMTOTFNetCDFFile(filename)
    T = numpy.zeros(naxes0*naxes1*naxes2)
    
    WT = numpy.zeros(naxes1*naxes2)
    
    TSYS = numpy.zeros(naxes1*naxes2)
    
    MAX_WT = numpy.ones(naxes1*naxes2) * (-1.0 * 10**30)
    INT_TIME = numpy.zeros(naxes1*naxes2)
    nchan = filename.hdu.header.nchan
    filename._reduce_data(biased)
    if (sigmaweight):
        filename.baseline()

    if (tsysweight):
        wt1 = numpy.zeros(filename.hdu.header.nhorns)
        for ih in range(filename.hdu.header.nhorns):
            if (filename.hdu.header.tsys[ih] != 0.0):
                wt1[ih] = 1.0 / (filename.hdu.header.tsys[ih] * filename.hdu.header.tsys[ih])
            else:
                wt1[ih] = 1.0
    else:
        wt1 = 1.0
    clight = 2.99792458e10 #cm/s

    #in the future, we can just check what telescope it was from
    #and then assign a diameter value based on that
    D = 1370.0 #FCRAO diameter in cm
    ic = int(filename.hdu.header.nchan / 2.0)
    if (filename.hdu.header.nchan > 128):
        ic = 505
    dx = RMAX / 256.0

    # for all uses of XPOS and YPOS, I'm going to bump up the idmp index
    # by 1 - this is to compensate for them being also having the reference
    # positions as their 0th and last indices as well as the actual data in between

    XMAX = 60.0 * crval2 # in arcseconds? (crvals must be in arcminutes?)
    YMIN = 60.0 * crval3 # in arcseconds?
    lambda_D = 206264.81 * (clight / filename.hdu.header.fsky / ((1.0e9))/ D)
    step = theta_n / lambda_D

    for idmp in range(filename.hdu.header.nsample - 2):
        print idmp, "of", filename.hdu.header.nsample - 2
        for ih in range(filename.hdu.header.nhorns):
            if (tsysweight):
                weight1 = wt1[ih]
            if (sigmaweight):
                weight1 = 1.0 / (filename.hdu.data.sigma[idmp,ih] * filename.hdu.data.sigma[idmp,ih])
            else:
                weight1 = wt1
            YBEG = 60.0 * (filename.hdu.data.YPOS[ih, idmp + 1] - YMIN) / lambda_D
            iybeg = int((YBEG - RMAX) / step  + 1)
            if (iybeg < 0):
                iybeg = 0
            iyend = int((YBEG + RMAX) / step + 1)
            if (iyend > naxes2):
                iyend = naxes2

            XBEG = 60.0 * (XMAX - filename.hdu.data.XPOS[ih, idmp + 1]) / lambda_D
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
                        IS = int(numpy.sqrt(rad2) / dx)
                        #print IS
                        weight0 = weights[IS]
                        wt = weight0 * weight1
                        #wt = 1.0 # temporary
                        ii = int(ix + iy * naxes1)
                        #print ii
                        WT[ii] += wt
                        TSYS[ii] += filename.hdu.header.tsys[ih] * wt
                        INT_TIME[ii] += wt * filename.hdu.header.tdump
                        if (weight0 > MAX_WT[ii]):
                            MAX_WT[ii] = weight0
                        jj = int(ii * naxes0)
                        #print jj
                        #print (ix, iy)
                        for k in range(naxes0):
                            if (filename.hdu.data.reduced[idmp, ih, k] < (10.0)**30):
                                T[jj + k] += (wt*filename.hdu.data.reduced[idmp, ih, k])
    return T, WT, TSYS

    
