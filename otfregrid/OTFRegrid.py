import numpy
from scipy.special import j1
from file_compatibility.LMTOTFFile import LMTOTFNetCDFFile
from collections import OrderedDict

from netCDF4 import Dataset, Variable #, _private_atts
from file_compatibility import pynetcdf4



class LMTOTFRegrid(object):
    def __init__(self, xmin, xmax, ymin, ymax, filelist, case = 1, theta_n = 25.0, RMAX = 3.0, biased=True, tsysweight=False, sigmaweight=True):
        setattr(self, 'xmin', xmin)
        setattr(self, 'xmax', xmax)
        setattr(self, 'ymin', ymin)
        setattr(self, 'ymax', ymax)
        setattr(self, 'theta_n', theta_n)
        setattr(self, 'RMAX', RMAX)
        setattr(self, 'case', case)
        setattr(self, 'biased', biased)
        setattr(self, 'filelist', filelist)
        self.get_parameters()
        self.make_weights()
        #self.get_reduc_prefs()
        self.make_filelist()
        self.make_arrays()
        #self.make_grid()
        #self.create_netcdf()
    def make_filelist(self):
        # this finds the inputted filetype and then assigns some initial variables
        
        # this needs to be a literal list of string names of files for this
        # check to work. more durable solution coming soon.

        #setattr(self, 'filelist', filelist)
        file_0 = self.filelist[0]
        file_type = (file_0.split('.'))[1]
        if (file_type == 'fits'):
            self.filetype = 'fits'
        if (file_type == 'nc'):
            self.filetype = 'netCDF'
            print "File 0  is ", self.filelist[0]
            temp_otfscan = LMTOTFNetCDFFile(self.filelist[0])
            self.nchan = temp_otfscan.hdu.header.nchan)
    
    def get_parameters(self):
        
        i1 = abs(self.xmax) / (self.theta_n / 60.0)
        crval2 = 0.0
        if (self.xmax != 0.0):
            crval2 = i1 * (self.theta_n / 3600.0) * (abs(self.xmax) / self.xmax)
        i2 = abs(self.ymin) / (self.theta_n / 60.0)
        crval3 = 0.0
        if (self.ymin != 0.0):
            crval3 = i2 * (self.theta_n / 3600.0) * (abs(self.ymin) / self.ymin)
        setattr(self, 'crval2', crval2)
        setattr(self, 'crval3', crval3)
    #def get_reduc_prefs(self):
        # the class method this is used for
        # expects the following:
        #setattr(self, 'biased', biased) # true / false
        #setattr(self, 'tdep', tdep) # true / false
    
    #def get_convolve_weights(self, rmsweight, noiseweight):
    #    # both of these should be 1 or 0
    #    # I suspect they should be the same
    #    # for all data that goes into the final map
    #    setattr(self, 'rmsweight', rmsweight)
    #    setattr(self, 'noiseweight', rmsweight)
        
    def jinc(self, x):
        # a 'jinc' function
        if (len(x) == 1):
            if (x == 0.0):
                return 0.5
            else:
                return j1(x)/x
        else:
            y = numpy.zeros(len(x))
            index_0 = numpy.where(x == 0.0)
            index_not_0 = numpy.where(x != 0.0)
            y[index_0] = 0.5
            y[index_not_0] = j1(x[index_not_0])/x[index_not_0]
        return y

            
        
    def make_weights(self):
        #case is a stand in for the case/switch thing
        #from C - I am not sure how the different cases
        #are triggered 
        weights = numpy.zeros(256)
        dx = self.RMAX / 256.0
        i = numpy.arange(256)
        x = dx * i
        if (self.case == 1):
            # case 1
            weights[0] = numpy.exp(-1.0 * (2.0*x[0]/4.75)**2)
            weights[1:] = 2.0*self.jinc(6.28318531*(x[1:]) / 1.1) * numpy.exp(-1.0 * (2.0 * (x[1:]) / 4.75)**2) * 2.0*self.jinc(3.831706*(x[1:])/(self.RMAX))
        if (self.case == 2):
            # case 2
            wt1 = numpy.zeros(256)
            wt1[0] = 1.0
            wt1[1:] = numpy.sin(6.28318531*(x[1:])/1.1)/(6.28318531*(x[1:])/1.1)
            weights = wt1 * numpy.exp(-1.0 * (2.0*x/4.75)**2)
        if (self.case == 3):
            #case 3
            weights[0] = 1.0
            weights[1:] = numpy.sin(6.28318531*(x[1:])/1.1)/(6.28318531*(x[1:])/1.1)
        if (self.case == 4):
            #case 4
            weights = numpy.exp(-2.77258872*(x**2))
        setattr(self, 'weights', weights)

    def make_arrays(self):
        naxes0 = self.nchan
        naxes1 = (self.xmax - self.xmin) / (self.theta_n / 60.0) + 1
        naxes2 = (self.ymax - self.ymin) / (self.theta_n / 60.0) + 1
        setattr(self, 'naxes0', int(naxes0))
        setattr(self, 'naxes1', int(naxes1))
        setattr(self, 'naxes2', int(naxes2))
        nvox = (self.naxes0 * self.naxes1 * self.naxes2)
        nelements = self.naxes1 * self.naxes2
        setattr(self, 'T', numpy.zeros(nvox) )
        setattr(self, 'WT', numpy.zeros(nelements) )
        setattr(self, 'TSYS', numpy.zeros(nelements) )
        setattr(self, 'MAX_WT', numpy.ones(nelements) * (-1.0 * 10**30))
        setattr(self, 'INT_TIME', numpy.zeros(nelements))

    #def check_mask(self):
    #    pass
    #def get_dump(self):
    #    pass
    #def WeightREF(self):
    #    pass
    #def get_index(self, a, b, c):
    #    return (a + b*c)
    def convolve(self, filename):
        # this function as currently written (adapted from convolve.c)
        # is designed to work inside of two for loops
        # one to loop through each horn in the focal plane array (ih)
        # and one to loop through each data dump in each OTF run (idmp)
        if (self.tsysweight):
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
        dx = self.RMAX / 256.0
        
        # for all uses of XPOS and YPOS, I'm going to bump up the idmp index
        # by 1 - this is to compensate for them being also having the reference
        # positions as their 0th and last indices as well as the actual data in between
        
        XMAX = 60.0 * self.crval2 # in arcseconds? (crvals must be in arcminutes?)
        YMIN = 60.0 * self.crval3 # in arcseconds?
        lambda_D = 206264.81 * (clight / filename.hdu.header.fsky / ((1.0e9))/ D)
        step = self.theta_n / lambda_D

        for idmp in range(filename.hdu.header.nsample - 2):
            print idmp, "of", filename.hdu.header.nsample - 2
            for ih in range(filename.hdu.header.nhorns):
                if (self.tsysweight):
                    weight1 = wt1[ih]
                if (self.sigmaweight):
                    weight1 = 1.0 / (filename.hdu.data.sigma[idmp,ih] * filename.hdu.data.sigma[idmp,ih])
                else:
                    weight1 = wt1
                YBEG = 60.0 * (filename.hdu.data.YPOS[ih, idmp + 1] - YMIN) / lambda_D
                iybeg = int((YBEG - self.RMAX) / step  + 1)
                if (iybeg < 0):
                    iybeg = 0
                iyend = int((YBEG + self.RMAX) / step + 1)
                if (iyend > self.naxes2):
                    iyend = self.naxes2

                XBEG = 60.0 * (XMAX - filename.hdu.data.XPOS[ih, idmp + 1]) / lambda_D
                ixbeg = int((XBEG - self.RMAX) / step + 1)
                if (ixbeg < 0):
                    ixbeg = 0
                ixend = int((XBEG + self.RMAX) / step + 1)
                if (ixend > self.naxes1):
                    ixend = self.naxes1
                
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
                        if (rad2 < (self.RMAX * self.RMAX)):
                            IS = int(numpy.sqrt(rad2) / dx)
                            #print IS
                            weight0 = self.weights[IS]
                            wt = weight0 * weight1
                            #wt = 1.0 # temporary
                            ii = int(ix + iy * self.naxes1)
                            #print ii
                            self.WT[ii] += wt
                            self.TSYS[ii] += filename.hdu.header.tsys[ih] * wt
                            self.INT_TIME[ii] += wt * filename.hdu.header.tdump
                            if (weight0 > self.MAX_WT[ii]):
                                self.MAX_WT[ii] = weight0
                            jj = int(ii * self.naxes0)
                            #print jj
                            #print (ix, iy)
                            for k in range(self.naxes0):
                                if (filename.hdu.data.reduced[idmp, ih, k] < (10.0)**30):
                                    self.T[jj + k] += (wt*filename.hdu.data.reduced[idmp, ih, k])
                                    #if (wt*filename.hdu.data.reduced[idmp, ih, k] == 0.0):
                                        #print "index", k, "is 0.0 ."
                        #else:
                            #print "rad2 > rmax^2"
                
                        #iy = iy + 1
                    #ix = ix + 1

        #setattr(self, 'T', T)
    def normalize_grid(self):
        for iy in range(int(self.naxes2)):
            for ix in range(int(self.naxes1)):
                ii = ix + iy * self.naxes1
                # MAXWEIGHT is defined as .25, so I replaced it here
                if (self.WT[ii] != 0.0) and (self.MAX_WT[ii] > .25):
                    self.TSYS[ii] = self.TSYS[ii] / self.WT[ii]
                    jj = ii * self.naxes0
                    for k in range(int(self.naxes0)):
                        if (self.T[jj + k] < 10**30):
                            self.T[jj + k] = self.T[jj + k] / self.WT[ii]
                        else:
                            self.T[jj + k] = 0.0
                            print k
    def make_grid(self):
        for filename in range(len(self.filelist)):
            otffile = LMTOTFNetCDFFile(self.filelist[filename])
            otffile._reduce_data(self.biased)
            if (self.sigmaweight):
                otffile.baseline()

            #print Z
            print self.filelist[filename]
            #for idmp in range(otffile.hdu.header.nsample - 2):
                #print idmp, "of", otffile.hdu.header.nsample - 2
                # the C file checks if XOFF[idmp] < 21,600 here
                # skipping that for now
                #for ih in range(otffile.hdu.header.nhorns):
                    #for k in range(otffile.hdu.header.nchan):
                        # in the C code, here it runs 'checkMASK'
                        # skipping that too, because we don't have
                        # the masks (they're in a .x file?) [spooky]
                    
                        # then it runs 'WeightREF' - this just
                        # makes the weights for reducing
                        # the data dump
                        # however: this is now a method of the
                        # LMTOTFNetCDFFile class
                        # (self._reduce_data())
            self.convolve(otffile)
        self.normalize_grid()
        self.T = self.T.reshape((self.naxes2, self.naxes1, self.naxes0))
        self.create_netcdf()

    def create_netcdf(self):
        #I'm going to name the output grid's
        #netcdf file after the first file in the grid
        #and add the total filelist as an attribute
        
        newfilename = (self.filelist[0]).strip("'")
        newfilename = newfilename.strip(".loa")
        newfilename += ("_grid.nc")
        if (len(self.filelist) > 8):
            newfilename = "large_grid.nc"
        
        newfile = Dataset(newfilename, mode='w', clobber=True)
        newfile.createDimension('naxes0', self.naxes0)
        newfile.createDimension('naxes1', self.naxes1)
        newfile.createDimension('naxes2', self.naxes2)
        
        var = newfile.createVariable('otfmap', numpy.dtype(numpy.float32), (('naxes2', 'naxes1', 'naxes0')))
        
        var[:] = self.T
        
        var.__setattr__('filenames', self.filelist)
        var.__setattr__('xmax', self.xmax)
        var.__setattr__('ymax', self.ymax)
        var.__setattr__('xmin', self.xmin)
        var.__setattr__('ymin', self.ymin)
        
        newfile.close()




if __name__ == "__main__":
    files =  ['35065305.nc','35065304.nc','35065306.nc','35065300.nc','35065301.nc','35065302.nc','35065303.nc','35065307.nc','35065308.nc','35065309.nc']
    otf = LMTOTFRegrid(525.0, 535.0, 645.0, 655.0, files)
