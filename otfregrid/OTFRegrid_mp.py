import numpy
from scipy.special import j1
from file_compatibility.LMTOTFFile import LMTOTFNetCDFFile
from collections import OrderedDict
from convolve_mp import convolve_mp
from netCDF4 import Dataset, Variable #, _private_atts
from file_compatibility import pynetcdf4

#new!
from multiprocessing import Pool

class LMTOTFRegrid_mp(object):
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
        setattr(self, 'tsysweight', tsysweight)
        setattr(self, 'sigmaweight', sigmaweight)
        self.get_parameters()
        self.make_weights()
        #self.get_reduc_prefs()
        self.make_filelist()
        self.make_arrays()
        self.make_grid()
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
            self.nchan = temp_otfscan.hdu.header.nchan
    
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
        for i in range(len(self.filelist)):
            convolve_mp(self.filelist[i], self.biased, self.sigmaweight, self.tsysweight, self.RMAX, self.crval2, self.crval3, self.weights, self.naxes0, self.naxes1, self.naxes2, self.theta_n)
            
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
