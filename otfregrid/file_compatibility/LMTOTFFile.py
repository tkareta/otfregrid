### writen by Theodore Kareta, with major help from Dr. Gopal Narayanan
import lmtnetcdffile
from dreampy.lmtheader import LMTHeader
import os
from netCDF4 import _private_atts
from dreampy.utils import DreampyGeneralError
from LMTOTFHeader import LMTOTFHeader_old
from LMTOTFData import LMTOTFData_old
import numpy

class LMTOTFNetCDFFile(LMTNetCDFFile):
    _private_atts.extend(['old'])
    def __init__(self, filename, mode='r',
                 old=True):
        self.old = old
        if mode == 'w':
            super(LMTOTFNetCDFFile, self).__init__(filename)
        else:
            if os.path.exists(filename):
                super(LMTOTFNetCDFFile, self).__init__(filename)
            else:
                raise DreampyGeneralError("Does Not Exist", "File %s does not exist" % filename)
        self._populate_headers(self.variables, self.dimensions)
        #self._reduce_data(weightmode=equal)

    def _populate_headers(self, variables, dimensions):
        """
        redefine the conventional _populate_headers method
        for a holography specific header
        """
        if self.old:
            return LMTOTFHeader_old(ncvariables=variables,
                                    dimensions=dimensions
                                    )
        else:
            return LMTHeader(ncvariables=variables,
                             dimensions=dimensions
                             )

    def _reduce_data(self, biased):
        #this subtracts the bandpass characteristic of the receiever
        #by removing the references position observations
        #currently set to only work with old data because it can
        #only handle 2 references driven by a ref array

        # biased or equal

        # first, we make a blank array the size of the data dump array
        # minus the references
        newshape = ((self.hdu.header.nsample - 2,self.hdu.header.nhorns, self.hdu.header.nchan))
        reduced = numpy.zeros(newshape)
        ref1 = self.hdu.data.refs[0,:,:]
        ref1.reshape((1, self.hdu.header.nhorns, self.hdu.header.nchan))
        ref2 = self.hdu.data.refs[1,:,:]
        ref2.reshape((1, self.hdu.header.nhorns, self.hdu.header.nchan))
        if (biased):
            bias = numpy.arange(self.hdu.header.nsample - 2) / (float(self.hdu.header.nsample - 2))
            wt1 = 1.0 - bias
            wt2 = bias
            wt1 = wt1.reshape((self.hdu.header.nsample - 2, 1, 1))
            wt2 = wt2.reshape((self.hdu.header.nsample - 2, 1, 1))
        else:
            wt1 = .5
            wt2 = .5
        #print wt1 + wt2
        for i in range(self.hdu.header.nsample - 2):
            reduced[i,:,:] = self.hdu.data.dump[i+1,:,:]
        #print self.hdu.data.dump.shape
        #print reduced.shape
        #print reduced.mean()
        #print reduced.max()
        TCAL = self.hdu.header.tsys.reshape(1,32,1)
        reduced = (reduced - wt1 * ref1 - wt2 * ref2) / (wt1 * ref1 + wt2 * ref2) * TCAL
        #note: reduced now has units of Kelvin

        #print reduced.mean()
        #print reduced.max()
         
        setattr(self.hdu.data, 'reduced', reduced)
        
    def baseline(self, order = 0, subtract = False):
        #this will only work if the data has already been
        #reduced using the _reduce_data method
        if not hasattr(self.hdu.data, 'reduced'):
            raise ValueError("You have no 'reduced' array!")


        defaultwindows = [(0,450),(550,1023)]
        setattr(self.hdu.data, 'order', order)
        #windows over which to do the fit
        #ideally excluding any large lines
        #for now, just using channel numbers
        
 
        windows = []
        sigma = numpy.zeros((self.hdu.header.nsample - 2, self.hdu.header.nhorns))
        
        if hasattr(self, 'compwindows'):
            pass
        else:
            for twople in range(len(defaultwindows)):
                c1, c2 = defaultwindows[twople]
                #c1, c2 = sorted((c1, c2))
                windows.append((c1,c2))

        x = numpy.arange(self.hdu.header.nchan)
        c_loop = 0
        for win in windows:  
            if (len(win) != 2):
                print "Each element in the window list must be a 2-tuple"
            c1, c2 = win
            c1, c2 = sorted((c1, c2))
            ind = numpy.logical_and(x>=c1, x<=c2)
            if (c_loop == 0):
                final_ind = numpy.logical_or(ind,ind)
            else:
                final_ind = numpy.logical_or(final_ind,ind)
                
        x_windows = x[final_ind]
        for idmp in range(self.hdu.header.nsample - 2):
            for ih in range(self.hdu.header.nhorns):
                spectra = self.hdu.data.reduced[idmp,ih,:]
                spec_windows = spectra[final_ind]
                    
                    
                p = numpy.polyfit(x_windows,spec_windows,self.hdu.data.order)
                spec_windows -= numpy.polyval(p,len(x_windows))
                
                sigma[idmp, ih] = spec_windows.std()
                if (subtract):
                    self.hdu.data.reduced[idmp,ih,:] -= numpy.polyval(p,len(x))
                    
        setattr(self.hdu.data, "sigma", sigma)
