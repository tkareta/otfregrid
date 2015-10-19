"""This module implements the NetCDF-4 read/write capabilities.
It is a wrapper around python-netCDF4 which in turn depends on
libhdf5 and libnetcdf4
Gopal Narayanan <gopal@astro.umass.edu>
"""

from netCDF4 import Dataset, Variable, _private_atts
import numpy


class NetCDFFile(Dataset, object):
    _private_atts.extend(['filename', 'hdu', 'hdus'])
    def __init__(self, filename, mode='r', 
                 clobber=True, format='NETCDF4'):
        super(NetCDFFile, self).__init__(filename, mode=mode,
                                         clobber=clobber, format=format)
        #self.__dict__['filename'] = filename
        #self.filename = filename
        

def Variable_repr(self):
    if self.dtype == numpy.dtype('c'):
        return "%s" % self.tostring().strip()
    else:
        return "%s" % self[:]

NetCDFVariable = Variable

#Variable.__str__ = Variable_repr

