"""
LMT Specific NetCDF File Reader/Writer Module.

Gopal Narayanan <gopal@astro.umass.edu>

This is a thin wrapper around the pynetcdf4 library
or the pupynere.py library to interface to NetCDF v4/v3
files.
"""

try:
    #Try importing the NetCDF4 version first
    from pynetcdf4 import NetCDFFile
except ImportError:
    #from pynetcdf import NetCDFFile, NetCDFVariable
    from pynetcdf import NetCDFFile
    
#
from lmtheader import LMTHeader

from lmtholodata import LMTHoloData
from ordered_dict import OrderedDict

class LMThdu(object):
    """
    The LMThdu is a class that contains the basic
    Header-Data Unit. This is meant to be similar to pyfits
    HDU. The two main attributes are the header and data.
    In addition, there can be a filename attribute that contains
    the name of the file from which the data was read
    """
    def __init__(self, data=None, header=None,
                 filename=None):
        """
        The initialization is usually performed by passing in
        three parameters.

        @param data: A data instance derived from L{LMTData} class. 
        @type data: instance of L{LMTData} class
        @param header: A header instance derived from L{LMTHeader} class.
        @type header: instance of L{LMTHeader} class
        @param filename: name of the input netCDF file
        @type filename: string
        @return: an instance of the LMThdu class
        @rtype: L{LMThdu} instance
        """
        self.data = data
        self.header = header
        self.filename = filename

    def __repr__(self):
        if self.filename:
            return "HDU for file: %s" % self.filename
        else:
            return "HDU"
            
class LMTNetCDFFile(NetCDFFile):
    """
    A LMTNetCDFFile is an object that opens up a LMT netcdf v3 file
    and reads the header variables and data variables.

    It can also be used to create LMT netcdf v3 files given header and data
    objects.
    """
    def __init__(self, filename, mode='r'):
        NetCDFFile.__init__(self, filename, mode)
        self.filename = filename
        #if self.file_format=='NETCDF3_CLASSIC' and mode in ('r', 'a'):
        if len(self.groups) == 0 and mode in ('r', 'a'):
            self._make_hdu(self.variables, self.dimensions)
        #self.populate_data()
        self.hdus = OrderedDict()
        
    def _make_hdu(self, variables, dimensions):
        """
        A HDU is a Header-Data Unit. This is a way to keep the
        header (meta-data) and data attributes of a dataset together
        First get the data and then the headers. The order is important
        for some type of files like holography files.

        @param variables: A netCDF variables object which is an OrderedDict
        of netCDF variables
        @type variables: An OrderedDict of netCDF variables
        @param dimensions: A netCDF dimensionss object which is an OrderedDict
        of netCDF dimensions
        @type dimensions: An OrderedDict of netCDF dimensions
        @return: A LMThdu instance
        """
        data = self._populate_data(variables)
        #print "Make headers"
        header = self._populate_headers(variables, dimensions)
        self.hdu = LMThdu(data=data, header=header,
                          filename=self.filename)
        
    def _populate_headers(self, variables, dimensions):
        return LMTHeader(ncvariables=variables,
                         dimensions=dimensions)

    def _populate_data(self, variables):
        """
        here we will actually distinguish what type of data
        type we are dealing with - viz. Holography, AzTEC, redshift, etc.
        and then use the appropriate class to instantiate"""
        #for now use the LMTHoloData class by default
        return LMTHoloData(variables)
        

    
