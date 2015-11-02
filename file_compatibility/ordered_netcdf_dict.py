from ordered_dict import OrderedDict
import numpy
try:
    from netCDF4 import Variable as NetCDFVariable
except ImportError:
    from pynetcdf import NetCDFVariable
from otfregrid_error import OTFRegridGeneralError, OTFRegridArgumentError #replaced the dreampy errors with these
import types

def netcdf_repr(key, val):
    if isinstance(val, NetCDFVariable):
        if val.dtype == numpy.dtype('c'):
            variab = val[:].tostring().strip().strip('\x00')
        else:
            variab = val[:]    
    else:
        variab = val
    return variab

def netcdf_getValue(val):
    """Given a NetCDFVariable type of input,
    this function will return either the array type if the size
    of the variable > 1 or if size == 1, the scalar value
    typecast to the correct numeric type
    """
    if isinstance(val, NetCDFVariable):
        if val.dtype == numpy.dtype('c'):
            variab = val[:].tostring().strip().strip('\x00')
            return variab
        else:
            if val.size > 1:
                return val[:]
            elif val.size == 1:
                dtype = val.dtype
                return dtype.type(val[:])
            else:
                return val
    else:
        return val

class OrderedNetCDFDict(OrderedDict):
    """
    Inherits from Ordered Dictionary object type.
    The netCDF variables are stored as OrederedNetCDFDict instants.
    This class allows the setting and getting of values
    """
    def _reprlist(self):
        mystr_list = ['{ ']
        for i, key in enumerate(self._sequence):
            mystr_list.append("  '%s': %s" % (key, netcdf_repr(key, self[key])))
        mystr_list.append('}')
        return mystr_list
    
    def __repr__(self):
        """
        Used for __repr__ and __str__
        """
        #return '{%s}' % ('\n'. join(['%r: %r' % (key, netcdf_repr(self[key])) for key in self._sequence]))
        mystr_list = ['{ ']
        for i, key in enumerate(self._sequence):
            if isinstance(self[key], OrderedNetCDFDict):
                mystr_list.append("  '%s':" % key)
                mystr_list.extend([' '*(len(key)+4)+s for s in self[key]._reprlist()])
            else:
                mystr_list.append("  '%s': %s" % (key, netcdf_repr(key, self[key])))
        mystr_list.append('}')
        #if isinstance(
        return '\n'.join(mystr_list)

    def get(self, item, default=None, debug=False):
        """
        Given an item such as 'Dcs.Receiver', returns the value
        of the NetCDF variable corresponding to self['Dcs']['Receiver']
        If you just specify, 'Dcs', it will look for a NetCDF variable
        with that name.
        """
        if item.find('.'):
            header_tree = item.split('.')
            val = self.copy()
            for h in header_tree:
                if val.has_key(h):
                    val = val.__getitem__(h)
                else:
                    if debug:
                        print "Key %s not found in header" % h
                    return default
            return netcdf_getValue(val)
        else:
            #no dots return the header item directly
            if val.has_key(h):
                val = self.__getitem__(item)
            else:
                if debug:
                    print "Key %s not found in header" % item
                return default
            return netcdf_getValue(val)

    def get_scalar(self, item, default=None):
        """
        Given a scalar item such as 'Dcs.Receiver', returns the value
        of the NetCDF variable corresponding to self['Dcs']['Receiver']
        If you just specify, 'Dcs', it will look for a NetCDF variable
        with that name.
        """
        var = self.get(item, default=default)
        if var is None:
            return default
        else:
            if type(var) == types.StringType:
                return var
            elif isinstance(var, numpy.ndarray):
                if var.size == 1:
                    return var.dtype.type(var[0])
                else:
                    print "Not scalar"
                    return var
            else:
                return var
        
    def set(self, item, value, nc=None,
            dtype=None, dimensions=()):
        """
        Given an item such as Dcs.Receiver and an appropriate
        value, an existing NetCDF variable is modified by the
        new value, or else a new NetCDF variable is created.
        """
        #first find out if item exists
        create = False
        if item.find('.'):
            header_tree = item.split('.')
            var = self  #.copy()
            for i, h in enumerate(header_tree):
                if not var.has_key(h):
                    #raise DreampyGeneralError("No Such item", "Item %s is not found in header" % h)
                    print "Creating item %s in header" % h
                    if i < len(header_tree)-1:
                        var[h] = OrderedNetCDFDict()
                        var = var.__getitem__(h)
                    else:
                        var[h] = None
                    create = True
                else:
                    var = var.__getitem__(h)
        else:
            #no dots
            var = self  #.copy()
            if not var.has_key(item):
                print "Creating item %s in header" % item
                create = True
                h = item
                var[h] = None
                #var[item] = OrderedNetCDFDict()
            else:
                var = self.__getitem__(item)
        if not create:
            #existing variable is to be modified
            #shape = var.shape
            dimensions = var.dimensions
            dtype = var.dtype
            if isinstance(value,(type(None),str,int,float,bool)):
                #this is a scalar value
                if var.size == 1:
                    #var is also a scalar
                    var[:] = value
                else:
                    raise OTFRegridArgumentError("Wrong Size", "Variable %s has size %d, not the same as value given %s" % (item, len(var), value))
            else:
                if len(var) != len(value):
                    raise OTFRegridArgumentError("Wrong Size", "Variable %s has size %d, not the same as size of value given %s" % (item, len(var), value))
                if var.shape != value.shape:
                    raise OTFRegridArgumentError("Wrong Shape", "Variable %s has shape %s, not the same as shape of value given by %s" % (item, var.shape, value.shape))
                var[:] = value
            print "Item %s updated to %s" % (item, self.get(item))
        else:
            #create a brand new variable
            if nc is None:
                raise OTFRegridGeneralError("Create Error", "New NetCDF variable is to be created. Please make sure to pass the netCDF dataset instance in nc")
            new_var = nc.createVariable("Header.%s" % item, dtype,
                                        dimensions=dimensions)
            new_var[:] = value
            var[h] = new_var
            

                                        
