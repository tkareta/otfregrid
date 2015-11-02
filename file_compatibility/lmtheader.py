import numpy
try:
    from netCDF4 import Variable as NetCDFVariable
except ImportError:
    from pynetcdf4 import NetCDFVariable
from ordered_netcdf_dict import OrderedNetCDFDict

class LMTHeader(OrderedNetCDFDict):
    """The most generic LMT Header class. This class is not usually
    instantiated by the user. It is part of the initialization of
    opening an LMT NetCDF file."""
    def __init__(self, ncvariables=None, dimensions=None,
                 fromold=False):
        """
        if fromold is True, this is a header initialized
        from old data formats"""
        OrderedNetCDFDict.__init__(self, init_val=(), strict=False)
        #self.ncvariables = ncvariables
        if not fromold:
            self.make_header_keys(ncvariables)
        self.dimensions = dimensions

    def make_header_keys(self, ncvariables):
        heads = [name for name in ncvariables.keys() if name.find('Header') != -1]
        for head in heads:
            htype = head.split('.')[1]
            self[htype] = OrderedNetCDFDict()
        for key in self.keys():
            for subhead in [n for n in ncvariables.keys() if n.startswith('Header.%s.' % key)]:
                shead = subhead.split('.')[-1]
                self[key][shead] = ncvariables[subhead]

    def _pprint_subhead(self, k, v):
        hstr = '{'
        for i, (key, val) in enumerate(v.items()):
            if isinstance(val, NetCDFVariable):
                if val.dtype == numpy.dtype('c'):
                    variab = val[:].tostring().strip()
                else:
                    variab = val[:]
                if i == 0:
                    hstr += "'%s': %s,\n" % (key, variab)
                else:
                    hstr += ' '*(len(k)+4)+"'%s': %s,\n" % (key, variab)
        hstr += '}\n'
        return hstr

    def utdate(self):
        """Returns python datetime equivalent for TimePlace.UTDate"""
        import datetime
        utd = self.get_scalar('TimePlace.UTDate')
        year = int(utd)
        if year % 4 == 0:
            days_in_year = 366.
        else:
            days_in_year = 365.
        secs = 24.*3600.*days_in_year*(utd-year)
        dt = datetime.datetime(year,1,1)+datetime.timedelta(seconds=secs)
        return dt



    
#     def __repr__(self):
#         hstr = ''
#         for k, v in self.items():
#             if type(v) == types.DictType:
#                 hstr += "{\n%s : %s\n}" % (k, self._pprint_subhead(k, v))
#         return hstr
