
class LMTData(dict):
    """The most generic LMT Data class. This class is not usually
    instantiated by the user. It is part of the initialization of
    opening an LMT NetCDF file."""
    def __init__(self, ncvariables):
        #self.ncvariables = ncvariables
        self.make_data(ncvariables)

    def _make_data_keys(self, ncvariables):
        dnames = [name for name in ncvariables.keys() if name.find('Data') != -1]
        for dname in dnames:
            dtype =dname.split('.')[1]
            self[dtype] = {}
        for key in self.keys():
            for subdata in [n for n in ncvariables.keys() if n.startswith('Data.%s.' % key)]:
                sdata = subdata.split('.')[-1]
                self[key][sdata] = ncvariables[subdata]

    def make_data(self, ncvariables):
        self._make_data_keys(ncvariables)
