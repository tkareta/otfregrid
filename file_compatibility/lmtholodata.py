from lmtdata import LMTData
import numpy

class LMTHoloData(LMTData):
    """
    This is the main class for the LMT Holography
    Data type object"""
    def __init__(self, ncvariables):
        #LMTData.__init__(self, ncvariables)
        self.make_data(ncvariables)

    def make_data(self, ncvariables):
        """Make the actual holography data"""
        datas = [name for name in ncvariables.keys() if name.find('Data') == 0]
        for d in datas:
            data_type = d.split('.')[1]
            self[data_type] = []
        for data_type in self.keys():
            keys = [name.split('.')[-1] for name in ncvariables.keys() if name.find('Data.%s.' % data_type) == 0]
            for key in keys:
                try:
                    #self.__setattr__(key, ncvariables['Data.%s.%s' % (data_type, key)].get())
                    if hasattr(self, key):
                        # there  is already a data attribute of same name so let's
                        # change how we write this
                        #print "Key %s already present" % key
                        self.__setattr__("%s%s" % (data_type, key), ncvariables['Data.%s.%s' % (data_type, key)][:])
                    else:
                        #print "Adding %s to data" % key
                        self.__setattr__(key, ncvariables['Data.%s.%s' % (data_type, key)][:])
                    self[data_type].append(key)
                except:
                    self.__setattr__(key, None)

    def update_data(self, nc):
        """
        Updates data changed to the netcdf Dataset instance
        nc
        """
        datas = [name for name in nc.variables.keys() if name.find('Data') == 0]
        datavar_name = []
        for d in datas:
            datavar_name.append(d.split('.')[-1])
        for i, ncvar in enumerate(datas):
            nc.variables[ncvar][:] = getattr(self, datavar_name[i])
        nc.sync()

    def TotalPower(self, off=1):
        x = self.AbAmp*numpy.cos(self.AbPhase*numpy.pi/180.)
        y = self.AbAmp*numpy.sin(self.AbPhase*numpy.pi/180.)
        idx = numpy.where(self.BufPos == 0)   # Find the ONpositions
        xmean1 = x[idx].mean()
        xstd1 = x[idx].std()
        ymean1 = y[idx].mean()
        ystd1 = y[idx].std()
        amp1 = numpy.sqrt(xmean1**2+ymean1**2)/self.BPower[idx].mean()
        samp1 = numpy.sqrt((xmean1*xstd1)**2+(ymean1*ystd1)**2)/amp1
        if off == 0:
            jdx = numpy.where(self.BufPos == 1)    # Find Offpositions
            xmean0 = x[jdx].mean()
            xstd0 = x[jdx].std()
            ymean0 = y[jdx].mean()
            ystd0 = y[jdx].std()
            amp0 = numpy.sqrt(xmean0**2+ymean0**2)/self.BPower[jdx].mean()
            samp0 = numpy.sqrt((xmean0*xstd0)**2+(ymean0*ystd0)**2)/amp0
            amp = amp1-amp0
            samp = numpy.sqrt(samp1**2+samp0**2)
        else:
            amp = amp1
            samp = samp1
        return amp, samp
      
