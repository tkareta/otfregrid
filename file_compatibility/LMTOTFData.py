from lmtnetcdffile import LMTNetCDFFile
from collections import OrderedDict

class LMTOTFData_old(LMTNetCDFFile):
    def __init__(self):
        self.make_data()
    
    def make_data(self):
        setattr(self.hdu.data, 'dump', self.get('OldOTF.dump'))
        setattr(self.hdu.data, 'refs', self.get('OldOTF.refs'))

