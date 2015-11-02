from lmtheader import LMTHeader
from collections import OrderedDict


variablelist = OrderedDict((
        # oldname, newname, is_a_scalar
        ('otflab', (1, False)),
        ('name', ('sourcename', False)),
        ('bunit', (1, False)),
        ('vlabel', (1, False)),
        ('xlabel', (1, False)),
        ('ylabel', (1, False)),
        ('observer', (1, False)),
        
        ('utd', (1, True)),
        ('f0', (1, True)),
        ('fsky', (1, True)),
        ('azs', (1, True)),
        ('els', (1, True)),
        ('XCEN', (1, True)),
        ('YCEN', (1, True)),
        ('xvel1', (1, True)),
        ('xfr1', (1, True)),
        ('dveldin', (1, True)),
        ('dfrdin', (1, True)),
        ('fr1', (1, True)),
        ('fr2', (1, True)),
        
        ('tsid', (1, True)),
        ('utc', (1, True)),
        ('targetrms', (1, True)),
        ('tdump', (1, True)),
        ('tref', (1, True)),
        ('telrate', (1, True)),
        ('epoch', (1, True)),
        ('bw', (1, True)),
        ('vlsr', (1, True)),
        ('vclsr', (1, True)),
        ('vcsun', (1, True)),
        ('NOS', (1, True)),
        ('tau', (1, True)),
        ('h2omm', (1, True)),
        ('pamb', (1, True)),
        ('tamb', (1, True)),
        ('rotang', (1, True)),
        ('XSIZE', (1, True)),
        ('YSIZE', (1, True)),
        ('XOFFMAP', (1, True)),
        ('YOFFMAP', (1, True)),
        ('XRAMP', (1, True)),
        ('YRAMP', (1, True)),
        
        ('iscan', (1, True)),
        ('nhorns', (1, True)),
        ('nchan', (1, True)),
        ('nsample', (1, True)),
        ('nblock', (1, True)),
        ('blksize', (1, True)),
        ('telcoord', (1, True)),
        ('kx', (1, True)),
        ('ky', (1, True)),
        ('klo', (1, True)),
        ('khm', (1, True)),
        ('harm', (1, True)),
        ('krecv', (1, True)),
        ('khorn', (1, True)),
        ('caltype', (1, True)),
        ('calscan', (1, True)),
        ('cfact', (1, True)),
        ('etab', (1, True)),
        ('kmap', (1, True)),
        ('raoff', (1, True)),
        ('decoff', (1, True)),
        ('fsrate', (1, True)),
        ('stype', (1, True)),
        ('iconf', (1, True)),
        ('tsys', (1, False)),
        ('gain', (1, False)),
        ('DXHORNS', (1, False)),
        ('DYHORNS', (1, False)),
        ('hornID', (1, False)),
        
        ('ANGLE', (1, False)),
        ('XOFF', (1, False)),
        ('YOFF', (1, False))
        ))

class LMTOTFHeader_old(LMTHeader):
    def __init__(self, ncvariables=None, dimensions=None, fromold=False):
        super(LMTOTFHeader_old, self).__init__(ncvariables=ncvariables,
                                               dimensions=dimensions)
        self.make_header()

    def make_header(self):
        #this takes values from the LMTNetCDFFile
        #and arranges them in an OTF-oriented
        #and simplified way
        #
        #starting with variables, values
        for oldname, (newname, is_a_scalar) in variablelist.items():
            if (newname == 1):
                newname = oldname
            oldname = "OldOTF." + oldname
            if (is_a_scalar):
                setattr(self, newname, self.get_scalar(oldname))
            else:
                setattr(self, newname, self.get(oldname))
