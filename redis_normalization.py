import numpy
from netCDF4 import Dataset, Variable

def redis_normalization(filelist, writeloc):
    print "grabbing .npz files to reconstruct the OTF map"
    for i in range(len(filelist)):
        tempfile = filelist[i]
        tempfile = tempfile.strip(".nc")
        tempfile = tempfile+".npz"
        filename = writeloc + tempfile
        array_load = numpy.load(filename)
        if (i==0):
            MAX_WT = array_load[3]
        else:
            MAX_WT = numpy.maximum(MAX_WT, array_load[3])
        T += array_load[0]
        WT += array_load[1]
        TSYS += array_load[2]
        naxes = array_load[3]
        naxes2 = naxes[0]
        naxes1 = naxes[1]
        naxes0 = naxes[2]

    print "files all loaded, starting grid normalization"
    for iy in range(int(naxes2)):
        for ix in range(int(naxes1)):
            ii = ix + iy * naxes1
            if (WT[ii] != 0.0) and (MAX_WT[ii] > .25):
                TSYS[ii] = TSYS[ii] / WT[ii]
                jj = ii * naxes0
                if (T[jj+k] < 10**30):
                    for k in range(int(naxes0)):
                        T[jj+k] = T[jj+k] / WT[ii]
                else:
                    self.T[jj+k] = 0.0
                    print k
    T = T.reshape((naxes2, naxes1, naxes0))
    create_netcdf(filelist, T, WT, TSYS, naxes0, naxes1, naxes2)

def create_netcdf(filelist, T, WT, TSYS, naxes0, naxes1, naxes2):
    newfilename = (filelist[0].strip(".nc"))
    if (len(filelist)<7):
        newfilename = newfilename+"_smallgrid.nc"
    else:
        newfilename = newfilename+"_largegrid.nc"
    newfile = Dataset(newfilename, mode='w', clobber=True)
    newfile.createdimension('naxes0', naxes0)
    newfile.createdimension('naxes1', naxes1)
    newfile.createdimension('naxes2', naxes2)

    var1 = newfile.createVariable('otfmap', numpy.dtype(numpy.float32), (('naxes2', 'naxes1', 'naxes0')))
    var1[:] = T
    var1.__setattr__('filelist', filelist)

    var2 = newfile.createVariable('WT', numpy.dtype(numpy.float32), (('naxes2', 'naxes1')))
    var2[:] = WT

    var3 = newfile.createVariable('TSYS', numpy.dtype(numpy.float32), (('naxes2', 'naxes1')))
    var3[:] = TSYS
    
    var4 = newfile.createVariable('MAX_WT', numpy.dtype(numpy.float32), (('naxes2', 'naxes1')))
    
    newfile.close()
    
