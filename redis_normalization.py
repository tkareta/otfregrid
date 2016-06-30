import numpy
from netCDF4 import Dataset, Variable
import os

def redis_normalization(filelist, writeloc):
    print "grabbing .npz files to reconstruct the OTF map"
    for i in range(len(filelist)):
        tempfile = filelist[i]
        tempfile = tempfile.strip(".nc")
        
        tempfile = tempfile+".npz"
        
        filename = writeloc + tempfile
        
        array_load = numpy.load(filename)
        T = 0
        WT = 0
        TSYS = 0
        
        #print len(array_load['arr_0'])
        if (i==0):
            MAX_WT = array_load['arr_3']
        else:
            MAX_WT = numpy.maximum(MAX_WT, array_load['arr_3'])
        T += array_load['arr_0']
        WT += array_load['arr_1']
        TSYS += array_load['arr_2']
        naxes = array_load['arr_4']
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
                for k in range(int(naxes0)):
                    if (T[jj+k] < 10**30):
                        T[jj+k] = T[jj+k] / WT[ii]
                    else:
                        T[jj+k] = 0.0
                        print k
        
    T = T.reshape((naxes2, naxes1, naxes0))
    create_netcdf(filelist, T, WT, TSYS, naxes0, naxes1, naxes2, writeloc)

def create_netcdf(filelist, T, WT, TSYS, naxes0, naxes1, naxes2, writeloc):
    newfilename = (filelist[0].strip(".nc"))
    if (len(filelist)<7):
        newfilename = newfilename+"_smallgrid.nc"
    else:
        newfilename = newfilename+"_largegrid.nc"
    newfilename = writeloc + newfilename
    if (os.path.exists(newfilename)):
        os.unlink(newfilename)
        print "Deleting pre-existing grid."

    newfile = Dataset(newfilename, mode='w', clobber=True)
    newfile.createDimension('naxes0', naxes0)
    newfile.createDimension('naxes1', naxes1)
    newfile.createDimension('naxes2', naxes2)

    var1 = newfile.createVariable('otfmap', numpy.dtype(numpy.float32), (('naxes2', 'naxes1', 'naxes0')))
    var1[:] = T
    var1.__setattr__('filelist', filelist)

    var2 = newfile.createVariable('WT', numpy.dtype(numpy.float32), (('naxes2', 'naxes1')))
    var2[:] = WT

    var3 = newfile.createVariable('TSYS', numpy.dtype(numpy.float32), (('naxes2', 'naxes1')))
    var3[:] = TSYS
    
    var4 = newfile.createVariable('MAX_WT', numpy.dtype(numpy.float32), (('naxes2', 'naxes1')))
    
    newfile.close()
    

if __name__=="__main__":
    filelist = ['35065304.nc','35065306.nc','35065300.nc','35065301.nc','35065302.nc','35065303.nc','35065307.nc','35065308.nc','35065305.nc', '35065309.nc']
    #redis_normalization(filelist,'/archives/fcrao/otfdataout/')
    redis_normalization(filelist, '/otftmp/')
