print "starting up a worker to process OTF data"
from redisworkq import WorkQueue
from gridmaker_dumps import gridmaker_dumps

def redisworkerbee(user, dataloc='null', writeloc='null'):
    wq = WorkQueue('otfworkq', redishost='cln.astro.umass.edu')
    if (user==1):
        dataloc = '/archives/fcrao/otfdata/'
        writeloc = '/archives/fcrao/otfdataout/'
    if (user==2):
        dataloc = '/otfraw/'
        writeloc = '/otftmp/'
    else:
        print "non-standard user specified, writing will probably fail"
    while (wq.queue_length() > 0):
        filename = wq.get_item()
        filename = filename
        filelist = [filename]

        gridmaker_dumps(525.0,535.0, 645.0, 655.0, filelist, dataloc=dataloc, writeloc=writeloc, normalize=False)

        print "if this prints, ", filename," was probably completed"
        wq.ack_item()

if __name__=="__main__":
    redisworkerbee(user=2)
