print "starting up a worker to process OTF data"
from redisworkq import WorkQueue
from gridmaker_dumps import gridmaker_dumps

wq = WorkQueue('otfworkq', redishost='cln.astro.umass.edu')

while (wq.queue_length() > 0):
    filename = wq.get_item()
    filelist = [filename]

    gridmaker_dumps(525.0,535.0, 645.0, 655.0, filelist, normalize=False)

    print "if this prints, ", filename," was probably completed"
    wq.ack_item()
