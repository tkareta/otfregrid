from redisworkq import WorkQueueManager

def otfredishost(filelist):
    print "starting a REDIS workqueue to process OTF data"
    print "this is the host machine (presumably CLN)"
    wq = WorkQueueManager('otfworkq', redishost='cln.astro.umass.edu')
    for i in range(len(filelist)):
        wq.add_item(filelist[i])
    print "Number of Files: ", len(filelist)
    print "Initial Length of WorkQueue", WorkQueueManager.queue_length()
    


if __name__=="__main__":
    filelist = ['35065304.nc','35065306.nc','35065300.nc','35065301.nc','35065302.nc','35065303.nc','35065307.nc','35065308.nc','35065305.nc', '35065309.nc']
    otfredishost(filelist)

