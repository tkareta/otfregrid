from redisworkq import WorkQueueManager
from time import sleep

def redisqueenbee(filelist):
    print "starting a REDIS workqueue to process OTF data"
    print "this is the host machine (presumably CLN)"
    wq = WorkQueueManager('otfworkq', redishost='cln.astro.umass.edu')
    wq.delete_queue() #this clears any pre-existing jobs
    for i in range(len(filelist)):
        wq.add_item(filelist[i])
    print "Number of Files: ", len(filelist)
    print "Initial Length of WorkQueue", wq.queue_length()
    continuing = True
    while (continuing):
        sleep(5*60) #sleep for 60 seconds
        if (wq.queue_length() > 0):
            continuing = True
            print "files still available to be grabbed"
        if (wq.queue_length() == 0):
            for item in wq.working_item_iterator():
                wq.fail_item_with_requeue()
            if(wq.queue_length() == 0):
                print "queue finished"
            else:
                "unacknowledged items re-added"
            
    


if __name__=="__main__":
    filelist = ['35065304.nc','35065306.nc','35065300.nc','35065301.nc','35065302.nc','35065303.nc','35065307.nc','35065308.nc','35065305.nc', '35065309.nc']
    redisqueenbee(filelist)

