import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import time

filenames = sys.argv

#print(filenames)

f=h5py.File(filenames[1],"r")
nevents=len(f.keys()) #number of events stored
print(nevents,"events")
events=list(f.items())
def mysort(val): return int(val[0])
events.sort(key=mysort)

xdata=[]
ydata=[]
starttime=time.time()
for evtindex in range(0,nevents):
    e=events[evtindex] #get the event
    if evtindex == 0: print("first event:", e)
    evtnum=int(e[0])
    data=e[1]
    nsamples=data.shape[-1]
    lendatashape=len(data.shape)
    etime=data.attrs.get("time")
    ysubdata=[]
    for chan in range(0,2):
        if lendatashape == 3:
            xdata = data[chan][0]
            ysubdata = data[chan][1]
            print(evtindex,"event",evtnum,"time",etime,"nsamp",nsamples,"chan",chan,xdata,ydata)
        elif lendatashape == 2:
            ysubdata = data[chan]
            #print(evtindex,"event",evtnum,"time",etime,"nsamp",nsamples,"chan",chan,xdata,ydata)
        else: print("unknown data shape length",lendatashape,"!")
        ydata.append(ysubdata)
    plt.figure(figsize=(7,5))
    plt.plot(xdata, ydata[-2], label="chan 0")
    plt.plot(xdata, ydata[-1], label="chan 1")
    ax = plt.gca()
    ax.set_xlabel('time in ns')
    ax.set_ylabel('signal strength in V')
    plt.grid()
    #plt.show()    
    plt.savefig(filenames[1]+'_event_'+str(evtnum)+'.jpg')
    plt.close()
