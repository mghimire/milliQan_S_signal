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

xdata		= []
bigdata		= []
smalldata	= []

starttime=time.time()
for evtindex in range(0,nevents):
	e=events[evtindex] #get the event
	#if evtindex == 0: print("first event:", e)
	evtnum=int(e[0])
	data=e[1]
	nsamples=data.shape[-1]
	lendatashape=len(data.shape)
	etime=data.attrs.get("time")
	for chan in range(0,4):
		if lendatashape == 3:
			xdata = data[chan][0]
			if chan == 0: smalldata.append(data[chan][1])
			if chan == 1: bigdata.append(data[chan][1])
			#print(evtindex,"event",evtnum,"time",etime,"nsamp",nsamples,"chan",chan,xdata,data[chan][1])
		elif lendatashape == 2:
			if chan == 0: smalldata.append(data[chan])
			if chan == 1: bigdata.append(data[chan])
			#print(evtindex,"event",evtnum,"time",etime,"nsamp",nsamples,"chan",chan,xdata,data[chan])
		else: print("unknown data shape length",lendatashape,"!")

xdata 		= np.asarray(xdata)
bigdata 	= np.asarray(bigdata)
smalldata 	= np.asarray(smalldata)
timeres		= xdata[1]-xdata[0]

print("temporal resolution is: ", timeres, "ns")

baseline 	= -3
peakthres	= -1


bigrms		= np.sqrt(np.mean((bigdata[0]- baseline)**2))
print('big rms is', bigrms)
bigthres	= bigrms*3 + baseline
print('big threshold is', bigthres)

smallrms	= np.sqrt(np.mean((smalldata[0]- baseline)**2))
print('small rms is', smallrms)
smallthres 	= smallrms*3 + baseline
print('small threshold is', smallthres)


monosgnl   	= []	#signals with a single peak
multisgnl   	= []	#signals with multiple peaks

for i in range(nevents):
	peaks = np.where(bigdata[i] > peakthres)
	if np.any(np.ediff1d(peaks) > 1300):
		multisgnl.append(i)
	else:
		monosgnl.append(i)

monosgnl 	= np.asarray(monosgnl)
print('There are', np.shape(monosgnl)[0], 'single signal events')
multisgnl 	= np.asarray(multisgnl)
print('There are', np.shape(multisgnl)[0], 'multiple signal events')

"""#examples of mono signal, multi signal and anomaly
plt.figure(figsize=(7,5))
plt.plot(xdata, bigdata[monosgnl[0]], 'y', xdata, smalldata[monosgnl[0]], 'b')
#plt.legend()
plt.grid()
#plt.show()
plt.savefig('mono.jpg')
plt.close()

for i in range(np.shape(multisgnl)[0]):
    plt.figure(figsize=(7,5))
    plt.plot(xdata, bigdata[multisgnl[i]], 'y', xdata, smalldata[multisgnl[i]], 'b')
    #plt.legend()
    plt.grid()
    #plt.show()
    plt.savefig('multi_'+str(i)+'.jpg')
    plt.close()"""

#single signal analysis
#----------------------
#idea is to compare the integrals under the peaks of the signal event and the expanded signal event

mononum		= np.shape(monosgnl)[0]
bigint		= np.zeros(mononum)
smallint	= np.zeros(mononum)

bigmax		= np.zeros(mononum)
smallmax	= np.zeros(mononum)

bigloc		= np.zeros(mononum)
smallloc 	= np.zeros(mononum)

j = 0 #counter variable for integral arrays
for i in monosgnl:
	bigpeaks 	= np.where(bigdata[i]>bigthres)
	bigint[j]	= np.sum(bigdata[i][bigpeaks])*timeres
	bigmax[j]	= np.amax(bigdata[i])
	bigloc[j]	= xdata[np.argmax(bigdata[i])]

	smallint[j]	= 0	
	smallpeaks	= np.where(smalldata[i]>smallthres)
	smallint[j]	= np.sum(smalldata[i][smallpeaks])*timeres
	smallmax[j]	= np.amax(smalldata[i])
	smallloc[j]	= xdata[np.argmax(smalldata[i])]

	j += 1


minpeakval 	= np.amin(bigmax)
maxpeakval	= np.amax(bigmax)

print('high voltage data saturates at', maxpeakval)

newpeaks	= np.where(bigmax > minpeakval + bigthres)
newmonosgnl	= monosgnl[newpeaks]
newbigint 	= bigint[newpeaks]
newsmallint	= smallint[newpeaks]
newbigmax	= bigmax[newpeaks]
newsmallmax	= smallmax[newpeaks]
newbigloc	= bigloc[newpeaks]
newsmallloc	= smallloc[newpeaks]

print('There are', np.shape(newpeaks)[1], 'filtered events')

plt.figure(figsize=(7,5))
plt.scatter(smallint, bigint, marker='.')
ax = plt.gca()
ax.set_xlabel('integral around small signal peak')
ax.set_ylabel('integral around big signal peak')
plt.grid()
#plt.show()
plt.savefig('peakint_correlation.jpg')
plt.close()

plt.figure(figsize=(7,5))
plt.scatter(smallint, bigint, marker='.')
plt.plot(np.unique(smallint), np.poly1d(np.polyfit(smallint, bigint, 1))(np.unique(smallint)), 'r')
ax = plt.gca()
ax.set_xlabel('integral around small signal peak')
ax.set_ylabel('integral around big signal peak')
plt.grid()
plt.text(0.7, 0.1,'m = ' + '%.5f' % np.polyfit(smallint, bigint, 1)[0] + ' b = ' + '%.5f' % np.polyfit(smallint, bigint, 1)[1],
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)
#plt.show()
plt.savefig('peakint_lobf.jpg')
plt.close()

plt.figure(figsize=(7,5))
plt.scatter(smallint, bigint, marker='.')
plt.plot(np.unique(smallint), np.poly1d(np.polyfit(smallint, bigint, 2))(np.unique(smallint)), 'r')
ax = plt.gca()
ax.set_xlabel('integral around small signal peak')
ax.set_ylabel('integral around big signal peak')
plt.grid()
plt.text(0.7, 0.1,'a = ' + '%.5f' % np.polyfit(smallint, bigint, 2)[0] + ' b = ' + '%.5f' % np.polyfit(smallint, bigint, 2)[1] + ' c = ' + '%.5f' % np.polyfit(smallint, bigint, 2)[2], horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
#plt.show()
plt.savefig('peakint_pobf.jpg')
plt.close()

plt.figure(figsize=(7,5))
plt.scatter(newsmallint, newbigint, marker='.')
plt.plot(np.unique(newsmallint), np.poly1d(np.polyfit(newsmallint, newbigint, 2))(np.unique(newsmallint)), 'r')
ax = plt.gca()
ax.set_xlabel('integral around small signal peak')
ax.set_ylabel('integral around big signal peak')
plt.grid()
plt.text(0.7, 0.1,'a = ' + '%.5f' % np.polyfit(newsmallint, newbigint, 2)[0] + ' b = ' + '%.5f' % np.polyfit(newsmallint, newbigint, 2)[1] + ' c = ' + '%.5f' % np.polyfit(newsmallint, newbigint, 2)[2], horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
#plt.show()
plt.savefig('newpeakint_pobf.jpg')
plt.close()

plt.figure(figsize=(7,5))
plt.scatter(smallmax, bigmax, marker='.')
ax = plt.gca()
ax.set_xlabel('highest value of small signal')
ax.set_ylabel('highest value of big signal')
plt.grid()
#plt.show()
plt.savefig('peaks_correlation.jpg')
plt.close()

plt.figure(figsize=(7,5))
plt.scatter(newsmallmax, newbigmax, marker='.')
ax = plt.gca()
ax.set_xlabel('highest value of small signal')
ax.set_ylabel('highest value of big signal')
plt.grid()
#plt.show()
plt.savefig('newpeaks_correlation.jpg')
plt.close()

plt.figure(figsize=(7,5))
plt.scatter(smallloc, bigloc, marker='.')
ax = plt.gca()
ax.set_xlabel('location of small signal peak')
ax.set_ylabel('location of big signal peak')

j = 0


plt.grid()
# #plt.show()
plt.savefig('peakloc_correlation.jpg')
plt.close()

plt.figure(figsize=(7,5))
plt.scatter(newsmallloc, newbigloc, marker='.')
plt.plot(np.unique(newsmallloc), np.poly1d(np.polyfit(newsmallloc, newbigloc, 1))(np.unique(newsmallloc)), 'r')
ax = plt.gca()
ax.set_xlabel('location of small signal peak')
ax.set_ylabel('location of big signal peak')

j = 0

plt.grid()
plt.text(0.7, 0.1,'m = ' + '%.5f' % np.polyfit(newsmallloc, newbigloc, 1)[0] + ' b = ' + '%.5f' % np.polyfit(newsmallloc, newbigloc, 1)[1],
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)
#plt.show()
plt.savefig('newpeakloc_lobf.jpg')
plt.close()

"""trnsgrsrs = []

def error(i):
	return newbigloc[i] - np.polyfit(newsmallloc, newbigloc, 1)[0]*newsmallloc[i] - np.polyfit(newsmallloc, newbigloc, 1)[1]

for i in range(np.shape(newpeaks)[1]):
	err = error(i)
	if np.abs(err) > 2:
		trnsgrsrs.append([i,err])
		plt.figure(figsize=(7,5))
		plt.plot(time, expndata[monosgnl[i]], 'y', time, sgnldata[monosgnl[i]], 'b')
 		#plt.legend()
		plt.grid()
 		#plt.show()
		plt.savefig(filenames[newmonosgnl[i]] + '.jpg')
		plt.close()

trnsgrsrs = np.asarray(trnsgrsrs).astype(int)

#print(trnsgrsrs[:,0])

cleanbigint = np.delete(newbigint, trnsgrsrs[:,0], 0)
cleansmallint = np.delete(newsmallint, trnsgrsrs[:,0], 0)
cleanbigloc = np.delete(newbigloc, trnsgrsrs[:,0], 0)
cleansmallloc = np.delete(newsmallloc, trnsgrsrs[:,0], 0)
cleanmonosgnl = np.delete(newmonosgnl, trnsgrsrs[:,0], 0)

plt.figure(figsize=(7,5))
plt.scatter(cleansmallint, cleanbigint, marker='.')
ax = plt.gca()
ax.set_xlabel('integral around small signal peak without anomalies')
ax.set_ylabel('integral around big signal peak without anomalies')
plt.grid()
#plt.show()
plt.savefig('cleanint_correlation.jpg')
plt.close()

plt.figure(figsize=(7,5))
plt.scatter(cleansmallint, cleanbigint, marker='.')
ax = plt.gca()
ax.set_xlabel('integral around small signal peak without anomalies')
ax.set_ylabel('integral around big signal peak without anomalies')

j = 0
for i in cleanmonosgnl:
	ax.annotate(Original[i][0], (cleansmallint[j], cleanbigint[j]))
	j += 1

plt.grid()
#plt.show()
plt.savefig('labeled_cleanint_correlation.jpg')
plt.close()

plt.figure(figsize=(7,5))
plt.scatter(cleansmallint, cleanbigint, marker='.')
plt.plot(np.unique(cleansmallint), np.poly1d(np.polyfit(cleansmallint, cleanbigint, 1))(np.unique(cleansmallint)), 'r')
ax = plt.gca()
ax.set_xlabel('integral around small signal peak')
ax.set_ylabel('integral around big signal peak')
plt.grid()
plt.text(0.7, 0.1,'m = ' + '%.5f' % np.polyfit(cleansmallint, cleanbigint, 1)[0] + ' b = ' + '%.5f' % np.polyfit(cleansmallint, cleanbigint, 1)[1],
      horizontalalignment='center',
      verticalalignment='center',
      transform = ax.transAxes)
#plt.show()
plt.savefig('cleanint_lobf.jpg')
plt.close()

# plt.figure(figsize=(7,5))
# plt.scatter(cleansmallint, cleanbigint, marker='.')
# plt.plot(np.unique(cleansmallint), np.poly1d(np.polyfit(cleansmallint, cleanbigint, 2))(np.unique(cleansmallint)), 'r')
# ax = plt.gca()
# ax.set_xlabel('integral around small signal peak')
# ax.set_ylabel('integral around big signal peak')
# plt.grid()
# plt.text(0.7, 0.1,'a = ' + '%.5f' % np.polyfit(cleansmallint, cleanbigint, 2)[0] + ' b = ' + '%.5f' % np.polyfit(cleansmallint, cleanbigint, 2)[1] + ' c = ' + '%.5f' % np.polyfit(cleansmallint, cleanbigint, 2)[2], horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
# #plt.show()
# plt.savefig('cleanint_pobf.jpg')
# plt.close()

# plt.figure(figsize=(7,5))
# plt.scatter(cleansmallloc, cleanbigloc, marker='.')
# plt.plot(np.unique(cleansmallloc), np.poly1d(np.polyfit(cleansmallloc, cleanbigloc, 1))(np.unique(cleansmallloc)), 'r')
# ax = plt.gca()
# ax.set_xlabel('location of small signal peak')
# ax.set_ylabel('location of big signal peak')

# j = 0
# for i in cleanmonosgnl:
# 	ax.annotate(Original[i][0], (cleansmallloc[j], cleanbigloc[j]))
# 	j += 1

# plt.grid()
# plt.text(0.7, 0.1,'m = ' + '%.5f' % np.polyfit(cleansmallloc, cleanbigloc, 1)[0] + ' b = ' + '%.5f' % np.polyfit(cleansmallloc, cleanbigloc, 1)[1],
#      horizontalalignment='center',
#      verticalalignment='center',
#      transform = ax.transAxes)
# #plt.show()
# plt.savefig('cleanpeakloc_lobf.jpg')
# plt.close()

locdiff = cleansmallloc - cleanbigloc

# plt.figure(figsize=(7,5))
# plt.scatter(cleansmallint, locdiff, marker='.')
# plt.plot(np.unique(cleansmallint), np.poly1d(np.polyfit(cleansmallint, locdiff, 1))(np.unique(cleansmallint)), 'r')
# ax = plt.gca()
# ax.set_xlabel('integral around small signal peak')
# ax.set_ylabel('small signal peak location - big signal peak location')

# #j = 0
# #for i in monosgnl:
# #	ax.annotate(Original[i][0], (smallloc[j], bigloc[j]))
# #	j += 1

# plt.grid()
# plt.text(0.7, 0.1,'m = ' + '%.5f' % np.polyfit(cleansmallint, locdiff, 1)[0] + ' b = ' + '%.5f' % np.polyfit(cleansmallint, locdiff, 1)[1],
#      horizontalalignment='center',
#      verticalalignment='center',
#      transform = ax.transAxes)
# #plt.show()
# plt.savefig('cleansmallint_peaklocdiff_correlation.jpg')
# plt.close()

plt.figure(figsize=(7,5))
plt.hist(locdiff, bins=100, label = 'error distribution')
ax = plt.gca()
ax.set_xlabel('difference between peak locations (small peak location - big peak location)')
plt.grid()
plt.text(0.7, 0.1,'mean = ' + '%.5f' % np.mean(locdiff) + ' std dev = ' + '%.5f' % np.std(locdiff),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)
#plt.show()
plt.savefig('timeres_hist.jpg')
plt.close()

cleanerr = sgnlint - (expnint + 226.62)/11.43

plt.figure(figsize=(7,5))
plt.hist(cleanerr, bins=100, label = 'error distribution')
ax = plt.gca()
ax.set_xlabel('difference between predicted and actual area under signal')
plt.grid()
plt.text(0.7, 0.1,'mean = ' + '%.5f' % np.mean(cleanerr) + ' std dev = ' + '%.5f' % np.std(cleanerr),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)
#plt.show()
plt.savefig('err_hist.jpg')
plt.close()

# sgnlpksval	= np.zeros(datanum)
# expnpksval	= np.zeros(datanum)
# sgnlpksloc	= np.zeros(datanum)
# expnpksloc	= np.zeros(datanum)

# for i in range(datanum):
# 	sgnlpksval[i] = np.amax(data[i][:,2])
# 	expnpksval[i] = np.amax(data[i][:,1])
# 	sgnlpksloc[i] = data[i][np.argmax(data[i][:,2]),0]
# 	expnpksloc[i] = data[i][np.argmax(data[i][:,1]),0]
     


# sgnlpksloc = sgnlpksloc*1e7
# expnpksloc = expnpksloc*1e7

# plt.figure(figsize=(7,5))
# plt.scatter(sgnlpksval, expnpksval)
# ax = plt.gca()
# ax.set_xlabel('value of original signal peak')
# ax.set_ylabel('value of stretched signal peak')
# plt.grid()
# #plt.show()
# plt.savefig('peakval_correlation.pdf')
# plt.close()

# plt.figure(figsize=(7,5))
# plt.scatter(sgnlpksval, expnpksval)
# plt.plot(np.unique(sgnlpksval), np.poly1d(np.polyfit(sgnlpksval, expnpksval, 1))(np.unique(sgnlpksval)), 'r')
# ax = plt.gca()
# ax.set_xlabel('value of original signal peak')
# ax.set_ylabel('value of stretched signal peak')
# plt.grid()
# plt.text(0.7, 0.1,'m = ' + '%.5f' % np.polyfit(sgnlpksval, expnpksval, 1)[0] + ' b = ' + '%.5f' % np.polyfit(sgnlpksval, expnpksval, 1)[1],
#      horizontalalignment='center',
#      verticalalignment='center',
#      transform = ax.transAxes)
# #plt.show()
# plt.savefig('peakval_lobf.pdf')
# plt.close()

# plt.figure(figsize=(7,5))
# plt.scatter(sgnlpksval, expnpksval)
# plt.plot(np.unique(sgnlpksval), np.poly1d(np.polyfit(sgnlpksval, expnpksval, 2))(np.unique(sgnlpksval)), 'r')
# ax = plt.gca()
# ax.set_xlabel('value of original signal peak')
# ax.set_ylabel('value of stretched signal peak')
# plt.grid()
# plt.text(0.7, 0.1,'a = ' + '%.5f' % np.polyfit(sgnlpksval, expnpksval, 2)[0] + ' b = ' + '%.5f' % np.polyfit(sgnlpksval, expnpksval, 2)[1] + ' c = ' + '%.5f' % np.polyfit(sgnlpksval, expnpksval, 2)[2],
#      horizontalalignment='center',
#      verticalalignment='center',
#      transform = ax.transAxes)
# #plt.show()
# plt.savefig('peakval_pobf.pdf')
# plt.close()

# plt.figure(figsize=(7,5))
# plt.scatter(sgnlpksloc, expnpksloc)
# ax = plt.gca()
# ax.set_xlabel('location of original signal peak')
# ax.set_ylabel('location of stretched signal peak')
# plt.grid()
# #plt.show()
# plt.savefig('peakloc_correlation.pdf')
# plt.close()

# outliers = np.argwhere(sgnlpksloc > -2)

# sgnlpksloc = np.delete(sgnlpksloc, outliers)
# expnpksloc = np.delete(expnpksloc, outliers)

# plt.figure(figsize=(7,5))
# plt.scatter(sgnlpksloc, expnpksloc)
# plt.plot(np.unique(sgnlpksloc), np.poly1d(np.polyfit(sgnlpksloc, expnpksloc, 1))(np.unique(sgnlpksloc)), 'r')
# ax = plt.gca()
# ax.set_xlabel('location of original signal peak')
# ax.set_ylabel('location of stretched signal peak')
# plt.text(0.7, 0.1,'m = ' + '%.5f' % np.polyfit(sgnlpksloc, expnpksloc, 1)[0] + ' b = ' + '%.5f' % np.polyfit(sgnlpksloc, expnpksloc, 1)[1],
#      horizontalalignment='center',
#      verticalalignment='center',
#      transform = ax.transAxes)
# plt.grid()
# #plt.show()
# plt.savefig('peakloc_lobf_no_outliers.pdf')
# plt.close()"""
