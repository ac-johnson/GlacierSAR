#Melt utilities

import numpy as np
from datetime import datetime,timedelta
import asf_search as asf
import fnmatch
from matplotlib import pyplot as plt
from osgeo import gdal
import random



def timediff(time1,time2,form='seconds'):
    '''time difference between time2 and time1 expressed as either days or seconds'''
    # dd = 1
    # if time1>time2:
    #     dd=1
    diff = time2-time1
    if form=='seconds':
        diff = diff.seconds+diff.days*86400
    elif form=='days':
        diff = diff.seconds/86400+diff.days
    return diff

def getasfresults(WTKstr,year,trimresults=True,winterN=4,showplots=False,showstats=False,listpaths=False,maxsize=None):
    """trim results removes images if there are less than winterN number of winter scenes to establish the baseline"""
    results = asf.geo_search(platform=[asf.PLATFORM.SENTINEL1], intersectsWith=WTKstr, maxResults=1000,start=f'{str(year)}-01-01',
                         end=f'{str(year+1)}-01-01',processingLevel='SLC',beamMode='IW')

    datelist = []
    pflist = [] #path-frame list
    for res in results:
        res = res.properties
        # print(res['startTime'])
        rdate = datetime.strptime(res['startTime'],'%Y-%m-%dT%H:%M:%S.%fZ')
        rpf = str(res['pathNumber'])+'-'+str(res['frameNumber'])
        # print(rdate)
        res['datetime']=rdate
        res['pfval'] = rpf
        datelist.append(rdate)
        pflist.append(rpf)

    pfunique,count = np.unique(pflist,return_counts=True)
    pfvals = pflist
    pflist = [np.where(pfunique==i)[0][0] for i in pflist]

    
    if showplots==True:
        plt.figure()
        plt.plot(datelist,pflist,'.',markersize=12)
        
    trem = 0
    pfremove = []
    if trimresults==True:
        wstart = datetime(year,1,1,0,0,0)
        wend = datetime(year,4,1,1,0,0,0)
        pfremove=[]
        trem = 0
        #identify path-frames that do not meet requirements
        for pf in pfunique:
            eles = [i for i in results if i.properties['pfval']==pf]
            datelist = np.array([i.properties['datetime'] for i in eles])
            dl2 = datelist[datelist > wstart]
            dl2 = dl2[dl2 < wend]
            if len(dl2)<winterN:
                pfremove.append(pf)
                trem += len(datelist)

        for res in results:
            if res.properties['pfval'] in pfremove:
                results.remove(res)
    
    if maxsize:
        pfremoveold = []
        pfremove = []
        tremstart = trem
        pfuniquesort,countsort = np.array(pfunique),np.array(count)

        countsort,pfuniquesort = zip(*sorted(zip(countsort, pfuniquesort)))

        pfuniquesort = np.flip(np.array(pfuniquesort))
        countsort = np.flip(np.array(countsort))


        if len(pfremoveold)>0:
            for pf in pfuniquesort:
                countsort = np.remove(countsort,pfuniquesort==pf)
                pfuniquesort = np.remove(pfuniquesort,pfuniquesort==pf)
                
        j = 0
        currentsize = 0
        while currentsize+countsort[j]<=maxsize:
            currentsize+=countsort[j]
            j+=1
        for i in range(j+1,len(pfuniquesort)):
            # print(i)
            pfremove.append(pfuniquesort[i])
            trem+=countsort[i]
        
        if trem>tremstart:
            maxcount = currentsize

            oldrlen = len(results)
            print(len(pfunique)-len(pfremove))
            while oldrlen>maxcount:
                k=0
                while k < len(results)-1:
                    for k,res in enumerate(results[k:],start=k):
                        if res.properties['pfval'] in pfremove:
                            results.remove(res)
                            break
                oldrlen = len(results)
                if len(results)==oldrlen:
                    break
            print(len(results))

                  
    if trem>0:        
        print(f'total removed: {trem}')
        # print(pfremove)
        
        #now plot again:
        datelist = [i.properties['datetime'] for i in results]
        pflist = [i.properties['pfval'] for i in results]
        pfunique = np.unique(pflist)
        pflist = [np.where(pfunique==i)[0][0] for i in pflist]

        if showplots==True:
            plt.figure()
            plt.plot(datelist,pflist,'.',markersize=12)
        
    if showstats==True:
        dlsort = [i for i in datelist]
        dlsort.sort()
        dldelta = [timediff(dlsort[i],dlsort[i+1],form='days') for i in range(len(dlsort)-1)]
        dldelta = np.array(dldelta)

        print('Acquisition stats')
        print(f'Number of acquisitions:       {len(dldelta)+1}')
        print(f'Number of images >1min apart: {len(dldelta[dldelta>1/1440])+1}\n')
        print(f'Mean days apart:   {np.mean(dldelta):.3f}')
        print(f'Mean days >1min:   {np.mean(dldelta[dldelta>1/1440]):.3f}')
        print(f'Median days apart: {np.median(dldelta):.3f}')
        print(f'Median days >1min: {np.median(dldelta[dldelta>1/1440]):.3f}')
        print(f'Stdev days apart:  {np.std(dldelta):.3f}')
        print(f'Longest gap:       {np.max(dldelta):.3f}')
        print('Note: these are over the whole area, not necessarily by pixel')
    
    if listpaths==True:
        for pf in pfunique:
            print(f'Path-frame: {pf}, count: {len([i for i in pfvals if i==pf])}')
    
    return results



def getSARdate(file,returndatestring=False):
    '''Obtains the date of the start of a SAR acquisition, returns as a datetime
    or a string.'''

    file = str(file)
    query = "_????????T??????_"
    seqquery = str(query)
    seqquery.replace('?','[0-9]')

    substrings = fnmatch.filter((file[i:i+len(query)] for i in range(len(file) - len(query))), seqquery)
    dtstr = substrings[0][1:-1].split('T')
    
    dtint = [dtstr[0][0:4],dtstr[0][4:6],dtstr[0][6:8],dtstr[1][0:2],dtstr[1][2:4],dtstr[1][4:6]]
    dtint = [int(i)for i in dtint]
    dtimg = datetime(dtint[0],dtint[1],dtint[2],dtint[3],dtint[4],dtint[5])

    if returndatestring:
        dtstr = dtstr[0]+dtstr[1]
        return dtimg,dtstr
    else:
        return dtimg

def mapbin(dem,ibin):
    return np.logical_and(dem>=ibin[0],dem<ibin[1])

def get_rgb(imgfilecp,imgfilexp,usept=None):
    """Creates an RGB image from co- and cross-polarized images"""
    imgcp = gdal.Open(str(imgfilecp)).ReadAsArray()
    imgxp = gdal.Open(str(imgfilexp)).ReadAsArray()

    
    if usept != None:
        ipt,jpt = usept
        imgcp = np.array([[imgcp[ipt,jpt]]])
        imgxp = np.array([[imgxp[ipt,jpt]]])
    
    rgb = hyp3_rgb(imgcp, imgxp)
    rgb[rgb>1]=1.0

    imgcp = None
    imgxp = None
    
    return rgb

def hyp3_rgb(imgcp,imgxp,k=-24):
    """Runs the algorithm to create an RGB image after files have been opened and prepared"""
    imgcp = 10**(imgcp/10)
    imgxp = 10**(imgxp/10)
    k = 10**(k/10)
    
    maskb = imgxp < k
    maskr = imgxp >= k
    maskx = imgxp > 0
    
    ps = imgcp-3*imgxp
    pr = np.array(ps)
    pb = -np.array(ps)
    pr[pr<0]=0
    pb[pb<0]=0

    imgdiff = imgcp-imgxp
    imgdiff[imgdiff<0]=0
    z = 2/np.pi*maskb*np.arctan(imgdiff**0.5)
    
    rgbscale = 1.0 #for 1 to 255, set these to 254 and 1
    rgbn = 0

    ar = rgbscale * maskx * (2 * maskr * pr**0.5 + z) + rgbn
    ag = rgbscale * maskx * (3 * maskr * imgxp**0.5 + 2*z) + rgbn
    ab = rgbscale * maskx * (2 * pb**0.5 + 5*z) + rgbn
    # ab[~maskb]=0
    # ag[~maskr]=0
    # ar[~maskr]=0
    
    rgb = np.stack((ar,ag,ab),axis=2)
    return rgb


def meltseason(meltvec,datelist):
    """will return the start date and end date of each melt season, as well as the uncertainty of each"""
    if type(datelist)==tuple:
        datelist = np.array(datelist)
    
    useind = ~np.isnan(meltvec)
    # print(useind)
    meltvec = meltvec[useind]
    datelist = datelist[useind]
    
    mdiff = meltvec[1:]-meltvec[:-1]
    dstart = []
    dstartstd = []
    dend = []
    dendstd = []
    
    mstartind = np.where(mdiff==1)

    mstartind = [i+1 for i in mstartind[0]]
    # print(mstartind)
    mlen = len(meltvec)
    
    
    while len(mstartind)>0:
        mstart = mstartind[0]
        # print(mstartind)
        mstartind.remove(mstart)

        d1 = datelist[mstart-1]
        # print(mstart)
        # print(len(datelist))
        d2 = datelist[mstart]
        dstart.append( d1+(d2-d1)/2)
        # print(d2-d1)
        # print(d1)
        dstartstd.append((d2-d1).days*0.2886751345) #This is 1/sqrt(12)        
    
        cind = mstart
        lastone = mstart #the last one that was seen. Not the last zero or last two.
        cdstart = dstart[-1]
        
        while cind < mlen-1:
            cind+=1
            # print(cind)
            if meltvec[cind]==1.0:
                lastone = cind
                # print('still melt')
                if cind in mstartind:
                    mstartind.remove(cind)
            if meltvec[cind] == 0:
                ddiff = (datelist[cind]-datelist[lastone]).days

                if ddiff >= 12:
                    d1 = datelist[lastone]
                    d2 = datelist[lastone+1]
                    d12diff = (d2-d1).days+(d1-d2).seconds/86400
                    dend.append( d1+timedelta(days=(d12diff)/2))
                    dendstd.append(d12diff*0.2886751345) #This is 1/sqrt(12)
                    break

            if cind==mlen-1:
                if meltvec[cind]==1:
                    dendstd.append(np.nan)
                    dend.append(datetime(datelist[-1].year+1,1,1))

                else:

                    d1 = datelist[lastone]
                    d2 = datelist[lastone+1]
                    d12diff = (d2-d1).days+(d1-d2).seconds/86400
                    dend.append( d1+timedelta(days=(d12diff)/2))
                    dendstd.append(d12diff*0.2886751345) #This is 1/sqrt(12)


    if len(dstart)>len(dend):
        dend.append(datetime(datelist[-1].year+1,1,1))
        dendstd.append(np.nan)

    mslen = 0
    msstd = 0
    for i in range(len(dstart)):
        imslen = timediff(dstart[i],dend[i],form='days')
        imsstd = np.sqrt(dstartstd[i]**2+dendstd[i]**2)
        if imslen>mslen:
            mslen=imslen
            msstd=imsstd
            
            
    return mslen,msstd

def stacksample(onsetstack,gmask,n=1e4):
    """sample onset timeseries, for the PCA decomposion.
    (this output can be fed into np.cov() and then np.eig()"""
    nanmask = ~np.isnan(np.sum(onsetstack,axis=0))
    ipts,jpts = np.where((gmask==1) & (nanmask))
    ksamples = np.array(random.sample(range(len(ipts)),int(n)))
    onsetsamp = onsetstack[:,ipts[ksamples],jpts[ksamples]]
    return onsetsamp,ipts[ksamples],jpts[ksamples]


# from osgeo import gdal
def tifxy(tiffile):
    """this returns and array of the x and y values for each point of a geotiff"""
    ds = gdal.Open(tiffile)
    c, a, b, f, d, e = ds.GetGeoTransform()
    arr = ds.GetRasterBand(1).ReadAsArray()
    imgshp = np.shape(arr)
    
    X = np.zeros(imgshp)
    Y = np.zeros(imgshp)
    
    
    for row in range(imgshp[0]):
        for col in range(imgshp[1]):
            # unravel GDAL affine transform parameters
            X[row,col] = a * col + b * row + a * 0.5 + b * 0.5 + c
            Y[row,col] = d * col + e * row + d * 0.5 + e * 0.5 + f
    
    return X,Y

