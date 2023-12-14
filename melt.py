#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 14:10:30 2023

@author: acjohnson16
"""

import numpy as np
from pathlib import Path
from datetime import datetime
import mu
import pandas as pd
from osgeo import gdal,osr


class glacierscene:
    def __init__(self,name,year,data_dir='default',pol='VV',createmeta='auto',
                 WTKstr=None):
        """loads in the meta data. sets the folders and stuff
        if creating the metadata, then WTKstr must also be provided"""
        
        self.name = name
        self.year = year
        self.pol = pol
        
        if data_dir=='default':
            data_dir = Path.cwd() / 'data' / self.name
        
        self.cwd = Path.cwd()
        
        self.data_dir = data_dir
        self.year_dir = self.data_dir / str(year)
        self.rawdata_dir = self.year_dir / 'rtc_clipped' / self.pol
        
        if createmeta=='auto':
            createmeta=True
            if (self.year_dir / f'imgmeta_{self.pol}.pkl').is_file():
                createmeta=False
        if createmeta==True:
            print('creating metadata')
            if WTKstr==None:
                raise ValueError('WTKstr is empty!!')
            results = mu.getasfresults(WTKstr, self.year,showstats=True)
            print(len(results))
            self.createmeta(results)       
        self.meta = pd.read_pickle(self.year_dir / f'imgmeta_{self.pol}.pkl')
        self.mpfu = np.unique(self.meta['mpfulist'])
        self.imgtot = len(self.meta)
        dst = gdal.Open(self.meta['pnamelist'].iloc[0])
        self.imgshp = np.shape(dst.ReadAsArray())
        del dst
        
        self.startdate = datetime(self.year,1,1)
        doy = []
        for d in self.meta['mdatelist'][:]:
            doy.append(mu.timediff(self.startdate,d,form='days'))            
        self.doy = np.array(doy)
        


    def createmeta(self,results):
        # for usepol in pols:
        
        mnamelist = []
        mdatelist = []
        mpflist = []
        mpfulist = []
    
        files = (self.rawdata_dir).glob('*.tif*')
        filelist = []
        pdatelist = []
        pnamelist = []
   
        for f in files:
            pdatelist.append(mu.getSARdate(f))
            filelist.append(f)
    
        pdatelist=np.array(pdatelist)
        filelist = np.array(filelist)
        
        for res in results:
            res = res.properties
            if self.pol in res['polarization']:
                rdate = datetime.strptime(res['startTime'],'%Y-%m-%dT%H:%M:%S.%fZ')
                pind = np.where(rdate==pdatelist)

                if len(pind[0])>=1:
                    pind = pind[0][0]
                    pnamelist.append(str(filelist[pind]))
                    rpf = str(res['pathNumber'])+'-'+str(res['frameNumber'])
                    mpflist.append(rpf)
                    mnamelist.append(res['sceneName'])
                    mdatelist.append(rdate)
                else:
                    print('WARNING: Missing product')
                    print(res['sceneName'])
                    
        pfunique = np.unique(mpflist)
        mpfulist = [np.where(pfunique==i)[0][0] for i in mpflist]

        pnamelist = np.array(pnamelist)
        mnamelist = np.array(mnamelist)
        mdatelist = np.array(mdatelist)
        mpflist = np.array(mpflist)
        mpfulist = np.array(mpfulist)
    
        mdatelist,pnamelist,mnamelist,mpflist,mpfulist = zip(*sorted(zip(mdatelist,pnamelist,mnamelist,mpflist,mpfulist)))

        meta = pd.DataFrame()
        meta['pnamelist'] = np.array(pnamelist)
        meta['mnamelist'] = np.array(mnamelist)
        meta['mdatelist'] = np.array(mdatelist)
        meta['mpflist']   = np.array(mpflist)
        meta['mpfulist']  = np.array(mpfulist)
        
        meta.to_pickle(self.year_dir/f'imgmeta_{self.pol}.pkl')        


    def createmeltth(self,meltth_val=-3.0):
        print('creating melt threshold')
        #step 1, add a list of 'is it winter' to the list
        #step 2, go through mpfu and get the mean value of each winter scene
        #step 3, save that as useful info somewhere
        
        #this is hardcoded now but could be variables:
        winterstart = datetime(self.year,1,1)
        winterend = datetime(self.year,4,1)
        meltth = np.zeros((len(self.mpfu),self.imgshp[0],self.imgshp[1]))
        self.meltth_val=meltth_val
        
        for u in self.mpfu:
            wmeta = self.meta.loc[self.meta['mpfulist']==u]
            # print(len(wmeta))
            wmeta = wmeta.loc[(wmeta['mdatelist']>=winterstart)&(wmeta['mdatelist']<winterend)]

            wimgs = np.array(wmeta['pnamelist'])
            wlen = len(wimgs)

            wintermap = np.zeros(self.imgshp)
            wintermap = meltth[u]

            for f in wimgs:
                img = gdal.Open(f).ReadAsArray()
                img[img==-60]=np.nan
                wintermap += img
                del img
            
            # print(wlen)
            if wlen>0:
                wintermap = 1/wlen * wintermap
            else:
                print('Warning: something very wrong has occured. There are not winter scenes for one mpfu')
                wintermap = wintermap-27
            wintermap = wintermap + meltth_val
            meltth[u] = wintermap
        np.save(self.year_dir/f'meltth_{self.pol}.npy',meltth)        

    def avgmonthmap(self,month,meltth_val=-3.0): 
        '''similar to above but for an arbitrary time period.
        month should be int
        
        something not properly accounted for: what if there is a nan within an
        mpfu, rather than without. No, I can fix that'''
        startdate = datetime(self.year,month,1)
        if month==12:
            stopdate = datetime(self.year+1,1,1)
        else:
            stopdate = datetime(self.year,month+1,1)
        self.load_meltth()
        avgmap = np.zeros(self.imgshp)

        weights = np.zeros((len(self.mpfu),self.imgshp[0],self.imgshp[1]))
        
        
        #iterate once to generate weights
        for u in self.mpfu:
            mmeta = self.meta.loc[self.meta['mpfulist']==u]
            mmeta = mmeta.loc[(mmeta['mdatelist']>=startdate)&(mmeta['mdatelist']<stopdate)]
            mimgs = np.array(mmeta['pnamelist'])
            weightmap = np.zeros_like(avgmap)
            
            for f in mimgs:
                img = gdal.Open(f).ReadAsArray()
                img[img==-60]=np.nan
                weightmap += ~np.isnan(img)
                del img

            weights[u] = weightmap

        weights = weights/np.sum(weights,axis=0)
        
        #iter twice to get data
        for u in self.mpfu:
            mmeta = self.meta.loc[self.meta['mpfulist']==u]
            mmeta = mmeta.loc[(mmeta['mdatelist']>=startdate)&(mmeta['mdatelist']<stopdate)]

            mimgs = np.array(mmeta['pnamelist'])

            mmap = np.zeros(self.imgshp)
            
            #iter twice to get data
            for f in mimgs:
                img = gdal.Open(f).ReadAsArray()
                img[img==-60]=np.nan
                mmask = ~np.isnan(img)
                mmap[mmask] += img[mmask]
                del img

            iterupdate = weights[u] * (mmap-self.meltth[u]-meltth_val)
            mmask = ~np.isnan(iterupdate)
            avgmap[mmask] += iterupdate[mmask]
               
        self.load_meltth(unload=True)
        return avgmap

    def load_meltth(self,unload=False):
        if unload == True:
            del self.meltth
        else:
            self.meltth = np.load(self.year_dir/f'meltth_{self.pol}.npy')
        
        
    def melt_imgnum(self,imgnum,loadmeltth=True):
        u = self.meta['mpfulist'].loc[imgnum]
        meltmap = np.ma.masked_array(np.zeros(self.imgshp,dtype=bool))
        img = gdal.Open(self.meta['pnamelist'].loc[imgnum]).ReadAsArray()

        if loadmeltth:
           self.meltth = np.load(self.year_dir/f'meltth_{self.pol}.npy')

        meltmap[img<=self.meltth[u]]=True
        meltmap.mask = img<=-60.0
        
        if loadmeltth:
            del self.meltth

        del img
        return meltmap
    
    def meltonset(self,adduncertainty=False):
        meltonsetmap = np.zeros(self.imgshp)
        self.load_meltth()
        doydiff = (self.doy[1:]+self.doy[:-1])*0.5
        doydiff = np.concatenate(([self.doy[0]],doydiff))
        
        if adduncertainty:
            doystd = np.array((self.doy[1:]-self.doy[:-1]))*0.2886751345 #This is 1/sqrt(12)
            doystd = np.concatenate(([0],doystd))
            doyunc = np.array([[doydiff[i]-x,doydiff[i]+x] for i,x in enumerate(doystd)]) #doy uncertainty    
    
        for i in range(len(self.meta)):
            cmelt = self.melt_imgnum(i,loadmeltth=False)
            upmask = (cmelt==True) & (meltonsetmap==0)
            meltonsetmap[upmask] = doydiff[i]
            if adduncertainty:
                crng = doyunc[i][1]-doyunc[i][0]
                uncvec = crng * np.random.random_sample(np.sum(upmask))-crng/2#+doyunc[i][0]
                meltonsetmap[upmask] = meltonsetmap[upmask] + uncvec
                # print(uncvec)
            
        self.load_meltth(unload=True)
    
        meltonsetmap[meltonsetmap==0]=np.nan
        self.meltonsetmap = meltonsetmap
    
    
    def loadmask(self,outlinefile=None):
        if outlinefile==None:
            outlinefile = self.data_dir / 'RGIraster.tiff'
        self.gmask = gdal.Open(str(outlinefile)).ReadAsArray()
        self.gmask = np.array(self.gmask,dtype=bool)


    def loaddem(self,demfile=None,usegmask=True,makebins=True,window=200):
        if demfile==None:
            demfile = self.data_dir / 'dem_align.npy'
        if str(demfile)[-4:]=='.tif':
            self.dem = gdal.Open(str(demfile)).ReadAsArray()
        elif str(demfile)[-4:]=='.npy':
            self.dem = np.load(demfile) # this is bad and it should feel bad
        if np.shape(self.dem)[0]!=self.imgshp[0]:
            print('adding line to dem')
            self.dem = np.concatenate((self.dem,np.array([self.dem[:,-1]]).transpose()),axis=1)
            self.dem = np.concatenate((self.dem,np.array([self.dem[-1,:]])),axis=0)
        self.dem = np.ma.masked_array(self.dem)
        
        
        if usegmask:
            self.loadmask()
            self.dem.mask=True
            self.dem.mask = ~self.gmask

        #ok now do the binning
        #too much binning
        
        if makebins:
            # window = 200 #m
            bin_start = np.floor(np.nanmin(self.dem)/window)*window
            bin_stop = np.ceil(np.nanmax(self.dem)/window)*window
            binn = int((bin_stop-bin_start)/window)
            
            # print(bin_start)
            # print(bin_stop)
            # print(binn)
            
            self.elebins = [[window*i+bin_start,window*(i+1)+bin_start] for i in range(binn)]
            self.binmids = np.array([window*i+bin_start+window/2 for i in range(binn)])
            
    def loadaspect(self,aspectfile=None):
        '''this is made in QGIS (I think)'''
        if aspectfile==None:
            aspectfile = self.data_dir / 'dem_aspect.tif'

        if str(aspectfile)[-4:]=='.tif':
            self.aspect = gdal.Open(str(aspectfile)).ReadAsArray()
        elif str(aspectfile)[-4:]=='.npy':
            self.aspect = np.load(aspectfile) # this is bad and it should feel bad
        print(self.imgshp)
        print(np.shape(self.aspect))
        self.aspect = np.concatenate((self.aspect,np.array([self.aspect[:,-1]]).transpose()),axis=1)
        self.aspect = np.concatenate((self.aspect,np.array([self.aspect[-1,:]])),axis=0)
        self.aspect = np.ma.masked_array(self.aspect)
        print(np.shape(self.aspect))

    def timeseries1pt(self,ipt,jpt,makeplot=False,saveplot=False,meltth_val=-3.0,correctvals=True):
        #iterate through each image. This is slow, but you're only doing plots for one

        self.load_meltth()
        self.datelist = self.meta['mdatelist']
        ptvec = []
        
        for i,img in self.meta.iterrows():
            cdata = gdal.Open(img['pnamelist']).ReadAsArray()
            adj_val = 0
            if correctvals:
                adj_val = self.meltth[img['mpfulist'],ipt,jpt]-meltth_val
            ptvec.append(cdata[ipt,jpt]-adj_val)
            del cdata
        self.load_meltth(unload=True)
        
        
        ptvec = np.array(ptvec)
        if makeplot:
            plt.figure()
            plt.plot(self.datelist,ptvec,'.')
            plt.axhline(0,linestyle='--',color='k')#,color='C01')
            plt.axhline(-3,linestyle='--',color='C01')
        

        return ptvec        
        
    def decomp(self,imgnum,decomp='hyp3_rgb'):
        img = self.meta.iloc[imgnum]
        
        
        if decomp=='hyp3_rgb':
            return mu.hyp3_rgb(imgcp,imgxp)
        
        
    def rgbsnow(self,rgb):
        arr = rgb[:,:,0]
        gee = rgb[:,:,1]
        bee = rgb[:,:,2]
        
        geemap = gee > arr+bee
        beemap = bee > arr+gee
        
        snowmap = geemap+2*beemap
        return snowmap
    
    def icefallmask(self,maxslope=0.5,dem_spacing=30,saveslope=False):
        '''takes numpy gradient of the dem'''
        self.loaddem(makebins=False)
        gradi,gradj = np.gradient(self.dem)
        gradi = gradi/dem_spacing
        gradj = gradj/dem_spacing
        slope = np.sqrt(gradi**2+gradj**2)
        slope[slope>50]=0
        
        self.ifmask = np.ones(self.imgshp,dtype='bool')
        self.ifmask[slope>maxslope]=False
        if saveslope==True:
            self.slope=slope

    def firnaq(self,mon1=3,mon2=10,firnaqth = -12):
        
        map1 = self.avgmonthmap(mon1)
        map2 = self.avgmonthmap(mon2)
        
        self.firnaqmap = np.array(map2-map1<firnaqth)
        
        
        
# def mapbin(dem,ibin):
#     return np.logical_and(dem>=ibin[0],dem<ibin[1])

