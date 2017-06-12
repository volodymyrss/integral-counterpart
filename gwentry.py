import plot

import json

import numpy as np

import astropy.io.fits as pyfits
from astropy.coordinates import SkyCoord
from astropy import units as u

import healpy
import healtics

from integralclient import converttime
import ephem

import os
import sys
from counterpart import Counterpart

from integralvisibility import Visibility

rootd="/Integral3/throng/savchenk/projects/ligo/predictions/o2/"


class GWEntry:
    def __init__(self,l):
        self.read_em_str(l)
        if not self: return

        self.dir=rootd+"/integral/sim_20160929/"+self.utc

        try:
            os.makedirs(self.dir)
        except OSError:
            pass
        
    def read_em_str(self,entry_str):
        if entry_str[0]=="#":
            return
        
        
        fields="GPS UTC mass1 mass2 dist SNR RA dec inclination skymap".split()        
        for f,v in zip(fields,entry_str.split(";")):
           # print f,v
            setattr(self,f,v.replace("\"",""))

        self.ra=float(self.RA)/pi*180
        self.dec=float(self.dec)/pi*180
            
        self.utc=self.UTC.replace(" ","T")
        self.ijd=float(converttime("UTC",self.utc,"IJD"))
        
        self.skymap=self.skymap[:-1]
            
    def __nonzero__(self):
        try:
            a=str(self.GPS)
            return True
        except:
            return False
            
    def __repr__(self):
        try:
            return str(self.GPStime)
        except:
            return "undef"
        
    def compute_visibility(self):
 #       coord=SkyCoord(ra.flatten(),dec.flatten(),unit=u.deg)
 #       m=Visibility().for_time("2016-10-09T17:31:06",coord=coord)

  #      figure(figsize=(20,10))
  #      scatter(ra,dec,s=100,c=m,lw=0,alpha=0.4)

        visibility=Visibility()

        gwm=healpy.read_map(rootd+"/"+self.skymap)
        nsides=healpy.npix2nside(gwm.shape[0])
        vmap=visibility.for_time(self.utc,nsides=nsides)

        healpy.mollview(gwm*(vmap*2-1),title="visibility for INTEGRAL due to sun constrains\n"+self.utc)

        healpy.projscatter(self.ra,self.dec,lonlat=True)
        healpy.graticule()
        plot.plot(self.dir+"/visibility.png")

        healpy.write_map(self.dir+"/visibility.fits",vmap)

        source_theta=(90-self.dec)/180*pi
        source_phi=(self.ra)/180*pi

        visibility=dict(
                    probability_visible=sum(gwm*vmap),
                    source_visible=vmap[healpy.ang2pix(nsides,source_theta,source_phi)]
                )

        json.dump(visibility,open(self.dir+"/visibility.json","w"))
        
    def compute_extend(self):
        # how many points
        pass

        
    def compute_counterpart(self):
        self.mutc=self.utc.replace("2016","2015")
        self.mutc=self.mutc.replace("2017","2016")

        os.chdir(self.dir)
        Counterpart(use_target_map_fn=rootd+"/"+self.skymap,use_utc=self.mutc).get()

def iterate_entries():
    gw_table=[]

    for l in open(rootd+"/sim_20160929/summary_sim_20160929.txt"):
        if len(gw_table)>20: continue
               
        gw=GWEntry(l)
        if not gw: continue

        gw.compute_visibility()
        gw.compute_counterpart()
            
        gw_table.append(gw)
        
        
        

