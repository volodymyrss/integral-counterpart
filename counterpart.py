from __future__ import print_function

import sys
import os
import string
import sys
import gzip
import time

from scipy.stats import norm

import healpy

import plot
import healtics
import astropy.io.fits as pyfits
import json

from astropy.coordinates import SkyCoord, PhysicsSphericalRepresentation
from astropy import units as u

import integralvisibility

import healpy

import integralclient as ic

from dataanalysis import core as da
from dataanalysis import hashtools
from dataanalysis.hashtools import shhash
import  dataanalysis

import datetime

import numpy as np
import integralclient

import healpy
import healtics
from matplotlib import pylab as p


def transform_rmap(rmap):
    nside = healpy.npix2nside(rmap.shape[0])
    npx = np.arange(healpy.nside2npix(nside))
    theta, phi = healpy.pix2ang(nside, npx)
    # SkyCoord(phi, theta, 1, unit=(u.rad, u.rad), representation="physicsspherical")
    ntheta = np.pi - theta
    nphi = np.pi + phi  # + or -???
    nphi[nphi > np.pi * 2] -= np.pi * 2
    return healpy.get_interp_val(rmap, ntheta, nphi)


def healpix_fk5_to_galactic(mp):
    nside = healpy.npix2nside(mp.shape[0])
    theta, phi = healpy.pix2ang(nside, np.arange(healpy.nside2npix(nside)))
    ra = phi / np.pi * 180
    dec = 90 - theta / np.pi * 180
    coord_galactic_gmap = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="galactic")
    coord_galactic_map = coord_galactic_gmap.transform_to("fk5")
    return healpy.get_interp_val(mp, coord_galactic_map.ra.degree, coord_galactic_map.dec.degree, lonlat=True)

def healpix_galactic_to_fk5(mp):
    nside = healpy.npix2nside(mp.shape[0])
    theta, phi = healpy.pix2ang(nside, np.arange(healpy.nside2npix(nside)))
    ra = phi / np.pi * 180
    dec = 90 - theta / np.pi * 180
    coord_galactic_gmap = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="fk5")
    coord_galactic_map = coord_galactic_gmap.transform_to("galactic")
    return healpy.get_interp_val(mp, coord_galactic_map.l.degree, coord_galactic_map.b.degree, lonlat=True)


def grid_for_healpix_map(mp):
    return healpy.pix2ang(healpy.npix2nside(mp.shape[0]),np.arange(mp.shape[0]),lonlat=True)

class CacheCounterparts(dataanalysis.caches.cache_core.CacheNoIndex):
    def get_burst(self, hashe):
        if isinstance(hashe, tuple):
            if hashe[0] == "analysis":
                if hashe[2].startswith('IceCubeEvent') or hashe[2].startswith('LIGOEvent'):
                    return hashe[2]
                return self.get_burst(hashe[1])
            if hashe[0] == "list":  # more universaly
                for k in hashe[1:]:
                    r = self.get_burst(k)
                    if r is not None:
                        return r
                return None
            raise Exception("unknown tuple in the hash:" + str(hashe))
        if hashe is None:
            return None  # 'Any'
        if isinstance(hashe, str):
            return None
        raise Exception("unknown class in the hash:" + str(hashe))

    def get_rev(self, hashe):
        # if dataanalysis.printhook.global_log_enabled: print("search for rev in",hashe)
        if isinstance(hashe, tuple):
            if hashe[
                0] == "analysis":  # more universaly
                if hashe[2].startswith('Revolution'):
                    return hashe[1]
                return self.get_rev(hashe[1])
            if hashe[
                0] == "list":  # more universaly
                for k in hashe[1:]:
                    r = self.get_rev(k)
                    if r is not None:
                        return r
                return None
            raise Exception("unknown tuple in the hash:" + str(hashe))
        if hashe is None:
            return None
        if isinstance(hashe, str):
            return None
        raise Exception("unknown class in the hash:" + str(hashe))

    def hashe2signature(self, hashe):
        scw = self.get_burst(hashe)
        if scw is not None:
            if isinstance(hashe, tuple):
                if hashe[0] == "analysis":
                    return hashe[2] + ":" + scw + ":" + shhash(hashe)[:8]
            return shhash(hashe)[:8]

        return hashe[2] + ":" + shhash(hashe)[:16]

    def construct_cached_file_path(self, hashe, obj=None):
        # print("will construct INTEGRAL cached file path for",hashe)

        burst = self.get_burst(hashe)

        def hash_to_path2(hashe):
            return shhash(repr(hashe[1]))[:8]

        if not isinstance(burst, str):
            burst = None

        if burst == "Any":
            burst = None

        rev = self.get_rev(hashe)
        print("Rev in hashe:",rev)

        if burst is None or len(burst.split(".",2))!=3:
            if dataanalysis.printhook.global_log_enabled: print("not burst-grouped cache")

            if rev is None:
                r = self.filecacheroot + "/global/" + hashe[2] + "/"  + hash_to_path2(hashe) + "/"
            else:
                if dataanalysis.printhook.global_log_enabled: print("is rev-grouped cache")
                hashe = hashtools.hashe_replace_object(hashe, rev, "any")
                r = self.filecacheroot + "/byrev/" + rev + "/" + hashe[2] + "/" + hash_to_path2(hashe) + "/"  # choose to avoid overlapp
        else:
            #utc=burst.split(".")[2]
            event_name = burst.split(".",2)[2]
            #utc_date=datetime.datetime.strptime(utc, "%Y-%m-%dT%H:%M:%S.%f")

            hashe = hashtools.hashe_replace_object(hashe, burst, "any")

            # str("reduced hashe:",hashe,hash_to_path2(hashe))
            # if dataanalysis.printhook.global_log_enabled: print("reduced hashe:",hashe,hash_to_path2(hashe))
            #open("reduced_hashe.txt", "w").write(hash_to_path2(hashe) + "\n\n" + pprint.pformat(hashe) + "\n")
            #print(burst, hashe[2], burst)

            r = self.filecacheroot + "/byevent/"+event_name
            #r += "/%.4i" % utc_date.year
            #r += "/%.4i-%.2i" % (utc_date.year,utc_date.month)
            #r += "/%.4i-%.2i-%.2i" % (utc_date.year,utc_date.month,utc_date.day)
            #r += "/%s" % utc.replace(":","-")

            r += "/" + hashe[2] + "/" + hash_to_path2(hashe) + "/"

        # if dataanalysis.printhook.global_log_enabled: print("cached path:",r)

        print(self, "cached path:", r)

        return r  # choose to avoid overlapp

IntegralCacheRoots=os.environ['INTEGRAL_DDCACHE_ROOT']

CacheStack=[]

for IntegralCacheRoot in IntegralCacheRoots.split(":"):
    ro_flag=False
    if IntegralCacheRoot.startswith("ro="):
        ro_flag=True
        IntegralCacheRoot=IntegralCacheRoot.replace("ro=","")

    mcgfb=CacheCounterparts(IntegralCacheRoot)
    mcgfb.readonly_cache=ro_flag
    if CacheStack==[]:
        CacheStack=[mcgfb]
    else:
        CacheStack[-1].parent=mcgfb
        CacheStack.append(mcgfb)

cache_counterpart=CacheStack[-1]

class DataAnalysis(da.DataAnalysis):
    store_preview_yaml=True

    cache=cache_counterpart

    def get_burst(self):
        if self._da_locally_complete is not None:
            try:
                return "(completeburst:%s)"%self.cache.get_scw(self._da_locally_complete)
            except:
                return "(complete)"

        for a in self.assumptions:
            if isinstance(a,LIGOEvent) or isinstance(a,IceCubeEvent):
                return "(burst:%s)"%str(a.t0_utc)

        return ""

class LIGOEvent(DataAnalysis):

    @property
    def loc_map(self):
        if not hasattr(self,'_loc_map'):
            self._loc_map=healpy.read_map(self.loc_map_path)
        return self._loc_map

    @property
    def nside(self):
        return healpy.npix2nside(self.loc_map.shape[0])

    @property
    def loc_region(self):
        indices = np.argsort(-self.loc_map)
        target_cum_map = np.empty(self.loc_map.shape)
        target_cum_map[indices] = 100 * np.cumsum(self.loc_map[indices])

        return target_cum_map

    def plot_map(self):
        healpy.mollview(self.loc_region,title="LIGO/Virgo localization of "+self.gname)
        healpy.graticule()
        for tr in np.linspace(0, 360 - 30, (360 / 30)):
            healpy.projtext(tr - 3, 3, "%i" % tr, lonlat=True)


class IceCubeEvent(DataAnalysis):
    nside=128

    ic_name="unnamed"

    @property
    def gname(self):
        return self.ic_name

    @property
    def trigger_time(self):
        return self.event_time['utc']

    def get_version(self):
        return self.get_signature()+"."+self.version+"."+self.ic_name #+"."+self.trigger_time

    def get_ra_dec(self):
        npx = np.arange(healpy.nside2npix(self.nside))
        return healpy.pix2ang(self.nside, npx, lonlat=True)

    @property
    def loc_map(self):
        if not hasattr(self,'_loc_map'):
            #self._loc_map=np.zeros(healpy.nside2npix(self.nside))

            ra,dec = self.get_ra_dec()

            p50_to_1sigma=1/norm.isf((1 - .50) / 2)/2


            dist_ra=(ra - self.loc_parameters['ra'])
            dist_ra[abs(dist_ra)>180]-=360

            self._loc_map=np.exp(
                -0.5*(dist_ra/(p50_to_1sigma*np.mean(map(abs,self.loc_parameters['dra_p50']))))**2 \
                -0.5*((dec - self.loc_parameters['dec']) / (p50_to_1sigma*np.mean(map(abs, self.loc_parameters['ddec_p50']))))**2
            )

            self._loc_map/=sum(self._loc_map)

        return self._loc_map

    @property
    def loc_region(self):
        indices = np.argsort(-self.loc_map)
        target_cum_map = np.empty(self.loc_map.shape)
        target_cum_map[indices] = 100 * np.cumsum(self.loc_map[indices])

        return target_cum_map

    def plot_map(self):
        healpy.mollview(self.loc_region,title="IceCube localization of "+self.gname)
        healpy.graticule()
        for tr in np.linspace(0, 360 - 30, (360 / 30)):
            healpy.projtext(tr - 3, 3, "%i" % tr, lonlat=True)

    def describe_loc_regions(self):
        ra, dec= self.get_ra_dec()

        i_peak=self.loc_region.argmin()

        print("peak",ra[i_peak],dec[i_peak])

        for limit in 50, 90:
            m = self.loc_region < limit
            print(limit,"% containment")
            print("  RA",round(ra[i_peak]-ra[m].min(),2),'..',round(ra[m].max()-ra[i_peak],2),)
            print("  RA", round(ra[m].min(), 2), '..', round(ra[m].max(), 2), )
            print("  Dec", round(dec[i_peak] - dec[m].min(), 2), '..', round(dec[m].max() - dec[i_peak], 2), )
            print("  Dec", round(dec[m].min(), 2), '..', round(dec[m].max(), 2), )


class Event(da.DataAnalysis):
    run_for_hashe=True

    event_kind=LIGOEvent

    def main(self):
        return self.event_kind

class INTEGRALVisibility(DataAnalysis):
    input_target=Event

    minsolarangle=40
    follow_up_delay=0

    _da_settings=['minsolarangle','follow_up_delay']

    def main(self):
        ijd0=float(ic.converttime("UTC", self.input_target.trigger_time, "IJD"))
        ijd=ijd0+self.follow_up_delay

        utc= ic.converttime("IJD", "%.20lg"%ijd, "UTC")

        self.visibility = integralvisibility.Visibility()
        self.visibility.minsolarangle = self.minsolarangle
        self.visibility_map = self.visibility.for_time(utc, nsides=self.input_target.nside)

    @property
    def total_visible(self):
        return sum(self.input_target.loc_map[self.visibility_map > 0])

    def peak_target_visible(self):
        ra,dec=grid_for_healpix_map(self.visibility_map)

        i=np.argmax(self.input_target.loc_map)
        return ra[i],dec[i]


    def plot(self):
        healpy.mollview(self.visibility_map)
        healpy.graticule()

class SourceAssumptions(DataAnalysis):

    version="v2"

    def main(self):
        data = {}
        data.update(dict(
            byinstrument={},
            instruments=["JEM-X", "ISGRI", "PICsIT", "SPI", "SPI-ACS", "IBIS/Veto"],
            zones=["fov", "z1", "z2", "z3"],
            typical_spectra=dict(
                hard=dict(model="compton", alpha=-0.5, epeak=600, beta=-9, ampl=1),
                soft=dict(model="band", alpha=-1, epeak=300, beta=-2.5, ampl=1),
                soft1s=dict(model="band", alpha=-1, epeak=300, beta=-2.5, ampl=1),
                crabby8s=dict(model="powerlaw", alpha=-2, epeak=300, beta=-2.5, ampl=1),
                crabby1s=dict(model="powerlaw", alpha=-2, epeak=300, beta=-2.5, ampl=1),
            )
        ))

        data['byinstrument']['JEM-X'] = dict(emin=3, emax=30, fov_area=130, angres="3'", zones=["fov", "onaxis"])
        data['byinstrument']['ISGRI'] = dict(emin=20, emax=200, fov_area=900, angres="12'", zones=["fov", "onaxis", "z1"])
        data['byinstrument']['PICsIT'] = dict(emin=260, emax=8000, fov_area=900, angres="30'",
                                              zones=["fov", "onaxis", "z1"])
        data['byinstrument']['SPI'] = dict(emin=25, emax=8000, fov_area=1225, angres="2.5$^\circ$", zones=["fov", "onaxis"])
        data['byinstrument']['SPI-ACS'] = dict(emin=75, emax=100000, fov_area=None, angres=None,
                                               zones=["fov", "onaxis", "z1", "z2", "z3"])
        data['byinstrument']['IBIS/Veto'] = dict(emin=100, emax=100000, fov_area=None, angres=None,
                                                 zones=["fov", "onaxis", "z1", "z2", "z3"])
        # data['instruments']['Compton']=dict(emin=260,emax=800)

        data['byinstrument']['JEM-X']['erange_ag'] = dict(emin=3, emax=30)
        data['byinstrument']['ISGRI']['erange_ag'] = dict(emin=20, emax=200)
        data['byinstrument']['PICsIT']['erange_ag'] = dict(emin=260, emax=2600)
        data['byinstrument']['SPI']['erange_ag'] = dict(emin=25, emax=1000)

        data['byinstrument']['SPI-ACS']['erange_grb'] = dict(emin=75, emax=2000)
        data['byinstrument']['IBIS/Veto']['erange_grb'] = dict(emin=100, emax=2000)
        data['byinstrument']['JEM-X']['erange_grb'] = dict(emin=3, emax=30)
        data['byinstrument']['ISGRI']['erange_grb'] = dict(emin=20, emax=200)
        data['byinstrument']['PICsIT']['erange_grb'] = dict(emin=260, emax=2600)
        data['byinstrument']['SPI']['erange_grb'] = dict(emin=25, emax=8000)

        self.data=data
        self.duration_by_kind = dict(hard=1, soft=8, soft1s=1, crabby1s=1, crabby8s=8)

class OperationStatus(DataAnalysis):
    cached=True

    @property
    def targets(self):
        targets = []
        if self.ibis_on:
            targets += ["ISGRI", "IBIS/Veto"]
        if self.spi_on:
            targets += ["SPI-ACS", ]
        return targets

    def main(self):
        ibis_on = True
        spi_on = True

class OperationsReport(DataAnalysis):
    input_opsstatus=OperationStatus

    def main(self):
        if self.input_opsstatus.ibis_on and self.input_opsstatus.spi_on:
            self.text="INTEGRAL was operating in nominal mode"
        elif self.input_opsstatus.spi_on and not self.input_opsstatus.ibis_on:
            self.text="INTEGRAL SPI-ACS was operating in nominal mode, while IBIS was not taking data"
        elif self.input_opsstatus.ibis_on and not self.input_opsstatus.spi_on:
            self.text="INTEGRAL IBIS was operating in nominal mode, while SPI-ACS was not taking data"
        else:
            self.text = "INTEGRAL was not operational"



class CountLimits(DataAnalysis):
    input_target=Event
    input_assumptions=SourceAssumptions
    input_operation_status=OperationStatus

    cached=True

    version="v2"

    @property
    def targets(self):
        targets=[]
        if self.input_operation_status.ibis_on:
            targets+=["ISGRI", "IBIS/Veto"]
        if self.input_operation_status.spi_on:
            targets+=["SPI-ACS",]
        return targets

    def get_count_limit(self,target, scale, nsig=None):

        n_attempts=10
        while n_attempts>=0:
            try:
                hk = ic.get_hk(
                                target=target,
                                utc=self.input_target.trigger_time,
                                span=600, t1=0, t2=0, ra=0, dec=0, rebin=scale,
                                vetofiltermargin=0.03
                            )['lc']
            except ic.ServiceException as e:
                print("waiting...",e)
                time.sleep(1)
                n_attempts-=1

            else:
                break
        if n_attempts<0:
            raise Exception("unable to reach the server in %i attempts!"%n_attempts)

        if not hasattr(self,'hk'):
            self.hk={}

        if target not in self.hk:
            self.hk[target]={}

        self.hk[target][scale]=hk

        print(target, ":", scale, hk)
        return hk['count limit 3 sigma']

    def main(self):

        hkname = {'ISGRI': 'ISGRI', "PICsIT": "SPTI234", "IBIS/Veto": "IBIS_VETO", "SPI-ACS": "ACS"}

        count_limits = {}
        for kind, model in self.input_assumptions.data['typical_spectra'].items():
            if kind in self.input_assumptions.duration_by_kind:
                count_limits[kind] = {}
                for target in self.targets:
                    count_limits[kind][target] = self.get_count_limit(hkname[target], scale=self.input_assumptions.duration_by_kind[kind])

        self.count_limits=count_limits

class BackgroundStabilityAssertion(DataAnalysis):
    input_countlimits=CountLimits

    def main(self):
        excvar=[]
        instabilities=[]
        for target,t in self.input_countlimits.hk.items():
            for scale,d in t.items():
                print(d)

                excvar.append(d['excvar'])
                if d['maxsig']>6.:
                    instabilities.append(d['excvar'])

        if max(excvar)<1.3:
            self.text="very stable"
        elif max(excvar) < 1.7:
            self.text = "stable"
        elif max(excvar) < 3:
            self.text = "somewhat unstable"
        else:
            self.text = "very unstable"

        if len(instabilities)>0:
            self.text+=" with some isolated features"




class ScSystem(DataAnalysis):
    input_target=Event

    version="v2"

    def main(self):
        self.sc = ic.get_sc(self.input_target.trigger_time,
                            ra=self.input_target.loc_parameters['ra'],
                            dec=self.input_target.loc_parameters['dec'])

class Responses(DataAnalysis):
    input_assumptions=SourceAssumptions
    input_target=Event

    cached = True

    def main(self):
        CT = Counterpart()
        CT.nside = self.input_target.nside
        CT.utc = self.input_target.trigger_time
        CT.sc = ic.get_sc(CT.utc)
        CT.ra = CT.sc['scx']['ra']
        CT.dec = CT.sc['scx']['dec']
        CT.compute_transform_grids()

        responses = {}

        rname = {'ISGRI': 'ISGRI', "PICsIT": "PICsIT", "IBIS/Veto": "VETO", "SPI-ACS": "ACS"}
        lt_byt = dict(ISGRI=30, PICsIT=260, ACS="map2", VETO=100)

        for kind, model in self.input_assumptions.data['typical_spectra'].items():
            responses[kind] = {}
            print(model)

            for target in ["SPI-ACS", "ISGRI", "IBIS/Veto"]:

                n_attempts=10
                while True:
                    try:
                        mp = transform_rmap(np.array(ic.get_response_map(target=rname[target],
                                                                                  lt=lt_byt[rname[target]],
                                                                                  model=model['model'],
                                                                                  epeak=model['epeak'],
                                                                                  alpha=model['alpha'],
                                                                                  emin=40,
                                                                                  emax=40000,
                                                                                  )))
                    except Exception as e:
                        print("waiting",e)
                        time.sleep(1)
                        n_attempts-=1
                    else:
                        break

                skymp = CT.sc_map_in_sky(np.array(mp))
                responses[kind][target] = skymp

        self.responses=responses

    def plot(self):
        for kind,v in self.responses.items():
            for instrument,mp in v.items():
                healpy.mollview(mp,title=kind+" "+instrument)


class OrientationComment(DataAnalysis):
    input_responses = Responses
    input_event = Event
    input_opsstatus = OperationStatus

    version="v2"

    def main(self):
        ra, dec = self.input_event.get_ra_dec()

        assert len(self.input_responses.responses['hard']['SPI-ACS'])==len(ra)

        i_peak = self.input_event.loc_region.argmin()
        self.response_onpeak = dict([
            [k, self.input_responses.responses['hard'][k][i_peak]] for k in self.input_opsstatus.targets
        ])

        self.response_best = dict([
            [k, self.input_responses.responses['hard'][k].min()] for k in self.input_opsstatus.targets
        ])

        self.response_quality={}

        for k in self.response_onpeak:
            if self.response_onpeak[k]/self.response_best[k]<1.5:
                self.response_quality[k]="near-optimal"
            if self.response_onpeak[k]/self.response_best[k]<2.5:
                self.response_quality[k]="somewhat suppressed"
            else:
                self.response_quality[k]="strongly suppressed"


        print("peak", ra[i_peak], dec[i_peak], "response", self.response_onpeak)


        self.text="This orientation implies "

        if len(self.response_onpeak)==1:
            k=self.response_onpeak.keys()[0]
            self.text+=self.response_quality[k]+" response of "+k

        #near-optimal
        #    response of SPI-ACS, and this instrument provides best sensitivity to
        #    both short and long GRBs."""


class Sensitivities(DataAnalysis):
    input_assumptions=SourceAssumptions
    input_target=Event
    input_countlimits=CountLimits
    input_responses=Responses

    cached = True

    def main(self):
        sens_maps = {}
        sens_maps_gal = {}

        theta_lvt, phi_lvt=grid_for_healpix_map(self.input_target.loc_map)

        for kind, model in self.input_assumptions.data['typical_spectra'].items():
            if kind in self.input_assumptions.duration_by_kind:
                sens_maps[kind] = {}
                sens_maps_gal[kind] = {}
                for target in self.input_countlimits.targets:
                    r = self.input_responses.responses[kind][target] * self.input_countlimits.count_limits[kind][target]
                    #r[CT.sky_coord.separation(self.input_bodiesbody_coord).degree < bd["body_size"]] = 1e9

                    sens_maps[kind][target] = healpy.get_interp_val(r, theta_lvt, phi_lvt, lonlat=True)
                    sens_maps_gal[kind][target] = healpix_fk5_to_galactic(sens_maps[kind][target])

                sens_maps[kind]['best'] = np.array([m for k, m in sens_maps[kind].items()]).min(0)
                sens_maps_gal[kind]['best'] = np.array([m for k, m in sens_maps_gal[kind].items()]).min(0)

        self.sens_maps = sens_maps
        self.sens_maps_gal = sens_maps_gal

        self.summarize()

    def plot(self):
        for kind, v in self.sens_maps.items():
            for instrument, mp in v.items():
                healpy.mollview(mp, title=kind + " " + instrument)



    def plot_fancy(self):
        p.figure()

        rangemin=dict(hard=1.5,soft=5)
        rangemax = dict(hard=6, soft=15)

        for kind in ['hard','soft']:
            fig = healtics.plot_with_ticks(self.sens_maps[kind]['best'] / 1e-7,
                                           title="INTEGRAL 3-sigma upper limit on " + self.input_target.gname+ "; model "+kind,
                                       overplot=[[(self.input_target.loc_region, "k", 90), (self.input_target.loc_region, "k", 50),
                                                  #(mp, 'g', 0.5),
                                                  #   (fermi_shadow,'m',0.5),
                                                  # (region2,"r",90),(region2,"r",50),
                                                  # (sens_mcrab,"g",10),
                                                  #       (GC,"m",0.5),
                                                  ]],
                                       # overplot=[(region,"k",10),(healpy.smoothing(region2,0.2),"r",10)],
                                       # ilevels=array([50,10]),
                                       vmin=rangemin[kind], vmax=rangemax[kind], invra=False,
                                       unit="$10^{-7} \mathrm{erg^{}cm^{-2} s^{-1}}$")

    def summarize(self):
        lvt = self.input_target.loc_region
        lvt_prob = self.input_target.loc_map

        self.summary={}

        for kind in "hard","soft", "soft1s", "crabby1s", "crabby8s":
            mp=self.sens_maps[kind]['best']

            self.summary[kind]=dict(
                p90_min=mp[lvt < 90].min(),
                p90_max=mp[lvt < 90].max(),
                p50_min=mp[lvt < 50].min(),
                p50_max=mp[lvt < 50].max(),
                prob_within_1p5best=sum(lvt_prob[mp < 1.5 * mp[lvt < 90].min()]),
            )

            print("model", kind, "sensitivity in 90% containment extreme from",mp[lvt<90].min(),"to",mp[lvt<90].max())
            print("model", kind, "sensitivity in 50% containment extreme from", mp[lvt < 50].min(), "to",mp[lvt < 50].max())
            print("model", kind, "sensitivity within 50% of best contains",sum(lvt_prob[mp<1.5*mp[lvt<90].min()]))




class Counterpart(DataAnalysis):
    input_target=Event

    syst = 0.2

    t1 = 0
    t2 = 0

    do_burst_analysis = False

    def main(self):
        print(self.utc)

        self.sc = integralclient.get_sc(self.utc, ra=0, dec=0)

        if self.target_map_fn == "":
            self.target_map = np.zeros(healpy.nside2npix(16))
        else:
            self.target_map = healpy.read_map(self.target_map_fn)

        indices = np.argsort(-self.target_map)
        self.target_cum_map = np.empty(self.target_map.shape)
        self.target_cum_map[indices] = 100 * np.cumsum(self.target_map[indices])

        self.nside = healpy.npix2nside(self.target_map.shape[0])
        self.compute_transform_grids()

        if self.target_map_fn == "":
            cat = integralclient.get_cat(self.utc)
            locerr = cat['locerr']
            ra = cat['ra']
            dec = cat['dec']
            if locerr < 50:
                if locerr < 0.5:
                    locerr = 0.5
                self.target_map = np.array(
                    np.exp(-(self.sky_coord.separation(SkyCoord(ra, dec, unit=(u.deg, u.deg))) / u.deg) / locerr ** 2 / 2),
                    dtype=float)
                print(self.target_map.shape)

        if False:
            scales = [8]
            alpha = -1
            epeak = 300
            beta = -2.5
            model = "band"
            self.isgri_cl = 1.2
        else:
            scales = [1]
            alpha = -0.5
            epeak = 600
            beta = -2.5
            model = "compton"
            self.isgri_cl = 1.5

        self.response_mp_acs = transform_rmap(np.array(
            integralclient.get_response_map(target="ACS", lt='map2', alpha=alpha, epeak=epeak, beta=beta, model=model,
                                            kind="response")))
        self.response_mp_veto = transform_rmap(np.array(
            integralclient.get_response_map(target="VETO", lt=100, alpha=alpha, epeak=epeak, beta=beta, model=model,
                                            kind="response")))
        self.response_mp_isgri = transform_rmap(np.array(
            integralclient.get_response_map(target="ISGRI", lt=30, alpha=alpha, epeak=epeak, beta=beta, model=model,
                                            kind="response")))
        self.response_mp_picsit = transform_rmap(np.array(
            integralclient.get_response_map(target="PICsIT", lt=250, alpha=alpha, epeak=epeak, beta=beta, model=model,
                                            kind="response")))

        print("best response ACS, Veto, ISGRI", self.response_mp_acs.min(), self.response_mp_veto.min(), self.response_mp_isgri.min())

        def get_count_limit(target, scale):

            try:
                r = integralclient.get_hk(target=target, utc=self.utc, span=30.01, t1=0, t2=0, ra=0, dec=0, rebin=scale)
                hk = r.json()['lc']
            except Exception as e:
                print(e, r.content)

            print(target, ":", scale, hk)
            return hk['std bkg'] * 3 * hk['timebin']

        if self.do_burst_analysis:
            def get_burst_counts(target):
                span = (self.t2 - self.t1) * 2. + 100
                hk = integralclient.get_hk(target=target, utc=self.utc, span=span, t1=self.t1, t2=self.t2, ra=0, dec=0,
                                           rebin=0)['lc']
                print(target, ":", hk)
                return hk['burst counts'], hk['burst counts error'], hk['burst region']

            self.acs_counts = get_burst_counts("ACS")
            self.veto_counts = get_burst_counts("VETO")
            self.isgri_counts = get_burst_counts("ISGRI")
            self.picsit_counts = get_burst_counts("SPTI1234")

        # scales=[1,]
        # scales=[1,8,0.1]

        for scale in scales:
            acs_lim = get_count_limit("ACS", scale)
            veto_lim = get_count_limit("VETO", scale)
            isgri_lim = get_count_limit("ISGRI", scale)
            picsit_lim = isgri_lim * (3000 / 600) ** 0.5

            # acs_lim=300
            # veto_lim=300
            # isgri_lim=300


            print("ACS, Veto, ISGRI", acs_lim, veto_lim, isgri_lim)

            self.syst = 0.

            sens_map = self.response_mp_acs * acs_lim * (self.syst + 1)
            sens_map_acs = self.response_mp_acs * acs_lim * (self.syst + 1)
            sens_map_veto = self.response_mp_veto * veto_lim * (self.syst + 1)
            sens_map_isgri = self.response_mp_isgri * isgri_lim * (self.syst + 1) * self.isgri_cl
            sens_map_picsit = self.response_mp_picsit * picsit_lim * (self.syst + 1)
            sens_map[sens_map > sens_map_veto] = sens_map_veto[sens_map > sens_map_veto]
            sens_map[sens_map > sens_map_isgri] = sens_map_isgri[sens_map > sens_map_isgri]
            sens_map[sens_map > sens_map_picsit] = sens_map_picsit[sens_map > sens_map_picsit]

            na = sens_map[~np.isnan(sens_map) & (sens_map > 0)].min()

            na_e = int(np.log10(na)) - 1
            na_b = int(na * 10 / 10 ** na_e) / 10.

            ### 

            # na_b=3
            # na_e=-7
            na_b = 1.3

            ###

            na = na_b * 10 ** na_e

            print("best ACS", na)
            nv = sens_map_veto[~np.isnan(sens_map_veto) & (sens_map_veto > 0)].min()
            print("best VETO", nv)
            self.sens_scale = na_b * 10 ** na_e
            sens_map /= self.sens_scale
            sens_map_veto /= self.sens_scale
            sens_map_acs /= self.sens_scale
            sens_map_isgri /= self.sens_scale
            sens_map_picsit /= self.sens_scale

            self.sens_scale_e = na_e
            self.sens_scale_b = na_b

            self.sens_map = sens_map
            self.sens_map_acs = sens_map_acs
            self.sens_map_veto = sens_map_veto
            self.sens_map_isgri = sens_map_isgri
            self.sens_map_picsit = sens_map_picsit

            self.tag = "_%.5lgs" % scale

            if self.do_burst_analysis:
                self.localize()
                self.localize2()
                self.localization()
            self.compute_maps()

    def localize(self):

        allmaps = []
        alldet = [[(c, ce, br), m, n] for (c, ce, br), m, n in [(self.acs_counts, self.response_mp_acs, "ACS"),
                                                                (self.veto_counts, self.response_mp_veto, "VETO"),
                                                                (self.isgri_counts, self.response_mp_isgri, "ISGRI"),
                                                                (self.picsit_counts, self.response_mp_picsit, "PICsIT")]
                  if c / ce > -3 and
                  ce > 0 and
                  br > 0 and
                  c != 0
                  ]

        print("will localize with", [n for (c, ce, br), m, n in alldet])

        total_region = []
        c_vec = np.array([c for (c, ce, br), m, n in alldet])
        ce_vec = np.array([ce for (c, ce, br), m, n in alldet])
        # ce_vec=(ce_vec**2+(c_vec*0.05)**2)**0.5

        response_mt = np.array([m for ((c, ce, br), m, n) in alldet])
        print(response_mt.shape)

        nc_mt = np.outer(c_vec, np.ones(response_mt.shape[1])) * response_mt
        nce_mt = np.outer(ce_vec, np.ones(response_mt.shape[1])) * response_mt

        mean_map = sum(nc_mt / nce_mt ** 2, axis=0) / sum(1. / nce_mt ** 2, axis=0)
        err_map = 1 / sum(1. / nce_mt ** 2, axis=0) ** 0.5
        chi2_map = sum((nc_mt - np.outer(np.ones_like(c_vec), mean_map)) ** 2 / nce_mt ** 2, axis=0)

        min_px = chi2_map.argmin()
        print("minimum prediction", response_mt[:, min_px], mean_map[min_px] / response_mt[:, min_px], chi2_map[min_px])
        print("measure", c_vec)
        print("measure err", ce_vec)
        print("sig", (c - mean_map[min_px] * response_mt[:, min_px]) / ce_vec)

        healpy.mollview(chi2_map)
        plot.plot("chi2_map.png")

        self.locmap = chi2_map / chi2_map.min()

        print(self.locmap.min(), self.locmap.max())

        healpy.mollview(mean_map)
        plot.plot("mean_map.png")

    def localize2(self):

        allmaps = []
        alldet = [(self.acs_counts, self.response_mp_acs, "ACS"),
                  (self.veto_counts, self.response_mp_veto, "VETO"),
                  (self.picsit_counts, self.response_mp_picsit, "PICsIT")]

        #         (self.isgri_counts,self.response_mp_isgri,"ISGRI"),

        total_region = np.ones_like(self.response_mp_acs, dtype=bool)
        for i1, ((c1, ce1, br1), m1, n1) in enumerate(alldet):
            for i2, ((c2, ce2, br2), m2, n2) in enumerate(alldet):
                if c1 == 0 or c2 == 0: continue  # ???
                if ce1 == 0 or ce2 == 0: continue
                if br1 == 0 or br2 == 0: continue
                if i2 >= i1: continue

                ang0 = np.arctan2(m2, m1)  # inversed for responose
                h0 = np.histogram(ang0.flatten(), 100)

                ang1 = np.arctan2((c1 - ce1) * (1 - self.syst), (c2 + ce2) * (1 + self.syst))
                ang2 = np.arctan2((c1 + ce1) * (1 + self.syst), (c2 - ce2) * (1 - self.syst))

                print(n1, ":", c1, ce1, "; ", n2, c2, ce2, " => ", ang1, ang2, " while ", ang0.min(), ang0.max())

                region = (ang0 > ang1) & (ang0 < ang2)

                healpy.mollview(region)
                plot.plot("region_%s_%s.png" % (n1, n2))

                total_region &= region

        healpy.mollview(total_region)
        plot.plot("total_region.png")

    def localization(self):
        syst = 0.02

        fluence = self.sens_map_acs.min() * 10

        print("fluence for 10 sigma in ACS:", fluence)

        sig_map_acs = fluence / self.sens_map_acs
        sig_map_veto = fluence / self.sens_map_veto
        sig_map_isgri = fluence / self.sens_map_isgri
        sig_map_picsit = fluence / self.sens_map_picsit

        if False:
            figure()
            scatter(sig_map_acs, sig_map_picsit, c=coord.theta)
            colorbar()

            figure()
            scatter(sig_map_acs, sig_map_veto, c=coord.theta)
            colorbar()

            figure()
            scatter(sig_map_acs, sig_map_isgri, c=coord.theta)
            colorbar()

            figure()
            scatter(sig_map_veto, sig_map_picsit, c=coord.theta)
            colorbar()

        allmaps = []
        allsig = [sig_map_acs,
                  sig_map_veto,
                  sig_map_isgri,
                  sig_map_picsit]

        for m1 in allsig:
            for m2 in allsig:
                # m1=m1/sig_map_acs*sig_map_acs.max()
                # m2=m1/sig_map_acs*sig_map_acs.max()

                ang0 = np.arctan2(m1, m2)
                h0 = np.histogram(ang0.flatten(), 100)
                ang1 = np.arctan2((m1 - 1) * (1 - syst), (m2 + 1) * (1 + syst))
                ang2 = np.arctan2((m1 + 1) * (1 + syst), (m2 - 1) * (1 - syst))

                areas = [sum((ang0[i] > ang1) & (ang0[i] < ang2)) for i in range(ang0.shape[0])]

                area = (np.array(areas) * healpy.nside2pixarea(16))
                area /= 4 * np.pi

                allmaps.append(area)

        bestarea = np.array(allmaps).min(0)
        healpy.mollview(bestarea)
        self.bestarea = bestarea
        plot.plot("bestlocalization.png")

        plot.p.figure()
        plot.p.a = plot.p.hist(bestarea, np.logspace(-5, -0.5, 100), log=True)
        plot.p.semilogx()
        plot.plot("bestlocalization_hist.png")

    def get_grid(self, nside=None):
        nside = nside if nside is not None else self.nside
        npx = np.arange(healpy.nside2npix(nside))
        theta, phi = healpy.pix2ang(nside, npx)
        return SkyCoord(phi, theta, 1, unit=(u.rad, u.rad), representation="physicsspherical")

    def compute_transform_grids(self):
        sky_coord = self.get_grid()
        self.sky_coord = sky_coord

        x = sky_coord.cartesian.x
        y = sky_coord.cartesian.y
        z = sky_coord.cartesian.z

        self.scX = SkyCoord(self.sc['scx']['ra'], self.sc['scx']['dec'], frame='icrs', unit='deg')
        self.scY = SkyCoord(self.sc['scy']['ra'], self.sc['scy']['dec'], frame='icrs', unit='deg')
        self.scZ = SkyCoord(self.sc['scz']['ra'], self.sc['scz']['dec'], frame='icrs', unit='deg')
        self.scz = self.scX
        self.scy = self.scY
        self.scx = self.scZ

        # maps of sky coordinates in detector coordinate grid
        x_sky_in_sc = sky_coord.cartesian.x * self.scx.cartesian.x + sky_coord.cartesian.y * self.scy.cartesian.x + sky_coord.cartesian.z * self.scz.cartesian.x
        y_sky_in_sc = sky_coord.cartesian.x * self.scx.cartesian.y + sky_coord.cartesian.y * self.scy.cartesian.y + sky_coord.cartesian.z * self.scz.cartesian.y
        z_sky_in_sc = sky_coord.cartesian.x * self.scx.cartesian.z + sky_coord.cartesian.y * self.scy.cartesian.z + sky_coord.cartesian.z * self.scz.cartesian.z

        sky_in_sc = SkyCoord(x=x_sky_in_sc, y=y_sky_in_sc, z=z_sky_in_sc, representation="cartesian")

        on_scz = self.scz.cartesian.x * sky_coord.cartesian.x + self.scz.cartesian.y * sky_coord.cartesian.y + self.scz.cartesian.z * sky_coord.cartesian.z
        x_inxy = x - on_scz * self.scz.cartesian.x
        y_inxy = y - on_scz * self.scz.cartesian.y
        z_inxy = z - on_scz * self.scz.cartesian.z

        r_inxy = (x_inxy * x_inxy + y_inxy * y_inxy + z_inxy * z_inxy) ** 0.5
        x_inxy /= r_inxy
        y_inxy /= r_inxy
        z_inxy /= r_inxy

        sc_in_sky = SkyCoord(
            np.arctan2(x_inxy * self.scy.cartesian.x + y_inxy * self.scy.cartesian.y + z_inxy * self.scy.cartesian.z,
                    x_inxy * self.scx.cartesian.x + y_inxy * self.scx.cartesian.y + z_inxy * self.scx.cartesian.z),
            np.arccos(on_scz),
            1,
            unit=(u.rad, u.rad),
            representation="physicsspherical")

        self.sky_in_sc = sky_in_sc
        self.sc_in_sky = sc_in_sky

        return sky_in_sc, sc_in_sky

    def sky_map_in_sc(self, sky_map):
        return healpy.get_interp_val(sky_map,
                                     theta=self.sky_in_sc.represent_as("physicsspherical").theta.rad,
                                     phi=self.sky_in_sc.represent_as("physicsspherical").phi.rad)

    def sc_map_in_sky(self, sc_map):
        return healpy.get_interp_val(sc_map,
                                     theta=self.sc_in_sky.represent_as("physicsspherical").theta.rad,
                                     phi=self.sc_in_sky.represent_as("physicsspherical").phi.rad)

    def plot_sky_diagram(self):
        target_map_in_sc = self.sky_map_in_sc(self.target_map)
        target_map_in_sc[
            abs(self.sky_coord.represent_as("physicsspherical").theta.deg) < 20] += target_map_in_sc.max() / 5.
        target_map_in_sc[
            abs(self.sky_coord.represent_as("physicsspherical").theta.deg) > 140] += target_map_in_sc.max() / 5.
        target_map_in_sc[(self.sky_coord.represent_as("physicsspherical").phi.deg < 30) | (
            self.sky_coord.represent_as("physicsspherical").phi.deg > 360 - 30)] += target_map_in_sc.max() / 10.
        healpy.mollview(target_map_in_sc, cmap='YlOrBr')
        healpy.projscatter(np.pi / 2, 0)
        healpy.graticule()
        plot.plot()

    def compute_maps(self):
        healpy.mollview(self.sens_map_acs, cmap="YlOrBr")
        healpy.graticule()

        plot.plot()

        sens_map_sky = healpy.sphtfunc.smoothing(self.sc_map_in_sky(self.sens_map), 5. / 180. * np.pi)
        sens_map_sky_acs = healpy.sphtfunc.smoothing(self.sc_map_in_sky(self.sens_map_acs), 5. / 180. * np.pi)
        sens_map_sky_veto = healpy.sphtfunc.smoothing(self.sc_map_in_sky(self.sens_map_veto), 5. / 180. * np.pi)
        sens_map_sky_isgri = healpy.sphtfunc.smoothing(self.sc_map_in_sky(self.sens_map_isgri), 5. / 180. * np.pi)
        sens_map_sky_picsit = healpy.sphtfunc.smoothing(self.sc_map_in_sky(self.sens_map_picsit), 5. / 180. * np.pi)

        if self.do_burst_analysis:
            bestarea_sky = healpy.sphtfunc.smoothing(self.sc_map_in_sky(self.bestarea), 5. / 180. * np.pi)
            locmap_sky = healpy.sphtfunc.smoothing(self.sc_map_in_sky(self.locmap), 5. / 180. * np.pi)

        good_mask = lambda x: sens_map_sky < sens_map_sky.min() * x
        print("good for", [sum(self.target_map[good_mask(x)]) for x in [1.01, 1.1, 1.2, 1.5, 2.]])

        print("very good", sens_map_sky.min())
        print("typical bad", sum(self.target_map[~good_mask(1.2)] * sens_map_sky[~good_mask(1.2)]) / sum(),
            self.target_map[~good_mask(1.2)])
        print("typical good", sum(self.target_map[good_mask(1.2)] * sens_map_sky[good_mask(1.2)]) / sum(),
            self.target_map[good_mask(1.2)])

        # map_sc=healpy.sphtfunc.smoothing(map_sc,5./180.*pi)
        # target_cum_map_sm=healpy.sphtfunc.smoothing(self.target_cum_map,2./180.*pi)

        # overplot=[]
        try:
            o_isgrimap = np.loadtxt(gzip.open("isgri_sens.txt.gz"))
            o_isgrimap[np.isnan(o_isgrimap) | np.isinf(o_isgrimap)] = 0

            isgrimap = healpy.get_interp_val(o_isgrimap,
                                             self.sky_coord.represent_as("physicsspherical").theta.rad,
                                             self.sky_coord.represent_as("physicsspherical").phi.rad,
                                             )

            o_jemxmap = np.loadtxt(gzip.open("jemx_sens.txt.gz"))
            o_jemxmap[np.isnan(o_jemxmap) | np.isinf(o_jemxmap)] = 0
            jemxmap = healpy.get_interp_val(o_jemxmap,
                                            self.sky_coord.represent_as("physicsspherical").theta.rad,
                                            self.sky_coord.represent_as("physicsspherical").phi.rad,
                                            )

            o_spimap = np.loadtxt(gzip.open("spi_sens.txt.gz"))
            o_spimap[np.isnan(o_spimap) | np.isinf(o_spimap)] = 0
            spimap = healpy.get_interp_val(o_spimap,
                                           self.sky_coord.represent_as("physicsspherical").theta.rad,
                                           self.sky_coord.represent_as("physicsspherical").phi.rad,
                                           )

            cover_info = {}

            for detname, detmap, o_detmap in [("isgri", isgrimap, o_isgrimap),
                                              ("jemx", jemxmap, o_jemxmap),
                                              ("spi", spimap ** 2, o_spimap ** 2)]:
                bestsens = min(o_detmap[o_detmap > 5.9e-6 ** 2])
                print("min", bestsens, bestsens ** 0.5)
                cover = (detmap < bestsens * 20 ** 2) & (detmap > 0)
                print("contained in", detname, "area", sum(cover) / float(cover.shape[0]) * 4 * np.pi * (180 / np.pi) ** 2, self.target_map.sum(), \
                    self.target_map[cover].sum(), sum(cover) * 1. / cover.shape[0], sum(cover) * 1. / cover.shape[
                    0] * 4 * np.pi * (180 / np.pi) ** 2)

                cover_info[detname] = dict(
                    area_deg2=sum(cover) / float(cover.shape[0]) * 4 * np.pi * (180 / np.pi) ** 2,
                    target_coverage=self.target_map[cover].sum()
                )

            json.dump(cover_info, open("coverage_info.json", "w"))

            overplot = [(self.target_cum_map, "k", [50, 90]),
                        (spimap, "g", spimap[spimap > 0].min() * 20),
                        (jemxmap, "b", jemxmap[jemxmap > 0].min() * 20 ** 2),
                        (isgrimap, "r", isgrimap[isgrimap > 0].min() * 20 ** 2)],
        except Exception as e:
            raise
            if sum(target_cum_map_sm > 0) > 0:
                overplot = [(target_map_sm, "k", None)]
                # overplot.append((locmap_sky, "gist_gray", None))



                # cover=theta_sc_rad>120./180*pi
                # print("contained in >120",map_px[cover].sum(),sum(cover)*1./cover.shape[0],sum(cover)*1./cover.shape[0]*4*pi*(180/pi)**2)

                # cover=(theta_sc_rad>80./180*pi) & (theta_sc_rad<120./180*pi)
                # print("contained in 80-120",map_px[cover].sum(),sum(cover)*1./cover.shape[0],sum(cover)*1./cover.shape[0]*4*pi*(180/pi)**2)

                # cover=(theta_sc_rad<80./180*pi)
                # print("contained in <80",map_px[cover].sum(),sum(cover)*1./cover.shape[0],sum(cover)*1./cover.shape[0]*4*pi*(180/pi)**2)
                #       sens_mp_sky[]*=2

        for body_name in "earth", "moon", "sun":
            bd = self.sc['bodies'][body_name]
            body_coord_sc = SkyCoord(bd['body_in_sc'][1], bd['body_in_sc'][0], 1, unit=(u.deg, u.deg),
                                     representation="physicsspherical")
            body_coord = SkyCoord(bd['body_ra'], bd['body_dec'], unit=(u.deg, u.deg))
            print("body:", body_name, bd)
            print("body coordinates:", bd['body_ra'], bd['body_dec'], body_coord)
            sens_map_sky[self.sky_coord.separation(body_coord).degree < bd["body_size"]] = 1e9

        healpy.write_map("sens_map_sky_" + self.tag + ".fits", sens_map_sky)

        if self.do_burst_analysis:
            p = healtics.plot_with_ticks(locmap_sky, cmap="jet", title="",
                                         unit="",
                                         vmin=1, vmax=9,
                                         overplot=overplot)
            plot.plot("sky_locmap.png", format='png', dpi=100)

            p = healtics.plot_with_ticks(locmap_sky, cmap="jet", title="",
                                         unit="",
                                         vmin=1, vmax=locmap_sky.max(),
                                         overplot=overplot)
            plot.plot("sky_locmap_full.png", format='png', dpi=100)

        def saveplot(prefix):
            plot.plot(prefix + self.tag + ".svg", format='svg', dpi=200)

        p = healtics.plot_with_ticks(sens_map_sky * self.sens_scale_b, cmap="YlOrBr", title="",
                                     unit="$10^{%i} \mathrm{erg^{ }cm^{-2}}$" % self.sens_scale_e,
                                     overplot=overplot,
                                     vmin=self.sens_scale_b, vmax=10 * self.sens_scale_b)
        saveplot("sky_sens_")

        p = healtics.plot_with_ticks(sens_map_sky_acs * self.sens_scale_b, cmap="YlOrBr", title="",
                                     overplot=overplot,
                                     unit="$10^{%i} \mathrm{erg^{ }cm^{-2}}$" % self.sens_scale_e,
                                     vmin=self.sens_scale_b, vmax=10 * self.sens_scale_b)
        saveplot("sky_sens_acs_")

        p = healtics.plot_with_ticks(sens_map_sky_veto * self.sens_scale_b, cmap="YlOrBr", title="",
                                     overplot=overplot,
                                     unit="$10^{%i} \mathrm{erg^{ }cm^{-2}}$" % self.sens_scale_e,
                                     vmin=self.sens_scale_b, vmax=10 * self.sens_scale_b)
        saveplot("sky_sens_veto_")

        p = healtics.plot_with_ticks(sens_map_sky_isgri * self.sens_scale_b, cmap="YlOrBr", title="",
                                     overplot=overplot,
                                     unit="$10^{%i} \mathrm{erg^{ }cm^{-2}}$" % self.sens_scale_e,
                                     vmin=self.sens_scale_b, vmax=10 * self.sens_scale_b)
        saveplot("sky_sens_isgri_")

        p = healtics.plot_with_ticks(sens_map_sky_picsit * self.sens_scale_b, cmap="YlOrBr", title="",
                                     overplot=overplot,
                                     unit="$10^{%i} \mathrm{erg^{ }cm^{-2}}$" % self.sens_scale_e,
                                     vmin=self.sens_scale_b, vmax=10 * self.sens_scale_b)
        saveplot("sky_sens_picsit_")

        if self.do_burst_analysis:
            p = healtics.plot_with_ticks(bestarea_sky * 100, cmap="YlOrBr", title="",
                                         overplot=overplot,
                                         unit="% of the sky",
                                         vmin=1, vmax=100)
            plot.plot("sky_sens_locarea_" + self.tag + ".png", format='png', dpi=100)

            # p=healtics.plot_with_ticks(sens_mp_sky,cmap="jet",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=[(map_px,"gist_gray",None),(spimap,"summer",20),(jemxmap,"winter",20**2),(isgrimap,"autumn",20**2)],vmin=0,vmax=sens_mp_sky.max())
            # p=healtics.plot_with_ticks(sens_mp_sky,cmap="YlOrBr",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=[(map_px,"gist_gray",None),(spimap,"summer",20),(jemxmap,"winter",20**2),(isgrimap,"autumn",20**2)],vmin=0,vmax=sens_mp_sky.max())
            # p=healtics.plot_with_ticks(sens_mp_sky,cmap="YlOrBr",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=healpy.sphtfunc.smoothing(map_px,5./180.*pi),vmin=1.5,vmax=15)

            # plot.plot("sky_sens.svg", format='svg', dpi=200)

    def plot_raw_sky(self):
        healpy.mollview(self.target_map, cmap='YlOrBr')
        healpy.projtext(self.scx.represent_as("physicsspherical").theta.rad,
                        self.scx.represent_as("physicsspherical").phi.rad,
                        "scx: SCZ")
        healpy.projscatter(self.scx.represent_as("physicsspherical").theta.rad,
                           self.scx.represent_as("physicsspherical").phi.rad)
        healpy.projtext(self.scy.represent_as("physicsspherical").theta.rad,
                        self.scy.represent_as("physicsspherical").phi.rad,
                        "scy: SCY")
        healpy.projscatter(self.scy.represent_as("physicsspherical").theta.rad,
                           self.scy.represent_as("physicsspherical").phi.rad)

        healpy.projtext(self.scz.represent_as("physicsspherical").theta.rad,
                        self.scz.represent_as("physicsspherical").phi.rad,
                        "scz: SCX")
        healpy.projscatter(self.scz.represent_as("physicsspherical").theta.rad,
                           self.scz.represent_as("physicsspherical").phi.rad)
        healpy.graticule()

        plot.plot("rawsky.png")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='make counterpart')
    parser.add_argument('--target-map', dest='targetmap', action='store', default="",
                        help='map of the target')
    parser.add_argument('--target-position', dest='targetposition', action='store', default="",
                        help='location of the target')
    parser.add_argument('--utc', dest='utc', action='store',
                        help='utc')
    parser.add_argument('--t1', dest='t1', action='store', default="0",
                        help='t1')
    parser.add_argument('--t2', dest='t2', action='store', default="0",
                        help='t2')

    args = parser.parse_args()

    Counterpart(use_target_map_fn=args.targetmap, use_utc=args.utc, use_t1=float(args.t1), use_t2=float(args.t2)).get()

