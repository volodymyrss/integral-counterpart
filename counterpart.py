import healpy
import dataanalysis as da
from numpy import *
import plot
import healtics
import pyfits
            
from astropy.coordinates import SkyCoord,PhysicsSphericalRepresentation
from astropy import units as u

#from numpy import *
import numpy as np
import string
import sys
import gzip

import integralclient

def transform_rmap(rmap):
    nside = healpy.npix2nside(rmap.shape[0])
    npx = arange(healpy.nside2npix(nside))
    theta, phi = healpy.pix2ang(nside, npx)
    #SkyCoord(phi, theta, 1, unit=(u.rad, u.rad), representation="physicsspherical")
    ntheta=pi-theta
    nphi=pi+phi #+ or -???
    nphi[nphi>pi*2]-=pi*2
    return healpy.get_interp_val(rmap,ntheta,nphi)


class Counterpart(da.DataAnalysis):
    target_map_fn=None
    utc=None


    def main(self):
        self.sc=integralclient.get_sc(self.utc)

        self.target_map = healpy.read_map(self.target_map_fn)
        self.nside = healpy.npix2nside(self.target_map.shape[0])
        self.compute_transform_grids()

        response_mp_acs = transform_rmap(array(integralclient.get_response_map(target="ACS", lt=100, alpha=-0.5, epeak=500, kind="response")))
        response_mp_veto = transform_rmap(array(integralclient.get_response_map(target="VETO", lt=80, alpha=-0.5, epeak=500, kind="response")))
        response_mp_isgri = transform_rmap(array(integralclient.get_response_map(target="ISGRI", lt=30, alpha=-0.5, epeak=500, kind="response")))

        def get_count_limit(target,scale):
            hk=integralclient.get_hk(target=target,utc=self.utc,span=30.01,t1=0,t2=0,ra=0,dec=0,rebin=scale)['lc']
            print target,":",scale,hk
            return hk['std bkg']*(3+hk['maxsig'])*hk['timebin']

        #scales=[1,]
        scales=[8]
        
        for scale in scales:
            acs_lim=get_count_limit("ACS",scale)
            veto_lim=get_count_limit("VETO",scale)
            isgri_lim=get_count_limit("ISGRI",scale)
            #acs_lim=300
            #veto_lim=300
            #isgri_lim=300


            print "ACS, Veto, ISGRI",acs_lim,veto_lim,isgri_lim
       
            sens_map_acs=response_mp_acs*acs_lim
            sens_map_veto=response_mp_veto*veto_lim
            sens_map_isgri=response_mp_isgri*isgri_lim
            sens_map=response_mp_acs*acs_lim
            sens_map[sens_map>sens_map_veto]=sens_map_veto[sens_map>sens_map_veto]

            na=sens_map[~isnan(sens_map) & (sens_map>0)].min()

            na_e=int(log10(na))
            na_b=int(na*10/10**na_e)/10.

            na=na_b*10**na_e

            print "best ACS",na
            nv=sens_map_veto[~isnan(sens_map_veto) & (sens_map_veto>0)].min()
            print "best VETO",nv
            sens_map/=na
            sens_map_veto/=na

            self.sens_scale=na
            self.sens_scale_e=na_e
            self.sens_scale_b=na_b

            self.sens_map=sens_map
            self.sens_map_acs = sens_map_acs
            self.sens_map_veto = sens_map_veto
            self.sens_map_isgri=sens_map_isgri
            self.compute_map()

    def get_grid(self,nside=None):
        nside=nside if nside is not None else self.nside
        npx=arange(healpy.nside2npix(nside))
        theta,phi=healpy.pix2ang(nside,npx)
        return SkyCoord(phi,theta, 1,unit=(u.rad, u.rad),representation="physicsspherical")


    def compute_transform_grids(self):
        sky_coord = self.get_grid()
        self.sky_coord=sky_coord

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
            arctan2(x_inxy * self.scy.cartesian.x + y_inxy * self.scy.cartesian.y + z_inxy * self.scy.cartesian.z,
                    x_inxy * self.scx.cartesian.x + y_inxy * self.scx.cartesian.y + z_inxy * self.scx.cartesian.z),
            arccos(on_scz),
            1,
            unit=(u.rad, u.rad),
            representation="physicsspherical")

        self.sky_in_sc=sky_in_sc
        self.sc_in_sky=sc_in_sky

        return sky_in_sc,sc_in_sky

    def sky_map_in_sc(self,sky_map):
        return healpy.get_interp_val(sky_map,
                        theta=self.sky_in_sc.represent_as("physicsspherical").theta.rad,
                        phi=self.sky_in_sc.represent_as("physicsspherical").phi.rad)

    def sc_map_in_sky(self,sc_map):
        return healpy.get_interp_val(sc_map,
                        theta=self.sc_in_sky.represent_as("physicsspherical").theta.rad,
                        phi=self.sc_in_sky.represent_as("physicsspherical").phi.rad)

    def plot_sky_diagram(self):
        target_map_in_sc = self.sky_map_in_sc(self.target_map)
        target_map_in_sc[abs(self.sky_coord.represent_as("physicsspherical").theta.deg) < 20] += target_map_in_sc.max() / 5.
        target_map_in_sc[abs(self.sky_coord.represent_as("physicsspherical").theta.deg) > 140] += target_map_in_sc.max() / 5.
        target_map_in_sc[(self.sky_coord.represent_as("physicsspherical").phi.deg < 30) | (
            self.sky_coord.represent_as("physicsspherical").phi.deg > 360 - 30)] += target_map_in_sc.max() / 10.
        healpy.mollview(target_map_in_sc, cmap='YlOrBr')
        healpy.projscatter(np.pi/2,0)
        healpy.graticule()
        plot.plot()

    def compute_map(self):
        healpy.mollview(self.sens_map_acs, cmap="YlOrBr")
        healpy.graticule()

        plot.plot()

        sens_map_sky = healpy.sphtfunc.smoothing(self.sc_map_in_sky(self.sens_map), 5. / 180. * pi)
        sens_map_sky_acs = healpy.sphtfunc.smoothing(self.sc_map_in_sky(self.sens_map_acs), 5. / 180. * pi)
        sens_map_sky_veto = healpy.sphtfunc.smoothing(self.sc_map_in_sky(self.sens_map_veto), 5. / 180. * pi)
        sens_map_sky_isgri = healpy.sphtfunc.smoothing(self.sc_map_in_sky(self.sens_map_isgri), 5. / 180. * pi)

        good_mask=lambda x:sens_map_sky<sens_map_sky.min()*x
        print "good for",[sum(self.target_map[good_mask(x)]) for x in [1.01,1.1,1.2,1.5,2.]]

        print "very good",sens_map_sky.min()
        print "typical bad",sum(self.target_map[~good_mask(1.2)]*sens_map_sky[~good_mask(1.2)])/sum(self.target_map[~good_mask(1.2)])
        print "typical good",sum(self.target_map[good_mask(1.2)]*sens_map_sky[good_mask(1.2)])/sum(self.target_map[good_mask(1.2)])

        #map_sc=healpy.sphtfunc.smoothing(map_sc,5./180.*pi)
        target_map_sm=healpy.sphtfunc.smoothing(self.target_map,2./180.*pi)

        if True:
            o_isgrimap=loadtxt(gzip.open("isgri_sens.txt.gz"))
            o_isgrimap[isnan(o_isgrimap) | isinf(o_isgrimap)]=0
            isgrimap=healpy.get_interp_val(o_isgrimap,
                                            self.sky_coord.represent_as("physicsspherical").theta.rad,
                                            self.sky_coord.represent_as("physicsspherical").phi.rad,
                                           )

            o_jemxmap=loadtxt(gzip.open("jemx_sens.txt.gz"))
            o_jemxmap[isnan(o_jemxmap) | isinf(o_jemxmap)]=0
            jemxmap=healpy.get_interp_val(o_jemxmap,
                                            self.sky_coord.represent_as("physicsspherical").theta.rad,
                                            self.sky_coord.represent_as("physicsspherical").phi.rad,
                                           )

            o_spimap=loadtxt(gzip.open("spi_sens.txt.gz"))
            o_spimap[isnan(o_spimap) | isinf(o_spimap)]=0
            spimap=healpy.get_interp_val(o_spimap,
                                            self.sky_coord.represent_as("physicsspherical").theta.rad,
                                            self.sky_coord.represent_as("physicsspherical").phi.rad,
                                           )

        for detname,detmap,o_detmap in [("isgri",isgrimap,o_isgrimap),
                               ("jemx",jemxmap,o_jemxmap),
                               ("spi",spimap**2,o_spimap**2)]:
            bestsens=min(o_detmap[o_detmap>5.9e-6**2])
            print "min",bestsens,bestsens**0.5
            cover=(detmap<bestsens*20**2) & (detmap>0)
            print "contained in",detname,self.target_map.sum(),self.target_map[cover].sum(),sum(cover)*1./cover.shape[0],sum(cover)*1./cover.shape[0]*4*pi*(180/pi)**2

        #cover=theta_sc_rad>120./180*pi
        #print "contained in >120",map_px[cover].sum(),sum(cover)*1./cover.shape[0],sum(cover)*1./cover.shape[0]*4*pi*(180/pi)**2
        
        #cover=(theta_sc_rad>80./180*pi) & (theta_sc_rad<120./180*pi)
        #print "contained in 80-120",map_px[cover].sum(),sum(cover)*1./cover.shape[0],sum(cover)*1./cover.shape[0]*4*pi*(180/pi)**2
        
        #cover=(theta_sc_rad<80./180*pi)
        #print "contained in <80",map_px[cover].sum(),sum(cover)*1./cover.shape[0],sum(cover)*1./cover.shape[0]*4*pi*(180/pi)**2
 #       sens_mp_sky[]*=2

        for body_name in "earth","moon","sun":
            bd=self.sc['bodies'][body_name]
            body_coord_sc=SkyCoord(bd['body_in_sc'][1],bd['body_in_sc'][0],1,unit=(u.deg,u.deg),representation="physicsspherical")
            body_coord=SkyCoord(bd['body_ra'],bd['body_dec'],unit=(u.deg,u.deg))
            print("body:",body_name,bd)
            print("body coordinates:",bd['body_ra'],bd['body_dec'],body_coord)
            sens_map_sky[self.sky_coord.separation(body_coord).degree<bd["body_size"]]=1e9

        p = healtics.plot_with_ticks(sens_map_sky * self.sens_scale_b, cmap="YlOrBr", title="",
                                     overplot=[(target_map_sm, "gist_gray", None), (spimap, "summer", 20),
                                               (jemxmap, "winter", 20 ** 2), (isgrimap, "autumn", 20 ** 2)],
                                     unit="$10^{%i} \mathrm{erg^{}cm^{-2} s^{-1}}$" % self.sens_scale_e,
                                     vmin=self.sens_scale_b, vmax=10 * self.sens_scale_b)
        plot.plot("sky_sens.png", format='png', dpi=100)

        p = healtics.plot_with_ticks(sens_map_sky_veto * self.sens_scale_b, cmap="YlOrBr", title="",
                                     overplot=[(target_map_sm, "gist_gray", None), (spimap, "summer", 20),
                                               (jemxmap, "winter", 20 ** 2), (isgrimap, "autumn", 20 ** 2)],
                                     unit="$10^{%i} \mathrm{erg^{}cm^{-2} s^{-1}}$" % self.sens_scale_e,
                                     vmin=self.sens_scale_b, vmax=10 * self.sens_scale_b)
        plot.plot("sky_sens_veto.png", format='png', dpi=100)

        p = healtics.plot_with_ticks(sens_map_sky_isgri * self.sens_scale_b, cmap="YlOrBr", title="",
                                     overplot=[(target_map_sm, "gist_gray", None), (spimap, "summer", 20),
                                               (jemxmap, "winter", 20 ** 2), (isgrimap, "autumn", 20 ** 2)],
                                     unit="$10^{%i} \mathrm{erg^{}cm^{-2} s^{-1}}$" % self.sens_scale_e,
                                     vmin=self.sens_scale_b, vmax=10 * self.sens_scale_b)
        plot.plot("sky_sens_isgri.png", format='png', dpi=100)

        #p=healtics.plot_with_ticks(sens_mp_sky,cmap="jet",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=[(map_px,"gist_gray",None),(spimap,"summer",20),(jemxmap,"winter",20**2),(isgrimap,"autumn",20**2)],vmin=0,vmax=sens_mp_sky.max())
        #p=healtics.plot_with_ticks(sens_mp_sky,cmap="YlOrBr",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=[(map_px,"gist_gray",None),(spimap,"summer",20),(jemxmap,"winter",20**2),(isgrimap,"autumn",20**2)],vmin=0,vmax=sens_mp_sky.max())
        #p=healtics.plot_with_ticks(sens_mp_sky,cmap="YlOrBr",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=healpy.sphtfunc.smoothing(map_px,5./180.*pi),vmin=1.5,vmax=15)

        #plot.plot("sky_sens.svg", format='svg', dpi=200)


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


if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='make counterpart')
    parser.add_argument('--target-map', dest='targetmap', action='store',
                                                    help='map of the target')
    parser.add_argument('--utc', dest='utc', action='store',
                                                    help='utc')

    args = parser.parse_args()

    Counterpart(use_target_map_fn=args.targetmap,use_utc=args.utc).get()
