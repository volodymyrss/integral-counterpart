import healpy
import dataanalysis as da
from numpy import *
import plot
import healtics
import pyfits
            
from astropy.coordinates import SkyCoord,PhysicsSphericalRepresentation
from astropy import units as u

from numpy import *
import string
import sys

import integralclient

class Counterpart(da.DataAnalysis):
    target_map_fn=None
    utc=None

    def main(self):

        self.sc=integralclient.get_sc(self.utc)

        response_mp=array(integralclient.get_response_map(target="ACS",lt=100,alpha=-0.5,epeak=1000,kind="response"))
        response_mp_veto=array(integralclient.get_response_map(target="VETO",lt=80,alpha=-0.5,epeak=1000,kind="response"))
        response_mp_picsit=array(integralclient.get_response_map(target="PICsIT",lt=80,alpha=-0.5,epeak=1000,kind="response"))

        def get_count_limit(target,scale):
            acshk=integralclient.get_hk(target=target,utc=self.utc,span=30.01,t1=0,t2=0,ra=0,dec=0,rebin=scale)['lc']
            print "ACS:",scale,acshk
            return acshk['std bkg']*3+acshk['maxsig']

        scales=[0.1,1,10]
        
        for scale in scales:
            acs_lim=get_count_limit("ACS",scale)
            veto_lim=get_count_limit("VETO",scale)
            picsit_lim=get_count_limit("SPTI1234",scale)
       
            sens_mp=response_mp*acs_lim
            sens_mp_veto=response_mp_veto*veto_lim
            sens_mp_picsit=response_mp_picsit*picsit_lim
            sens_mp[sens_mp>sens_mp_veto]=sens_mp_veto[sens_mp>sens_mp_veto]

            self.sens_mp=sens_mp
            self.sens_mp_picsit=sens_mp_picsit
            self.compute_map()

    def get_grid_for_map(self,mp):
        npx=arange(mp.shape[0])
        nside=healpy.npix2nside(npx.shape[0])
        theta_rad,ra_rad=healpy.pix2ang(nside,npx)
        ra=ra_rad/pi*180
        phi_rad=ra_rad
        phi=ra
        theta=theta_rad/pi*180
        dec=90-theta 
        dec_rad=dec/180*pi
        #return ra,dec
        return SkyCoord(ra,dec, unit=(u.deg, u.deg))

    def compute_map(self):
        sens_mp=self.sens_mp
        ####

        map_px=healpy.read_map(self.target_map_fn)
        
        nside=healpy.npix2nside(map_px.shape[0])

        npx=arange(healpy.nside2npix(nside))
        theta_rad,ra_rad=healpy.pix2ang(nside,npx)
        ra=ra_rad/pi*180
        phi_rad=ra_rad
        phi=ra
        theta=theta_rad/pi*180
        dec=90-theta #!!!
        #dec=theta-90 #!!!
        dec_rad=dec/180*pi
        
        ### source!
        sky_coord=self.get_grid_for_map(map_px)
        sens_coord=self.get_grid_for_map(sens_mp)

        for body_name in "earth","moon","sun":
            bd=self.sc['bodies'][body_name]
            body_coord=SkyCoord(bd['body_in_sc'][1],bd['body_in_sc'][0],1,unit=(u.deg,u.deg),representation="physicsspherical")
            print("body:",body_name,bd)
            sens_mp[sens_coord.separation(body_coord).degree<bd["body_size"]]=1e9
        #map_px=s
        #map_px[map_px>3.]=0
        ### /source

        healpy.mollview(map_px,cmap="YlOrBr")
        healpy.graticule()

        #fig=healtics.plot_with_ticks(map_px,cmap="YlOrBr")

        x=cos(ra_rad)*cos(dec_rad)
        y=sin(ra_rad)*cos(dec_rad)
        z=sin(dec_rad)
        
        sc=self.sc
        ra_scx  = sc['scx']['ra']
        dec_scx = sc['scx']['dec']
        ra_scz  = sc['scz']['ra']
        dec_scz = sc['scz']['dec']

        scx=SkyCoord(ra_scx, dec_scx, frame='icrs', unit='deg')
        scz=SkyCoord(ra_scz, dec_scz, frame='icrs', unit='deg')

        print ra_scx,dec_scx

        #ra_b_h,dec_b_h="+05 12 40.80", "-07 02 56.0"
        #b=SkyCoord(ra_b_h+" "+dec_b_h, frame='icrs', unit=(u.hourangle, u.deg))



        anti_scx=SkyCoord(180+scx.ra.deg, -scx.dec.deg, frame='icrs', unit='deg')


        ra_scx_rad=ra_scx/180*pi
        dec_scx_rad=dec_scx/180*pi
        theta_scx=90-dec_scx
        phi_scx=ra_scx
        
        
        x_scx=cos(ra_scx_rad)*cos(dec_scx_rad)
        y_scx=sin(ra_scx_rad)*cos(dec_scx_rad)
        z_scx=sin(dec_scx_rad)
        
        
        #ra_scz,dec_scz=211.450,-7.413
        theta_scz=90-dec_scz
        phi_scz=ra_scz

        
        ra_scz_rad=ra_scz/180*pi
        dec_scz_rad=dec_scz/180*pi
        x_scz=cos(ra_scz_rad)*cos(dec_scz_rad)
        y_scz=sin(ra_scz_rad)*cos(dec_scz_rad)
        z_scz=sin(dec_scz_rad)
    
        x_scy=y_scx*z_scz-z_scx*y_scz
        y_scy=z_scx*x_scz-x_scx*z_scz
        z_scy=x_scx*y_scz-y_scx*x_scz

        print "X dot Z",x_scx*x_scz+y_scx*y_scz+z_scx*z_scz
        print "Y dot Z",x_scy*x_scz+y_scy*y_scz+z_scy*z_scz
        print "Y dot X",x_scy*x_scx+y_scy*y_scx+z_scy*z_scx

        phi_scy=arctan2(y_scy,x_scy)/pi*180
        theta_scy=arccos(z_scy)/pi*180
        ra_scy=phi_scy
        dec_scy=90-theta_scy

        
        scx=SkyCoord(ra_scx, dec_scx, frame='icrs', unit='deg')
        scz=SkyCoord(ra_scz, dec_scz, frame='icrs', unit='deg')
        scy=SkyCoord(ra_scy, dec_scy, frame='icrs', unit='deg')
        print "scx",scx
        print "scy",scy
        print "scz",scz
        print "X-Z",scx.separation(scz).degree
        print "Y-Z",scz.separation(scy).degree
        print "X-Y",scx.separation(scy).degree
        #print "X-Y sep",angsep(ra_scx,dec_scx,ra_scy,dec_scy)
        #print "X-Z sep",angsep(ra_scx,dec_scx,ra_scz,dec_scz)
        #print "Y-Z sep",angsep(ra_scy,dec_scy,ra_scz,dec_scz)

        healpy.projtext(ra_scx,dec_scx,"SCX %.5lg,%.5lg"%(ra_scx,dec_scx),lonlat=True)
        healpy.projscatter(ra_scx,dec_scx,lonlat=True)

        healpy.projtext(ra_scx+180,-dec_scx,"anti-SCX %.5lg,%.5lg"%(180+ra_scx,-dec_scx),lonlat=True)
        healpy.projscatter(ra_scx+180,-dec_scx,lonlat=True)

        healpy.projtext(ra_scz,dec_scz,"SCZ",lonlat=True)#.set_z_order(20)
        healpy.projscatter(ra_scz,dec_scz,lonlat=True)
        healpy.projtext(ra_scy,dec_scy,"SCY",lonlat=True)#.set_z_order(20)
        healpy.projscatter(ra_scy,dec_scy,lonlat=True)
        
        healpy.projtext(255.562226,-47.681015,"candidate",lonlat=True)#.set_z_order(20)
        healpy.projscatter(255.562226,-47.681015,lonlat=True)

        healpy.projtext(60,30,"60,30",lonlat=True)#.set_z_order(20)
        healpy.projtext(60,-30,"60,-30",lonlat=True)#.set_z_order(20)
        healpy.projtext(-60,-30,"-60,-30",lonlat=True)#.set_z_order(20)
        healpy.projtext(-60,30,"-60,30",lonlat=True)#.set_z_order(20)


        print "scx ra,dec:",ra_scx,dec_scx
        print "scy ra,dec:",ra_scy,dec_scy
        print "scz ra,dec:",ra_scz,dec_scz

        #plot.plot("rawsky.png")


        print x_scx*x_scy+y_scx*y_scy+z_scx*z_scy
        print x_scz*x_scy+y_scz*y_scy+z_scz*z_scy

        b_on_scx=x_scx*x+y_scx*y+z_scx*z

        print b_on_scx.max()

        theta_rad=arccos(b_on_scx)
        theta=theta_rad/pi*180
    
        x_scpx=x-b_on_scx*x
        y_scpx=y-b_on_scx*y
        z_scpx=z-b_on_scx*z
        n=(x_scpx**2+y_scpx**2+z_scpx**2)**0.5
        x_scpx/=n*2**0.5
        y_scpx/=n*2**0.5
        z_scpx/=n*2**0.5

        print x_scpx.max(),y_scpx.max(),z_scpx.max()
        a=x_scpx*x_scz+y_scpx*y_scz+z_scpx*z_scz
        print a.max(),a.min()

        #phi_rad=arccos(a)
        #phi=phi_rad/pi*180

        #nprint phi_rad.max(),phi_rad.min()

        #m=x_scpx*x_scy+y_scpx*y_scy+z_scpx*z_scy<0
        #phi_rad[m]=360-phi_rad[m]

        Theta_sc_rad,Phi_sc_rad=healpy.pix2ang(nside,npx)
        Theta_sc=Theta_sc_rad/pi*180
        Phi_sc=Phi_sc_rad/pi*180

        #print "sat grid Theta,Phi",Theta_sc,Phi_sc
        
        Ra_sc_rad=Phi_sc_rad
        Dec_sc_rad=Theta_sc_rad-pi/2.
       
        #print "sat grid Ra,Dec",Ra_sc,Dec_sc
        
        Y_sc=sin(Ra_sc_rad)*sin(Theta_sc_rad)
        Z_sc=cos(Ra_sc_rad)*sin(Theta_sc_rad)
        X_sc=cos(Theta_sc_rad)

        i=argmin(Theta_sc_rad)
                

        #print "sat grid X,Y,Z",X_sc,Y_sc,Z_sc

        x_sc=X_sc*x_scx+Y_sc*x_scy+Z_sc*x_scz
        y_sc=X_sc*y_scx+Y_sc*y_scy+Z_sc*y_scz
        z_sc=X_sc*z_scx+Y_sc*z_scy+Z_sc*z_scz

        r_sc=(x_sc**2+y_sc**2+z_sc**2)**0.5

        #print "r range",r_sc.max(),r_sc.min()
        print "min Theta_sc_rad",x_sc[i],y_sc[i],z_sc[i]

        phi_sc_rad=arctan2(y_sc,x_sc)
        theta_sc_rad=arccos(z_sc/r_sc)
        theta_sc=theta_sc_rad*180/pi
        phi_sc=phi_sc_rad*180/pi

        map_px2=map_px[:]

        #map_px[abs(theta_sc)<3]=map_px.max()
        #map_px[(abs(theta_sc-90)<2) & (abs(phi_sc-0)<2)]=map_px.max()/2
        #map_px[(abs(theta_sc-30)<2)]=map_px.max()/2
        #map_px[(abs(theta_sc-60)<2)]=map_px.max()/2
        #map_px[(abs(theta_sc-120)<2)]=map_px.max()/2
        #map_px[(abs(theta_sc-150)<2)]=map_px.max()/2
        #map_px[abs(theta_sc-90)<2]=map_px.max()
        #map_px[abs(theta_sc-135)<2]=map_px.max()
        #map_px[abs(theta_sc-180)<2]=map_px.max()
        
        theta_rad,phi_rad=healpy.pix2ang(nside,npx)
        theta=theta_rad/pi*180
        phi=phi_rad/pi*180
    
        x=sin(theta_rad)*cos(phi_rad)
        y=sin(theta_rad)*sin(phi_rad)
        z=cos(theta_rad)

        on_scx=x*x_scx+y*y_scx+z*z_scx
        x_inyz=x-on_scx*x_scx
        y_inyz=y-on_scx*y_scx
        z_inyz=z-on_scx*z_scx
        r_inyz=(x_inyz*x_inyz+y_inyz*y_inyz+z_inyz*z_inyz)**0.5
        x_inyz/=r_inyz
        y_inyz/=r_inyz
        z_inyz/=r_inyz

        theta_insc_rad=arccos(on_scx)
        phi_insc_rad=arctan2(
                            x_inyz*x_scy+y_inyz*y_scy+z_inyz*z_scy,
                            x_inyz*x_scz+y_inyz*y_scz+z_inyz*z_scz,
                            )
        theta_insc=theta_insc_rad*180/pi
        phi_insc=phi_insc_rad*180/pi

        map_sc=healpy.get_interp_val(map_px,theta_sc_rad,phi_sc_rad)

        ##!!!
        map_sc[abs(theta)<20]+=map_sc.max()/5.
        map_sc[abs(theta)>140]+=map_sc.max()/5.
        map_sc[(phi<30) | (phi>360-30)]+=map_sc.max()/10.
        
        i=argmax(z_sc)
        print ra_scx,dec_scx,"vs",theta[i],phi[i],theta_sc[i],phi_sc[i],z_sc[i],r_sc[i],x_scx

        
        map_px2=healpy.get_interp_val(map_sc,theta_insc_rad,phi_insc_rad)
        map_px2[abs(theta_insc)<20]+=map_px.max()/5.
        map_px2[abs(phi_insc)<30]+=map_px.max()/10.
        map_px2[abs(theta_insc)>140]+=map_px.max()/5.
        

        #map_px=healpy.get_interp_val(map_sc,theta_rad,phi_rad)
       # map_px[theta_insc<20]+=map_px.max()/5.
        #map_px[abs(phi_sc)<30]+=map_px.max()/10.
        #map_px[abs(theta_sc)>140]+=map_px.max()/5.
        ##!!!

        p=healpy.mollview(map_px2,cmap='YlOrBr')
        healpy.projtext(phi_scx,dec_scx,"SCX",lonlat=True)
        healpy.projscatter(phi_scx,dec_scx,lonlat=True)
        healpy.projtext(phi_scx+180,-dec_scx,"anti-SCX",lonlat=True)
        healpy.projscatter(phi_scx+180,-theta_scx,lonlat=True)
        healpy.projtext(phi_scz,dec_scz,"SCZ",lonlat=True) #.set_z_order(20)
        healpy.projscatter(phi_scz,dec_scz,lonlat=True)
        healpy.projtext(phi_scy,dec_scy,"SCY",lonlat=True) #.set_z_order(20)
        healpy.projscatter(phi_scy,dec_scy,lonlat=True)
        healpy.graticule()
        #plot.plot("sky_sens.png")

        
        if False:
            for sel_ra,sel_dec in zip(selected_positions['RA'],selected_positions['DEC'])+[(ra_scx,dec_scx),(ra_scz,dec_scz),(ra_scy,dec_scy)]:
                print "selected",sel_ra,sel_dec
                sel_theta=sel_dec+90
                sel_phi=sel_ra-pi
                sel_d2=(sel_theta-theta)**2+(sel_phi-phi)**2
                sel_m=sel_d2<(0.2)**2
                print theta_sc[sel_m],phi_sc[sel_m]

        #sens_mp=loadtxt("/Integral/throng/savchenk/projects/spiacs/spiacsmodel/python/ligo/LVT151012/acs_map.txt")*1810/1e-7 * 1.15
        #sens_mp_veto=loadtxt("/Integral/throng/savchenk/projects/spiacs/spiacsmodel/python/ligo/LVT151012/veto_map.txt")*1010/1e-7*8**0.5 * 1.15

        #sens_mp=loadtxt("../../../mp_0.075-2_grbm_al_-0.5_bet_-2.5_e0_1.txt")*1410*1e-10/1e-7 * 1.15

        #sens_mp=loadtxt("compton_map.txt")
        #sens_mp=loadtxt("../mp_0.075-2_grbm_al_-0.5_bet_-1.5_e0_1.txt")*3
        #sens_mp=loadtxt("../mp_0.075-2_grbm_al_-1_bet_-2.5_e0_0.5.txt")*3
        #sens_mp=loadtxt("../mp_0.075-2_grbm_al_-1.5_bet_-2.5_e0_0.3.txt")*3
        #sens_mp=loadtxt("../mp.txt")*3
        sens_mp_sky=healpy.get_interp_val(sens_mp,pi-theta_insc_rad,phi_insc_rad-pi) ##
        sens_mp_sc=healpy.get_interp_val(sens_mp,pi-theta_rad,phi_rad-pi) ##
        #sens_mp_sky=healpy.get_interp_val(sens_mp,theta_insc_rad,phi_insc_rad-pi) ##
       # sens_mp_sc=healpy.get_interp_val(sens_mp,theta_rad,phi_rad-pi) ##


        sens_mp_picsit=self.sens_mp_picsit
        sens_mp_picsit=sens_mp_picsit.max()/sens_mp_picsit
        sens_mp_picsit[isnan(sens_mp_picsit)]=0
        sens_mp_picsit[isinf(sens_mp_picsit)]=0

        sens_mp_picsit_sky=healpy.get_interp_val(sens_mp_picsit,theta_sc_rad,phi_sc_rad-pi)
        #sens_mp_picsit_sky=healpy.sphtfunc.smoothing(sens_mp_picsit_sky,10./180.*pi)
        
        sens_mp_picsit_sc=healpy.get_interp_val(sens_mp_picsit,theta_rad,phi_rad-pi)

        print sens_mp_picsit_sky

       # healtics.plot_with_ticks(mp,cmap="YlOrBr",title="sky coordinates, Ra, Dec",overplot=mp)

        sens_mp_sky=healpy.sphtfunc.smoothing(sens_mp_sky,5./180.*pi)
        sens_mp_sky=sens_mp_sky

        good_mask=lambda x:sens_mp_sky<sens_mp_sky.min()*x
        print "good for",[sum(map_px[good_mask(x)]) for x in [1.01,1.1,1.2,1.5,2.]]

        print "very good",sens_mp_sky.min()
        print "typical bad",sum(map_px[~good_mask(1.2)]*sens_mp_sky[~good_mask(1.2)])/sum(map_px[~good_mask(1.2)])
        print "typical good",sum(map_px[good_mask(1.2)]*sens_mp_sky[good_mask(1.2)])/sum(map_px[good_mask(1.2)])

        map_sc=healpy.sphtfunc.smoothing(map_sc,5./180.*pi)
        map_px=healpy.sphtfunc.smoothing(map_px,2./180.*pi)

        if True:
            o_isgrimap=loadtxt("isgri_sens.txt")
            o_isgrimap[isnan(o_isgrimap) | isinf(o_isgrimap)]=0
            isgrimap=healpy.get_interp_val(o_isgrimap,theta_rad,phi_rad)

            o_jemxmap=loadtxt("jemx_sens.txt")
            o_jemxmap[isnan(o_jemxmap) | isinf(o_jemxmap)]=0
            jemxmap=healpy.get_interp_val(o_jemxmap,theta_rad,phi_rad)
            
            o_spimap=loadtxt("spi_sens.txt")
            o_spimap[isnan(o_spimap) | isinf(o_spimap)]=0
            spimap=healpy.get_interp_val(o_spimap,theta_rad,phi_rad)
        
        ## mess up
        theta_rad,phi_rad=healpy.pix2ang(nside,npx)
   #     phi_rad+=pi/2.*3
        phi_rad[phi_rad>2*pi]=phi_rad[phi_rad>2*pi]-2*pi
        sens_mp_sky=healpy.get_interp_val(sens_mp_sky,theta_rad,phi_rad)
        map_px=healpy.get_interp_val(map_px,theta_rad,phi_rad)

        #p=healtics.plot_with_ticks(sens_mp_picsit_sky,cmap="YlOrBr",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=map_px,vmin=1.5,vmax=15/2.)
        #p=healtics.plot_with_ticks(sens_mp_sky,cmap="YlOrBr",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=[(map_px,"gist_gray",None)],vmin=1.5,vmax=15)
        
        # containment 


        for detname,detmap,o_detmap in [("isgri",isgrimap,o_isgrimap),
                               ("jemx",jemxmap,o_jemxmap),
                               ("spi",spimap**2,o_spimap**2)]:
            bestsens=min(o_detmap[o_detmap>5.9e-6**2])
            print "min",bestsens,bestsens**0.5
            cover=(detmap<bestsens*20**2) & (detmap>0)
            print "contained in",detname,map_px.sum(),map_px[cover].sum(),sum(cover)*1./cover.shape[0],sum(cover)*1./cover.shape[0]*4*pi*(180/pi)**2

        # //

        cover=theta_sc_rad>120./180*pi
        print "contained in >120",map_px[cover].sum(),sum(cover)*1./cover.shape[0],sum(cover)*1./cover.shape[0]*4*pi*(180/pi)**2
        
        cover=(theta_sc_rad>80./180*pi) & (theta_sc_rad<120./180*pi)
        print "contained in 80-120",map_px[cover].sum(),sum(cover)*1./cover.shape[0],sum(cover)*1./cover.shape[0]*4*pi*(180/pi)**2
        
        cover=(theta_sc_rad<80./180*pi) 
        print "contained in <80",map_px[cover].sum(),sum(cover)*1./cover.shape[0],sum(cover)*1./cover.shape[0]*4*pi*(180/pi)**2
 #       sens_mp_sky[]*=2


        p=healtics.plot_with_ticks(sens_mp_sky,cmap="YlOrBr",title="",overplot=[(map_px,"gist_gray",None),(spimap,"summer",20),(jemxmap,"winter",20**2),(isgrimap,"autumn",20**2)],vmin=1,vmax=10)
        #p=healtics.plot_with_ticks(sens_mp_sky,cmap="jet",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=[(map_px,"gist_gray",None),(spimap,"summer",20),(jemxmap,"winter",20**2),(isgrimap,"autumn",20**2)],vmin=0,vmax=sens_mp_sky.max())
        #p=healtics.plot_with_ticks(sens_mp_sky,cmap="YlOrBr",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=[(map_px,"gist_gray",None),(spimap,"summer",20),(jemxmap,"winter",20**2),(isgrimap,"autumn",20**2)],vmin=0,vmax=sens_mp_sky.max())
        #p=healtics.plot_with_ticks(sens_mp_sky,cmap="YlOrBr",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=healpy.sphtfunc.smoothing(map_px,5./180.*pi),vmin=1.5,vmax=15)
        healpy.projtext(phi_scx,theta_scx,"SCX",lonlat=True)
        healpy.projscatter(phi_scx,theta_scx,lonlat=True)
        healpy.projtext(ra_scz,theta_scz,"SCZ",lonlat=True) #.set_z_order(20)
        healpy.projscatter(ra_scz,theta_scz,lonlat=True)
        healpy.projtext(ra_scy,dec_scy,"SCY",lonlat=True) #.set_z_order(20)
        healpy.projscatter(ra_scy,dec_scy,lonlat=True)
        plot.plot("sky_sens.svg", format='svg', dpi=200)
        #plot.plot("sky_sens.svg", format='svg', dpi=1000)
        #p=healtics.plot_with_ticks(sens_mp_sky,cmap="YlOrBr",title="INTEGRAL SPI-ACS 3 sigma upper limit in 1 second",overplot=healpy.sphtfunc.smoothing(map_px,5./180.*pi),vmin=1.5,vmax=15)
        
        #p=healpy.mollview(map_px,cmap='YlOrBr')
        p=healpy.mollview(map_px2,cmap='YlOrBr')
        #p=healpy.mollview(sens_mp_sky,cmap='YlOrBr')
        healpy.projtext(phi_scx,dec_scx,"SCX",lonlat=True)
        healpy.projscatter(phi_scx,dec_scx,lonlat=True)
        healpy.projtext(phi_scx+180,-dec_scx,"anti-SCX",lonlat=True)
        healpy.projscatter(phi_scx+180,-dec_scx,lonlat=True)
        healpy.projtext(ra_scz,dec_scz,"SCZ",lonlat=True) #.set_z_order(20)
        healpy.projscatter(ra_scz,dec_scz,lonlat=True)
        healpy.projtext(ra_scy,dec_scy,"SCY",lonlat=True) #.set_z_order(20)
        healpy.projscatter(ra_scy,dec_scy,lonlat=True)
        healpy.graticule()
#        p=healpy.mollview(map_px2,cmap='YlOrBr')
        #plot.plot("sky_sens.png", format='svg', dpi=1000)

        #healtics.plot_with_ticks(sens_mp_sc,cmap="YlOrBr",title="detector coordinates",overplot=map_sc)
        #plot.plot("sccoord.png")
        #healtics.plot_with_ticks(sens_mp_picsit_sc,cmap="YlOrBr",title="detector coordinates",overplot=map_sc)
        #plot.plot("sccoord2.png")
        #healtics.plot_with_ticks(map_sc,cmap="YlOrBr",title="detector coordinates",overplot=sens_mp)
        #plot.plot("sccoord3.png")


if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='make counterpart')
    parser.add_argument('--target-map', dest='targetmap', action='store',
                                                    help='map of the target')
    parser.add_argument('--utc', dest='utc', action='store',
                                                    help='utc')

    args = parser.parse_args()

    Counterpart(use_target_map_fn=args.targetmap,use_utc=args.utc).get()
