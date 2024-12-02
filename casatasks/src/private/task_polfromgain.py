# enable local tools
import os
from math import pi,floor,atan2,sin,cos,sqrt
import pylab as mypl

import casatools
from casatasks import casalog

mytb=casatools.table()
myme=casatools.measures()
mymd=casatools.msmetadata()

def polfromgain(vis,tablein,caltable,paoffset,minpacov):

    casalog.origin('polfromgain')

    casalog.post("Deriving calibrator linear polarization from gain ratios.")

    casalog.post("Requiring at least "+str(minpacov)+" deg of parallactic angle coverage for each antenna solution.")
    minpacovR=minpacov*pi/180.0

    try:

        if ((type(vis)==str) & (os.path.exists(vis))):
            mymd.open(vis)
        else:
            raise ValueError('Visibility data set not found - please verify the name')

        nant=mymd.nantennas()
        nfld=mymd.nfields()
        fldnames=mymd.fieldnames()
        nspw=mymd.nspw()

        rempol=False
        if ((type(tablein)==str) & (os.path.exists(tablein))):
            if type(caltable)==str and len(caltable)>0:

                if os.path.exists(caltable):
                    raise ValueError('Output caltable='+caltable+' exists.  Choose another name or delete it.')

                casalog.post("New caltable, "+caltable+", corrected for linear polarization, will be generated.")
                mytb.open(tablein)
                myout=mytb.copy(newtablename=caltable,deep=True)
                mytb.close()
                myout.close()
                rempol=True
            else:
                casalog.post("No new caltable will be generated")
                caltable=tablein
        else:
            raise ValueError('input calibration table not found - please verify the name')
            

        if paoffset!=0.0:
            casalog.post("NB: default band position angle will be offset by "+str(paoffset)+"deg.")

        # Field coords
        mytb.open(caltable+'/FIELD')
        dirs=mytb.getcol('DELAY_DIR')[:,0,:]
        mytb.close()
        
        # Must retrieve nominal feed angles from MS.FEED!
        mytb.open(vis+'/FEED')
        nfeed=mytb.nrows()
        fang=mytb.getcol('RECEPTOR_ANGLE')
        fspw=mytb.getcol('SPECTRAL_WINDOW_ID')
        fant=mytb.getcol('ANTENNA_ID')
        rang=mypl.zeros((nant,nspw));
        for ifeed in range(nfeed):
            rang[fant[ifeed],fspw[ifeed]]=fang[0,ifeed]
        mytb.close()

        R=mypl.zeros((nspw,nfld))
        Q=mypl.zeros((nspw,nfld))
        U=mypl.zeros((nspw,nfld))
        mask=mypl.zeros((nspw,nfld),dtype=bool)

        IQUV={}
        nomod=not rempol
        mytb.open(caltable,nomodify=nomod)
        uflds=mypl.unique(mytb.getcol('FIELD_ID'))
        uspws=mypl.unique(mytb.getcol('SPECTRAL_WINDOW_ID'))
        for ifld in uflds:
            rah=dirs[0,ifld]*12.0/pi
            decr=dirs[1,ifld]
            IQUV[fldnames[ifld]]={}
            for ispw in uspws:

                r=mypl.zeros(nant)
                q=mypl.zeros(nant)
                u=mypl.zeros(nant)
                antok=mypl.zeros(nant,dtype=bool)
                parangbyant={}  # we will remember parang outside main iant loop (for rempol on bad ants)
                
                casalog.post('Fld='+fldnames[ifld]+' Spw='+str(ispw)+':')

                for iant in range(nant):
                    qstring='FIELD_ID=='+str(ifld)+' && SPECTRAL_WINDOW_ID=='+str(ispw)+' && ANTENNA1=='+str(iant)
                    st=mytb.query(query=qstring)
                    nrows=st.nrows()
                    if nrows > 0:

                        times=st.getcol('TIME')
                        gains=st.getcol('CPARAM')
                        flags=st.getcol('FLAG')
                        flags=mypl.logical_or(flags[0,0,:],flags[1,0,:])  # 1D

                        # parang
                        parang=mypl.zeros(len(times))
                
                        apos=mymd.antennaposition(iant)
                        latr=myme.measure(apos,'WGS84')['m1']['value']
                        myme.doframe(apos)
                        har=mypl.zeros(nrows)
                        for itim in range(len(times)):
                            tm=myme.epoch('UTC',str(times[itim])+'s')
                            last=myme.measure(tm,'LAST')['m0']['value']
                            last-=floor(last)  # days
                            last*=24.0  # hours
                            ha=last-rah  # hours
                            har[itim]=ha*2.0*pi/24.0  # radians
                    
                        parang=mypl.arctan2( (cos(latr)*mypl.sin(har)),
                                             (sin(latr)*cos(decr)-cos(latr)*sin(decr)*mypl.cos(har)) )
                        # correct for cycle at +/-180.
                        #   (makes values crossing +/-180 deg all positive, and thus continuous)
                        #   (hmm, this will still fail at inferior circumpolar transit, but that is very rare)
                        if (latr<decr):
                            parang[parang<0.0]+=(2*pi)

                        parang+=rang[iant,ispw]
                        parang+=(paoffset*pi/180.)       # manual feed pa offset

                        # save parang values for this antenna
                        parangbyant[iant]=parang

                        # Escape if insufficient samples
                        nsamp=nrows-mypl.sum(flags)
                        if (nsamp<3):
                            antok[iant]=False
                            casalog.post('  Ant='+str(iant)+' has insuffiencient sampling: nsamp='+
                                         str(nsamp)+' < 3')
                            st.close()
                            continue

                        # Check parang coverage
                        dparang=abs(parang[~flags].max()-parang[~flags].min())  # rad
                        if dparang<minpacovR:
                            antok[iant]=False
                            casalog.post('  Ant='+str(iant)+' has insuffiencient parang cov: '+
                                         str(round(dparang*180/pi,2))+' < '+str(minpacov)+'deg')
                            continue



                        # indep var matrix
                        A=mypl.ones((nrows,3))
                        A[:,1]=mypl.cos(2*parang)
                        A[:,2]=mypl.sin(2*parang)
                        A[flags,:]=0.0  # zero flagged rows
                    
                        # squared gain amplitude ratio
                        amps=mypl.absolute(gains)
                        amps[amps==0.0]=1.0
                        gratio2=mypl.square(amps[0,0,:]/amps[1,0,:])
                        gratio2[flags]=0.0  # zero flagged samples
                
                        fit=mypl.lstsq(A,gratio2,rcond=None)

                        r2=fit[0][0]
                        if r2<0.0:
                            casalog.post('  Ant='+str(iant)+' yielded an unphysical solution; skipping.')
                            continue

                        # Reaching here, we have a nominally good solution
                        antok[iant]=True;

                        q[iant]=fit[0][1]/r2/2.0
                        u[iant]=fit[0][2]/r2/2.0
                        p=sqrt(q[iant]**2+u[iant]**2)
                        x=0.5*atan2(u[iant],q[iant])*180/pi

                        casalog.post('  Ant='+str(iant)+
                                     ' (PA offset='+str(round(rang[iant,ispw]*180/pi+paoffset,2))+'deg)'+
                                     ' q='+str(round(q[iant],4))+' u='+str(round(u[iant],4))+' p='+str(round(p,4))+' x='+str(round(x,3))+
                                     ' Gx/Gy='+str(round(sqrt(r2),4))+
                                     ' (parang cov='+str(round(dparang*180/pi,1))+'deg; nsamp='+str(nsamp))

                        if rempol:
                            if p<1.0:
                                Qpsi=q[iant]*mypl.cos(2*parang) + u[iant]*mypl.sin(2*parang)
                                gains[0,0,:]/=mypl.sqrt(1.0+Qpsi)
                                gains[1,0,:]/=mypl.sqrt(1.0-Qpsi)
                                st.putcol('CPARAM',gains)
                            else:
                                st.close()
                                raise RuntimeError('Spurious fractional polarization!')
                    st.close()

                nantok=mypl.sum(antok)
                if nantok==0:
                    casalog.post('Found no good polarization solutions for Fld='+fldnames[ifld]+' Spw='+str(ispw),'WARN')
                    mask[ispw,ifld]=False
                else:
                    Q[ispw,ifld]=mypl.sum(q)/nantok
                    U[ispw,ifld]=mypl.sum(u)/nantok
                    R[ispw,ifld]=mypl.sum(r)/nantok
                    mask[ispw,ifld]=True
                
                    P=sqrt(Q[ispw,ifld]**2+U[ispw,ifld]**2)
                    X=0.5*atan2(U[ispw,ifld],Q[ispw,ifld])*180/pi

                    casalog.post('  Ant=<*>  '+
                                 ' Q='+str(round(Q[ispw,ifld],4))+
                                 ' U='+str(round(U[ispw,ifld],4))+
                                 ' P='+str(round(P,4))+' X='+str(round(X,3)))

                    IQUV[fldnames[ifld]]['Spw'+str(ispw)]=[1.0,Q[ispw,ifld],U[ispw,ifld],0.0]

                    # if required, remove ant-averaged polarization from non-ok antennas (if any)
                    if rempol and nantok<nant:
                        badantlist=[i for i in range(nant) if not antok[i]]
                        casalog.post('  (Correcting undersampled antennas ('+str(badantlist)+
                                     ') with <*> solution.)')                                    
                        for iant in badantlist:
                            if iant in parangbyant.keys():
                                qstring='FIELD_ID=='+str(ifld)+' && SPECTRAL_WINDOW_ID=='+str(ispw)+' && ANTENNA1=='+str(iant)
                                st=mytb.query(query=qstring)
                                if st.nrows()>0 and P<1.0:
                                    gains=st.getcol('CPARAM')
                                    parang=parangbyant[iant]
                                    Qpsi=Q[ispw,ifld]*mypl.cos(2*parang) + U[ispw,ifld]*mypl.sin(2*parang)
                                    gains[0,0,:]/=mypl.sqrt(1.0+Qpsi)
                                    gains[1,0,:]/=mypl.sqrt(1.0-Qpsi)
                                    st.putcol('CPARAM',gains)
                                st.close()


            if sum(mask[:,ifld])>0:
                casalog.post('For field='+fldnames[ifld]+' there are '+str(sum(mask[:,ifld]))+' good spws.')
                Qm=mypl.mean(Q[mask[:,ifld],ifld])
                Um=mypl.mean(U[mask[:,ifld],ifld])
                IQUV[fldnames[ifld]]['SpwAve']=[1.0,Qm,Um,0.0]
                Qe=mypl.std(Q[mask[:,ifld],ifld])
                Ue=mypl.std(U[mask[:,ifld],ifld])
                Pm=sqrt(Qm**2+Um**2)
                Xm=0.5*atan2(Um,Qm)*180/pi
                casalog.post('Spw mean: Fld='+fldnames[ifld]+' Q='+str(round(Qm,4))+' U='+str(round(Um,4))+' P='+str(round(Pm,4))+' X='+str(round(Xm,3)))
            else:
                casalog.post('Found no good polarization solutions for Fld='+fldnames[ifld]+' in any spw.','WARN')

        mytb.close()

        casalog.post("NB: Returning dictionary containing fractional Stokes results.")
        return IQUV

    finally:
        mymd.close()
        myme.done()
        
