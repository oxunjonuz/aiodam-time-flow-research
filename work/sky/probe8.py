import numpy as np, math
from astropy.io import fits
from astropy_healpix import HEALPix
import astropy.units as u
h=fits.open('LALInference_skymap.fits')
prob=np.asarray(h[1].data['PROB']).reshape(-1).astype(float)
hp=HEALPix(nside=512, order='nested', frame='icrs')
pix=np.arange(len(prob))
lon,lat=hp.healpix_to_lonlat(pix)
ra=lon.to_value(u.rad); dec=lat.to_value(u.rad)
n=np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)],axis=1)
C=299792458.0
H1=np.array([-2.16141492636e6,-3.83469517889e6,4.60035022664e6])
L1=np.array([-7.42760447238e4,-5.49628371971e6,3.22425701744e6])
dt=(H1-L1)@n.T/C*1e3   # arrival_H1 - arrival_L1 in ms
b=np.arange(-10.25,10.25,0.5)
hist,_=np.histogram(dt,bins=b,weights=prob)
tot=prob.sum()
print("delay histogram (arrival_H1-arrival_L1), prob per 0.5ms bin:")
for i in range(len(hist)):
    if hist[i]>0.005:
        print("  %+6.2f ms : %.4f"%(b[i]+0.25, hist[i]/tot*1.0 if False else hist[i]))
print("sum",hist.sum())
# where is the peak of the ring?
k=np.argmax(prob)
print("ML pixel delay %.3f ms, RA %.2f Dec %.2f"%(dt[k], np.degrees(ra[k]), np.degrees(dec[k])))
