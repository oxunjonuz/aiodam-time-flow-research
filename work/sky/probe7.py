import numpy as np, math
from astropy.io import fits
from astropy_healpix import HEALPix
import astropy.units as u
h=fits.open('LALInference_skymap.fits')
prob=np.asarray(h[1].data['PROB']).reshape(-1).astype(float)
nside=512
hp=HEALPix(nside=nside, order='nested', frame='icrs')
pix=np.arange(len(prob))
lon,lat=hp.healpix_to_lonlat(pix)
ra=lon.to_value(u.rad); dec=lat.to_value(u.rad)
n=np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)],axis=1)
C=299792458.0
VERT={"H1":np.array([-2.16141492636e6,-3.83469517889e6,4.60035022664e6]),
      "L1":np.array([-7.42760447238e4,-5.49628371971e6,3.22425701744e6])}
dH=(VERT["H1"]@n.T)/C*1e3
dL=(VERT["L1"]@n.T)/C*1e3
print("weighted mean RA %.2f Dec %.2f"%(np.degrees((ra*prob).sum()/prob.sum()), np.degrees((dec*prob).sum()/prob.sum())))
for name,delay in (("H1-L1", (VERT["H1"]-VERT["L1"])@n.T/C*1e3),
                   ("L1-H1", (VERT["L1"]-VERT["H1"])@n.T/C*1e3)):
    b=np.arange(-10.5,10.5,0.5)
    hist,_=np.histogram(delay,bins=b,weights=prob)
    print(name,"max bin prob %.4f at %.2f ms"%(hist.max(), b[np.argmax(hist)]))
    print("  top5 bins:", [(round(b[i],1), round(hist[i],4)) for i in np.argsort(hist)[::-1][:5]])
