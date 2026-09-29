import numpy as np
from astropy.io import fits
from astropy_healpix import HEALPix
import astropy.units as u
h=fits.open('LALInference_skymap.fits')
raw=np.asarray(h[1].data['PROB']).astype(float)
C=299792458.0
H1=np.array([-2.16141492636e6,-3.83469517889e6,4.60035022664e6])
L1=np.array([-7.42760447238e4,-5.49628371971e6,3.22425701744e6])
pix=np.arange(3145728)
def delays(order):
    hp=HEALPix(nside=512, order=order, frame='icrs')
    lon,lat=hp.healpix_to_lonlat(pix)
    ra=lon.to_value(u.rad); dec=lat.to_value(u.rad)
    n=np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)],axis=1)
    return (H1-L1)@n.T/C*1e3
for oname in ("nested","ring"):
    dt=delays(oname)
    for name,prob in (("row",raw.reshape(-1)),("col",raw.T.reshape(-1))):
        order=np.argsort(prob)[::-1]
        cum=np.cumsum(prob[order]); n90=int(np.searchsorted(cum,0.9))+1
        sel=order[:n90]
        w=prob[sel]/prob[sel].sum()
        mu=(dt[sel]*w).sum(); sd=np.sqrt(((dt[sel]-mu)**2*w).sum())
        print(f"{oname:7s} {name:4s}  90%%region delay mean={mu:+.3f} std={sd:.3f} ms  ML={dt[order[0]]:+.3f}")
